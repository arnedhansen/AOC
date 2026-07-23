%% AOC EEG ERD/ERS Sternberg time course plus topos
% Trial-averaged TFR is baselined [-1.5 -0.5] dB per subject.
% ERSD band matches stats/feature extraction: condition-wise IAF-4/+2 Hz,
% with fallback to fixed [8 14] Hz when IAF is undefined. Band selection is
% applied per subject x condition before grand-averaging (do not GA first).

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end
subjects = setdiff(subjects, {'361'});

baseline_window = [-1.5 -0.5];
alphaRange = [8 14];
cfgb = [];
cfgb.baseline = baseline_window;
cfgb.baselinetype = 'db';

%% Load TFR, baseline, collapse to ERSD band per subject x condition
tfr2_erd_all = {};
tfr4_erd_all = {};
tfr6_erd_all = {};
tfr2_raw_alpha_all = {};
tfr4_raw_alpha_all = {};
tfr6_raw_alpha_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tfr = fullfile(dp, 'tfr_stern.mat');
    f_iaf = fullfile(dp, 'IAF_sternberg.mat');
    T = load(f_tfr, 'tfr2', 'tfr4', 'tfr6');
    I = load(f_iaf, 'IAF2', 'IAF4', 'IAF6');

    tfr2_bl = ft_freqbaseline(cfgb, T.tfr2);
    tfr4_bl = ft_freqbaseline(cfgb, T.tfr4);
    tfr6_bl = ft_freqbaseline(cfgb, T.tfr6);

    % ERSD: IAF band with [8 14] fallback (matches AOC_eeg_fex_sternberg)
    tfr2_erd_all{end + 1} = select_ersd_band(tfr2_bl, I.IAF2, alphaRange, true); %#ok<AGROW>
    tfr4_erd_all{end + 1} = select_ersd_band(tfr4_bl, I.IAF4, alphaRange, true); %#ok<AGROW>
    tfr6_erd_all{end + 1} = select_ersd_band(tfr6_bl, I.IAF6, alphaRange, true); %#ok<AGROW>

    % Baseline-window raw alpha: IAF band only (no fallback; skip if IAF missing)
    raw2 = select_ersd_band(T.tfr2, I.IAF2, alphaRange, false);
    raw4 = select_ersd_band(T.tfr4, I.IAF4, alphaRange, false);
    raw6 = select_ersd_band(T.tfr6, I.IAF6, alphaRange, false);
    if ~isempty(raw2), tfr2_raw_alpha_all{end + 1} = raw2; end %#ok<AGROW>
    if ~isempty(raw4), tfr4_raw_alpha_all{end + 1} = raw4; end %#ok<AGROW>
    if ~isempty(raw6), tfr6_raw_alpha_all{end + 1} = raw6; end %#ok<AGROW>

    fprintf('Loaded Sternberg TFR %d / %d\n', subj, numel(subjects));
end

nSubj = numel(tfr2_erd_all);
ref = tfr2_erd_all{1};
chPlot = {};
for i = 1:numel(ref.label)
    label = ref.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        chPlot{end + 1} = label;
    end
end

%% Grand average (keep individual subjects); frequency already collapsed
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2_erd = ft_freqgrandaverage(cfg, tfr2_erd_all{:});
ga_sb_4_erd = ft_freqgrandaverage(cfg, tfr4_erd_all{:});
ga_sb_6_erd = ft_freqgrandaverage(cfg, tfr6_erd_all{:});

ga_sb_2_raw_alpha = ft_freqgrandaverage(cfg, tfr2_raw_alpha_all{:});
ga_sb_4_raw_alpha = ft_freqgrandaverage(cfg, tfr4_raw_alpha_all{:});
ga_sb_6_raw_alpha = ft_freqgrandaverage(cfg, tfr6_raw_alpha_all{:});

%% Plot ERSD time course (occipital ROI, mean +/- SEM across subjects)
cfg = [];
cfg.channel = chPlot;
cfg.avgoverchan = 'yes';
cfg.latency = [-.5 3];
tlk2 = ft_selectdata(cfg, ga_sb_2_erd);
tlk4 = ft_selectdata(cfg, ga_sb_4_erd);
tlk6 = ft_selectdata(cfg, ga_sb_6_erd);
timeVec = tlk2.time(:)';

close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 50;
legendFontSize = fontSize * 0.666;
mask = timeVec >= -.5 & timeVec <= 2;
x = timeVec(mask);

tc2 = squeeze(tlk2.powspctrm);
tc4 = squeeze(tlk4.powspctrm);
tc6 = squeeze(tlk6.powspctrm);
if isvector(tc2)
    tc2 = tc2(:)';
    tc4 = tc4(:)';
    tc6 = tc6(:)';
end

y6 = mean(tc6(:, mask), 1, 'omitnan');
e6 = std(tc6(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y4 = mean(tc4(:, mask), 1, 'omitnan');
e4 = std(tc4(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y2 = mean(tc2(:, mask), 1, 'omitnan');
e2 = std(tc2(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);

hold on;
eb6 = shadedErrorBar(x, y6, e6, 'lineProps', {'-', 'Color', colors(3, :)});
eb4 = shadedErrorBar(x, y4, e4, 'lineProps', {'-', 'Color', colors(2, :)});
eb2 = shadedErrorBar(x, y2, e2, 'lineProps', {'-', 'Color', colors(1, :)});
ebs = [eb6, eb4, eb2];
cIdx = [3, 2, 1];
for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', 5, 'Color', colors(cIdx(k), :));
    set(ebs(k).patch, 'FaceColor', colors(cIdx(k), :), 'FaceAlpha', 0.25);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

set(gca, 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Power Change [dB]');
xlim([-0.5 2]);
ylim([-3 0.5]);
leg_p2 = patch(NaN, NaN, colors(1, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(1, :), 'LineWidth', 1.5);
leg_p4 = patch(NaN, NaN, colors(2, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(2, :), 'LineWidth', 1.5);
leg_p6 = patch(NaN, NaN, colors(3, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(3, :), 'LineWidth', 1.5);
legend([leg_p2, leg_p4, leg_p6], {' WM load 2', ' WM load 4', ' WM load 6'}, ...
    'Location', 'southeast', 'FontSize', legendFontSize, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_timecourse.png'));

%% Topoplots (late retention [1 2] s)
close all
xlimTopo = [1 2];
fontSize = 40;

topoData = {ga_sb_2_erd, ga_sb_4_erd, ga_sb_6_erd};
topoTitles = {'WM load 2', 'WM load 4', 'WM load 6'};

cfg_zlim = [];
cfg_zlim.latency = xlimTopo;
cfg_zlim.avgovertime = 'yes';
zlimMax = 0;
for k = 1:3
    tmp = ft_selectdata(cfg_zlim, topoData{k});
    dat = tmp.powspctrm;
    if size(dat, 1) > 1 && numel(tmp.label) > 1
        dat = mean(dat, 1, 'omitnan');
    end
    dat = dat(isfinite(dat(:)));
    zlimMax = max(zlimMax, max(abs(dat(:))));
end
zlimTopo = [-zlimMax zlimMax];

cmap = interp1(linspace(0, 1, 5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0, 1, 64));

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.zlim = zlimTopo;
cfg.xlim = xlimTopo;
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';
cfg.gridscale = 300;
cfg.colormap = cmap;

figure('Position', [0 0 1512 982*0.6], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
topoAxes = gobjects(1, 3);
for k = 1:3
    ax = nexttile(tl, k);
    topoAxes(k) = ax;
    cfg.figure = ax;
    ft_topoplotER(cfg, topoData{k});
    set(ax, 'FontSize', fontSize);
    title(ax, topoTitles{k}, 'Interpreter', 'none', 'FontSize', fontSize);
end
set(topoAxes, 'XLim', xlim(topoAxes(1)), 'YLim', ylim(topoAxes(1)));
for k = 1:3
    axis(topoAxes(k), 'equal');
    axis(topoAxes(k), 'off');
end
cb = colorbar(topoAxes(end));
cb.Layout.Tile = 'east';
cb.Label.String = 'Power Change [dB]';
cb.FontSize = fontSize;
cb.Label.FontSize = fontSize;
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_topos.png'));

%% Topoplots (baseline window)
close all
xlimTopo = baseline_window;

topoData = {ga_sb_2_raw_alpha, ga_sb_4_raw_alpha, ga_sb_6_raw_alpha};
topoTitles = {'WM load 2', 'WM load 4', 'WM load 6'};

cfg_sel = [];
cfg_sel.channel = chPlot;
cfg_sel.latency = xlimTopo;
cfg_sel.avgoverchan = 'yes';
cfg_sel.avgovertime = 'yes';
bl_vals = [];
for k = 1:3
    tmp = ft_selectdata(cfg_sel, topoData{k});
    bl_vals = [bl_vals; tmp.powspctrm(:)];
end
bl_vals = bl_vals(isfinite(bl_vals));
global_max = 5%prctile(bl_vals, 99);

cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.zlim = [global_max / 3 global_max];
cfg.xlim = xlimTopo;
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';
cfg.gridscale = 300;
cfg.colormap = cmap;

figure('Position', [0 0 1512 982*0.6], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
topoAxes = gobjects(1, 3);
for k = 1:3
    ax = nexttile(tl, k);
    topoAxes(k) = ax;
    cfg.figure = ax;
    ft_topoplotER(cfg, topoData{k});
    set(ax, 'FontSize', fontSize);
    title(ax, topoTitles{k}, 'Interpreter', 'none', 'FontSize', fontSize);
end
set(topoAxes, 'XLim', xlim(topoAxes(1)), 'YLim', ylim(topoAxes(1)));
for k = 1:3
    axis(topoAxes(k), 'equal');
    axis(topoAxes(k), 'off');
end
cb = colorbar(topoAxes(end));
cb.Layout.Tile = 'east';
cb.Label.String = 'Power [\muV^2/Hz]';
cb.FontSize = fontSize;
cb.Label.FontSize = fontSize;
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_topos_baseline_window.png'));

%% Local helpers (match 2_feature_extraction/eeg/AOC_eeg_fex_sternberg.m)
function out = select_ersd_band(tfr, iaf, alphaRange, allowFallback)
band = ersd_alpha_band(iaf, alphaRange, allowFallback);
if isempty(band)
    out = [];
    return
end
cfg = [];
cfg.frequency = band;
cfg.avgoverfreq = 'yes';
out = ft_selectdata(cfg, tfr);
end

function band = ersd_alpha_band(IAF, alphaRange, allowFallback)
% ERSD: (IAF-4, IAF+2) when IAF is valid; else fixed alphaRange if allowFallback.
% AlphaPower-style (allowFallback=false): empty band when IAF missing.
if ~isfinite(IAF)
    if allowFallback
        band = alphaRange;
    else
        band = [];
    end
    return
end
band = [IAF - 4, IAF + 2];
if any(~isfinite(band)) || band(1) >= band(2)
    if allowFallback
        band = alphaRange;
    else
        band = [];
    end
end
end
