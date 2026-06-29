%% AOC EEG ERD/ERS N Back time course plus topos
% Trial-averaged TFR is baselined [-1.5 -0.5] dB per subject, then grand-averaged

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end
subjects = setdiff(subjects, {'361'});

baseline_window = [-1.5 -0.5];
cfgb = [];
cfgb.baseline = baseline_window;
cfgb.baselinetype = 'db';

%% Load trial-averaged TFR, baseline per subject
tfr1_all = {};
tfr2_all = {};
tfr3_all = {};
tfr1_raw_all = {};
tfr2_raw_all = {};
tfr3_raw_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tfr = fullfile(dp, 'tfr_nback.mat');
    T = load(f_tfr, 'tfr1', 'tfr2', 'tfr3');
    tfr1_raw_all{end + 1} = T.tfr1;
    tfr2_raw_all{end + 1} = T.tfr2;
    tfr3_raw_all{end + 1} = T.tfr3;
    tfr1_all{end + 1} = ft_freqbaseline(cfgb, T.tfr1);
    tfr2_all{end + 1} = ft_freqbaseline(cfgb, T.tfr2);
    tfr3_all{end + 1} = ft_freqbaseline(cfgb, T.tfr3);
    fprintf('Loaded n back TFR %d / %d\n', subj, numel(subjects));
end

nSubj = numel(tfr1_all);
ref = tfr1_all{1};
chPlot = {};
for i = 1:numel(ref.label)
    label = ref.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        chPlot{end + 1} = label;
    end
end

%% Grand average (keep individual subjects) and select alpha band
cfg = [];
cfg.keepindividual = 'yes';
ga_nb_1 = ft_freqgrandaverage(cfg, tfr1_all{:});
ga_nb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_nb_3 = ft_freqgrandaverage(cfg, tfr3_all{:});

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_nb_1_erd = ft_selectdata(cfg, ga_nb_1);
ga_nb_2_erd = ft_selectdata(cfg, ga_nb_2);
ga_nb_3_erd = ft_selectdata(cfg, ga_nb_3);

cfg_ga = [];
cfg_ga.keepindividual = 'yes';
ga_nb_1_raw = ft_freqgrandaverage(cfg_ga, tfr1_raw_all{:});
ga_nb_2_raw = ft_freqgrandaverage(cfg_ga, tfr2_raw_all{:});
ga_nb_3_raw = ft_freqgrandaverage(cfg_ga, tfr3_raw_all{:});
ga_nb_1_raw_alpha = ft_selectdata(cfg, ga_nb_1_raw);
ga_nb_2_raw_alpha = ft_selectdata(cfg, ga_nb_2_raw);
ga_nb_3_raw_alpha = ft_selectdata(cfg, ga_nb_3_raw);

%% Plot ERSD time course (occipital ROI, mean +/- SEM across subjects)
cfg = [];
cfg.channel = chPlot;
cfg.avgoverchan = 'yes';
cfg.latency = [-.5 2];
tlk1 = ft_selectdata(cfg, ga_nb_1_erd);
tlk2 = ft_selectdata(cfg, ga_nb_2_erd);
tlk3 = ft_selectdata(cfg, ga_nb_3_erd);
timeVec = tlk1.time(:)';

close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 50;
legendFontSize = fontSize * 0.666;
mask = timeVec >= -0.5 & timeVec <= 2;
x = timeVec(mask);

tc1 = squeeze(tlk1.powspctrm);
tc2 = squeeze(tlk2.powspctrm);
tc3 = squeeze(tlk3.powspctrm);
if isvector(tc1)
    tc1 = tc1(:)';
    tc2 = tc2(:)';
    tc3 = tc3(:)';
end

y3 = mean(tc3(:, mask), 1, 'omitnan');
e3 = std(tc3(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y2 = mean(tc2(:, mask), 1, 'omitnan');
e2 = std(tc2(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y1 = mean(tc1(:, mask), 1, 'omitnan');
e1 = std(tc1(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);

hold on;
eb3 = shadedErrorBar(x, y3, e3, 'lineProps', {'-', 'Color', colors(3, :)});
eb2 = shadedErrorBar(x, y2, e2, 'lineProps', {'-', 'Color', colors(2, :)});
eb1 = shadedErrorBar(x, y1, e1, 'lineProps', {'-', 'Color', colors(1, :)});
ebs = [eb3, eb2, eb1];
cIdx = [3, 2, 1];
for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', 2, 'Color', colors(cIdx(k), :));
    set(ebs(k).patch, 'FaceColor', colors(cIdx(k), :), 'FaceAlpha', 0.20);
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
leg_p1 = patch(NaN, NaN, colors(1, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(1, :), 'LineWidth', 1.5);
leg_p2 = patch(NaN, NaN, colors(2, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(2, :), 'LineWidth', 1.5);
leg_p3 = patch(NaN, NaN, colors(3, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(3, :), 'LineWidth', 1.5);
legend([leg_p1, leg_p2, leg_p3], {'1-back', '2-back', '3-back'}, ...
    'Location', 'southeast', 'FontSize', legendFontSize, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_timecourse.png'));

%% Topoplots (retention [0 2] s)
xlimTopo = [0 2];
fontSize = 40;

topoData = {ga_nb_1_erd, ga_nb_2_erd, ga_nb_3_erd};
topoTitles = {'1-back', '2-back', '3-back'};

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
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos.png'));

%% Topoplots (baseline window)
close all
xlimTopo = baseline_window;

topoData = {ga_nb_1_raw_alpha, ga_nb_2_raw_alpha, ga_nb_3_raw_alpha};
topoTitles = {'1-back', '2-back', '3-back'};

cfg_sel = [];
cfg_sel.channel = chPlot;
cfg_sel.latency = xlimTopo;
cfg_sel.avgoverchan = 'yes';
cfg_sel.avgovertime = 'yes';
bl_vals = [];
for k = 1:3
    tmp = ft_selectdata(cfg_sel, topoData{k});
    bl_vals = [bl_vals; tmp.powspctrm(:)]; %#ok<AGROW>
end
bl_vals = bl_vals(isfinite(bl_vals));
global_max = 4%prctile(bl_vals, 99);

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
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos_baseline_window.png'));