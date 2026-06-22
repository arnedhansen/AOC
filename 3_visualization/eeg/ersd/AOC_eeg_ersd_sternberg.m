%% AOC EEG ERD/ERS Sternberg time course plus topos
% Trial-averaged TFR is baselined [-1.5 -0.5] dB per subject, then grand-averaged
% (matches supervisor occ_alpha_power_script.m).

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
tfr2_all = {};
tfr4_all = {};
tfr6_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tfr = fullfile(dp, 'tfr_stern.mat');
    T = load(f_tfr, 'tfr2', 'tfr4', 'tfr6');
    tfr2_all{end + 1} = ft_freqbaseline(cfgb, T.tfr2);
    tfr4_all{end + 1} = ft_freqbaseline(cfgb, T.tfr4);
    tfr6_all{end + 1} = ft_freqbaseline(cfgb, T.tfr6);
    fprintf('Loaded Sternberg TFR %d / %d\n', subj, numel(subjects));
end

nSubj = numel(tfr2_all);
ref = tfr2_all{1};
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
ga_sb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_sb_4 = ft_freqgrandaverage(cfg, tfr4_all{:});
ga_sb_6 = ft_freqgrandaverage(cfg, tfr6_all{:});

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_sb_2_erd = ft_selectdata(cfg, ga_sb_2);
ga_sb_4_erd = ft_selectdata(cfg, ga_sb_4);
ga_sb_6_erd = ft_selectdata(cfg, ga_sb_6);

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
fontSize = 30;
mask = timeVec >= -0.5 & timeVec <= 2;
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
    set(ebs(k).mainLine, 'LineWidth', 2, 'Color', colors(cIdx(k), :));
    set(ebs(k).patch, 'FaceColor', colors(cIdx(k), :), 'FaceAlpha', 0.20);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

set(gca, 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Power [dB]');
xlim([-0.5 2]);
ylim([-3 3]);
leg_p6 = patch(NaN, NaN, colors(3, :), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2, :), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(1, :), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'Location', 'southeast', 'FontSize', 25, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_timecourse.png'));

%% Topoplots (late retention [1 2] s)
xlimTopo = [1 2];

cmap = interp1(linspace(0, 1, 5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0, 1, 64));

topoData = {ga_sb_2_erd, ga_sb_4_erd, ga_sb_6_erd};
topoTitles = {'WM load 2', 'WM load 4', 'WM load 6'};

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.zlim = 'maxabs';
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
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
for k = 1:3
    ax = nexttile(k);
    cfg.figure = ax;
    ft_topoplotER(cfg, topoData{k});
    if k == 3
        cb = colorbar(ax, 'eastoutside');
        cb.Label.String = 'Power [dB]';
    end
    set(ax, 'FontSize', fontSize);
    title(ax, topoTitles{k}, 'Interpreter', 'none');
end
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_topos.png'));
