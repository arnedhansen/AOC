%% AOC EEG ERD/ERS N Back time course plus topos
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
tfr1_all = {};
tfr2_all = {};
tfr3_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tfr = fullfile(dp, 'tfr_nback.mat');
    T = load(f_tfr, 'tfr1', 'tfr2', 'tfr3');
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
fontSize = 30;
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
ylabel('Power [dB]');
xlim([-0.5 2]);
ylim([-3.25 0.65]);
leg_p3 = patch(NaN, NaN, colors(3, :), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(2, :), 'EdgeColor', 'none');
leg_p1 = patch(NaN, NaN, colors(1, :), 'EdgeColor', 'none');
legend([leg_p1, leg_p2, leg_p3], {'1-back', '2-back', '3-back'}, ...
    'Location', 'southeast', 'FontSize', 25, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_timecourse.png'));

%% Topoplots (retention [0 2] s)
xlimTopo = [0 2];

cmap = interp1(linspace(0, 1, 5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0, 1, 64));

topoData = {ga_nb_1_erd, ga_nb_2_erd, ga_nb_3_erd};
topoTitles = {'1-back', '2-back', '3-back'};

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.zlim = [-2 2];
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
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos.png'));
