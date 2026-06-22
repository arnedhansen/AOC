%% AOC EEG ERD/ERS N Back time course plus topos
% Time course is loaded from saved ERSD outputs produced during feature extraction.
% No ERSD recomputation is performed here.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end
subjects = setdiff(subjects, {'361'});

%% Load subject ERSD timecourses and cached TFR for topoplots
tc1 = [];
tc2 = [];
tc3 = [];
timeVec = [];
tfr1_all = {};
tfr2_all = {};
tfr3_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tc = fullfile(dp, 'ersd_nback_timecourse.mat');
    f_tfr = fullfile(dp, 'tfr_nback.mat');
    S = load(f_tc, 'ersd_timecourse');
    E = S.ersd_timecourse;
    conds = E.condition(:);
    timeVec = E.time(:)';
    tcMat = E.ersd_occ_8_14_db;
    tc1 = [tc1; tcMat(conds == 1, :)];
    tc2 = [tc2; tcMat(conds == 2, :)];
    tc3 = [tc3; tcMat(conds == 3, :)];

    T = load(f_tfr, 'tfr1_bl', 'tfr2_bl', 'tfr3_bl');
    tfr1_all{end + 1} = T.tfr1_bl;
    tfr2_all{end + 1} = T.tfr2_bl;
    tfr3_all{end + 1} = T.tfr3_bl;
    fprintf('Loaded n back ERSD timecourse %d / %d\n', subj, numel(subjects));
end

nSubj = size(tc1, 1);
ref = tfr1_all{1};
chPlot = {};
for i = 1:numel(ref.label)
    label = ref.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        chPlot{end + 1} = label;
    end
end

%% Plot ERSD Timecourse 
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 30;
mask = timeVec >= -0.5 & timeVec <= 2;
x = timeVec(mask);

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

%% Topoplots
close all
cfg = [];
cfg.keepindividual = 'yes';
ga_nb_1 = ft_freqgrandaverage(cfg, tfr1_all{:});
ga_nb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_nb_3 = ft_freqgrandaverage(cfg, tfr3_all{:});

xlimTopo = [0 2];

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_nb_1_erd = ft_selectdata(cfg, ga_nb_1);
ga_nb_2_erd = ft_selectdata(cfg, ga_nb_2);
ga_nb_3_erd = ft_selectdata(cfg, ga_nb_3);

cmap = interp1(linspace(0, 1, 5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0, 1, 64));

topoData = {ga_nb_1_erd, ga_nb_2_erd, ga_nb_3_erd};
topoTitles = {'1-back', '2-back', '3-back'};

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
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos.png'));
