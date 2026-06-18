%% AOC EEG ERD/ERS Sternberg time course plus topos
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
tc2 = [];
tc4 = [];
tc6 = [];
timeVec = [];
tfr2_all = {};
tfr4_all = {};
tfr6_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tc = fullfile(dp, 'ersd_sternberg_timecourse.mat');
    f_tfr = fullfile(dp, 'tfr_stern.mat');
    S = load(f_tc, 'ersd_timecourse');
    E = S.ersd_timecourse;
    conds = E.condition(:);
    if isempty(timeVec)
        timeVec = E.time(:)';
    elseif numel(timeVec) ~= numel(E.time) || any(abs(timeVec - E.time(:)') > 1e-12)
        warning('Time vector mismatch in subject %s', subjects{subj});
        continue
    end

    tcMat = E.ersd_occ_8_14_db;
    tc2 = [tc2; tcMat(conds == 2, :)];
    tc4 = [tc4; tcMat(conds == 4, :)];
    tc6 = [tc6; tcMat(conds == 6, :)];

    T = load(f_tfr, 'tfr2_bl', 'tfr4_bl', 'tfr6_bl');
    tfr2_all{end + 1} = T.tfr2_bl;
    tfr4_all{end + 1} = T.tfr4_bl;
    tfr6_all{end + 1} = T.tfr6_bl;
    fprintf('Loaded Sternberg ERSD timecourse %d / %d\n', subj, numel(subjects));
end

nSubj = size(tc2, 1);
ref = tfr2_all{1};
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
mask = timeVec >= -0.5 & timeVec <= 3;
x = timeVec(mask);

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
ylim([-3.25 0.65]);
leg_p6 = patch(NaN, NaN, colors(3, :), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2, :), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(1, :), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'Location', 'southeast', 'FontSize', 25, 'Box', 'off');
drawnow; pause(0.05);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_timecourse.png'));

%% Topoplots from cached baselined TFR
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_sb_4 = ft_freqgrandaverage(cfg, tfr4_all{:});
ga_sb_6 = ft_freqgrandaverage(cfg, tfr6_all{:});

xlimTopo = [1 2];

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_sb_2_erd = ft_selectdata(cfg, ga_sb_2);
ga_sb_4_erd = ft_selectdata(cfg, ga_sb_4);
ga_sb_6_erd = ft_selectdata(cfg, ga_sb_6);

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
