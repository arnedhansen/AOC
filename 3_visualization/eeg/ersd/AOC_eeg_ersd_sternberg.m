%% AOC EEG ERD/ERS — Sternberg (time course + topos)
% Loads baselined TFR per subject, grand-averages with individuals, plots fixed-band
% [8 14] Hz occipital time course (mean +- SEM) and topoplots averaged over [1 2]s.
% Occipital channels: labels containing O or I (same rule as AOC_eeg_fex_* IAF ROI).
%
% Key outputs:
%   Figures under paths.figures/eeg/ersd

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end

subjects = setdiff(subjects, {'361'});

%% Load TFR
tfr2_all = cell(1, numel(subjects));
tfr4_all = cell(1, numel(subjects));
tfr6_all = cell(1, numel(subjects));
for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    fmat = fullfile(dp, 'tfr_stern.mat');
    if ~isfile(fmat)
        warning('Missing %s', fmat);
        continue
    end
    cd(dp);
    load tfr_stern tfr2_bl tfr4_bl tfr6_bl
    tfr2_all{subj} = tfr2_bl;
    tfr4_all{subj} = tfr4_bl;
    tfr6_all{subj} = tfr6_bl;
    fprintf('Loaded Sternberg TFR %d / %d\n', subj, numel(subjects));
end

nonEmpty = ~cellfun(@isempty, tfr2_all);
tfr2_all = tfr2_all(nonEmpty);
tfr4_all = tfr4_all(nonEmpty);
tfr6_all = tfr6_all(nonEmpty);
nSubj = numel(tfr2_all);
if nSubj == 0
    error('No tfr_stern.mat found for any subject. Run AOC_eeg_fex_sternberg_TFR.m.');
end

ref = tfr2_all{1};
chPlot = {};
for i = 1:numel(ref.label)
    label = ref.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        chPlot{end+1} = label;
    end
end
if isempty(chPlot)
    error('No occipital (O/I) channels present in TFR labels.');
end

%% Grand average with individuals
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

%% Time course figure (select [-0.5 3]s, display xlim [-0.5 2])
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 22;

cfgT = [];
cfgT.channel = chPlot;
cfgT.avgoverchan = 'yes';
cfgT.latency = [-0.5 3];

tlk2 = ft_selectdata(cfgT, ga_sb_2_erd);
tlk4 = ft_selectdata(cfgT, ga_sb_4_erd);
tlk6 = ft_selectdata(cfgT, ga_sb_6_erd);

plot_erds_line_sem(tlk6, [0.97, 0.26, 0.26], nSubj);
hold on;
plot_erds_line_sem(tlk4, [0.30, 0.75, 0.93], nSubj);
plot_erds_line_sem(tlk2, [0 0 0], nSubj);

set(gca, 'FontSize', fontSize);
title('Sternberg task');
xlabel('Time [s]');
ylabel('Power change [dB]');
xlim([-0.5 2]);
ylim([-3 3]);
grid on;
box on;
legend({'WM 6', 'WM 4', 'WM 2'}, 'Location', 'northeast', 'FontSize', 18);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_timecourse.png'));

%% Topoplots [1 2]s
cfg = [];
cfg.layout = headmodel.layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [1 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';

figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1, 3, 1);
ft_topoplotER(cfg, ga_sb_2_erd);
subplot(1, 3, 2);
ft_topoplotER(cfg, ga_sb_4_erd);
subplot(1, 3, 3);
ft_topoplotER(cfg, ga_sb_6_erd);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_topos.png'));

function plot_erds_line_sem(tlk, col, nSubj)
    x = tlk.time(:);
    ps = squeeze(tlk.powspctrm);
    if isvector(ps)
        y = ps(:);
        e = zeros(size(y));
    else
        y = mean(ps, 1)';
        e = std(ps, 0, 1)' ./ sqrt(nSubj);
    end
    low = y - e;
    high = y + e;
    hp = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], col, 'HandleVisibility', 'off');
    hold on;
    hl = line(x, y);
    set(hp, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    set(hl, 'Color', col, 'LineWidth', 2);
end
