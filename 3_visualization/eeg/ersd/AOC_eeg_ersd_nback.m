%% AOC EEG ERD/ERS — N-Back (time course + topos)
% Loads baselined TFR per subject, grand-averages with individuals, plots fixed-band
% [8 14] Hz occipital time course (mean +- SEM) and topoplots averaged over [0 2]s.
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
tfr1_all = cell(1, numel(subjects));
tfr2_all = cell(1, numel(subjects));
tfr3_all = cell(1, numel(subjects));
for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    fmat = fullfile(dp, 'tfr_nback.mat');
    if ~isfile(fmat)
        warning('Missing %s', fmat);
        continue
    end
    cd(dp);
    load tfr_nback tfr1_bl tfr2_bl tfr3_bl
    tfr1_all{subj} = tfr1_bl;
    tfr2_all{subj} = tfr2_bl;
    tfr3_all{subj} = tfr3_bl;
    fprintf('Loaded n-back TFR %d / %d\n', subj, numel(subjects));
end

nonEmpty = ~cellfun(@isempty, tfr1_all);
tfr1_all = tfr1_all(nonEmpty);
tfr2_all = tfr2_all(nonEmpty);
tfr3_all = tfr3_all(nonEmpty);
nSubj = numel(tfr1_all);
if nSubj == 0
    error('No tfr_nback.mat found for any subject. Run AOC_eeg_fex_nback_TFR.m.');
end

ref = tfr1_all{1};
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
ga_nb_1 = ft_freqgrandaverage(cfg, tfr1_all{:});
ga_nb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_nb_3 = ft_freqgrandaverage(cfg, tfr3_all{:});

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_nb_1_erd = ft_selectdata(cfg, ga_nb_1);
ga_nb_2_erd = ft_selectdata(cfg, ga_nb_2);
ga_nb_3_erd = ft_selectdata(cfg, ga_nb_3);

%% Time course figure
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 22;

cfgT = [];
cfgT.channel = chPlot;
cfgT.avgoverchan = 'yes';
cfgT.latency = [-0.5 2];

tlk1 = ft_selectdata(cfgT, ga_nb_1_erd);
tlk2 = ft_selectdata(cfgT, ga_nb_2_erd);
tlk3 = ft_selectdata(cfgT, ga_nb_3_erd);

% 3-back (red), 2-back (cyan), 1-back (black)
plot_erds_line_sem(tlk3, [0.97, 0.26, 0.26], nSubj);
hold on;
plot_erds_line_sem(tlk2, [0.30, 0.75, 0.93], nSubj);
plot_erds_line_sem(tlk1, [0 0 0], nSubj);

set(gca, 'FontSize', fontSize);
title('N-back task');
xlabel('Time [s]');
ylabel('Power change [dB]');
xlim([-0.5 2]);
ylim([-3 3]);
grid on;
box on;
legend({'3-back', '2-back', '1-back'}, 'Location', 'northeast', 'FontSize', 18);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_timecourse.png'));

%% Topoplots [0 2]s
cfg = [];
cfg.layout = headmodel.layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [0 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';

figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1, 3, 1);
ft_topoplotER(cfg, ga_nb_1_erd);
subplot(1, 3, 2);
ft_topoplotER(cfg, ga_nb_2_erd);
subplot(1, 3, 3);
ft_topoplotER(cfg, ga_nb_3_erd);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos.png'));

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
