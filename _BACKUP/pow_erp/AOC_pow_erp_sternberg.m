%% Alpha Power as ERP plot for Sternberg data

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

alpha_pow2 = [];
alpha_pow4 = [];
alpha_pow6 = [];
alpha_pow8 = [];

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load data_sternberg
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind2=find(data.trialinfo==52);
    ind4=find(data.trialinfo==54);
    ind6=find(data.trialinfo==56);
    ind8=find(data.trialinfo==58);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    load2 = ft_freqanalysis(cfg,data);
    cfg.trials = ind4;
    load4 = ft_freqanalysis(cfg,data);
    cfg.trials = ind6;
    load6 = ft_freqanalysis(cfg,data);
    cfg.trials = ind8;
    load8 = ft_freqanalysis(cfg,data);

    %% Save data
    window_time = load2.time;
    cd(datapath)
    save erp_pow_sternberg load2 load4 load6 load8 window_time
    cd('/Volumes/methlab/Students/Arne/AOC/scripts/eeg')
    fprintf('Subject %.2d done \n', subj)
end

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    try
        load erp_pow_sternberg
    end
    try
        load erp_pow
    end
    cfg =[];
    cfg.baseline = [-Inf -.5];
    cfg.baselinetype ='db';
    load2 = ft_freqbaseline(cfg,load2);
    load4 = ft_freqbaseline(cfg,load4);
    load6 = ft_freqbaseline(cfg,load6);
    load8 = ft_freqbaseline(cfg,load8);
    l2{subj}= load2;
    l4{subj}= load4;
    l6{subj}= load6;
    l8{subj}= load8;
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr2 = ft_freqgrandaverage([],l2{:});
gatfr4 = ft_freqgrandaverage([],l4{:});
gatfr6 = ft_freqgrandaverage([],l6{:});
gatfr8 = ft_freqgrandaverage([],l8{:});
% gatfr.powspctrm =  freq data with 124 channels, 27 frequencybins and 121 timebins

%% Plot aggregated alpha power
addpath('/Users/Arne/Documents/matlabtools/shadedErrorBar')
close all

% Occipital channels
load('power_stern.mat');
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

% Alpha Definition
alpha_range = 8:14;
alpha_power_2 = [];
alpha_power_4 = [];
alpha_power_6 = [];
alpha_power_8 = [];
compute_mean_alpha = @(data, coi, alpha_range) ...
    mean(mean(data.powspctrm(ismember(data.label, coi), alpha_range, :), 2), 1);

for subj = 1:length(subjects)
    % Compute alpha power for each subject
    alpha_power_2(subj, :) = compute_mean_alpha(l2{subj}, channels, alpha_range);
    alpha_power_4(subj, :) = compute_mean_alpha(l4{subj}, channels, alpha_range);
    alpha_power_6(subj, :) = compute_mean_alpha(l6{subj}, channels, alpha_range);
    alpha_power_8(subj, :) = compute_mean_alpha(l8{subj}, channels, alpha_range);
end

% Compute grand average and SEM across subjects
mean_alpha_2 = mean(alpha_power_2, 1);
mean_alpha_4 = mean(alpha_power_4, 1);
mean_alpha_6 = mean(alpha_power_6, 1);
mean_alpha_8 = mean(alpha_power_8, 1);

% Compute SEM
compute_sem = @(data) std(data) / sqrt(size(data, 1));
sem_2 = compute_sem(squeeze(mean(gatfr2.powspctrm(ismember(gatfr2.label, channels), alpha_range, :), 2)));
sem_4 = compute_sem(squeeze(mean(gatfr4.powspctrm(ismember(gatfr4.label, channels), alpha_range, :), 2)));
sem_6 = compute_sem(squeeze(mean(gatfr6.powspctrm(ismember(gatfr6.label, channels), alpha_range, :), 2)));
sem_8 = compute_sem(squeeze(mean(gatfr8.powspctrm(ismember(gatfr8.label, channels), alpha_range, :), 2)));

figure;
set(gcf, "Position", [0, 0, 2000, 1200], "Color", 'w')
hold on;
shadedErrorBar(gatfr2.time, mean_alpha_2, sem_2, 'lineprops', {'r', 'LineWidth', 1});
shadedErrorBar(gatfr4.time, mean_alpha_4, sem_4, 'lineprops', {'g', 'LineWidth', 1});
shadedErrorBar(gatfr6.time, mean_alpha_6, sem_6, 'lineprops', {'b', 'LineWidth', 1});
shadedErrorBar(gatfr8.time, mean_alpha_8, sem_8, 'lineprops', {'k', 'LineWidth', 1});
hold off;

ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15; 
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Alpha Power [dB]', 'FontName', 'Arial', 'FontSize', 20);
legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 20);
title('Posterior Alpha Power during Sternberg', 'FontName', 'Arial', 'FontSize', 25);
% xlim([-0.5 2.5])

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/pow_erp/pow_erp_sternberg.png');

%% Plot for individual subjects 
% close all
% figure;
% set(gcf, "Position", [0, 0, 2000, 1200], "Color", 'w');
% 
% for subj = 1:length(subjects)
%     % Get data
%     load2 = l2{subj};
%     load4 = l4{subj};
%     load6 = l6{subj};
%     load8 = l8{subj};
% 
%     % Mean alpha power per subject
%     alpha_power_2 = squeeze(compute_mean_alpha(load2, channels, alpha_range))';
%     alpha_power_4 = squeeze(compute_mean_alpha(load4, channels, alpha_range))';
%     alpha_power_6 = squeeze(compute_mean_alpha(load6, channels, alpha_range))';
%     alpha_power_8 = squeeze(compute_mean_alpha(load8, channels, alpha_range))';
% 
%     % SEM for each subject
%     sem_2 = compute_sem(squeeze(mean(load2.powspctrm(ismember(load2.label, channels), alpha_range, :), 2)));
%     sem_4 = compute_sem(squeeze(mean(load4.powspctrm(ismember(load4.label, channels), alpha_range, :), 2)));
%     sem_6 = compute_sem(squeeze(mean(load6.powspctrm(ismember(load6.label, channels), alpha_range, :), 2)));
%     sem_8 = compute_sem(squeeze(mean(load8.powspctrm(ismember(load8.label, channels), alpha_range, :), 2)));
% 
%     % Plot
%     subplot(5, 2, subj);
%     hold on;
%     shadedErrorBar(load2.time, alpha_power_2, sem_2, 'lineprops', {'r', 'LineWidth', 1});
%     %shadedErrorBar(load4.time, alpha_power_4, sem_4, 'lineprops', {'g', 'LineWidth', 1});
%     shadedErrorBar(load6.time, alpha_power_6, sem_6, 'lineprops', {'b', 'LineWidth', 1});
%     %shadedErrorBar(load8.time, alpha_power_8, sem_8, 'lineprops', {'k', 'LineWidth', 1});
%     hold off;
% 
%     ax = gca;
%     ax.XTickLabel = ax.XTickLabel;
%     ax.YTickLabel = ax.YTickLabel;
%     ax.FontSize = 5;
%     xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 10);
%     ylabel('Average Alpha Power [dB]', 'FontName', 'Arial', 'FontSize', 10);
%     legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 5);
%     title(['Subject ' num2str(subj)], 'FontName', 'Arial', 'FontSize', 15);
% 
%     % Set x-axis limits
%     % xlim([-0.5 2.5]);
%     xlim([0.5 2]);
% end
% 
% % Save the figure
% saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/pow_erp/pow_erp_subj_sternberg.png');

