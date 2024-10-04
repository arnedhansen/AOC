%% Alpha Power as ERP plot for Nback data

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

alpha_pow1 = [];
alpha_pow2 = [];
alpha_pow3 = [];

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load dataEEG_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==1);
    ind2=find(data.trialinfo==2);
    ind3=find(data.trialinfo==3);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    load1 = ft_freqanalysis(cfg,data);
    cfg.trials = ind2;
    load2 = ft_freqanalysis(cfg,data);
    cfg.trials = ind3;
    load3 = ft_freqanalysis(cfg,data);

    %% Save data
    window_time = load1.time;
    cd(datapath)
    save erp_pow_nback load1 load2 load3 window_time
    cd('/Volumes/methlab/Students/Arne/AOC/scripts/eeg')
    fprintf('Subject %.2d done \n', subj)
end

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load erp_pow_nback
    % baseline
    cfg =[];
    cfg.baseline = [-Inf -.5];% avoids taking -.2 activity which is the ERP onset
    cfg.baselinetype ='db';
    load1 = ft_freqbaseline(cfg,load1);
    load2 = ft_freqbaseline(cfg,load2);
    load3 = ft_freqbaseline(cfg,load3);
    l1{subj}= load1;
    l2{subj}= load2;
    l3{subj}= load3;
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr1 = ft_freqgrandaverage([],l1{:});
gatfr2 = ft_freqgrandaverage([],l2{:});
gatfr3 = ft_freqgrandaverage([],l3{:});
% gatfr.powspctrm =  freq data with 124 channels, 27 frequencybins and 121 timebins

%% Plot aggregated alpha power
addpath('/Users/Arne/Documents/matlabtools/shadedErrorBar')
close all
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', ...
    'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', ...
    'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
alpha_range = 8:13;

% Compute mean alpha power for coi
compute_mean_alpha = @(data, coi, alpha_range) ...
    mean(mean(data.powspctrm(ismember(data.label, coi), alpha_range, :), 2), 1);
alpha_power_1 = squeeze(compute_mean_alpha(gatfr1, coi, alpha_range))';
alpha_power_2 = squeeze(compute_mean_alpha(gatfr2, coi, alpha_range))';
alpha_power_3 = squeeze(compute_mean_alpha(gatfr3, coi, alpha_range))';

% Compute SEM
compute_sem = @(data) std(data) / sqrt(size(data, 1));
sem_1 = compute_sem(squeeze(mean(gatfr1.powspctrm(ismember(gatfr1.label, coi), alpha_range, :), 2)));
sem_2 = compute_sem(squeeze(mean(gatfr2.powspctrm(ismember(gatfr2.label, coi), alpha_range, :), 2)));
sem_3 = compute_sem(squeeze(mean(gatfr3.powspctrm(ismember(gatfr3.label, coi), alpha_range, :), 2)));

figure;
set(gcf, "Position", [0, 0, 2000, 1200], "Color", 'w')
hold on;
shadedErrorBar(gatfr1.time, alpha_power_1, sem_1, 'lineprops', {'r', 'LineWidth', 1});
shadedErrorBar(gatfr2.time, alpha_power_2, sem_2, 'lineprops', {'g', 'LineWidth', 1});
shadedErrorBar(gatfr3.time, alpha_power_3, sem_3, 'lineprops', {'b', 'LineWidth', 1});
hold off;

ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Alpha Power [dB]', 'FontName', 'Arial', 'FontSize', 20);
legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20);
title('Posterior Alpha Power during N-back', 'FontName', 'Arial', 'FontSize', 25);
% xlim([-0.5 2.5])

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/pow_erp/pow_erp_nback.png');

%% Plot for individual subjects
close all
figure;
set(gcf, "Position", [0, 0, 2000, 1200], "Color", 'w');

for subj = 1:length(subjects)
    % Get data
    load1 = l1{subj};
    load2 = l2{subj};
    load3 = l3{subj};

    % Mean alpha power per subject
    alpha_power_1 = squeeze(compute_mean_alpha(load1, coi, alpha_range))';
    alpha_power_2 = squeeze(compute_mean_alpha(load2, coi, alpha_range))';
    alpha_power_3 = squeeze(compute_mean_alpha(load3, coi, alpha_range))';

    % SEM for each subject
    sem_1 = compute_sem(squeeze(mean(load1.powspctrm(ismember(load1.label, coi), alpha_range, :), 2)));
    sem_2 = compute_sem(squeeze(mean(load2.powspctrm(ismember(load2.label, coi), alpha_range, :), 2)));
    sem_3 = compute_sem(squeeze(mean(load3.powspctrm(ismember(load3.label, coi), alpha_range, :), 2)));

    % Plot
    subplot(5, 2, subj);
    hold on;
    shadedErrorBar(load1.time, alpha_power_1, sem_1, 'lineprops', {'r', 'LineWidth', 1});
    shadedErrorBar(load2.time, alpha_power_2, sem_2, 'lineprops', {'g', 'LineWidth', 1});
    shadedErrorBar(load3.time, alpha_power_3, sem_3, 'lineprops', {'b', 'LineWidth', 1});
    hold off;

    ax = gca;
    ax.XTickLabel = ax.XTickLabel;
    ax.YTickLabel = ax.YTickLabel;
    ax.FontSize = 5;
    xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('Average Alpha Power [dB]', 'FontName', 'Arial', 'FontSize', 10);
    legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20);
    title(['Subject ' num2str(subj)], 'FontName', 'Arial', 'FontSize', 15);

    % Set x-axis limits
     xlim([-0.5 2.5]);
    %xlim([0.5 2]);
end

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/pow_erp/pow_erp_subj.png');

