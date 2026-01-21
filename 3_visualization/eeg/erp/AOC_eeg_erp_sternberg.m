%% AOC ERP — Sternberg
% Loads dataEEG_TFR_sternberg, computes ERPs per condition (load 2/4/6), plots. Saves figures.
%
% Key outputs:
%   ERP figures (per condition)

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data and compute ERPs
for subj = 1:length(subjects)
    % Load data
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath);
    load('dataEEG_TFR_sternberg.mat');

    % Identify indices of trials belonging to conditions
    ind2 = find(dataTFR.trialinfo == 22); % WM load 2
    ind4 = find(dataTFR.trialinfo == 24); % WM load 4
    ind6 = find(dataTFR.trialinfo == 26); % WM load 6

    % Select data per condition
    cfg.trials = ind2;
    dat2 = ft_selectdata(cfg, dataTFR);
    cfg.trials = ind4;
    dat4 = ft_selectdata(cfg, dataTFR);
    cfg.trials = ind6;
    dat6 = ft_selectdata(cfg, dataTFR);

    % Baseline
    cfg = [];
    cfg.baseline = [-1.5 0.5];
    dat2bl = ft_timelockbaseline(cfg, dat2);
    dat4bl = ft_timelockbaseline(cfg, dat4);
    dat6bl = ft_timelockbaseline(cfg, dat6);

    % Compute ERP
    cfg = [];
    cfg.keepindividual = 'no';
    erp2{subj} = ft_timelockanalysis(cfg, dat2);
    erp4{subj} = ft_timelockanalysis(cfg, dat4);
    erp6{subj} = ft_timelockanalysis(cfg, dat6);

    % Save
    cd(datapath);
    save ERPs_sternberg erp2 erp4 erp6 

    fprintf('Processed ERP for subject %d/%d \n', subj, numel(subjects));
end

%% Compute grand avg
cfg = [];
gaERP2 = ft_timelockgrandaverage(cfg, erp2{:});
gaERP4 = ft_timelockgrandaverage(cfg, erp4{:});
gaERP6 = ft_timelockgrandaverage(cfg, erp6{:});

%% Visualize ERP
close all
figure;
set(gcf, 'Position', [0 0 1200 800], 'Color', 'W')
cfg = [];
cfg.channel = channels;
%cfg.channel = {'POz'};
cfg. xlim = [-.5 2];
cfg.linewidth = 2;
ft_singleplotER(cfg, gaERP2, gaERP4, gaERP6);
hold on
xline(0, '--')
legend({'WM load 2', 'WM load 4', 'WM load 6'}, 'FontSize', 15)
set(gca, "FontSize", 20)
title('Sternberg GA ERPs', 'FontSize', 25)
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/erp/AOC_eeg_erp_sternberg.png')
