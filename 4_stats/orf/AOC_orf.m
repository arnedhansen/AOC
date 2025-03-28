%% Temporal Response Function for AOC - Ocular Response Function (ORF)
% This script computes TRFs to associate oculomotor data with EEG data for each
% subject, saves the computed models, and creates figures for both individual subjects
% and the grand average.

%% Setup
startup;
[subjects, ~, colors, ant128lay] = setup('AOC');

% Define base paths based on the operating system
if ispc
    base_data_path     = 'W:\Students\Arne\AOC\data\merged\';
    base_features_path = 'W:\Students\Arne\AOC\data\features\';
    base_figures_path  = 'W:\Students\Arne\AOC\figures\orf\';
    addpath(genpath('W:\Students\Arne\toolboxes\mTRF-Toolbox-master'))
else
    base_data_path     = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
    base_features_path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
    base_figures_path  = '/Volumes/methlab/Students/Arne/AOC/figures/orf/';
end

% Preallocate arrays to store averaged TRFs for each ocular regressor across subjects
all_TRF_gazeX = [];
all_TRF_gazeY = [];
all_TRF_pupil = [];
subject_count = 0;  % count successfully processed subjects

% Loop over subjects
for s = 1:length(subjects)
    subjectID = subjects{s};
    try
        % Load subject data using the appropriate base path
        datapath = fullfile(base_data_path, subjectID, [subjectID '_EEG_ET_RestingEO_merged.mat']);
        load(datapath);  % Expects EEG variable in the file
        disp(['Data loaded for Subject ' subjectID]);
    catch ME
        warning(['Skipping Subject ' subjectID ': ' ME.message]);
        continue;  % Skip to the next subject if an error occurs
    end

    %% Identify Ocular Channels
    % 130: L-GAZE-X, 131: L-GAZE-Y, 132: L-AREA (pupil size)
    idx_gazeX = 130;
    idx_gazeY = 131;
    idx_pupil  = 132;

    %% Extract Ocular Signals
    stim_gazeX = EEG.data(idx_gazeX, :);
    stim_gazeY = EEG.data(idx_gazeY, :);
    stim_pupil = EEG.data(idx_pupil, :);

    %% Define EEG Channels (excluding ocular channels)
    EEG_channels = 1:(size(EEG.data,1)-3);
    EEG_brain = EEG.data(EEG_channels, :);

    %% Compute TRFs Using the mTRF Toolbox
    fs = EEG.srate;       % Sampling rate (Hz)
    tmin = -100;          % Start lag in ms
    tmax = 400;           % End lag in ms
    lambda = 1e2;         % Regularisation parameter for ridge regression
    direction = 1;        % Forward model: stimulus (ocular) -> EEG

    nChannels = size(EEG_brain,1);
    % Initialise empty arrays for storing the TRF models
    clear modelTRF_gazeX modelTRF_gazeY modelTRF_pupil

    for chan = 1:nChannels
        modelTRF_gazeX(chan) = mTRFtrain(stim_gazeX, EEG_brain(chan, :), fs, direction, tmin, tmax, lambda);
        modelTRF_gazeY(chan) = mTRFtrain(stim_gazeY, EEG_brain(chan, :), fs, direction, tmin, tmax, lambda);
        modelTRF_pupil(chan) = mTRFtrain(stim_pupil, EEG_brain(chan, :), fs, direction, tmin, tmax, lambda);
        disp(['Computed TRF for channel ' num2str(chan) ' of ' num2str(nChannels) ' for subject ' subjectID '.']);
    end

    %% Compute the average TRF across channels for each ocular regressor for this subject
    nLags = length(modelTRF_gazeX(1).w);

    % For Gaze X:
    all_w_gazeX = zeros(nLags, nChannels);
    for i = 1:nChannels
        all_w_gazeX(:, i) = modelTRF_gazeX(i).w;
    end
    TRF_gazeX_avg = mean(all_w_gazeX, 2);

    % For Gaze Y:
    all_w_gazeY = zeros(nLags, nChannels);
    for i = 1:nChannels
        all_w_gazeY(:, i) = modelTRF_gazeY(i).w;
    end
    TRF_gazeY_avg = mean(all_w_gazeY, 2);

    % For Pupil:
    all_w_pupil = zeros(nLags, nChannels);
    for i = 1:nChannels
        all_w_pupil(:, i) = modelTRF_pupil(i).w;
    end
    TRF_pupil_avg = mean(all_w_pupil, 2);

    %% Save the computed TRF models for this subject
    models_folder = fullfile(base_features_path, subjectID, 'orf');
    if ~exist(models_folder, 'dir')
        mkdir(models_folder);
    end
    save(fullfile(models_folder, 'TRF_models.mat'), 'modelTRF_gazeX', 'modelTRF_gazeY', 'modelTRF_pupil');

    %% Plot the averaged TRFs for this subject
    time_lags = modelTRF_gazeX(1).t;  % Time vector for the TRF lags
    figure;

    subplot(1,3,1);
    plot(time_lags, TRF_gazeX_avg, 'LineWidth', 2);
    title(['Temporal Response Function for Gaze X - Subject ' subjectID]);
    xlabel('Lag [ms]');
    ylabel('Amplitude');

    subplot(1,3,2);
    plot(time_lags, TRF_gazeY_avg, 'LineWidth', 2);
    title(['Temporal Response Function for Gaze Y - Subject ' subjectID]);
    xlabel('Lag [ms]');
    ylabel('Amplitude');

    subplot(1,3,3);
    plot(time_lags, TRF_pupil_avg, 'LineWidth', 2);
    title(['Temporal Response Function for Pupil Size (L-AREA) - Subject ' subjectID]);
    xlabel('Lag [ms]');
    ylabel('Amplitude');

    %% Save the figure for this subject
    if ~exist(base_figures_path, 'dir')
        mkdir(base_figures_path);
    end
    saveas(gcf, fullfile(base_figures_path, ['AOC_orf_' subjectID '.png']));
    close(gcf);

    %% Accumulate the subject's averaged TRFs for the grand average computation
    subject_count = subject_count + 1;
    all_TRF_gazeX(:, subject_count) = TRF_gazeX_avg;
    all_TRF_gazeY(:, subject_count) = TRF_gazeY_avg;
    all_TRF_pupil(:, subject_count) = TRF_pupil_avg;
end

%% Load subject TRFs
nChannels = 129;
for subs = 1:length(subjects)
    try
        subjectID = subjects{subs};
        models_folder = fullfile(base_features_path, subjectID, 'orf');
        datapath = fullfile(models_folder, 'TRF_models.mat');
        load(datapath);

        % Compute the average TRF across channels for each ocular regressor for this subject
        nLags = length(modelTRF_gazeX(1).w);
        all_w_gazeX = zeros(nLags, nChannels);
        for i = 1:nChannels
            all_w_gazeX(:, i) = modelTRF_gazeX(i).w;
        end
        TRF_gazeX_avg = mean(all_w_gazeX, 2);
        all_w_gazeY = zeros(nLags, nChannels);
        for i = 1:nChannels
            all_w_gazeY(:, i) = modelTRF_gazeY(i).w;
        end
        TRF_gazeY_avg = mean(all_w_gazeY, 2);
        all_w_pupil = zeros(nLags, nChannels);
        for i = 1:nChannels
            all_w_pupil(:, i) = modelTRF_pupil(i).w;
        end
        TRF_pupil_avg = mean(all_w_pupil, 2);

        % Accumulate subject data
        all_TRF_gazeX(:, subs) = TRF_gazeX_avg;
        all_TRF_gazeY(:, subs) = TRF_gazeY_avg;
        all_TRF_pupil(:, subs) = TRF_pupil_avg;

        disp(['TRF models loaded for Subject ' subjectID]);
    catch ME
        warning(['Skipping Subject ' subjectID ': ' ME.message]);
        continue;  % Skip to the next subject if an error occurs
    end
end

%% Compute and Plot the Grand Average TRF across Subjects
close all

% Get average for all models
grand_TRF_gazeX = mean(all_TRF_gazeX, 2, 'omitnan');
grand_TRF_gazeY = mean(all_TRF_gazeY, 2, 'omitnan');
grand_TRF_pupil = mean(all_TRF_pupil, 2, 'omitnan');

% Smooth using a Gaussian window of width 5 samples
grand_TRF_gazeX = smoothdata(grand_TRF_gazeX, 'gaussian', 50);
grand_TRF_gazeY = smoothdata(grand_TRF_gazeY, 'gaussian', 50);
grand_TRF_pupil = smoothdata(grand_TRF_pupil, 'gaussian', 25);

% Plot
figure;
set(gcf, 'Position', [0 0 2000 1000])
time_lags = modelTRF_gazeX(1).t;  % Time vector for the TRF lags

subplot(1,2,1);
plot(time_lags, grand_TRF_gazeX, 'LineWidth', 2);
title('Grand Average Temporal Response Function for Gaze X');
xlabel('Lag [ms]');
ylabel('Amplitude');
xline(0, 'r', 'LineWidth', 1);
yline(0, '--');
set(gca, "YLim", [-0.05 0.05])
set(gca, 'FontSize', 15)

subplot(1,2,2);
plot(time_lags, grand_TRF_gazeY, 'LineWidth', 2);
title('Grand Average Temporal Response Function for Gaze Y');
xlabel('Lag [ms]');
ylabel('Amplitude');
xline(0, 'r', 'LineWidth', 1);
yline(0, '--');
set(gca, "YLim", [-0.03 0.03])
set(gca, 'FontSize', 15)

% Save the grand average figure
saveas(gcf, fullfile(base_figures_path, 'AOC_orf_all.png'));
