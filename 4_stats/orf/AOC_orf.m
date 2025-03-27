%% Temporal Response Function for AOC
%  Ocular Reponse Function (ORF)
%  for associating oculomotor data with EEG data

%% Setup
startup
[subjects, ~, colors, ant128lay] = setup('AOC');

%% Load resting data
for subj = 1%:length(subjects)
    try
        % Load data
        datapath = strcat('/Volumes/methlab/Students/Arne/AOC/data/merged/', subjects{subj}, filesep, subjects{subj}, '_EEG_ET_RestingEO_merged.mat');
        load(datapath)
        disp(['Data loaded for Subject ', num2str(subjects{subj})])
    catch ME
        ME.message
        error(['ERROR in Subject ' num2str(subjects{subj}) '!'])
    end
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
% Define TRF parameters
fs = EEG.srate;       % Sampling rate (Hz)
tmin = -500;          % Start lag in ms
tmax = 500;           % End lag in ms
lambda = 1e2;         % Regularisation parameter for ridge regression
direction = 1;        % Forward model: stimulus (ocular) -> EEG

% Compute TRFs for each channel
for chans = 1:size(EEG_brain,1)
    modelTRF_gazeX(chans) = mTRFtrain(stim_gazeX, EEG_brain(chans, :), fs, direction, tmin, tmax, lambda);
    modelTRF_gazeY(chans) = mTRFtrain(stim_gazeY, EEG_brain(chans, :), fs, direction, tmin, tmax, lambda);
    modelTRF_pupil(chans) = mTRFtrain(stim_pupil, EEG_brain(chans, :), fs, direction, tmin, tmax, lambda);
    disp(['Computation of TRFs for channel ' num2str(chans) '/', num2str(size(EEG_brain,1)), ' done.'])
end

%% Average the TRF Across EEG Channels for Visualisation
modelTRF_gazeX_avg = mean(modelTRF_gazeX);
modelTRF_gazeY_avg = mean(modelTRF_gazeY);
modelTRF_pupil_avg = mean(modelTRF_pupil);

%% Plot the TRFs
figure;

% Create Time Vector for the TRF Lags
num_lags = size(w_gazeX, 1);
time_lags = linspace(tmin, tmax, num_lags);

subplot(3,1,1);
plot(time_lags, TRF_gazeX_avg, 'LineWidth', 2);
title('Temporal Response Function for Gaze X');
xlabel('Lag [ms]');
ylabel('Amplitude');

subplot(3,1,2);
plot(time_lags, TRF_gazeY_avg, 'LineWidth', 2);
title('Temporal Response Function for Gaze Y');
xlabel('Lag [ms]');
ylabel('Amplitude');

subplot(3,1,3);
plot(time_lags, TRF_pupil_avg, 'LineWidth', 2);
title('Temporal Response Function for Pupil Size (L-AREA)');
xlabel('Lag [ms]');
ylabel('Amplitude');
