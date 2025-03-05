%% AOC Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation (Euclidean distances)
%   Gaze standard deviation
%   Pupil size
%   Microsaccades
%
% Gaze metrics labelled by eye-tracker (saccades, blinks and
% fixations) are extracted already in AOC_preprocessing_sternberg.m

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
gaze_data_sternberg = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeStdX', {}, 'GazeStdY', {}, 'PupilSize', {}, 'MSRate', {}, ...
    'Blinks', {}, 'Fixations', {}, 'Saccades', {});

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET_sternberg'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = length(dataet.trialinfo);
    condition = [];
    gazeDev = [];
    gazeSDx = [];
    gazeSDy = [];
    pupilSize = [];
    microsaccadeRate = [];

    %% Get trial-by-trial gaze data
    for trl = 1:length(dataet.trialinfo)
        close all
        data = dataet.trial{trl};

        %% Choose data 1000 - 2000ms after stimulus presentation to exclude evoked activity
        % Should already be this way: 500 samples for each 2ms
        analysis_period = [1 2];
        time_vector = dataet.time{trl};
        analysis_idx = (time_vector >= analysis_period(1)) & (time_vector <= analysis_period(2));
        data = data(:, analysis_idx);

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices);
        data = valid_data;
        data(2, :) = 600 - data(2, :); % Invert Y-axis

        %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        %% Extract gaze data and pupil size
        gaze_x{subj, trl} = data(1, :);
        gaze_y{subj, trl} = data(2, :);
        pupil_size{subj, trl} = mean(data(3, :), 'omitnan') / 1000;
        pups = pupil_size{subj, trl};

        %% Compute gaze deviation as euclidean distances from the center
        x_coords = gaze_x{subj, trl};
        y_coords = gaze_y{subj, trl};

        % Calculate Euclidean distances
        dx = x_coords - 400 - nanmean(x_coords - 400); % Distance from mean (to get rid of impact of gaze shifts)
        dy = y_coords - 300 - nanmean(y_coords - 300); % Distance from mean (to get rid of impact of gaze shifts)
        gaze_euclidean_dev = sqrt(dx.^2 + dy.^2);

        % Calculate the mean Euclidean distance
        mean_euclidean_distance = nanmean(gaze_euclidean_dev);
        gaze_standard_deviation_x = nanstd(x_coords);
        gaze_standard_deviation_y = nanstd(y_coords);

        %% Compute microsaccades
        fsample = 500; % Sample rate of 500 Hz
        velData = [gaze_x{subj, trl}; gaze_y{subj, trl}]; % Concatenate x and y gaze coordinates to compute the velocity of eye movements in a 2D space
        trlLength = length(dataet.time{trl});
        [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength);
        ms_data(trl) = microsaccade_details;

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; dataet.trialinfo(trl)-50];
        gazeDev = [gazeDev; mean_euclidean_distance];
        gazeSDx = [gazeSDx; gaze_standard_deviation_x];
        gazeSDy = [gazeSDy; gaze_standard_deviation_y];
        pupilSize = [pupilSize; pups];
        microsaccadeRate = [microsaccadeRate; microsaccade_rate];
    end

    %% Check data by visualizing raw gaze data
    % close all
    % % Preallocate arrays for averaged gaze data
    % num_samples = 1000;
    % mean_gaze_x = nan(num_samples, 1);
    % mean_gaze_y = nan(num_samples, 1);
    %
    % % Stack trials into matrices for averaging
    % gaze_x_matrix = cell2mat(cellfun(@(x) x(:), gaze_x(subj, :), 'UniformOutput', false)');
    % gaze_y_matrix = cell2mat(cellfun(@(y) y(:), gaze_y(subj, :), 'UniformOutput', false)');
    %
    % % Calculate the mean over trials for each sample
    % mean_gaze_x = nanmean(gaze_x_matrix, 2);
    % mean_gaze_y = nanmean(gaze_y_matrix, 2);
    %
    % % Plot the averaged gaze data
    % figure;
    % set(gcf, "Position", [200, 200, 1000, 600]);
    % plot(mean_gaze_x, mean_gaze_y, 'o');
    % hold on;
    % plot(400, 300, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Centre point
    % title('Averaged Gaze Data Distribution Across Samples');
    % xlabel('X Coordinates');
    % ylabel('Y Coordinates');
    % xlim([0 800]);
    % ylim([0 600]);

    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), ...
        'Trial', num2cell(trial_num), ...
        'Condition', num2cell(condition), ...
        'GazeDeviation', num2cell(gazeDev), ...
        'GazeStdX', num2cell(gazeSDx), ...
        'GazeStdY', num2cell(gazeSDy), ...
        'PupilSize', num2cell(pupilSize), ...
        'MSRate', num2cell(microsaccadeRate));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2gSDx = mean([l2.GazeStdX], 'omitnan');
    l2gSDy = mean([l2.GazeStdY], 'omitnan');
    l2pups = mean([l2.PupilSize], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');

    l4 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 4);
    l4gdev = mean([l4.GazeDeviation], 'omitnan');
    l4gSDx = mean([l4.GazeStdX], 'omitnan');
    l4gSDy = mean([l4.GazeStdY], 'omitnan');
    l4pups = mean([l4.PupilSize], 'omitnan');
    l4msrate = mean([l4.MSRate], 'omitnan');

    l6 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 6);
    l6gdev = mean([l6.GazeDeviation], 'omitnan');
    l6gSDx = mean([l2.GazeStdX], 'omitnan');
    l6gSDy = mean([l2.GazeStdY], 'omitnan');
    l6pups = mean([l6.PupilSize], 'omitnan');
    l6msrate = mean([l6.MSRate], 'omitnan');

    %% Load gaze metrics (extracted in GCP_preprocessing.m)
    load([datapath, filesep, 'gaze_metrics_sternberg'])

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:3)), ...
        'Condition', num2cell([2; 4; 6]), ...
        'GazeDeviation', num2cell([l2gdev; l4gdev; l6gdev]), ...
        'GazeStdX', num2cell([l2gSDx; l4gSDx; l6gSDx]), ...
        'GazeStdY', num2cell([l2gSDy; l4gSDy; l6gSDy]), ...
        'PupilSize', num2cell([l2pups; l4pups; l6pups]), ...
        'MSRate', num2cell([l2msrate; l4msrate; l6msrate]), ...
        'Blinks', num2cell([blinks_l2; blinks_l4; blinks_l6]), ...
        'Fixations', num2cell([fixations_l2; fixations_l4; fixations_l6]), ...
        'Saccades', num2cell([saccades_l2; saccades_l4; saccades_l6]));

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_sternberg_trial subj_data_gaze_trial
    save gaze_matrix_sternberg subj_data_gaze
    save gaze_dev_sternberg l2gdev l4gdev l6gdev
    save gaze_std_nback l2gSDx l4gSDx l6gSDx l2gSDy l4gSDy l6gSDy
    save pupil_size_sternberg l2pups l4pups l6pups
    save ms_rate_sternberg l2msrate l4msrate l6msrate
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_sternberg = [gaze_data_sternberg; subj_data_gaze];
end
trialinfo = dataet.trialinfo';
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_sternberg gaze_x gaze_y trialinfo
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg gaze_data_sternberg
