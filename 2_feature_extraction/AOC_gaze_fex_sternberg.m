%% AOC Gaze Feature Extraction Sternberg
%
% Extracted features:
%   Gaze deviation
%   Pupil size

%% Setup 
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
gaze_data_sternberg = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, 'PupilSize', {});

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
    pupilSize = [];

    %% Get trial-by-trial gaze data
    for trl = 1:length(dataet.trialinfo)
        close all
        data = dataet.trial{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices); % Excluding pupil size data

        %% Remove blinks with a window of 100ms (= 50 timepoints)
        data = remove_blink_window(data, 50);

        %% Extract gaze data and pupil size
        gaze_x{subj, trl} = data(1, :);
        gaze_y{subj, trl} = data(2, :);
        pupil_size = data(3, :);

        %% Compute gaze deviation as euclidean distances from the center
        x_coords = gaze_x{subj, trl};
        y_coords = gaze_y{subj, trl};
        gaze_euclidean_dev = zeros(1, length(x_coords) - 1);
        % Calculate Euclidean distances
        for samps = 1:length(x_coords)
            dx = x_coords(samps) - 400; % Distance from middle of x-axis (total 800 px)
            dy = y_coords(samps) - 300; % Distance from middle of y-axis (total 600 px)
            gaze_euclidean_dev(samps) = sqrt(dx^2 + dy^2);
        end
        % Calculate the mean Euclidean distance
        mean_euclidean_distance = mean(gaze_euclidean_dev, 'omitnan');

        % Sanity check
        % plot(gaze_euclidean_dev)

        %% Compute miscrosaccades

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; dataet.trialinfo(trl)-50];
        gazeDev = [gazeDev; mean_euclidean_distance];
        pupilSize = [pupilSize; mean(pupil_size, 'omitnan')/1000];
    end
    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'PupilSize', num2cell(pupilSize));

    %% Calculate subject-specific GazeDev by condition
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l4 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 4);
    l4gdev = mean([l4.GazeDeviation], 'omitnan');
    l6 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 6);
    l6gdev = mean([l6.GazeDeviation], 'omitnan');

    %% Calculate subject-specific PupilSize by condition
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2pups = mean([l2.PupilSize], 'omitnan');
    l4 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 4);
    l4pups = mean([l4.PupilSize], 'omitnan');
    l6 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 6);
    l6pups = mean([l6.PupilSize], 'omitnan');

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:3)), 'Condition', num2cell([2; 4; 6]), 'GazeDeviation', num2cell([l2gdev; l4gdev; l6gdev]), 'PupilSize', num2cell([l2pups; l4pups; l6pups]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_sternberg_trial subj_data_gaze_trial
    save gaze_matrix_sternberg subj_data_gaze
    save gaze_dev_sternberg l2gdev l4gdev l6gdev
    save pupil_size_sternberg l2pups l4pups l6pups
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_sternberg = [gaze_data_sternberg; subj_data_gaze];
end
trialinfo = dataet.trialinfo';
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_sternberg gaze_x gaze_y trialinfo
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg gaze_data_sternberg

%% Define function for blink removal
function cleaned_data = remove_blink_window(data, window_size)
blink_indices = find(all(data(1:2, :) == 0, 1));
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end
data(:, removal_indices) = NaN;
cleaned_data = data;
end
