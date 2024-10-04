%% AOC Gaze Feature Extraction N-back
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
gaze_data_nback = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, 'PupilSize', {});

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET_nback'])

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
        condition = [condition; dataet.trialinfo(trl)];
        gazeDev = [gazeDev; mean_euclidean_distance];
        pupilSize = [pupilSize; mean(pupil_size, 'omitnan')/1000];
    end
    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'PupilSize', num2cell(pupilSize));

    %% Calculate subject-specific GazeDev by condition
    l1 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 1);
    l1gdev = mean([l1.GazeDeviation], 'omitnan');
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l3 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 3);
    l3gdev = mean([l3.GazeDeviation], 'omitnan');

    %% Calculate subject-specific PupilSize by condition
    l1 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 1);
    l1pups = mean([l1.PupilSize], 'omitnan');
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2pups = mean([l2.PupilSize], 'omitnan');
    l3 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 3);
    l3pups = mean([l3.PupilSize], 'omitnan');

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:3)), 'Condition', num2cell([1; 2; 3]), 'GazeDeviation', num2cell([l1gdev; l2gdev; l3gdev]), 'PupilSize', num2cell([l1pups; l2pups; l3pups]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_nback_trial subj_data_gaze_trial
    save gaze_matrix_nback subj_data_gaze
    save gaze_dev_nback l1gdev l2gdev l3gdev
    save pupil_size_nback l1pups l2pups l3pups
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_nback = [gaze_data_nback; subj_data_gaze];
end
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_nback gaze_x gaze_y
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_nback gaze_data_nback

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
