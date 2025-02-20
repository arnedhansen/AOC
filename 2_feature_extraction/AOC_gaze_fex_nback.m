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
gaze_data_nback = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, 'PupilSize', {}, 'MSRate', {});

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
    microsaccadeRate = [];

    %% Get trial-by-trial gaze data
    for trl = 1:length(dataet.trialinfo)
        close all
        data = dataet.trial{trl};

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
        gaze_euclidean_dev = zeros(1, length(x_coords) - 1);
        % Calculate Euclidean distances
        for samps = 1:length(x_coords)
            dx = x_coords(samps) - 400; % Distance from middle of x-axis (total 800 px)
            dy = y_coords(samps) - 300; % Distance from middle of y-axis (total 600 px)
            gaze_euclidean_dev(samps) = sqrt(dx^2 + dy^2);
        end
        % Calculate the mean Euclidean distance
        mean_euclidean_distance = mean(gaze_euclidean_dev, 'omitnan');

        %% Compute microsaccades
        fsample = 500; % Sample rate of 500 Hz
        velData = [gaze_x{subj, trl}; gaze_y{subj, trl}]; % Concatenate x and y gaze coordinates to compute the velocity of eye movements in a 2D space
        trlLength = length(dataet.time{trl});
        microsaccade_rate = detect_microsaccades(fsample, velData, trlLength);

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; dataet.trialinfo(trl)];
        gazeDev = [gazeDev; mean_euclidean_distance];
        pupilSize = [pupilSize; pups];
        microsaccadeRate = [microsaccadeRate; microsaccade_rate];
    end
    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'GazeDeviation', num2cell(gazeDev), 'PupilSize', num2cell(pupilSize), 'MSRate', num2cell(microsaccadeRate));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    l1 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 21);
    l1gdev = mean([l1.GazeDeviation], 'omitnan');
    l1pups = mean([l1.PupilSize], 'omitnan');
    l1msrate = mean([l1.MSRate], 'omitnan');

    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 22);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2pups = mean([l2.PupilSize], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');

    l3 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 23);
    l3gdev = mean([l3.GazeDeviation], 'omitnan');
    l3pups = mean([l3.PupilSize], 'omitnan');
    l3msrate = mean([l3.MSRate], 'omitnan');

    %% Create across condition structure
    subj_data_gaze = struct('ID', num2cell(subject_id(1:3)), 'Condition', num2cell([1; 2; 3]), 'GazeDeviation', num2cell([l1gdev; l2gdev; l3gdev]), 'PupilSize', num2cell([l1pups; l2pups; l3pups]), 'MSRate', num2cell([l1msrate; l2msrate; l3msrate]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_nback_trial subj_data_gaze_trial
    save gaze_matrix_nback subj_data_gaze
    save gaze_dev_nback l1gdev l2gdev l3gdev
    save pupil_size_nback l1pups l2pups l3pups
    save ms_rate_nback l1msrate l2msrate l3msrate
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_nback = [gaze_data_nback; subj_data_gaze];
end
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_nback gaze_x gaze_y
save /Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_nback gaze_data_nback
