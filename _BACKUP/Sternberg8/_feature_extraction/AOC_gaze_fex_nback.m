%% AOC Gaze Feature Extraction N-back
%
% Extracted features:
%   Gaze deviation
%   Pupil size

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    load([datapath, filesep, 'dataET_nback'])
    ind1 = find(dataet.trialinfo==1);
    ind2 = find(dataet.trialinfo==2);
    ind3 = find(dataet.trialinfo==3);
    
    %% FT data structure
    cfg = [];
    cfg.latency = [0 2];
    cfg.trials = ind1;
    dataetL1 = ft_selectdata(cfg,dataet);
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind3;
    dataetL3 = ft_selectdata(cfg,dataet);

    %% Append data
    for loads = 1:3
        if loads == 1
            data = dataetL1;
            data = horzcat(dataetL1.trial{:});
        elseif loads == 2
            data = dataetL2;
            data = horzcat(dataetL2.trial{:});
        elseif loads == 3
            data = dataetL3;
            data = horzcat(dataetL3.trial{:});
        end

        %% Clean data
        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(:, valid_data_indices);

        % Remove data points that contain zeros
        window_size = 50;
        cleaned_data = remove_blink_window(data, window_size);
        data = cleaned_data;

        %% Extract gaze data
        gaze_x{subj,loads} = data(1, :);
        gaze_y{subj,loads} = data(2, :);
    end
    
    %% Compute gaze deviation
    conditions = [];
    gazeX = cell(1, length(subjects));
    gazeY = cell(1, length(subjects));
    timePoints = [];

    % Accumulate the gaze data for each subject
    numTrials = length(dataet.trial);
    lGazeX = cell(1, numTrials);
    lGazeY = cell(1, numTrials);
    for i = 1:numTrials
        lGazeX{i} = dataet.trial{1, i}(1, :);
        lGazeY{i} = dataet.trial{1, i}(2, :);
    end

    % Store unique conditions and time points
    if isempty(conditions)
        conditions = unique(dataet.trialinfo);
        timePoints = length(dataet.time{1});
    end

    % Average gaze data for current subject
    subjectAverageGazeX = cell(1, length(conditions));
    subjectAverageGazeY = cell(1, length(conditions));
    for condIdx = 1:length(conditions)
        cond = conditions(condIdx);
        % Get indices of trials corresponding to the current condition
        condTrialsIdx = find(dataet.trialinfo == cond);

        % Preallocate matrix to hold gaze data for current condition
        condGazeX = zeros(length(condTrialsIdx), timePoints);
        condGazeY = zeros(length(condTrialsIdx), timePoints);

        % Extract gaze data for the current condition
        for j = 1:length(condTrialsIdx)
            condGazeX(j, :) = lGazeX{condTrialsIdx(j)};
            condGazeY(j, :) = lGazeY{condTrialsIdx(j)};
        end

        % Compute average gaze position for each time point for the current subject
        subjectAverageGazeX{condIdx} = mean(condGazeX, 1);
        subjectAverageGazeY{condIdx} = mean(condGazeY, 1);
    end

    % Store the subject's average gaze data
    gazeX{subj} = subjectAverageGazeX;
    gazeY{subj} = subjectAverageGazeY;

    %% Compute pupil size
    

    %% Save gaze variables for individual participants
    cd(datapath)
    save nback_gaze_deviation.mat gaze X gazeY
    save nback_gaze_pupil_size.mat pupil_size
    clc
    fprintf('Subject %.3d/%.3d done \n', subj, length(subjects))
end

%% Compute overall gaze deviation
% grand average and standard error across all subjects for each condition
grandAverageGazeX = cell(1, length(conditions));
grandErrorGazeX = cell(1, length(conditions));
grandAverageGazeY = cell(1, length(conditions));
grandErrorGazeY = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    % Concatenate all subjects' average gaze data for the current condition
    allSubjectGazeX = zeros(length(subjects), timePoints);
    allSubjectGazeY = zeros(length(subjects), timePoints);

    for subj = 1:length(subjects)
        allSubjectGazeX(subj, :) = gazeX{subj}{condIdx};
        allSubjectGazeY(subj, :) = gazeY{subj}{condIdx};
    end

    % Compute the grand average and standard error for the current condition
    grandAverageGazeX{condIdx} = mean(allSubjectGazeX, 1);
    grandErrorGazeX{condIdx} = std(allSubjectGazeX, 0, 1) / sqrt(length(subjects));
    grandAverageGazeY{condIdx} = mean(allSubjectGazeY, 1);
    grandErrorGazeY{condIdx} = std(allSubjectGazeY, 0, 1) / sqrt(length(subjects));
end

%% Save aggregated gaze variables for all participants
save /Volumes/methlab/Students/Arne/AOC/data/features/nback_gaze_deviation gaze_deviation
save /Volumes/methlab/Students/Arne/AOC/data/features/nback_gaze_pupil_size.mat pupil_size

%% Define function for blink removal (N-back Task)
function cleaned_data = remove_blink_window(data, window_size)
blink_indices = find(all(data(1:2, :) == 0, 1));
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end
data(:, removal_indices) = [];
cleaned_data = data;
end