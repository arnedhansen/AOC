%% Master script for the AOC (Alpha Oculomotr Control) Study
%
% - Resting EEG
% - Sternberg
%   - Training (10 trials)
%   - Task (6 blocks x 50 trials)
% - N-back
%   - Training (10 trials)
%   - N-back task (6 blocks x 75 trials)

%% General settings, screens and paths

% Set up MATLAB workspace
clear;
close all;
clc;

% Clear parallel port
ppdev_mex('CloseAll');

% Define paths
PPDEV_PATH = '/home/methlab/Documents/MATLAB/ppdev-mex-master'; % For sending EEG triggers
DATA_PATH = '/home/methlab/Desktop/AOC_data'; % Folder to save data
FUNS_PATH = '/home/methlab/Desktop/AOC' ; % Folder with all functions
addpath(FUNS_PATH) % Add path to folder with functions
screenSettings % Manage screens

%% Collect ID and Age
dialogID;

%% Protect Matlab code from participant keyboard input
ListenChar(2);

%% Resting state EEG
% if ~isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_Resting.mat'])
%     Resting_EEG;
% end

%% Randomize order of Sternberg Task and NBack Task
% Check if the conversion was successful
if isnan(subject.ID)
    error('ERROR. Subject ID must be a numeric value.');
end
% Use subject ID for assignment of pseudorandom task Order (Sternberg & N-back)
if mod(subject.ID, 2) == 0
    taskOrder = 'NBackSternberg'; % Even subjectIDs do Nback first
elseif mod(subject.ID, 2) == 1
    taskOrder = 'SternbergNBack'; % Uneven subjectIDs do Sternberg first
else
    error('ERROR: TASK ORDER COULD NOT BE DETERMINED BY SUBJECT ID')
end
disp(['SubjectID AOC ', num2str(subjectID),'. Task order is: ', taskOrder])

% Shuffle n-back conditions (same rng for same subjectID in case of crashes) 
rng(subject.ID);
nback_condition = repmat(1:3, 1, 2);
nback_condition = nback_condition(randperm(length(nback_condition)));

%% Execute STERNBERG - NBACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(taskOrder, 'SternbergNBack')
    disp(['Subject ' num2str(subjectID) ': Running Sternberg followed by N-Back...']);

    %% Execute Sternberg Task
    % Training phase
    TASK = 'AOC_Sternberg';
    BLOCK = 0;
    TRAINING = 1;
    trainingFile = [num2str(subjectID), '_', TASK, '_block0_training.mat'];
    if isfile([DATA_PATH, '/', num2str(subjectID), '/', trainingFile])
        percentTotalCorrect = 60;
    else
        percentTotalCorrect = 0;
    end

    while percentTotalCorrect < 59
        disp([TASK, ' Training TASK...']);
        eval(TASK); % Run the task
        BLOCK = BLOCK + 1;
    end

    % Actual task
    TRAINING = 0;
    blockCount = 6;
    start = 1;
    for i = blockCount:-1:1
        if isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_', TASK, '_block', num2str(i), '_task.mat']])
            start = i + 1;
            break;
        end
    end

    for BLOCK = start:blockCount
        disp([TASK, ' STARTING...']);
        eval(TASK); % Run the task
    end

    %% Mandatory Break of at least 5 seconds
    disp('Waiting 5 seconds between tasks...');
    WaitSecs(5);

    %% Execute N-back Task
    % Training phase
    TASK = 'AOC_NBack';
    BLOCK = 0;
    TRAINING = 1;
    trainingFile = [num2str(subjectID), '_', TASK, '_block0_training.mat'];
    if isfile([DATA_PATH, '/', num2str(subjectID), '/', trainingFile])
        percentTotalCorrect = 60;
    else
        percentTotalCorrect = 0;
    end

    while percentTotalCorrect < 59
        disp([TASK, ' Training TASK...']);
        COND = nback_condition(1);
        eval(TASK); % Run the task
        BLOCK = BLOCK + 1;
    end

    % Actual task
    TRAINING = 0;
    blockCount = 6;
    start = 1;
    for i = blockCount:-1:1
        for condition = 1:3
            if isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_', TASK, '_block', num2str(i), '_', num2str(condition), 'back_task.mat']])
                start = i + 1;
                break;
            end
        end
        if start ~= 1
            break;
        end
    end

    for BLOCK = start:blockCount
        disp([TASK, ' STARTING...']);
        COND = nback_condition(BLOCK);
        eval(TASK); % Run the task
    end

    %% Execute NBACK - STERNBERG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp(['Subject ' num2str(subjectID) ': Running N-back followed by Sternberg...']);

    %% Execute N-back Task
    % Training phase
    TASK = 'AOC_NBack';
    BLOCK = 0;
    TRAINING = 1;
    trainingFile = [num2str(subjectID), '_', TASK, '_block0_training.mat'];
    if isfile([DATA_PATH, '/', num2str(subjectID), '/', trainingFile])
        percentTotalCorrect = 60;
    else
        percentTotalCorrect = 0;
    end

    while percentTotalCorrect < 59
        disp([TASK, ' Training TASK...']);
        COND = nback_condition(1);
        eval(TASK); % Run the task
        BLOCK = BLOCK + 1;
    end

    % Actual task
    TRAINING = 0;
    blockCount = 6;
    start = 1;
    for i = blockCount:-1:1
        for condition = 1:3
            if isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_', TASK, '_block', num2str(i), '_', num2str(condition), 'back_task.mat']])
                start = i + 1;
                break;
            end
        end
        if start ~= 1
            break;
        end
    end

    for BLOCK = start:blockCount
        disp([TASK, ' STARTING...']);
        COND = nback_condition(BLOCK);
        eval(TASK); % Run the task
    end

    %% Mandatory Break of at least 5 seconds
    disp('Waiting 5 seconds between tasks...');
    WaitSecs(5);

    %% Execute Sternberg Task
    % Training phase
    TASK = 'AOC_Sternberg';
    BLOCK = 0;
    TRAINING = 1;
    trainingFile = [num2str(subjectID), '_', TASK, '_block0_training.mat'];
    if isfile([DATA_PATH, '/', num2str(subjectID), '/', trainingFile])
        percentTotalCorrect = 60;
    else
        percentTotalCorrect = 0;
    end

    while percentTotalCorrect < 59
        disp([TASK, ' Training TASK...']);
        eval(TASK); % Run the task
        BLOCK = BLOCK + 1;
    end

    % Actual task
    TRAINING = 0;
    blockCount = 6;
    start = 1;
    for i = blockCount:-1:1
        if isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_', TASK, '_block', num2str(i), '_task.mat']])
            start = i + 1;
            break;
        end
    end

    for BLOCK = start:blockCount
        disp([TASK, ' STARTING...']);
        eval(TASK); % Run the task
    end
end

%% Allow keyboard input into Matlab code
ListenChar(0);

%% Load and display accuracy
try
AOC_display_accuracy;
catch 
end

%% Display payment plan
try
AOC_display_cash;
catch
end
