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
clear all;
close all;
clc;
rootFilepath = pwd; % Retrieve the present working directory

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
if ~isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_Resting.mat'])
    Resting_EEG;
end

%% Randomize order of Sternberg Task and NBack Task
% Use subject ID for assignment of pseudorandom task Order (Sternberg & N-back)
if mod(str2double(subjectID), 2) == 0
    taskOrder = 'NBackSternberg'; % Even subjectIDs do Nback first
else
    taskOrder = 'SternbergNBack'; % Uneven subjectIDs do Sternberg first
end

% Shuffle n-back conditions
rng();
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
cd(DATA_PATH)
for block = 1:6
    % Load Sternberg accuracy data
    tsk = 'AOC_Sternberg';
    fileName = [subjectID '_', tsk, '_block' num2str(block) '_task.mat'];
    dataSternberg = load(['/home/methlab/Desktop/AOC_data/', subjectID, '/', fileName]);
    percCorr(block, 1) = dataSternberg.saves.data.percentTotalCorrect;

    % Load Nback accuracy data
    tsk = 'AOC_NBack';
    dataDir = ['/home/methlab/Desktop/AOC_data/', subjectID, '/'];
    filePattern = [subjectID, '_', tsk, '_block', num2str(block), '_*_task.mat'];
    files = dir(fullfile(dataDir, filePattern));
    fileName = files.name;
    dataNback = load(['/home/methlab/Desktop/AOC_data/', subjectID, '/', fileName]);
    percCorr(block, 2) = dataNback.saves.data.percentTotalCorrect;
end

% Create a UI figure
uiFig = uifigure('Name', 'Task Accuracy', 'Position', [500, 400, 350, 300]);

% Define the tasks
tasks = {'Sternberg', 'NBack'};

% Round the percentage accuracy and convert to character arrays with a percentage sign
roundedPercCorr = round(percCorr); % Round to the nearest integer
percStrings = strcat(string(roundedPercCorr), '%'); % Append percentage sign
percStrings = cellstr(percStrings); % Convert to cell array of character arrays

% Create a vertical spacing offset
yOffset = 250;
rowHeight = 30;

% Add column headers for tasks
headerLabel1 = uilabel(uiFig, ...
    'Text', tasks{1}, ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'Position', [100, yOffset, 100, 20]);

headerLabel2 = uilabel(uiFig, ...
    'Text', tasks{2}, ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'Position', [200, yOffset, 100, 20]);

% Loop through each block to display the data
for block = 1:6
    % Block number
    blockLabel = uilabel(uiFig, ...
        'Text', ['Block ' num2str(block)], ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'left', ...
        'Position', [20, yOffset - block * rowHeight, 100, 20]);
    
    % NBack accuracy
    nBackLabel = uilabel(uiFig, ...
        'Text', percStrings{block, 1}, ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'center', ...
        'Position', [120, yOffset - block * rowHeight, 80, 20]);
    
    % Sternberg accuracy
    sternbergLabel = uilabel(uiFig, ...
        'Text', percStrings{block, 2}, ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'center', ...
        'Position', [220, yOffset - block * rowHeight, 80, 20]);
end
