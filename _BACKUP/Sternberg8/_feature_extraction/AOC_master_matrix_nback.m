%% AOC MASTER Matrix Nback

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
load('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_nback.mat');
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_nback.mat');
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_nback.mat');

%% Sort behavioral structure
conds = [behav_data_nback.Condition];
[~, sortedIndices] = sort(conds);
behav = behav_data_nback(sortedIndices);

%% Merge structures
% Initialize the merged structure array with the same size as the original structures
merged_data_nback = struct('ID', {behav.ID}, ...
                               'TrialTRUE', {behav.Trial}, ...
                               'Trial', {gaze_data_nback.Trial}, ...
                               'Condition', {behav.Condition}, ...
                               'Accuracy', {behav.Accuracy}, ...
                               'ReactionTime', {behav.ReactionTime}, ...
                               'Stimuli', {behav.Stimuli}, ...
                               'Match', {behav.Match}, ...
                               'GazeDeviation', {gaze_data_nback.GazeDeviation}, ...
                               'PupilSize', {gaze_data_nback.PupilSize}, ...
                               'AlphaPower', {eeg_data_nback.AlphaPower}, ...
                               'IAF', {eeg_data_nback.IAF});

% Save
save /Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.mat merged_data_nback