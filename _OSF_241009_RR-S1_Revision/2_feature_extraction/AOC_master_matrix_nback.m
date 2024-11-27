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
                               'Condition', {behav.Condition}, ...
                               'Accuracy', {behav.Accuracy}, ...
                               'ReactionTime', {behav.ReactionTime}, ...
                               'GazeDeviation', {gaze_data_nback.GazeDeviation}, ...
                               'PupilSize', {gaze_data_nback.PupilSize}, ...
                               'AlphaPower', {eeg_data_nback.AlphaPower}, ...
                               'IAF', {eeg_data_nback.IAF});

%% Save as .mat
save /Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.mat merged_data_nback

%% Save as .csv
merged_table = struct2table(merged_data_nback);
csv_filename = '/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv';
writetable(merged_table, csv_filename);