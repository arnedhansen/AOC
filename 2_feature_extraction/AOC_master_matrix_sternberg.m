%% AOC Master Matrix — Sternberg (Subject-Level)
% Loads behavioral, EEG, gaze matrices and demographics, merges by ID. Produces merged_data_sternberg.mat. Runs: load and merge.
%
% Key outputs:
%   merged_data_sternberg.mat (struct: ID, demographics, Condition, Accuracy, RT, gaze and EEG features)

%% Setup
clear
clc
close all
path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/automagic';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
% Demographics from methlab_vp
demog_data_sternberg = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog_data_sternberg = demog_data_sternberg(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_sternberg = table2struct(demog_data_sternberg(1:120, :));

% Behavioral
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/behavioral_matrix_sternberg.mat');

% EEG
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_sternberg.mat');

% Gaze
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/gaze_matrix_sternberg.mat');

%% Merge structures
demoIDs = [demog_data_sternberg.ID];
for i = 1:numel(behav_data_sternberg)
    idx = find(demoIDs == behav_data_sternberg(i).ID, 1);
    behav_data_sternberg(i).Gender           = demog_data_sternberg(idx).Gender;
    behav_data_sternberg(i).Alter            = demog_data_sternberg(idx).Alter;
    behav_data_sternberg(i).H_ndigkeit       = demog_data_sternberg(idx).H_ndigkeit;
    behav_data_sternberg(i).OcularDominance  = demog_data_sternberg(idx).OcularDominance;
end

% Merge all
merged_data_sternberg = struct( ...
    'ID', {behav_data_sternberg.ID}, ...
    'Gender', {behav_data_sternberg.Gender}, ...
    'Age', {behav_data_sternberg.Alter}, ...
    'Handedness', {behav_data_sternberg.H_ndigkeit}, ...
    'OcularDominance', {behav_data_sternberg.OcularDominance}, ...
    'Condition', {behav_data_sternberg.Condition}, ...
    'Accuracy', {behav_data_sternberg.Accuracy}, ...
    'ReactionTime', {behav_data_sternberg.ReactionTime}, ...
    'GazeDeviation', {gaze_data_sternberg.GazeDeviation}, ...
    'GazeStdX', {gaze_data_sternberg.GazeStdX}, ...
    'GazeStdY', {gaze_data_sternberg.GazeStdY}, ...
    'PupilSize', {gaze_data_sternberg.PupilSize}, ...
    'MSRate', {gaze_data_sternberg.MSRate}, ...
    'Blinks', {gaze_data_sternberg.Blinks}, ...
    'Fixations', {gaze_data_sternberg.Fixations}, ...
    'Saccades', {gaze_data_sternberg.Saccades}, ...
    'ScanPathLength', {gaze_data_sternberg.ScanPathLength}, ...
    'ConvexHullArea', {gaze_data_sternberg.ConvexHullArea}, ...
    'AlphaPower', {eeg_data_sternberg.AlphaPower}, ...
    'IAF', {eeg_data_sternberg.IAF},...
    'Lateralization', {eeg_data_sternberg.Lateralization});

%% Save as .mat
save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg.mat merged_data_sternberg

%% Save as .csv
merged_table = struct2table(merged_data_sternberg);
csv_filename = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg.csv';
writetable(merged_table, csv_filename);