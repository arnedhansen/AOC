%% AOC MASTER Matrix Nback

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/automagic_nohp';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
% Demographics from methlab_vp
demog_data_nback = readtable('/Volumes/methlab_vp/OCC/AOC/AOC_VPs.xlsx');
demog_data_nback = demog_data_nback(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_nback = table2struct(demog_data_nback(1:120, :));

% Behavioral
load('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_nback.mat');

% EEG
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_nback.mat');

% Gaze
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_nback.mat');

%% Merge structures
demoIDs = [demog_data_nback.ID];
for i = 1:numel(behav_data_nback)
    idx = find(demoIDs == behav_data_nback(i).ID, 1);
    behav_data_nback(i).Gender           = demog_data_nback(idx).Gender;
    behav_data_nback(i).Alter            = demog_data_nback(idx).Alter;
    behav_data_nback(i).H_ndigkeit       = demog_data_nback(idx).H_ndigkeit;
    behav_data_nback(i).OcularDominance  = demog_data_nback(idx).OcularDominance;
end

% Merge all
merged_data_nback = struct( ...
    'ID', {behav_data_nback.ID}, ...
    'Gender', {behav_data_nback.Gender}, ...
    'Age', {behav_data_nback.Alter}, ...
    'Handedness', {behav_data_nback.H_ndigkeit}, ...
    'OcularDominance', {behav_data_nback.OcularDominance}, ...
    'Condition', {behav_data_nback.Condition}, ...
    'Accuracy', {behav_data_nback.Accuracy}, ...
    'ReactionTime', {behav_data_nback.ReactionTime}, ...
    'GazeDeviation', {gaze_data_nback.GazeDeviation}, ...
    'GazeStdX', {gaze_data_nback.GazeStdX}, ...
    'GazeStdY', {gaze_data_nback.GazeStdY}, ...
    'PupilSize', {gaze_data_nback.PupilSize}, ...
    'MSRate', {gaze_data_nback.MSRate}, ...
    'Blinks', {gaze_data_nback.Blinks}, ...
    'Fixations', {gaze_data_nback.Fixations}, ...
    'Saccades', {gaze_data_nback.Saccades}, ...
    'AlphaPower', {eeg_data_nback.AlphaPower}, ...
    'IAF', {eeg_data_nback.IAF},...
    'Lateralization', {eeg_data_nback.Lateralization});

%% Save as .mat
save /Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.mat merged_data_nback

%% Save as .csv
merged_table = struct2table(merged_data_nback);
csv_filename = '/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv';
writetable(merged_table, csv_filename);