%% AOC MASTER Matrix Sternberg

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
demog_data_sternberg = readtable('/Volumes/methlab_vp/OCC/AOC/AOC_VPs.xlsx');
demog_data_sternberg = demog_data_sternberg(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_sternberg = table2struct(demog_data_sternberg(1:120, :));

% Behavioral
load('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_sternberg.mat');

% EEG
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg.mat');

% Gaze
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg.mat');


%% MERGE ON ACTUAL GLOBAL TRIAL IDS


safas

% convert structs to tables and inner-join on the canonical key
T_eeg  = struct2table(subj_data_eeg_trials);
T_gaze = struct2table(subj_data_gaze_trial);

T_join = innerjoin(T_eeg, T_gaze, ...
    'Keys', {'ID','GlobalTrial','Condition'}, ...
    'RightVariables', setdiff(T_gaze.Properties.VariableNames, {'ID','GlobalTrial','Condition','Trial'}));

% final guard
assert(height(T_join) == height(T_eeg), 'EEG–ET join dropped trials; investigate mismatched keys.');




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
    'AlphaPower', {eeg_data_sternberg.AlphaPower}, ...
    'IAF', {eeg_data_sternberg.IAF},...
    'Lateralization', {eeg_data_sternberg.Lateralization});

%% Save as .mat
save /Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.mat merged_data_sternberg

%% Save as .csv
merged_table = struct2table(merged_data_sternberg);
csv_filename = '/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv';
writetable(merged_table, csv_filename);