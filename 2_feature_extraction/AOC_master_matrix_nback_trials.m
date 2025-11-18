%% AOC MASTER Matrix N-back

%% Setup
clear
clc
close all
path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/automagic_nohp';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
% Demographics from methlab_vp
demog_data_nback = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog_data_nback = demog_data_nback(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_nback = table2struct(demog_data_nback(1:120, :));

% Behavioral
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/behavioral_matrix_nback_trials.mat');

% Gaze
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/gaze_matrix_nback_trials.mat');

% EEG
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_nback_trials.mat');

%% Merge structures
%  based on global trial IDs

% Add demographic infos to behavioral structure
demoIDs = [demog_data_nback.ID];
for i = 1:numel(behav_data_nback_trials)
    idx = find(demoIDs == behav_data_nback_trials(i).ID, 1);
    behav_data_nback_trials(i).Gender           = demog_data_nback(idx).Gender;
    behav_data_nback_trials(i).Alter            = demog_data_nback(idx).Alter;
    behav_data_nback_trials(i).H_ndigkeit       = demog_data_nback(idx).H_ndigkeit;
    behav_data_nback_trials(i).OcularDominance  = demog_data_nback(idx).OcularDominance;
end

% convert structs to tables
T_behav = struct2table(behav_data_nback_trials);
T_gaze = struct2table(gaze_data_nback_trials);
T_eeg  = struct2table(eeg_data_nback_trials);

% inner-join on key variables ID, Trial and Condition
mergeEEGxBehav = innerjoin(T_eeg, T_behav, 'Keys', {'ID','Trial','Condition'});
merged_data_nback_trials = innerjoin(mergeEEGxBehav, T_gaze, 'Keys', {'ID','Trial','Condition'});

%% Rename variables
merged_data_nback_trials.Properties.VariableNames{'Alter'} = 'Age';
merged_data_nback_trials.Properties.VariableNames{'H_ndigkeit'} = 'Handedness';

%% Re-arrange table
newOrder = [ ...
    {'ID', 'Trial', 'Condition'}, ...
    {'Gender', 'Age', 'Handedness', 'OcularDominance'}, ...
    {'Accuracy', 'ReactionTime', 'Stimuli', 'Probe', 'Match'}, ...
    {'GazeDeviationEarly', 'GazeDeviationEarlyBL', ...
    'GazeDeviationLate', 'GazeDeviationLateBL', ...
    'GazeDeviationFull', 'GazeDeviationFullBL', ...
    'ScanPathLengthEarly', 'ScanPathLengthEarlyBL', ...
    'ScanPathLengthLate', 'ScanPathLengthLateBL', ...
    'ScanPathLengthFull', 'ScanPathLengthFullBL', ...
    'PupilSizeEarly', 'PupilSizeEarlyBL', ...
    'PupilSizeLate', 'PupilSizeLateBL', ...
    'PupilSizeFull', 'PupilSizeFullBL', ...
    'MSRateEarly', 'MSRateEarlyBL', ...
    'MSRateLate', 'MSRateLateBL', ...
    'MSRateFull', 'MSRateFullBL'}, ...
    {'AlphaPowerEarly', 'AlphaPowerEarlyBL', ...
    'AlphaPowerLate', 'AlphaPowerLateBL', ...
    'AlphaPowerFull', 'AlphaPowerFullBL', ...
    'IAF', 'Lateralization'}]; 

merged_data_nback_trials = merged_data_nback_trials(:, newOrder)

%% Save as .mat
save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback_trials.mat merged_data_nback_trials

%% Save as .csv
csv_filename = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback_trials.csv';
writetable(merged_data_nback_trials, csv_filename);