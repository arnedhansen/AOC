%% AOC MASTER Matrix Sternberg

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
demog_data_sternberg = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog_data_sternberg = demog_data_sternberg(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_sternberg = table2struct(demog_data_sternberg(1:120, :));

% Behavioral
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/behavioral_matrix_sternberg_trials.mat');

% Gaze
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/gaze_matrix_sternberg_trials.mat');
gaze_data_sternberg_trials = rmfield(gaze_data_sternberg_trials, 'ScanPathSeriesT');
gaze_data_sternberg_trials = rmfield(gaze_data_sternberg_trials, 'ScanPathSeries');

% EEG
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat');

%% Merge structures
%  based on global trial IDs

% Add demographic infos to behavioral structure
demoIDs = [demog_data_sternberg.ID];
for i = 1:numel(behav_data_sternberg_trials)
    idx = find(demoIDs == behav_data_sternberg_trials(i).ID, 1);
    behav_data_sternberg_trials(i).Gender           = demog_data_sternberg(idx).Gender;
    behav_data_sternberg_trials(i).Alter            = demog_data_sternberg(idx).Alter;
    behav_data_sternberg_trials(i).H_ndigkeit       = demog_data_sternberg(idx).H_ndigkeit;
    behav_data_sternberg_trials(i).OcularDominance  = demog_data_sternberg(idx).OcularDominance;
end

% convert structs to tables
T_behav = struct2table(behav_data_sternberg_trials);
T_gaze = struct2table(gaze_data_sternberg_trials);
T_eeg  = struct2table(eeg_data_sternberg_trials);

% inner-join on key variables ID, Trial and Condition
mergeEEGxBehav = innerjoin(T_eeg, T_behav, 'Keys', {'ID','Trial','Condition'});
merged_data_sternberg_trials = innerjoin(mergeEEGxBehav, T_gaze, 'Keys', {'ID','Trial','Condition'});

%% Re-arrange table
newOrder = [ ...
    {'ID', 'Trial', 'Condition'}, ...
    {'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'}, ...
    {'Accuracy', 'ReactionTime', 'Stimuli', 'Probe', 'Match'}, ...
    {'GazeDeviationEarly', 'GazeDeviationEarlyBL', ...
    'GazeDeviationLate', 'GazeDeviationLateBL', ...
    'ScanPathLengthEarly', 'ScanPathLengthEarlyBL', ...
    'ScanPathLengthLate', 'ScanPathLengthLateBL', ...
    'PupilSizeEarly', 'PupilSizeEarlyBL', ...
    'PupilSizeLate', 'PupilSizeLateBL', ...
    'MSRateEarly', 'MSRateEarlyBL', ...
    'MSRateLate', 'MSRateLateBL'}, ...
    {'AlphaPowerEarly', 'AlphaPowerEarlyBL', ...
    'AlphaPowerLate', 'AlphaPowerLateBL', ...
    'IAF', 'Lateralization'}];

% Optionally overwrite
merged_data_sternberg_trials = merged_data_sternberg_trials(:, newOrder)

%% Rename variables
merged_data_sternberg_trials.Properties.VariableNames{'Alter'} = 'Age';
merged_data_sternberg_trials.Properties.VariableNames{'H_ndigkeit'} = 'Handedness';

%% Save as .mat
save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat merged_data_sternberg_trials

%% Save as .csv
csv_filename = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.csv';
writetable(merged_data_sternberg_trials, csv_filename);