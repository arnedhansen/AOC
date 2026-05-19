%% AOC Master Matrix — Sternberg (Trial-Level)
% Loads behavioral, EEG, gaze trial matrices and demographics, inner-joins
% on ID/Trial/Condition. Appends subject-level FOOOF columns from
% `AOC_eeg_matrix_sternberg_FOOOF.mat` (repeated per trial by ID/Condition).
% Produces merged_data_sternberg_trials.mat.
%
% Key outputs:
%   merged_data_sternberg_trials.mat (table: trial-wise behav, EEG, gaze, demographics, FOOOF alpha)

%% Setup
clear
clc
close all
startup
[~, paths, ~, ~] = setup('AOC', 0);
featPath = paths.features;

%% Load data
% Demographics from methlab_vp
demog_data_sternberg = readtable(paths.vp_table);
demog_data_sternberg = demog_data_sternberg(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_sternberg = table2struct(demog_data_sternberg(1:120, :));

% Behavioral
load(fullfile(featPath, 'AOC_behavioral_matrix_sternberg_trials.mat'));

% Gaze
load(fullfile(featPath, 'AOC_gaze_matrix_sternberg_trials.mat'));

% EEG
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg_trials.mat'));

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

%% Rename variables
merged_data_sternberg_trials.Properties.VariableNames{'Alter'} = 'Age';
merged_data_sternberg_trials.Properties.VariableNames{'H_ndigkeit'} = 'Handedness';

%% Add FOOOF alpha power (subject-level, repeated per trial)
% Uses split FOOOF EEG matrix from TFR pipeline.
nTrials = height(merged_data_sternberg_trials);
AlphaPower_FOOOF          = nan(nTrials, 1);
AlphaPower_FOOOF_bl       = nan(nTrials, 1);
AlphaPower_FOOOF_bl_early = nan(nTrials, 1);
AlphaPower_FOOOF_bl_late  = nan(nTrials, 1);
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg_FOOOF.mat')); % eeg_data_sternberg_FOOOF
T_fooof = struct2table(eeg_data_sternberg_FOOOF);
assert_unique_keys(T_fooof, {'ID', 'Condition'}, 'eeg_data_sternberg_FOOOF');

for i = 1:height(T_fooof)
    rowIdx = merged_data_sternberg_trials.ID == T_fooof.ID(i) & ...
        merged_data_sternberg_trials.Condition == T_fooof.Condition(i);
    if ~any(rowIdx), continue; end
    AlphaPower_FOOOF(rowIdx) = T_fooof.AlphaPower_FOOOF(i);
    AlphaPower_FOOOF_bl(rowIdx) = T_fooof.AlphaPower_FOOOF_bl(i);
    AlphaPower_FOOOF_bl_early(rowIdx) = T_fooof.AlphaPower_FOOOF_bl_early(i);
    AlphaPower_FOOOF_bl_late(rowIdx) = T_fooof.AlphaPower_FOOOF_bl_late(i);
end

merged_data_sternberg_trials.AlphaPower_FOOOF          = AlphaPower_FOOOF;
merged_data_sternberg_trials.AlphaPower_FOOOF_bl       = AlphaPower_FOOOF_bl;
merged_data_sternberg_trials.AlphaPower_FOOOF_bl_early = AlphaPower_FOOOF_bl_early;
merged_data_sternberg_trials.AlphaPower_FOOOF_bl_late  = AlphaPower_FOOOF_bl_late;

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
    'MSRateFull', 'MSRateFullBL', ...
    'BCEAEarly', 'BCEAEarlyBL', ...
    'BCEALate', 'BCEALateBL', ...
    'BCEAFull', 'BCEAFullBL', ...
    'BCEALatEarly', 'BCEALatEarlyBL', ...
    'BCEALatLate', 'BCEALatLateBL', ...
    'BCEALatFull', 'BCEALatFullBL'}, ...
    {'AlphaPowerEarly', 'AlphaPowerEarlyBL', ...
    'AlphaPowerLate', 'AlphaPowerLateBL', ...
    'AlphaPowerFull', 'AlphaPowerFullBL', ...
    'IAF', 'Lateralization'}, ...
    {'AlphaPower_FOOOF', 'AlphaPower_FOOOF_bl', ...
    'AlphaPower_FOOOF_bl_early', 'AlphaPower_FOOOF_bl_late'}];

merged_data_sternberg_trials = merged_data_sternberg_trials(:, newOrder);

%% Save as .mat
save(fullfile(featPath, 'AOC_merged_data_sternberg_trials.mat'), 'merged_data_sternberg_trials');

%% Save as .csv
writetable(merged_data_sternberg_trials, fullfile(featPath, 'AOC_merged_data_sternberg_trials.csv'));

function assert_unique_keys(T, keyVars, tableName)
[~, ia] = unique(T(:, keyVars), 'rows', 'stable');
if numel(ia) ~= height(T)
    error('Duplicate key rows detected in %s for keys: %s', tableName, strjoin(keyVars, ', '));
end
end