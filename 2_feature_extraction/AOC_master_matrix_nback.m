%% AOC Master Matrix — N-Back
% Merges four sources by `ID`,`Condition`:
%   (a) behavioral
%   (b) gaze
%   (c) EEG non-FOOOF (`AOC_eeg_matrix_nback.mat`)
%   (d) EEG FOOOF-only (`AOC_eeg_matrix_nback_FOOOF.mat`)
% Includes key-uniqueness checks before merge to prevent row inflation.
%
% Key outputs:
%   merged_data_nback.mat
%   merged_data_nback.csv
%
% Columns include:
%   EEG — Raw alpha power (IAF band from multitaper IAF, occ channels):
%     AlphaPower_raw_early / _late / _full
%     AlphaPower_bl_early / _late / _full   (dB baseline [-0.5 -0.25]s)
%     IAF (concat [0 2]s retention, mtmfft+DPSS peak); IAF_specParam (FOOOF alpha CF, median occ)
%
%   EEG — FOOOF alpha (IAF band, occ channels):
%     AlphaPower_FOOOF_full          [0 2]s, no baseline
%     AlphaPower_FOOOF_bl_full       [0 2]s, baselined absolute
%     AlphaPower_FOOOF_bl_early      [0 1]s, baselined
%     AlphaPower_FOOOF_bl_late       [1 2]s, baselined
%
%   Gaze — baselined (BL window [-0.5 -0.25]s; % change for all gaze metrics):
%     GazeDeviationFullBL / EarlyBL / LateBL
%     ScanPathLengthFullBL / EarlyBL / LateBL
%     PupilSizeFullBL / EarlyBL / LateBL          (% change)
%     MSRateFullBL / EarlyBL / LateBL
%     BCEAFullBL / BCEAEarlyBL / BCEALateBL
%     BCEALatFullBL / BCEALatEarlyBL / BCEALatLateBL

%% Setup
startup
[~, paths, ~, ~] = setup('AOC', 0);
featPath = paths.features;

%% Load data
% Demographics
demog = readtable(paths.vp_table);
demog = demog(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog = renamevars(demog, {'Alter', 'H_ndigkeit'}, {'Age', 'Handedness'});
demog = table2struct(demog(1:120, :));

% Behavioral
load(fullfile(featPath, 'AOC_behavioral_matrix_nback.mat'));  % behav_data_nback

% EEG (split products)
load(fullfile(featPath, 'AOC_eeg_matrix_nback.mat'));         % eeg_data_nback
load(fullfile(featPath, 'AOC_eeg_matrix_nback_FOOOF.mat'));   % eeg_data_nback_FOOOF

% Gaze
load(fullfile(featPath, 'AOC_gaze_matrix_nback.mat'));        % gaze_data_nback

%% Merge demographics into behavioral struct
demoIDs = [demog.ID];
for i = 1:numel(behav_data_nback)
    idx = find(demoIDs == behav_data_nback(i).ID, 1);
    behav_data_nback(i).Gender          = demog(idx).Gender;
    behav_data_nback(i).Age             = demog(idx).Age;
    behav_data_nback(i).Handedness      = demog(idx).Handedness;
    behav_data_nback(i).OcularDominance = demog(idx).OcularDominance;
end

%% Build merged struct
behav_table = struct2table(behav_data_nback);
gaze_table  = struct2table(gaze_data_nback);
eeg_table   = struct2table(eeg_data_nback);
eeg_fooof_table = struct2table(eeg_data_nback_FOOOF);

assert_unique_keys(eeg_table, {'ID', 'Condition'}, 'eeg_data_nback');
assert_unique_keys(eeg_fooof_table, {'ID', 'Condition'}, 'eeg_data_nback_FOOOF');

% Merge by ID and Condition; keep all existing columns from input matrices.
merged_table = outerjoin(behav_table, gaze_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_table = outerjoin(merged_table, eeg_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_table = outerjoin(merged_table, eeg_fooof_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_data_nback = table2struct(merged_table);

%% Save as .mat
save(fullfile(featPath, 'AOC_merged_data_nback.mat'), 'merged_data_nback');
fprintf('Saved merged_data_nback.mat\n');

%% Save as .csv
writetable(merged_table, fullfile(featPath, 'AOC_merged_data_nback.csv'));
fprintf('Saved merged_data_nback.csv\n');

%% Summary
disp(merged_data_nback)

function assert_unique_keys(T, keyVars, tableName)
[~, ia] = unique(T(:, keyVars), 'rows', 'stable');
if numel(ia) ~= height(T)
    error('Duplicate key rows detected in %s for keys: %s', tableName, strjoin(keyVars, ', '));
end
end