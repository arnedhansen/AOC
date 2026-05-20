%% AOC Master Matrix — Sternberg
% Merges four sources by `ID`,`Condition`:
%   (a) behavioral
%   (b) gaze
%   (c) EEG non-FOOOF (`AOC_eeg_matrix_sternberg.mat`)
%   (d) EEG FOOOF-only (`AOC_eeg_matrix_sternberg_FOOOF.mat`)
% Includes key-uniqueness checks before merge to prevent row inflation.
%
% Key outputs:
%   merged_data_sternberg.mat 
%   merged_data_sternberg.csv
%
% Columns include:
%   EEG — Raw alpha power (IAF band from multitaper IAF, occ channels):
%     AlphaPower_raw_early / _late / _full
%     AlphaPower_bl_early / _late / _full   (dB baseline [-0.5 -0.25]s)
%     IAF (concat [1 2]s retention, mtmfft+DPSS peak); IAF_specParam (FOOOF alpha CF, median occ)
%
%   EEG — FOOOF alpha (IAF band, occ channels):
%     AlphaPower_FOOOF_full          [0 2]s, no baseline
%     AlphaPower_FOOOF_bl_full       [0 2]s, baselined absolute
%     AlphaPower_FOOOF_bl_early      [0 1]s, baselined
%     AlphaPower_FOOOF_bl_late       [1 2]s, baselined
%
%   Gaze — baselined (BL window [-0.5 -0.25]s; % for GD/SPL/MS/BCEA,pupil):
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
load(fullfile(featPath, 'AOC_behavioral_matrix_sternberg.mat'));  % behav_data_sternberg

% EEG (split products)
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg.mat'));         % eeg_data_sternberg
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg_FOOOF.mat'));   % eeg_data_sternberg_FOOOF

% Gaze
load(fullfile(featPath, 'AOC_gaze_matrix_sternberg.mat'));        % gaze_data_sternberg

%% Merge demographics into behavioral struct
demoIDs = [demog.ID];
for i = 1:numel(behav_data_sternberg)
    idx = find(demoIDs == behav_data_sternberg(i).ID, 1);
    behav_data_sternberg(i).Gender          = demog(idx).Gender;
    behav_data_sternberg(i).Age             = demog(idx).Age;
    behav_data_sternberg(i).Handedness      = demog(idx).Handedness;
    behav_data_sternberg(i).OcularDominance = demog(idx).OcularDominance;
end

%% Build merged struct
behav_table = struct2table(behav_data_sternberg);
gaze_table  = struct2table(gaze_data_sternberg);
eeg_table   = struct2table(eeg_data_sternberg);
eeg_fooof_table = struct2table(eeg_data_sternberg_FOOOF);

assert_unique_keys(eeg_table, {'ID', 'Condition'}, 'eeg_data_sternberg');
if height(eeg_fooof_table) > 0
    assert_unique_keys(eeg_fooof_table, {'ID', 'Condition'}, 'eeg_data_sternberg_FOOOF');
end

% Merge by ID and Condition; keep all existing columns from input matrices.
merged_table = outerjoin(behav_table, gaze_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_table = outerjoin(merged_table, eeg_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
if height(eeg_fooof_table) > 0
    merged_table = outerjoin(merged_table, eeg_fooof_table, ...
        'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
else
    warning('AOC_master_matrix_sternberg:EmptyFOOOF', ...
        ['eeg_data_sternberg_FOOOF is empty; FOOOF columns are set to NaN. ', ...
        'Re-run AOC_eeg_fex_sternberg_TFR.m after non-FOOOF EEG (IAF CSV).']);
    merged_table.AlphaPower_FOOOF_full = nan(height(merged_table), 1);
    merged_table.AlphaPower_FOOOF_bl_full = nan(height(merged_table), 1);
    merged_table.AlphaPower_FOOOF_bl_early = nan(height(merged_table), 1);
    merged_table.AlphaPower_FOOOF_bl_late = nan(height(merged_table), 1);
end
merged_data_sternberg = table2struct(merged_table);

%% Save as .mat
save(fullfile(featPath, 'AOC_merged_data_sternberg.mat'), 'merged_data_sternberg');
fprintf('Saved merged_data_sternberg.mat\n');

%% Save as .csv
writetable(merged_table, fullfile(featPath, 'AOC_merged_data_sternberg.csv'));
fprintf('Saved merged_data_sternberg.csv\n');

%% Summary
disp(merged_data_sternberg)

function assert_unique_keys(T, keyVars, tableName)
[~, ia] = unique(T(:, keyVars), 'rows', 'stable');
if numel(ia) ~= height(T)
    error('Duplicate key rows detected in %s for keys: %s', tableName, strjoin(keyVars, ', '));
end
end