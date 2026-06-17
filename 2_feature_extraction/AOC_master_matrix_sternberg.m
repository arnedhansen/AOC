%% AOC Master Matrix — Sternberg
% Merges three sources by `ID`,`Condition`:
%   (a) behavioral
%   (b) gaze
%   (c) EEG (`AOC_eeg_matrix_sternberg.mat`)
% Includes key-uniqueness checks before merge to prevent row inflation.
%
% Key outputs:
%   merged_data_sternberg.mat 
%   merged_data_sternberg.csv
%
% Columns include:
%   EEG — Raw alpha power (IAF band from multitaper IAF, occ channels):
%     AlphaPower_raw_early / _late / _full
%     AlphaPower_bl_early / _late / _full   (dB baseline [-1.5 -0.5]s)
%     ERSD_early / ERSD_late / ERSD_full — fixed [8 14] Hz dB on tfr*_bl, O/I occ channels, windows [0 1] / [1 2] / [0 2]s
%         (AOC_eeg_fex_sternberg_TFR.m → AOC_eeg_matrix_sternberg_ERSD.mat)
%     IAF (concat [1 2]s retention, mtmfft+DPSS peak)
%
%   Gaze — baselined (BL window [-1.5 -0.5]s; % for GD/SPL/MS/BCEA,pupil):
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

% EEG
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg.mat'));         % eeg_data_sternberg

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
eeg_table   = join_ersd_table(eeg_table, featPath, 'AOC_eeg_matrix_sternberg_ERSD.mat', 'eeg_data_sternberg_ERSD');

assert_unique_keys(eeg_table, {'ID', 'Condition'}, 'eeg_data_sternberg');

% Merge by ID and Condition; keep all existing columns from input matrices.
merged_table = outerjoin(behav_table, gaze_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_table = outerjoin(merged_table, eeg_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
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

function T = join_ersd_table(T, featPath, ersdMatName, ersdVarName)
ersdMatPath = fullfile(featPath, ersdMatName);
if ~isfile(ersdMatPath)
    return
end
S = load(ersdMatPath, ersdVarName);
ersd_data = S.(ersdVarName);
if isempty(ersd_data)
    return
end
Te = struct2table(ersd_data);
if height(Te) == 0
    return
end
assert_unique_keys(Te, {'ID', 'Condition'}, ersdVarName);
for cn = {'ERSD_early', 'ERSD_late', 'ERSD_full'}
    if ismember(cn{1}, T.Properties.VariableNames)
        T.(cn{1}) = [];
    end
end
T = outerjoin(T, Te, 'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
end