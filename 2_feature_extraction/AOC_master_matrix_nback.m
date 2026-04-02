%% AOC Master Matrix — N-Back
% Builds merged_data_nback with:
%   (a) FOOOF-based alpha power (non-baselined, baselined full/early/late)
%   (b) Baselined gaze metrics (GazeDeviation, ScanPathLength, PupilSize,
%       MSRate, BCEA, BCEALateralization — averaged across trials per condition, dB or % change)
%   (c) Baselined raw alpha power (full/early/late, dB vs [-0.5 -0.25]s)
%
% All values are loaded from existing files — no recomputation needed.
%   FOOOF alpha  → per-subject power_nback_fooof.mat
%   Baselined gaze and alpha → AOC_gaze_matrix_nback.mat / AOC_eeg_matrix_nback.mat
%
% Key outputs:
%   merged_data_nback.mat
%   merged_data_nback.csv
%
% Columns include:
%   EEG — FOOOF alpha [8-14 Hz], occ channels:
%     AlphaPower_FOOOF           [0 2]s, no baseline
%     AlphaPower_FOOOF_bl        [0 2]s, baselined [-0.5 -0.25]s absolute
%     AlphaPower_FOOOF_bl_early  [0 1]s, baselined
%     AlphaPower_FOOOF_bl_late   [1 2]s, baselined
%
%   EEG — Raw alpha power (baselined, dB vs [-0.5 -0.25]s):
%     AlphaPower_bl              [0 2]s, baselined
%     AlphaPower_bl_early        [0 1]s, baselined
%     AlphaPower_bl_late         [1 2]s, baselined
%
%   Gaze — baselined (BL window [-0.5 -0.25]s; % change for all gaze metrics):
%     GazeDeviationFullBL / EarlyBL / LateBL
%     ScanPathLengthFullBL / EarlyBL / LateBL
%     PupilSizeFullBL / EarlyBL / LateBL          (% change)
%     MSRateFullBL / EarlyBL / LateBL
%     BCEAFullBL / BCEAEarlyBL / BCEALateBL
%     BCEALatFullBL / BCEALatEarlyBL / BCEALatLateBL

%% Setup
clear
clc
close all
startup
[~, paths, ~, ~] = setup('AOC', 0);
featPath = paths.features;

%% Load data
% Demographics
demog = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog = demog(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog = table2struct(demog(1:120, :));

% Behavioral
load(fullfile(featPath, 'AOC_behavioral_matrix_nback.mat'));  % behav_data_nback

% EEG
load(fullfile(featPath, 'AOC_eeg_matrix_nback.mat'));         % eeg_data_nback

% Gaze
load(fullfile(featPath, 'AOC_gaze_matrix_nback.mat'));        % gaze_data_nback

%% Merge demographics into behavioral struct
demoIDs = [demog.ID];
for i = 1:numel(behav_data_nback)
    idx = find(demoIDs == behav_data_nback(i).ID, 1);
    behav_data_nback(i).Gender          = demog(idx).Gender;
    behav_data_nback(i).Alter           = demog(idx).Alter;
    behav_data_nback(i).H_ndigkeit      = demog(idx).H_ndigkeit;
    behav_data_nback(i).OcularDominance = demog(idx).OcularDominance;
end

%% Build merged struct
behav_table = struct2table(behav_data_nback);
gaze_table  = struct2table(gaze_data_nback);
eeg_table   = struct2table(eeg_data_nback);

% Merge by ID and Condition; keep all existing columns from input matrices.
merged_table = outerjoin(behav_table, gaze_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_table = outerjoin(merged_table, eeg_table, ...
    'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
merged_data_nback = table2struct(merged_table);

%% Save as .mat
save(fullfile(featPath, 'AOC_merged_data_nback.mat'), 'merged_data_nback');
fprintf('Saved merged_data_nback.mat\n');

%% Save as .csv
writetable(merged_table, fullfile(featPath, 'AOC_merged_data_nback.csv'));
fprintf('Saved merged_data_nback.csv\n');

%% Summary
nCols = numel(fieldnames(merged_data_nback));
fprintf('\n=== N-Back merge complete ===\n');
fprintf('Rows: %d  |  Columns: %d\n', numel(merged_data_nback), nCols);
fprintf('Total columns:  %d\n', nCols);
