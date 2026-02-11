%% AOC Master Matrix — N-Back (Trial-Level)
% Loads behavioral, EEG, gaze trial matrices and demographics, inner-joins
% on ID/Trial/Condition. Adds subject-level FOOOF alpha power (repeated per
% trial). Produces merged_data_nback_trials.mat.
%
% Key outputs:
%   merged_data_nback_trials.mat (table: trial-wise behav, EEG, gaze, demographics, FOOOF alpha)

%% Setup
clear
clc
close all
featPath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';

%% Load data
% Demographics from methlab_vp
demog_data_nback = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog_data_nback = demog_data_nback(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog_data_nback = table2struct(demog_data_nback(1:120, :));

% Behavioral
load(fullfile(featPath, 'behavioral_matrix_nback_trials.mat'));

% Gaze
load(fullfile(featPath, 'gaze_matrix_nback_trials.mat'));

% EEG
load(fullfile(featPath, 'eeg_matrix_nback_trials.mat'));

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

%% Add FOOOF alpha power (subject-level, repeated per trial)
% Loads per-subject power_nback_fooof.mat; extracts scalar alpha [8-14 Hz]
% averaged over occipital channels for each condition.
alphaRange = [8 14];
nTrials = height(merged_data_nback_trials);
AlphaPower_FOOOF          = nan(nTrials, 1);
AlphaPower_FOOOF_bl       = nan(nTrials, 1);
AlphaPower_FOOOF_bl_early = nan(nTrials, 1);
AlphaPower_FOOOF_bl_late  = nan(nTrials, 1);

uIDs = unique(merged_data_nback_trials.ID);
for s = 1:numel(uIDs)
    subjID  = uIDs(s);
    subjStr = num2str(subjID);

    fooof_file = fullfile(featPath, subjStr, 'eeg', 'power_nback_fooof.mat');
    if ~isfile(fooof_file)
        warning('Missing power_nback_fooof.mat for subject %s — skipping.', subjStr);
        continue
    end
    fooof = load(fooof_file);

    % Occipital channels (labels containing 'O' or 'I')
    labels  = fooof.pow1_fooof.label;
    occ_idx = cellfun(@(l) contains(l, 'O') || contains(l, 'I'), labels);

    % Alpha frequency indices
    freqs     = fooof.pow1_fooof.freq;
    alpha_idx = freqs >= alphaRange(1) & freqs <= alphaRange(2);

    condPows          = {fooof.pow1_fooof,          fooof.pow2_fooof,          fooof.pow3_fooof};
    condPows_bl       = {fooof.pow1_fooof_bl,       fooof.pow2_fooof_bl,       fooof.pow3_fooof_bl};
    condPows_bl_early = {fooof.pow1_fooof_bl_early, fooof.pow2_fooof_bl_early, fooof.pow3_fooof_bl_early};
    condPows_bl_late  = {fooof.pow1_fooof_bl_late,  fooof.pow2_fooof_bl_late,  fooof.pow3_fooof_bl_late};
    condVals          = [1, 2, 3];

    for c = 1:3
        rowIdx = merged_data_nback_trials.ID == subjID & ...
                 merged_data_nback_trials.Condition == condVals(c);
        if ~any(rowIdx), continue; end

        AlphaPower_FOOOF(rowIdx)          = mean(mean(condPows{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl(rowIdx)       = mean(mean(condPows_bl{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl_early(rowIdx) = mean(mean(condPows_bl_early{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl_late(rowIdx)  = mean(mean(condPows_bl_late{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
    end
end

merged_data_nback_trials.AlphaPower_FOOOF          = AlphaPower_FOOOF;
merged_data_nback_trials.AlphaPower_FOOOF_bl       = AlphaPower_FOOOF_bl;
merged_data_nback_trials.AlphaPower_FOOOF_bl_early = AlphaPower_FOOOF_bl_early;
merged_data_nback_trials.AlphaPower_FOOOF_bl_late  = AlphaPower_FOOOF_bl_late;

%% Re-arrange table
newOrder = [ ...
    {'ID', 'Trial', 'Condition'}, ...
    {'Gender', 'Age', 'Handedness', 'OcularDominance'}, ...
    {'Accuracy', 'ReactionTime', 'Stimuli', 'Match'}, ...
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
    'IAF', 'Lateralization'}, ...
    {'AlphaPower_FOOOF', 'AlphaPower_FOOOF_bl', ...
    'AlphaPower_FOOOF_bl_early', 'AlphaPower_FOOOF_bl_late'}];

merged_data_nback_trials = merged_data_nback_trials(:, newOrder);

%% Save as .mat
save(fullfile(featPath, 'merged_data_nback_trials.mat'), 'merged_data_nback_trials');

%% Save as .csv
writetable(merged_data_nback_trials, fullfile(featPath, 'merged_data_nback_trials.csv'));