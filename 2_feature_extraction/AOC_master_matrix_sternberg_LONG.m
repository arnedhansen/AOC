%% AOC Master Matrix LONG — Sternberg (Subject-Level)
% Extends the standard merged_data_sternberg with:
%   (a) FOOOF-based alpha power (non-baselined, baselined full/early/late)
%   (b) Baselined gaze metrics (GazeDeviation, ScanPathLength, PupilSize,
%       MSRate — averaged across trials per condition, dB or % change)
%
% All values are loaded from existing files — no recomputation needed.
%   FOOOF alpha  → per-subject power_stern_fooof.mat
%   Baselined gaze → aggregate gaze_matrix_sternberg_trials.mat
%
% Key outputs:
%   merged_data_sternberg_LONG.mat  (struct with original + new columns)
%   merged_data_sternberg_LONG.csv  (same, as table)
%
% New columns (25 → 41):
%   EEG — FOOOF alpha [8-14 Hz], occ channels:
%     AlphaPower_FOOOF           [0 2]s, no baseline
%     AlphaPower_FOOOF_bl        [0 2]s, baselined [-0.5 -0.25]s absolute
%     AlphaPower_FOOOF_bl_early  [0 1]s, baselined
%     AlphaPower_FOOOF_bl_late   [1 2]s, baselined
%
%   Gaze — baselined (BL window [-0.5 -0.25]s; dB for GD/SPL/MS, % for pupil):
%     GazeDeviationFullBL / EarlyBL / LateBL
%     ScanPathLengthFullBL / EarlyBL / LateBL
%     PupilSizeFullBL / EarlyBL / LateBL          (% change)
%     MSRateFullBL / EarlyBL / LateBL

%% Setup
clear
clc
close all
featPath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';

%% Load data
% Demographics
demog = readtable('/Volumes/g_psyplafor_methlab$/VP/OCC/AOC/AOC_VPs.xlsx');
demog = demog(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog = table2struct(demog(1:120, :));

% Behavioral
load(fullfile(featPath, 'behavioral_matrix_sternberg.mat'));  % behav_data_sternberg

% EEG (original: AlphaPower, IAF, Lateralization)
load(fullfile(featPath, 'eeg_matrix_sternberg.mat'));         % eeg_data_sternberg

% Gaze (original: subject-level, non-baselined)
load(fullfile(featPath, 'gaze_matrix_sternberg.mat'));        % gaze_data_sternberg

% Gaze (trial-level with baselined metrics)
load(fullfile(featPath, 'gaze_matrix_sternberg_trials.mat')); % gaze_data_sternberg_trials

%% Merge demographics into behavioral struct
demoIDs = [demog.ID];
for i = 1:numel(behav_data_sternberg)
    idx = find(demoIDs == behav_data_sternberg(i).ID, 1);
    behav_data_sternberg(i).Gender          = demog(idx).Gender;
    behav_data_sternberg(i).Alter           = demog(idx).Alter;
    behav_data_sternberg(i).H_ndigkeit      = demog(idx).H_ndigkeit;
    behav_data_sternberg(i).OcularDominance = demog(idx).OcularDominance;
end

%% Indexing helpers
allIDs   = [behav_data_sternberg.ID];
allConds = [behav_data_sternberg.Condition];
uIDs     = unique(allIDs);
nRows    = numel(behav_data_sternberg);
alphaRange = [8 14];

%% ===== (A) FOOOF alpha power per subject =====
AlphaPower_FOOOF          = nan(nRows, 1);
AlphaPower_FOOOF_bl       = nan(nRows, 1);
AlphaPower_FOOOF_bl_early = nan(nRows, 1);
AlphaPower_FOOOF_bl_late  = nan(nRows, 1);

for s = 1:numel(uIDs)
    subjID  = uIDs(s);
    subjStr = num2str(subjID);

    fooof_file = fullfile(featPath, subjStr, 'eeg', 'power_stern_fooof.mat');
    if ~isfile(fooof_file)
        warning('Missing power_stern_fooof.mat for subject %s — skipping.', subjStr);
        continue
    end
    fooof = load(fooof_file);

    % Occipital channels (labels containing 'O' or 'I')
    labels  = fooof.pow2_fooof.label;
    occ_idx = cellfun(@(l) contains(l, 'O') || contains(l, 'I'), labels);

    % Alpha frequency indices
    freqs     = fooof.pow2_fooof.freq;
    alpha_idx = freqs >= alphaRange(1) & freqs <= alphaRange(2);

    condPows          = {fooof.pow2_fooof,          fooof.pow4_fooof,          fooof.pow6_fooof};
    condPows_bl       = {fooof.pow2_fooof_bl,       fooof.pow4_fooof_bl,       fooof.pow6_fooof_bl};
    condPows_bl_early = {fooof.pow2_fooof_bl_early, fooof.pow4_fooof_bl_early, fooof.pow6_fooof_bl_early};
    condPows_bl_late  = {fooof.pow2_fooof_bl_late,  fooof.pow4_fooof_bl_late,  fooof.pow6_fooof_bl_late};
    condVals          = [2, 4, 6];

    for c = 1:3
        rowIdx = find(allIDs == subjID & allConds == condVals(c));
        if isempty(rowIdx), continue; end

        AlphaPower_FOOOF(rowIdx)          = mean(mean(condPows{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl(rowIdx)       = mean(mean(condPows_bl{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl_early(rowIdx) = mean(mean(condPows_bl_early{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
        AlphaPower_FOOOF_bl_late(rowIdx)  = mean(mean(condPows_bl_late{c}.powspctrm(occ_idx, alpha_idx), 2, 'omitnan'), 1, 'omitnan');
    end

    clc
    fprintf('FOOOF: Subject %s (%d/%d) loaded.\n', subjStr, s, numel(uIDs));
end

%% ===== (B) Baselined gaze metrics (trial → subject-level) =====
% Average trial-level baselined values per subject x condition.
% Baseline correction (already applied in trial-level extraction):
%   GazeDeviation, ScanPathLength, MSRate: dB  = 10*log10(metric/baseline)
%   PupilSize:                             %Δ  = 100*(metric-baseline)/baseline
% Baseline window: [-0.5 -0.25]s

gazeFields_bl = { ...
    'GazeDeviationFullBL',    'GazeDeviationEarlyBL',    'GazeDeviationLateBL', ...
    'ScanPathLengthFullBL',   'ScanPathLengthEarlyBL',   'ScanPathLengthLateBL', ...
    'PupilSizeFullBL',        'PupilSizeEarlyBL',        'PupilSizeLateBL', ...
    'MSRateFullBL',           'MSRateEarlyBL',           'MSRateLateBL'};

nGazeFields = numel(gazeFields_bl);
gazeBL = nan(nRows, nGazeFields);

% Build look-up arrays from trial data
trlIDs   = [gaze_data_sternberg_trials.ID];
trlConds = [gaze_data_sternberg_trials.Condition];

for s = 1:numel(uIDs)
    subjID = uIDs(s);
    condVals = [2, 4, 6];

    for c = 1:3
        % Row in the subject-level output
        rowIdx = find(allIDs == subjID & allConds == condVals(c));
        if isempty(rowIdx), continue; end

        % Matching trials
        trlMask = trlIDs == subjID & trlConds == condVals(c);
        if ~any(trlMask), continue; end

        trlSubset = gaze_data_sternberg_trials(trlMask);

        for f = 1:nGazeFields
            vals = [trlSubset.(gazeFields_bl{f})];
            gazeBL(rowIdx, f) = mean(vals, 'omitnan');
        end
    end
end

fprintf('Baselined gaze metrics aggregated from trial data.\n');

%% Build LONG merged struct
merged_data_sternberg_LONG = struct( ...
    ... % --- Original columns (21) ---
    'ID',                      {behav_data_sternberg.ID}, ...
    'Gender',                  {behav_data_sternberg.Gender}, ...
    'Age',                     {behav_data_sternberg.Alter}, ...
    'Handedness',              {behav_data_sternberg.H_ndigkeit}, ...
    'OcularDominance',         {behav_data_sternberg.OcularDominance}, ...
    'Condition',               {behav_data_sternberg.Condition}, ...
    'Accuracy',                {behav_data_sternberg.Accuracy}, ...
    'ReactionTime',            {behav_data_sternberg.ReactionTime}, ...
    'GazeDeviation',           {gaze_data_sternberg.GazeDeviation}, ...
    'GazeStdX',                {gaze_data_sternberg.GazeStdX}, ...
    'GazeStdY',                {gaze_data_sternberg.GazeStdY}, ...
    'PupilSize',               {gaze_data_sternberg.PupilSize}, ...
    'MSRate',                  {gaze_data_sternberg.MSRate}, ...
    'Blinks',                  {gaze_data_sternberg.Blinks}, ...
    'Fixations',               {gaze_data_sternberg.Fixations}, ...
    'Saccades',                {gaze_data_sternberg.Saccades}, ...
    'ScanPathLength',          {gaze_data_sternberg.ScanPathLength}, ...
    'AlphaPower',              {eeg_data_sternberg.AlphaPower}, ...
    'IAF',                     {eeg_data_sternberg.IAF}, ...
    'Lateralization',          {eeg_data_sternberg.Lateralization}, ...
    ... % --- FOOOF alpha (4) ---
    'AlphaPower_FOOOF',          num2cell(AlphaPower_FOOOF)', ...
    'AlphaPower_FOOOF_bl',       num2cell(AlphaPower_FOOOF_bl)', ...
    'AlphaPower_FOOOF_bl_early', num2cell(AlphaPower_FOOOF_bl_early)', ...
    'AlphaPower_FOOOF_bl_late',  num2cell(AlphaPower_FOOOF_bl_late)', ...
    ... % --- Baselined gaze (12) ---
    'GazeDeviationFullBL',     num2cell(gazeBL(:, 1))', ...
    'GazeDeviationEarlyBL',    num2cell(gazeBL(:, 2))', ...
    'GazeDeviationLateBL',     num2cell(gazeBL(:, 3))', ...
    'ScanPathLengthFullBL',    num2cell(gazeBL(:, 4))', ...
    'ScanPathLengthEarlyBL',   num2cell(gazeBL(:, 5))', ...
    'ScanPathLengthLateBL',    num2cell(gazeBL(:, 6))', ...
    'PupilSizeFullBL',         num2cell(gazeBL(:, 7))', ...
    'PupilSizeEarlyBL',        num2cell(gazeBL(:, 8))', ...
    'PupilSizeLateBL',         num2cell(gazeBL(:, 9))', ...
    'MSRateFullBL',            num2cell(gazeBL(:,10))', ...
    'MSRateEarlyBL',           num2cell(gazeBL(:,11))', ...
    'MSRateLateBL',            num2cell(gazeBL(:,12))');

%% Save as .mat
save(fullfile(featPath, 'merged_data_sternberg_LONG.mat'), 'merged_data_sternberg_LONG');
fprintf('Saved merged_data_sternberg_LONG.mat\n');

%% Save as .csv
merged_table = struct2table(merged_data_sternberg_LONG);
writetable(merged_table, fullfile(featPath, 'merged_data_sternberg_LONG.csv'));
fprintf('Saved merged_data_sternberg_LONG.csv\n');

%% Summary
nCols = numel(fieldnames(merged_data_sternberg_LONG));
fprintf('\n=== Sternberg LONG merge complete ===\n');
fprintf('Rows: %d  |  Columns: %d\n', numel(merged_data_sternberg_LONG), nCols);
fprintf('Original columns (21):  ID … Lateralization\n');
fprintf('FOOOF alpha     ( 4):  AlphaPower_FOOOF, _bl, _bl_early, _bl_late\n');
fprintf('Baselined gaze  (12):  GazeDeviation/ScanPathLength/PupilSize/MSRate × Full/Early/Late BL\n');
fprintf('    GD/SPL/MS baseline correction: dB = 10*log10(metric/baseline)\n');
fprintf('    PupilSize baseline correction: %%change = 100*(metric-baseline)/baseline\n');
fprintf('    Baseline window: [-0.5 -0.25]s\n');
fprintf('Total columns:  %d\n', nCols);
