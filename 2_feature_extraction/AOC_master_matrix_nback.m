%% AOC Master Matrix — N-Back
% Merges three sources by `ID`,`Condition`:
%   (a) behavioral
%   (b) gaze
%   (c) EEG (`AOC_eeg_matrix_nback.mat`)
% Includes key-uniqueness checks before merge to prevent row inflation.
%
% Key outputs:
%   merged_data_nback.mat
%   merged_data_nback.csv
%
% Columns include:
%   EEG — Raw alpha power (IAF band from multitaper IAF, occ channels):
%     AlphaPower_raw_early / _late / _full
%     AlphaPower_bl_early / _late / _full   (dB baseline [-1.5 -0.5]s)
%     ERSD_early / ERSD_late / ERSD_full — (IAF-4, IAF+2) Hz dB on baselined tfr*_bl, O/I occ channels, windows [0 1] / [1 2] / [0 2]s (NaN if no valid IAF)
%     IAF
%     AlphaPower_baselineWindow — baselineWindow FFT [-1.5 -0.5] s, raw IAF-band power [μV²/Hz]
%
%   Gaze — baselined (BL window [-1.5 -0.5]s; % change for all gaze metrics):
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

% EEG
load(fullfile(featPath, 'AOC_eeg_matrix_nback.mat'));         % eeg_data_nback

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
eeg_table   = join_ersd_table(eeg_table, featPath, 'AOC_eeg_matrix_nback_ERSD.mat', 'eeg_data_nback_ERSD');
eeg_table   = join_baselineWindow_table(eeg_table, featPath, ...
    'AOC_eeg_matrix_nback_baselineWindow.mat', 'eeg_data_nback_baselineWindow');

assert_unique_keys(eeg_table, {'ID', 'Condition'}, 'eeg_data_nback');

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
disp(merged_data_nback)

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

function T = join_baselineWindow_table(T, featPath, blMatName, blVarName)
blMatPath = fullfile(featPath, blMatName);
if ~isfile(blMatPath)
    return
end
S = load(blMatPath, blVarName);
bl_data = S.(blVarName);
if isempty(bl_data)
    return
end
Tb = struct2table(bl_data);
if height(Tb) == 0
    return
end
assert_unique_keys(Tb, {'ID', 'Condition'}, blVarName);
% Keep AlphaPower_baselineWindow only; drop IAF_baselineWindow from merge.
if ismember('IAF_baselineWindow', Tb.Properties.VariableNames)
    Tb.IAF_baselineWindow = [];
end
for cn = {'AlphaPower_baselineWindow'}
    if ismember(cn{1}, T.Properties.VariableNames)
        T.(cn{1}) = [];
    end
end
if ismember('IAF_baselineWindow', T.Properties.VariableNames)
    T.IAF_baselineWindow = [];
end
T = outerjoin(T, Tb, 'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');
end