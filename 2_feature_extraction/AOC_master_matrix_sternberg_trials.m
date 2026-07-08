%% AOC Master Matrix Sternberg Trial Level
% Merge behavioral, gaze, and EEG trial matrices on ID Trial Condition.
% EEG contract is ERSD only for trial level outputs.
startup
[~, paths, ~, ~] = setup('AOC', 0);
featPath = paths.features;

demog = readtable(paths.vp_table);
demog = demog(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'OcularDominance'});
demog = table2struct(demog(1:120, :));

load(fullfile(featPath, 'AOC_behavioral_matrix_sternberg_trials.mat'));   % behav_data_sternberg_trials
load(fullfile(featPath, 'AOC_gaze_matrix_sternberg_trials.mat'));         % gaze_data_sternberg_trials
load(fullfile(featPath, 'AOC_eeg_matrix_sternberg_trials.mat'));          % eeg_data_sternberg_trials

demoIDs = [demog.ID];
for i = 1:numel(behav_data_sternberg_trials)
    idx = find(demoIDs == behav_data_sternberg_trials(i).ID, 1);
    behav_data_sternberg_trials(i).Gender = demog(idx).Gender;
    behav_data_sternberg_trials(i).Alter = demog(idx).Alter;
    behav_data_sternberg_trials(i).H_ndigkeit = demog(idx).H_ndigkeit;
    behav_data_sternberg_trials(i).OcularDominance = demog(idx).OcularDominance;
end

T_behav = struct2table(behav_data_sternberg_trials);
T_gaze = struct2table(gaze_data_sternberg_trials);
T_eeg = struct2table(eeg_data_sternberg_trials);

assert_unique_keys(T_behav, {'ID', 'Trial', 'Condition'}, 'behav_data_sternberg_trials');
assert_unique_keys(T_gaze, {'ID', 'Trial', 'Condition'}, 'gaze_data_sternberg_trials');
assert_unique_keys(T_eeg, {'ID', 'Trial', 'Condition'}, 'eeg_data_sternberg_trials');

assert_condition_set(T_behav, [2 4 6], 'behav_data_sternberg_trials');
assert_condition_set(T_gaze, [2 4 6], 'gaze_data_sternberg_trials');
assert_condition_set(T_eeg, [2 4 6], 'eeg_data_sternberg_trials');

T = innerjoin(T_behav, T_eeg, 'Keys', {'ID', 'Trial', 'Condition'});
T = innerjoin(T, T_gaze, 'Keys', {'ID', 'Trial', 'Condition'});

assert_has_vars(T, {'ERSD_early', 'ERSD_late', 'ERSD_full'});

T.Properties.VariableNames{'Alter'} = 'Age';
T.Properties.VariableNames{'H_ndigkeit'} = 'Handedness';

orderedVars = { ...
    'ID', 'Trial', 'Condition', ...
    'Gender', 'Age', 'Handedness', 'OcularDominance', ...
    'Accuracy', 'ReactionTime', 'Stimuli', 'Probe', 'Match', ...
    'ERSD_early', 'ERSD_late', 'ERSD_full', ...
    'GazeDeviationEarly', 'GazeDeviationEarlyBL', ...
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
    'BCEALatFull', 'BCEALatFullBL'};

orderedVars = orderedVars(ismember(orderedVars, T.Properties.VariableNames));
merged_data_sternberg_trials = T(:, orderedVars);

save(fullfile(featPath, 'AOC_merged_data_sternberg_trials.mat'), 'merged_data_sternberg_trials');
writetable(merged_data_sternberg_trials, fullfile(featPath, 'AOC_merged_data_sternberg_trials.csv'));

function assert_unique_keys(T, keyVars, tableName)
[~, ia] = unique(T(:, keyVars), 'rows', 'stable');
if numel(ia) ~= height(T)
    error('Duplicate key rows detected in %s for keys: %s', tableName, strjoin(keyVars, ', '));
end
end

function assert_has_vars(T, varsRequired)
missingVars = varsRequired(~ismember(varsRequired, T.Properties.VariableNames));
if ~isempty(missingVars)
    error('Missing required variable(s): %s', strjoin(missingVars, ', '));
end
end

function assert_condition_set(T, allowedSet, tableName)
vals = unique(T.Condition(:))';
if ~all(ismember(vals, allowedSet))
    error('Unexpected condition value(s) in %s: %s', tableName, mat2str(vals));
end
end
