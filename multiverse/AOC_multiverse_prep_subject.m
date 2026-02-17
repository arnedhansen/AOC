%% AOC Multiverse — Subject-Level Data Preparation
% Reads the trial-level multiverse CSVs and aggregates to subject means.
% Writes multiverse_sternberg_subject.csv and multiverse_nback_subject.csv.
%
% Much smaller output (~160k rows vs ~31M) for faster R processing.
% Subject-level FOOOF alpha = mean of per-trial FOOOF alpha values.
%
% Run AFTER AOC_multiverse_prep.m has generated the trial-level CSVs.

disp(upper('=== AOC MULTIVERSE SUBJECT-LEVEL PREP START ==='))

%% Paths
if ispc
    base_features = 'W:\Students\Arne\AOC\data\features';
else
    base_features = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features';
end
disp(upper(['Data directory: ' base_features]))

%% Process each task
tasks = {'sternberg', 'nback'};

for t = 1:length(tasks)
    task = tasks{t};
    csv_in  = fullfile(base_features, ['multiverse_' task '.csv']);
    csv_out = fullfile(base_features, ['multiverse_' task '_subject.csv']);

    if ~isfile(csv_in)
        disp(upper(['SKIPPING ' task ': trial-level CSV not found: ' csv_in]))
        continue
    end

    disp(upper(['READING TRIAL-LEVEL CSV: ' task]))
    tic
    T = readtable(csv_in, 'TextType', 'string');
    disp(upper(sprintf('  Loaded %d rows in %.1f seconds.', height(T), toc)))

    %% Grouping variables (everything except Trial, alpha, gaze_value)
    groupVars = {'task', 'universe_id', 'electrodes', 'fooof', 'latency_ms', ...
                 'alpha_type', 'gaze_measure', 'baseline_eeg', 'baseline_gaze', ...
                 'subjectID', 'Condition'};

    disp(upper('  AGGREGATING TO SUBJECT MEANS...'))
    tic
    T_subj = groupsummary(T, groupVars, 'mean', {'alpha', 'gaze_value'});
    disp(upper(sprintf('  Aggregated to %d rows in %.1f seconds.', height(T_subj), toc)))

    %% Clean up column names
    T_subj.Properties.VariableNames{'mean_alpha'}      = 'alpha';
    T_subj.Properties.VariableNames{'mean_gaze_value'} = 'gaze_value';
    T_subj.GroupCount = [];

    %% Write output
    disp(upper(['  WRITING: ' csv_out]))
    tic
    writetable(T_subj, csv_out);
    disp(upper(sprintf('  Written %d rows in %.1f seconds.', height(T_subj), toc)))

    %% Summary
    n_subjects  = length(unique(T_subj.subjectID));
    n_universes = length(unique(T_subj.universe_id));
    disp(upper(sprintf('  %s: %d subjects, %d universes, %d rows.', ...
        task, n_subjects, n_universes, height(T_subj))))
end

disp(upper('=== AOC MULTIVERSE SUBJECT-LEVEL PREP DONE ==='))
