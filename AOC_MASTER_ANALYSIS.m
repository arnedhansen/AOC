%% AOC Master Analysis Script
%
% Executes all analysis steps for the AOC study except cut and Automagic
%
% Executed analysis steps:
%
%   Preprocessing:
%       1_preprocessing/3_merge/AOC_mergeData.m
%       1_preprocessing/4_preprocessing/AOC_preprocessing_nback.m
%       1_preprocessing/4_preprocessing/AOC_preprocessing_sternberg.m
%
%   Feature Extraction:
%       2_feature_extraction/behavioral/AOC_behavioral_fex_nback.m
%       2_feature_extraction/behavioral/AOC_behavioral_fex_sternberg.m
%       2_feature_extraction/eeg/AOC_eeg_fex_nback.m
%       2_feature_extraction/eeg/AOC_eeg_fex_sternberg.m
%       2_feature_extraction/gaze/AOC_gaze_fex_nback.m
%       2_feature_extraction/gaze/AOC_gaze_fex_sternberg.m
%       2_feature_extraction/AOC_master_matrix_nback.m
%       2_feature_extraction/AOC_master_matrix_sternberg.m
%
%   Visualizations:   
%       3_visualization/behavioral/AOC_behav_dev_nback.m
%       3_visualization/behavioral/AOC_behav_dev_sternberg.m
%       3_visualization/eeg/alpha_over_trials/AOC_eeg_alpha_power_over_trials_nback.m
%       3_visualization/eeg/powspctrm/AOC_eeg_alpha_power_nback.m
%       3_visualization/eeg/powspctrm/AOC_eeg_alpha_power_sternberg.m
%       3_visualization/eeg/powspctrm/by_trial/AOC_eeg_alpha_power_by_trial_nback.m
%       3_visualization/eeg/powspctrm/by_trial/AOC_eeg_alpha_power_by_trial_sternberg.m
%       3_visualization/eeg/tfr/AOC_tfr_nback.m
%       3_visualization/eeg/tfr/AOC_tfr_sternberg.m
%       3_visualization/eeg/topos/AOC_eeg_topos_nback.m
%       3_visualization/eeg/topos/AOC_eeg_topos_sternberg.m
%       3_visualization/gaze/deviation/AOC_gaze_dev_nback.m
%       3_visualization/gaze/deviation/AOC_gaze_dev_sternberg.m
%       3_visualization/gaze/heatmap/AOC_gaze_heatmap_nback.m
%       3_visualization/gaze/heatmap/AOC_gaze_heatmap_sternberg.m
%       3_visualization/gaze/microsaccades/AOC_gaze_microsaccades_nback.m
%       3_visualization/gaze/microsaccades/AOC_gaze_microsaccades_sternberg.m
%
%   Stats:   
%       4_stats/orf/AOC_orf.m
%       4_stats/overview/AOC_stats_overview_nback.m
%       4_stats/overview/AOC_stats_overview_sternberg.m
%
%   Data Checks:   
%       controls/AOC_cal_val_ET.m
%       controls/AOC_fixCheck_exclusion_trials.m
%       controls/AOC_missing_data.m
%       controls/AOC_paradigm_durations.m
%       controls/AOC_recording_order.m

%% Setup
startup
clc

%% Run all
% Set base path depending on platform
if ispc
    basePath = 'C:\Users\dummy\Documents\GitHub\AOC';
else
    basePath = '/Users/Arne/Documents/GitHub/AOC';
end

% Get all subfolders
allDirs = genpath(basePath);
dirList = strsplit(allDirs, pathsep);

% Loop through each subdirectory
errorLog = {};
for i = 1:length(dirList)
    currentDir = dirList{i};

    % Skip empty or excluded directories
    [~, folderName] = fileparts(currentDir);
    if isempty(currentDir) || contains(currentDir, '_OSF') || contains(currentDir, '_BACKUP') || contains(currentDir, 'paradigms') || contains(currentDir, '1_cut') || contains(currentDir, 'MASTER')
        continue;
    end

    % Get all .m files in this folder
    mFiles = dir(fullfile(currentDir, '*.m'));

    % Sort files alphabetically
    mFileNames = sort({mFiles.name});
    if strcmp(mFileNames, 'AOC_MASTER_ANALYSIS.m')
        continue;
    end

    for j = 1:length(mFileNames)
        scriptPath = fullfile(currentDir, mFileNames{j});
        [~, scriptName, ~] = fileparts(scriptPath);

        % Display which script is running
        fprintf('Running script: %s\n', scriptPath);

        % Change to script directory
        oldPath = pwd;
        cd(currentDir);

        try
            runInFunction(scriptName);  % Execute the script
        catch ME
            warning('Error while running %s: %s', scriptName, ME.message);
            errorLog{end+1, 1} = scriptPath;
            errorLog{end, 2} = ME.message;
        end
        cd(oldPath);  % Return to previous directory
    end
end
if ispc
    run('C:\Users\dummy\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');
    run('C:\Users\dummy\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');
else
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/AOC_master_matrix_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/AOC_master_matrix_sternberg.m');
end
disp(upper('All AOC analysis scripts executed!'));

%% Display error summary
if ~isempty(errorLog)
    fprintf('\n===== ERROR SUMMARY =====\n');
    for i = 1:size(errorLog, 1)
        fprintf('Script: %s\nError: %s\n\n', errorLog{i,1}, errorLog{i,2});
    end
else
    disp('No errors occurred during execution.');
end
