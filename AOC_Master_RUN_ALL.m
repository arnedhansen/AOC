%% AOC Master Analysis
%
% Runs the AOC pipeline. Running preprocessing, feature extraction, visualizations,
% and creation of CSVs for stats. Change the paths to your file locations first.
%
% IMPORTANT: Scripts call startup/clear, which wipe the workspace. Use
% absolute paths for the 'run' commands.
%
% For raw data, DO THIS FIRST: in 1_preprocessing, run in order
%   (1) scripts in 1_cut to cut the raw data,
%   (2) Automagic (refer to paper for settings),
%   (3) scripts in 3_merge to get merged blocks data

%% Setup
startup

if ispc
    %% 1_preprocessing/4_preprocessing FieldTrip
    run('C:\Users\Administrator\Documents\GitHub\AOC\1_preprocessing\4_preprocessing\AOC_preprocessing_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\1_preprocessing\4_preprocessing\AOC_preprocessing_nback.m');

    %% 2 Subject-level feature extraction
    % Demographics
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_demographics.m');

    % Behavioral
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\behavioral\AOC_behavioral_fex_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\behavioral\AOC_behavioral_fex_nback.m');

    % Gaze
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_nback.m');

    % EEG
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg_FOOOF.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback_FOOOF.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\_controls\AOC_eeg_powspctrm_baseline_effects.m');

    % Master Matrices (CSVs)
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');

    %% 3 Visualization
    % Powspctrm (raw, baselined dB, FOOOF+baseline per task; baseline-window summary last)
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_nback_bl.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_nback_fooof_bl.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_sternberg_bl.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_sternberg_fooof_bl.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_baseline_window_fooof.m');

    % TFR
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback_fooof_abs.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg_fooof_abs.m');

    % Topos
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_sternberg.m');

    % Gaze
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_sternberg_split.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_nback_split.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\deviation\AOC_gaze_dev_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\deviation\AOC_gaze_dev_nback.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\microsaccades\AOC_gaze_microsaccades_sternberg.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\microsaccades\AOC_gaze_microsaccades_nback.m');

    %% 4 Split analyses
    run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed_GazeDev.m');
    run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed_MS.m');

    %% Stats
    % For the stats and raincloud plots, run the Python scripts.

else
    %% 1_preprocessing/4_preprocessing FieldTrip
    run('/Users/Arne/Documents/GitHub/AOC/1_preprocessing/4_preprocessing/AOC_preprocessing_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/1_preprocessing/4_preprocessing/AOC_preprocessing_nback.m');

    %% 2 Subject-level feature extraction
    % Demographics
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/AOC_demographics.m');

    % Behavioral
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/behavioral/AOC_behavioral_fex_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/behavioral/AOC_behavioral_fex_nback.m');

    % Gaze
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/gaze/AOC_gaze_fex_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/gaze/AOC_gaze_fex_nback.m');

    % EEG
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/eeg/AOC_eeg_fex_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/eeg/AOC_eeg_fex_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/eeg/AOC_eeg_fex_sternberg_FOOOF.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/eeg/AOC_eeg_fex_nback_FOOOF.m');
    run('/Users/Arne/Documents/GitHub/AOC/_controls/AOC_eeg_powspctrm_baseline_effects.m');

    % Master Matrices (CSVs)
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/AOC_master_matrix_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/2_feature_extraction/AOC_master_matrix_nback.m');

    %% 3 Visualization
    % Powspctrm (raw, baselined dB, FOOOF+baseline per task; baseline-window summary last)
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_nback_bl.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_nback_fooof_bl.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_sternberg_bl.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_sternberg_fooof_bl.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/powspctrm/AOC_eeg_powspctrm_baseline_window_fooof.m');

    % TFR
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/tfr/AOC_tfr_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/tfr/AOC_tfr_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/tfr/AOC_tfr_nback_fooof_abs.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/tfr/AOC_tfr_sternberg_fooof_abs.m');

    % Topos
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/topos/AOC_eeg_topos_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/eeg/topos/AOC_eeg_topos_sternberg.m');

    % Gaze
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/heatmap/AOC_gaze_heatmap_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/heatmap/AOC_gaze_heatmap_sternberg_split.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/heatmap/AOC_gaze_heatmap_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/heatmap/AOC_gaze_heatmap_nback_split.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/deviation/AOC_gaze_dev_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/deviation/AOC_gaze_dev_nback.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/microsaccades/AOC_gaze_microsaccades_sternberg.m');
    run('/Users/Arne/Documents/GitHub/AOC/3_visualization/gaze/microsaccades/AOC_gaze_microsaccades_nback.m');

    %% 4 Split analyses
    run('/Users/Arne/Documents/GitHub/AOC/splits/AOC_split_AlphaAmpRed_GazeDev.m');
    run('/Users/Arne/Documents/GitHub/AOC/splits/AOC_split_AlphaAmpRed_MS.m');

    %% Stats
    % For the stats and raincloud plots, run the Python scripts.
end

%%
disp(datestr(now))
