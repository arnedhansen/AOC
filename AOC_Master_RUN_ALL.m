%% AOC Master Analysis
%
% Runs the AOC pipeline. Running preprocessing, feature extraction and
% creation of CSVs for stats. Change the paths to your file locations first.
%
% IMPORTANT: Scripts call startup/clear, which wipe the workspace.

startup
disp(datestr(now))

% ----------------------------------------------------------------------------
% For raw data, DO THIS FIRST: in 1_preprocessing, run in order
%   (1) scripts in 1_cut to cut the raw data,
%   (2) Automagic (refer to paper for settings),
%   (3) scripts in 3_merge to get merged blocks data
% ----------------------------------------------------------------------------

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
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg_TFR.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback_TFR.m');

% Master Matrices (CSVs)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');

%% 3 Visualization
% Powspctrm
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_sternberg.m');

% TFR
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback_fooof_abs.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg_fooof_abs.m');

% Topos
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_sternberg.m');

%% 4 Split
run('C:\Users\Administrator\Documents\GitHub\AOC\AOC_split_AlphaAmpRed.m');

%% Stats
% For the stats and raincloud plots, run the Python scripts.

%%
disp(datestr(now))