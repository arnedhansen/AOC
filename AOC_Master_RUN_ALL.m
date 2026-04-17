%% AOC — master runner (startup-safe, subject-level only)
%
% This orchestration runs the MATLAB AOC pipeline in dependency order. Running
% it through the feature-extraction and master-matrix stages creates the CSV
% files required as inputs for the statistical analyses implemented in the
% Python scripts in this repository (alongside corresponding MAT outputs on the
% data share).
%
% IMPORTANT: Scripts call `startup` and/or `clear`, which wipe the workspace.
%

startup

% ----------------------------------------------------------------------------
% Manual prerequisite (not executed here): in 1_preprocessing, run in order
%   (1) scripts in 1_cut to cut the raw data,
%   (2) Automagic,
%   (3) scripts in 3_merge to build merged blocks,
%   then the AOC_preprocessing_* scripts below (4_preprocessing) apply.
% ----------------------------------------------------------------------------

disp('=== AOC Master RUN ALL (subject-level) start ===')
disp(datestr(now, 31))

%% 1_preprocessing/4_preprocessing — FieldTrip (merged blocks → epoched EEG)
disp('[RUN: AOC_preprocessing_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\1_preprocessing\4_preprocessing\AOC_preprocessing_sternberg.m');

disp('[RUN: AOC_preprocessing_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\1_preprocessing\4_preprocessing\AOC_preprocessing_nback.m');

%% 2 — Subject-level feature extraction
disp('[RUN: AOC_demographics.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_demographics.m');

disp('[RUN: AOC_behavioral_fex_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\behavioral\AOC_behavioral_fex_sternberg.m');

disp('[RUN: AOC_behavioral_fex_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\behavioral\AOC_behavioral_fex_nback.m');

disp('[RUN: AOC_gaze_fex_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_sternberg.m');

disp('[RUN: AOC_gaze_fex_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_nback.m');

disp('[RUN: AOC_eeg_fex_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg.m');

disp('[RUN: AOC_eeg_fex_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback.m');

disp('[RUN: AOC_eeg_fex_sternberg_TFR.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg_TFR.m');

disp('[RUN: AOC_eeg_fex_nback_TFR.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback_TFR.m');

disp('[RUN: AOC_master_matrix_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');

disp('[RUN: AOC_master_matrix_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');

%% 3 — Visualization (MATLAB figures)
disp('[RUN: AOC_eeg_powspctrm_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_nback.m');

disp('[RUN: AOC_eeg_powspctrm_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_sternberg.m');

disp('[RUN: AOC_tfr_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback.m');

disp('[RUN: AOC_tfr_nback_fooof_abs.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback_fooof_abs.m');

disp('[RUN: AOC_tfr_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg.m');

disp('[RUN: AOC_tfr_sternberg_fooof_abs.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg_fooof_abs.m');

disp('[RUN: AOC_eeg_topos_nback.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_nback.m');

disp('[RUN: AOC_eeg_topos_sternberg.m]')
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\topos\AOC_eeg_topos_sternberg.m');

%% 4 — SPLITS


%%
disp('=== AOC Master RUN ALL (subject-level) end ===')
disp(datestr(now, 31))
