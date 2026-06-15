%% AOC Master — rerun after EEG / gaze baseline window change
%
% Recomputes features and figures that depend on the baseline interval
% [-1.5 -0.5] s. Preprocessing is omitted (epochs already span that range).
%
% FOOOF: pipeline scripts use RUN_FOOOF = false (TFR, main EEG fex IAF_specParam,
% _controls/AOC_eeg_powspctrm_baseline_effects). Set those to true to re-enable.
% FOOOF-only run lines below are commented out so the driver matches that default.
%
% Change the paths if the repo is not under this location.
%
% IMPORTANT: Scripts call startup / clear, which wipe the workspace. Use
% absolute paths for the run commands.

%% Setup
startup

%% 2_feature_extraction — EEG (windowed power + TFR with baselines)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg_TFR.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback_TFR.m');

%% 2_feature_extraction — Gaze (percent-change baselines)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\AOC_gaze_fex_nback.m');

%% 2_feature_extraction — Master matrices (CSVs for stats)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');

%% 2_feature_extraction — Trial-level pipelines (optional; uncomment if needed)
% run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\trials\AOC_eeg_fex_sternberg_trials.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\trials\AOC_eeg_fex_nback_trials.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\trials\AOC_eeg_fex_sternberg_FOOOF_trials.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\trials\AOC_gaze_fex_sternberg_trials.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\gaze\trials\AOC_gaze_fex_nback_trials.m');

%% 3_visualization — EEG: baseline-window powspctrm, TFR (db / FOOOF abs), ERP
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_baseline_window.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\powspctrm\AOC_eeg_powspctrm_baseline_window_fooof.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg_bl.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback_bl.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_sternberg_fooof_abs.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\tfr\AOC_tfr_nback_fooof_abs.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\erp\AOC_eeg_erp_sternberg.m');

%% 3_visualization — Gaze
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_sternberg_split.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\heatmap\AOC_gaze_heatmap_nback_split.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\deviation\AOC_gaze_dev_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\deviation\AOC_gaze_dev_nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\velocity\AOC_gaze_velocity.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\gaze\microsaccades\AOC_gaze_microsaccades.m');

%% _controls — baseline vs retention powspctrm check
run('C:\Users\Administrator\Documents\GitHub\AOC\_controls\AOC_eeg_powspctrm_baseline_effects.m');

%% 4_splits (AlphaAmpRed* and AlphaLoads_TimeCourse need FOOOF-derived merged columns when RUN_FOOOF=false)
% run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed_Nback.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaLoads.m');
% run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaLoads_TimeCourse.m');

%% 4_stats
% After merged CSVs are updated, run from the repo or your stats environment, e.g.:
%   python 4_stats\AOC_stats_glmms_rainclouds.py
%   Rscript 4_stats\AOC_stats_fullModelTables.R

%%
disp(datestr(now))
