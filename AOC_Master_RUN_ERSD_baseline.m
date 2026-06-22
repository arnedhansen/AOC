%% AOC Master — ERSD baseline rerun (trial-avg then dB baseline)
%
% Minimal rerun after switching TFR baseline to trial-average-then-baseline.
% Rebuilds TFR products, ERSD scalars, merged CSVs, ERSD figures, and splits.
%
% IMPORTANT: Scripts call startup/clear, which wipe the workspace. Use
% absolute paths for the 'run' commands.
%
% After this script, rerun stats for updated rainclouds:
%   Rscript 4_stats/AOC_stats_glmm_nback.R
%   Rscript 4_stats/AOC_stats_glmm_sternberg.R
%   Rscript 4_stats/AOC_stats_glmm_export.R
%   python 4_stats/AOC_stats_glmms_rainclouds.py

%% Setup
startup

%% 2 Feature extraction — EEG TFR + ERSD (Sternberg + N-back)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\eeg\AOC_eeg_fex_nback.m');

%% 2 Master matrices (merged CSVs with updated ERSD_*)
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\2_feature_extraction\AOC_master_matrix_nback.m');

%% 3 ERSD figures (time course + topos)
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\ersd\AOC_eeg_ersd_sternberg.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\3_visualization\eeg\ersd\AOC_eeg_ersd_nback.m');

%% 4 Split analyses (ERSD thresholds, time courses, topos)
run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed.m');
run('C:\Users\Administrator\Documents\GitHub\AOC\splits\AOC_split_AlphaAmpRed_Nback.m');

%% Control — baselined retention spectra sanity check
run('C:\Users\Administrator\Documents\GitHub\AOC\_controls\AOC_eeg_powspctrm_baseline_effects.m');

%%
disp('AOC ERSD baseline rerun finished.')
disp(datestr(now))
