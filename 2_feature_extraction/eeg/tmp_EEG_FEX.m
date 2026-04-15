%% AOC EEG overnight tmp runner (startup-safe)
% IMPORTANT:
% - Scripts called below run `startup` internally and may clear workspace.
% - Therefore this file intentionally avoids loop/index state between calls.
% - Execution order follows the split pipeline dependencies.

startup

disp('=== AOC overnight pipeline start ===')
disp(datestr(now, 31))
repoRoot = 'C:\Users\Administrator\Documents\GitHub\AOC';

% 1) Non-FOOOF EEG (creates IAF CSV used by TFR FOOOF scripts)
disp('--- RUN: AOC_eeg_fex_sternberg.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'eeg', 'AOC_eeg_fex_sternberg.m'));

disp('--- RUN: AOC_eeg_fex_nback.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'eeg', 'AOC_eeg_fex_nback.m'));

% 2) TFR FOOOF EEG (reads IAF from non-FOOOF CSV)
disp('--- RUN: AOC_eeg_fex_sternberg_TFR.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'eeg', 'AOC_eeg_fex_sternberg_TFR.m'));

disp('--- RUN: AOC_eeg_fex_nback_TFR.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'eeg', 'AOC_eeg_fex_nback_TFR.m'));

% 3) Master matrices
disp('--- RUN: AOC_master_matrix_sternberg.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'AOC_master_matrix_sternberg.m'));

disp('--- RUN: AOC_master_matrix_nback.m ---')
run(fullfile(repoRoot, '2_feature_extraction', 'AOC_master_matrix_nback.m'));

disp('=== AOC overnight pipeline end ===')
disp(datestr(now, 31))