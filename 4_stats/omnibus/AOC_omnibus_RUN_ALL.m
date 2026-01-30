%% AOC Omnibus — Run All (plain + FOOOF)
% Runs in sequence:
%   1. AOC_omnibus_prep.m       → omnibus_data.mat (plain: raw TFR, dB baseline)
%   2. AOC_omnibus_prep_FOOOF.m → omnibus_data_FOOOF.mat, omnibus_data_NOBASELINE_FOOOF.mat
%   3. AOC_omnibus.m             → plain 2x2 figures + stat file
%   4. AOC_omnibus_FOOOF.m       → FOOOF 2x2 figures + stat file
%
% Run from project root or from 4_stats/omnibus. Ensures project path and runs all four scripts.

%% Setup path and startup
omnibus_dir = fullfile('/Users/Arne/Documents/GitHub/AOC/4_stats/omnibus');

%% 1. Plain omnibus data prep (raw TFR, dB baseline, 3–30 Hz)
disp(' ')
disp('========== 1/4 AOC_omnibus_prep.m (plain) ==========')
tic
run(fullfile(omnibus_dir, 'AOC_omnibus_prep.m'))
fprintf('Done in %.1f min.\n', toc/60)

%% 2. FOOOF omnibus data prep
disp(' ')
disp('========== 2/4 AOC_omnibus_prep_FOOOF.m ==========')
tic
run(fullfile(omnibus_dir, 'AOC_omnibus_prep_FOOOF.m'))
fprintf('Done in %.1f min.\n', toc/60)

%% 3. Plain omnibus figures and stats
disp(' ')
disp('========== 3/4 AOC_omnibus.m (plain figures) ==========')
tic
run(fullfile(omnibus_dir, 'AOC_omnibus.m'))
fprintf('Done in %.1f min.\n', toc/60)

%% 4. FOOOF omnibus figures and stats
disp(' ')
disp('========== 4/4 AOC_omnibus_FOOOF.m ==========')
tic
run(fullfile(omnibus_dir, 'AOC_omnibus_FOOOF.m'))
fprintf('Done in %.1f min.\n', toc/60)

%% Done
disp(' ')
disp('========== AOC omnibus RUN ALL finished ==========')
disp('Plain:  omnibus_data.mat + 2x2 figures + AOC_omnibus_FIGURE2x2_statFnb_statFsb.mat')
disp('FOOOF:  omnibus_data_FOOOF.mat (+ NOBASELINE) + 2x2 figures_FOOOF + AOC_omnibus_statFnb_statFsb_FOOOF.mat')
