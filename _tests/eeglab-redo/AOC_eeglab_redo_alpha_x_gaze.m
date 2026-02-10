%% AOC EEGLAB Independent Reanalysis — Posterior Alpha Power × Eye Movements
% =========================================================================
% Completely independent pipeline. Uses ONLY EEGLAB + MATLAB functions.
% Zero FieldTrip calls. Starts from the block-level merged EEG-ET data.
%
% PURPOSE
%   Explore whether posterior alpha power is linked to any eye-movement
%   measure across two working-memory tasks (Sternberg, N-back).
%
% ANALYSES
%   1. Trial-level within-subject Spearman correlations (alpha vs 10 gaze
%      metrics), Fisher-z → group t-test
%   2. Median split on alpha power → gaze metric comparison (paired t)
%   3. GLMMs: GazeMetric ~ AlphaZ * Condition + (1|ID)
%   4. Peri-microsaccadic alpha envelope (±1 s, cluster permutation)
%   5. Alpha envelope × eye velocity cross-correlation
%   6. Alpha–velocity magnitude-squared coherence (mscohere)
%   7. Pre-stimulus alpha → post-stimulus gaze prediction
%   8. Subject-level condition covariation (Δalpha vs Δgaze)
%   9. Time-resolved sliding-window alpha–velocity correlation
%
% GAZE METRICS (per trial, analysis window)
%   GazeDeviation, GazeStdX, GazeStdY, PupilSize, MSRate,
%   MeanVelocity, PeakVelocity, ScanPathLength
%
% OUTPUTS
%   Figures → figDir (all PNG)
%   Data   → dataDir (MAT)
%
% USAGE
%   Run the whole script on the science cloud (ispc == 1) or locally.
%   Requires: EEGLAB on the path, Signal Processing Toolbox, Statistics
%   and Machine Learning Toolbox.
% =========================================================================

%% ====================================================================
%% 0  CONFIGURATION
%% ====================================================================
clearvars; close all; clc;

% --- EEGLAB ---
if ispc
    addpath('W:\4marius_bdf\eeglab');
else
    addpath('/Volumes/g_psyplafor_methlab$/4marius_bdf/eeglab');
end
eeglab; close all; clc;

% --- Paths ---
if ispc
    mergedPath = 'W:\Students\Arne\AOC\data\merged';
    figDir     = 'W:\Students\Arne\AOC\figures\tests\eeglab-redo';
    dataDir    = 'W:\Students\Arne\AOC\data\features\eeglab-redo';
else
    mergedPath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/merged';
    figDir     = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/eeglab-redo';
    dataDir    = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeglab-redo';
end
if ~exist(figDir,  'dir'), mkdir(figDir);  end
if ~exist(dataDir, 'dir'), mkdir(dataDir); end

% --- Subject list ---
dirs    = dir(mergedPath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

% Hard exclusions (no data / extreme)
excludeList = {'333', '358', '361', '405'};
subjects = subjects(~ismember(subjects, excludeList));
nSubj = numel(subjects);
fprintf('Subjects to process: %d\n', nSubj);

% --- Task definitions ---
taskDefs = struct( ...
    'name',       {'sternberg',           'nback'}, ...
    'blockPat',   {'_EEG_ET_Sternberg_block%d_merged.mat', ...
                   '_EEG_ET_Nback_block%d_merged.mat'}, ...
    'triggers',   {{'22','24','26'},      {'21','22','23'}}, ...
    'epochWin',   {[-2 3.5],              [-1.5 2.5]}, ...
    'anaWin',     {[1 2],                 [0 2]}, ...
    'preStimWin', {[-1 0],                [-1 0]}, ...
    'condOffset', {20,                    20}, ...
    'condLabFun', {@(c) sprintf('Load %d', c), @(c) sprintf('%d-back', c)} ...
);
nTasks = numel(taskDefs);

% --- General parameters ---
fs = 500;               % sampling rate (Hz)
alphaRange = [8 14];     % fallback alpha band
nfft_welch = 512;        % FFT length for Welch
screenCenter = [400, 300]; % screen center in pixels
screenSize   = [800, 600]; % screen resolution

% Peri-microsaccadic envelope
periWin_s    = 1.0;                               % ±1000 ms
periSamp     = round(periWin_s * fs);             % 500 samples
nPeriSamp    = 2 * periSamp + 1;                  % 1001
periTime_ms  = linspace(-periWin_s, periWin_s, nPeriSamp) * 1000;

% Cross-correlation lags
maxLag_s   = 0.5;  % ±500 ms
maxLagSamp = round(maxLag_s * fs);

% Coherence settings
coh_nfft   = 256;
coh_novlap = 128;

% Cluster permutation
nPerm    = 2000;
alphaCBP = 0.05;

% Microsaccade detection (Engbert & Kliegl)
ms_lambda     = 6;
ms_minDurSec  = 0.006;
ms_minISISec  = 0.02;

% Sliding window for time-resolved analysis
slidWin_s  = 0.25;   % 250 ms window
slidStep_s = 0.05;   % 50 ms step

% Figure aesthetics
colors3 = [0.576 0.722 0.769;   % condition 1 (light blue)
           0.510 0.678 0.510;   % condition 2 (green)
           0.851 0.600 0.635];  % condition 3 (pink)
cmap_corr = [0.200 0.400 0.600]; % main analysis color
fontSize  = 22;
set(0, 'DefaultFigureColor', 'w');

% Gaze metric names (for tables & figures)
gazeMetricNames = {'GazeDeviation', 'GazeStdX', 'GazeStdY', 'PupilSize', ...
                   'MSRate', 'MeanVelocity', 'PeakVelocity', 'ScanPathLength'};
nGazeMetrics = numel(gazeMetricNames);

fprintf('Configuration complete.\n');

%% ====================================================================
%% 1  FEATURE EXTRACTION (subject × task loop)
%% ====================================================================
% Output: one table per task with trial-level features, plus per-subject
% time-series aggregates for peri-MS, cross-corr, coherence, time-resolved.

for ti = 1:nTasks
    taskName   = taskDefs(ti).name;
    blockPat   = taskDefs(ti).blockPat;
    triggers   = taskDefs(ti).triggers;
    epochWin   = taskDefs(ti).epochWin;
    anaWin     = taskDefs(ti).anaWin;
    preStimWin = taskDefs(ti).preStimWin;
    condOffset = taskDefs(ti).condOffset;
    nEpochSamp = round(diff(epochWin) * fs) + 1;

    fprintf('\n%s\n  FEATURE EXTRACTION: %s\n%s\n', ...
        repmat('=',1,60), upper(taskName), repmat('=',1,60));

    % --- Group-level containers ---
    % Trial table columns: SubjID, Condition, AlphaPower, PreStimAlpha, + gaze
    allTrialData = [];

    % Peri-microsaccadic alpha  [nSubj × nPeriSamp]
    periMS_subj   = nan(nSubj, nPeriSamp);
    periSurr_subj = nan(nSubj, nPeriSamp);
    nMS_subj      = zeros(nSubj, 1);

    % Cross-correlation  [nSubj × (2*maxLagSamp+1)]
    nLags  = 2 * maxLagSamp + 1;
    xcorr_subj = nan(nSubj, nLags);
    xcorr_lags = (-maxLagSamp:maxLagSamp) / fs * 1000; % ms

    % Coherence  [nSubj × nFreqBins] — filled after first subject
    coh_subj = [];
    coh_f    = [];

    % Time-resolved alpha–velocity correlation [nSubj × nWindows]
    % (window count depends on epoch, filled after first subject)
    timeCorr_subj = [];
    timeCorr_t    = [];

    % ------------------------------------------------------------------
    %                       SUBJECT LOOP
    % ------------------------------------------------------------------
    for si = 1:nSubj
        subjectID = subjects{si};
        fprintf('  %s | Subject %s (%d/%d) ... ', taskName, subjectID, si, nSubj);

        % ---- 1a. Load & epoch all blocks ----------------------------
        epoched_blocks = {};
        for block = 1:6
            blockFile = fullfile(mergedPath, subjectID, ...
                [subjectID, sprintf(blockPat, block)]);
            if ~exist(blockFile, 'file'), continue; end
            try
                tmp = load(blockFile, 'EEG');
                EEP = pop_epoch(tmp.EEG, triggers, epochWin, ...
                    'newname', sprintf('%s_%s_b%d', subjectID, taskName, block), ...
                    'epochinfo', 'yes');
                if EEP.trials > 0
                    epoched_blocks{end+1} = EEP; %#ok<SAGROW>
                end
                clear tmp EEP
            catch
            end
        end
        if isempty(epoched_blocks)
            fprintf('NO DATA — skipped.\n');
            continue
        end

        % ---- 1b. Identify channels (from first block) ---------------
        EEGtmp = epoched_blocks{1};
        allLabels = {EEGtmp.chanlocs.labels};

        % ET channels: contain GAZE or AREA (underscores or dashes)
        isET = cellfun(@(x) ~isempty(regexpi(x, 'GAZE|AREA')), allLabels);
        et_ch_idx = find(isET);
        % Find left-eye X, Y, Pupil among ET channels
        lx_idx = find(cellfun(@(x) ~isempty(regexpi(x, '^L.GAZE.X$')), allLabels), 1);
        ly_idx = find(cellfun(@(x) ~isempty(regexpi(x, '^L.GAZE.Y$')), allLabels), 1);
        la_idx = find(cellfun(@(x) ~isempty(regexpi(x, '^L.AREA$')),   allLabels), 1);
        if isempty(lx_idx) || isempty(ly_idx) || isempty(la_idx)
            fprintf('ET channels missing — skipped.\n');
            continue
        end
        et_xyz_idx = [lx_idx, ly_idx, la_idx]; % indices into full channel set

        % EOG channels
        eog_labels = {'HEOGR','HEOGL','VEOGU','VEOGL'};
        isEOG = ismember(upper(allLabels), upper(eog_labels));

        % Bad channels (B* prefix)
        isBad = cellfun(@(x) startsWith(x, 'B'), allLabels);

        % EEG channels = all - ET - EOG - bad
        eeg_ch_idx = find(~isET & ~isEOG & ~isBad);

        % Posterior ROI among EEG channels (contains O or I, or starts with P)
        eeg_labels = allLabels(eeg_ch_idx);
        isPosterior = cellfun(@(x) contains(x,'O') | contains(x,'I') | ...
                                    startsWith(x,'P'), eeg_labels);
        posterior_ch_idx = eeg_ch_idx(isPosterior);
        nPosterior = numel(posterior_ch_idx);
        if nPosterior == 0
            fprintf('No posterior channels — skipped.\n');
            continue
        end

        % Time vector (samples → seconds)
        timeVec = EEGtmp.xmin + (0:EEGtmp.pnts-1)/EEGtmp.srate;

        % Sample indices for analysis windows
        ana_samp  = timeVec >= anaWin(1) & timeVec <= anaWin(2);
        pre_samp  = timeVec >= preStimWin(1) & timeVec <= preStimWin(2);
        nAnaSamp  = sum(ana_samp);
        clear EEGtmp

        % ---- 1c. Collect all trials ----------------------------------
        % Store per-trial: posterior EEG [nPost × nPts], ET [3 × nPts], cond
        trial_eeg  = {};   % posterior EEG (full epoch)
        trial_et   = {};   % [X; Y; Pupil] (full epoch)
        trial_cond = [];

        for bi = 1:numel(epoched_blocks)
            EEP = epoched_blocks{bi};
            for tr = 1:EEP.trials
                % Condition code from epoch events
                evTypes = EEP.epoch(tr).eventtype;
                evLats  = EEP.epoch(tr).eventlatency;
                if iscell(evLats), evLats = cell2mat(evLats); end
                [~, zIdx] = min(abs(evLats));
                if iscell(evTypes)
                    condCode = str2double(evTypes{zIdx});
                else
                    condCode = evTypes(zIdx);
                    if ischar(condCode), condCode = str2double(condCode); end
                end

                trial_eeg{end+1}  = double(EEP.data(posterior_ch_idx, :, tr)); %#ok<SAGROW>
                trial_et{end+1}   = double(EEP.data(et_xyz_idx, :, tr));       %#ok<SAGROW>
                trial_cond(end+1) = condCode;                                  %#ok<SAGROW>
            end
        end
        clear epoched_blocks
        nTrials = numel(trial_cond);
        if nTrials < 5
            fprintf('%d trials — too few, skipped.\n', nTrials);
            continue
        end

        % ---- 1d. Compute IAF from pooled posterior PSD ---------------
        psd_pool = [];
        for tr = 1:nTrials
            seg = mean(trial_eeg{tr}(:, ana_samp), 1); % mean across posterior
            [pxx, f_welch] = pwelch(seg, hamming(min(nfft_welch, nAnaSamp)), ...
                round(min(nfft_welch, nAnaSamp)/2), ...
                nfft_welch, fs);
            psd_pool = [psd_pool; pxx']; %#ok<AGROW>
        end
        meanPSD = mean(psd_pool, 1);

        % Find IAF as the peak in 8-14 Hz
        alphaIdx = f_welch >= alphaRange(1) & f_welch <= alphaRange(2);
        alphaPSD = meanPSD(alphaIdx);
        alphaFreqs = f_welch(alphaIdx);
        [pks, locs] = findpeaks(alphaPSD);
        if isempty(pks)
            IAF = NaN;
        else
            [~, best] = max(pks);
            IAF = alphaFreqs(locs(best));
            % Reject edge peaks
            if IAF == alphaRange(1) || IAF == alphaRange(2)
                IAF = NaN;
            end
        end
        if isfinite(IAF)
            bandAlpha = [max(IAF-4, 1), IAF+2];
        else
            bandAlpha = alphaRange;
        end

        % ---- 1e. Design alpha bandpass filter (Butterworth, 4th order) --
        [bpB, bpA] = butter(4, bandAlpha / (fs/2), 'bandpass');

        % ---- 1f. Per-trial feature extraction -------------------------
        subj_features = nan(nTrials, 2 + nGazeMetrics); % [AlphaPower, PreStimAlpha, gaze×8]

        % Containers for time-series analyses (this subject)
        subj_periMS  = [];   % peri-MS alpha epochs [nEpochs × nPeriSamp]
        subj_periSur = [];   % surrogate
        subj_xcorr   = [];   % per-trial xcorr [nTrials × nLags] — average later
        subj_alpha_ts = [];  % alpha envelope in analysis window [nTrials × nAnaSamp]
        subj_vel_ts   = [];  % velocity magnitude in analysis window [nTrials × nAnaSamp]

        for tr = 1:nTrials
            eeg_tr = trial_eeg{tr};  % [nPosterior × nPts]
            et_tr  = trial_et{tr};   % [3 × nPts]  (X, Y, Pupil)

            % --- Alpha power (Welch, analysis window) ---
            seg_ana = mean(eeg_tr(:, ana_samp), 1);
            [pxx_tr, ~] = pwelch(seg_ana, hamming(min(nfft_welch, nAnaSamp)), ...
                round(min(nfft_welch, nAnaSamp)/2), nfft_welch, fs);
            bandIdx = f_welch >= bandAlpha(1) & f_welch <= bandAlpha(2);
            alphaPow = mean(pxx_tr(bandIdx));

            % --- Pre-stimulus alpha power ---
            nPreSamp = sum(pre_samp);
            if nPreSamp > 20
                seg_pre = mean(eeg_tr(:, pre_samp), 1);
                [pxx_pre, ~] = pwelch(seg_pre, hamming(min(nfft_welch, nPreSamp)), ...
                    round(min(nfft_welch, nPreSamp)/2), nfft_welch, fs);
                preAlpha = mean(pxx_pre(bandIdx));
            else
                preAlpha = NaN;
            end

            % --- Alpha envelope (full epoch, for time-series analyses) ---
            eeg_mean = mean(eeg_tr, 1); % average across posterior channels
            eeg_mean = fillmissing(eeg_mean, 'linear', 'EndValues', 'nearest');
            alpha_filt = filtfilt(bpB, bpA, eeg_mean);
            alpha_env  = abs(hilbert(alpha_filt));

            % --- Gaze: prepare data (analysis window) ---
            X_raw = et_tr(1, ana_samp);
            Y_raw = et_tr(2, ana_samp);
            A_raw = et_tr(3, ana_samp);

            % Remove off-screen & blinks
            offScreen = X_raw < 0 | X_raw > screenSize(1) | ...
                        Y_raw < 0 | Y_raw > screenSize(2);
            blink = ~isfinite(A_raw) | A_raw <= 0;
            bad_samps = offScreen | blink;
            X = X_raw; Y = Y_raw; A = A_raw;
            X(bad_samps) = NaN;
            Y(bad_samps) = NaN;

            % Invert Y for Cartesian convention
            Y = screenSize(2) - Y;

            % Interpolate NaN gaps
            X = fillmissing(X, 'linear', 'EndValues', 'nearest');
            Y = fillmissing(Y, 'linear', 'EndValues', 'nearest');

            % Skip trial if >50% data is bad
            if sum(bad_samps) / numel(bad_samps) > 0.5
                subj_features(tr, :) = NaN;
                continue
            end

            % --- Gaze metrics ---
            % 1. Gaze deviation (Euclidean from center)
            dx = X - screenCenter(1);
            dy = Y - screenCenter(2);
            gazeDev = mean(sqrt(dx.^2 + dy.^2), 'omitnan');

            % 2-3. Gaze variability
            gazeStdX = std(X, 'omitnan');
            gazeStdY = std(Y, 'omitnan');

            % 4. Pupil size
            pupil_vals = A_raw(~bad_samps);
            pupilSize = mean(pupil_vals, 'omitnan') / 1000;

            % 5. Microsaccade rate (Engbert & Kliegl)
            ms = ms_detect_engbert_local(X, Y, fs, ms_lambda, ms_minDurSec, ms_minISISec);
            T_dur = nAnaSamp / fs;
            msRate = numel(ms.onsets) / T_dur;

            % 6-7. Eye velocity (Savitzky-Golay)
            [vx, vy] = compute_velocity_sg_local(X, Y, fs, 3);
            velMag = sqrt(vx.^2 + vy.^2);
            meanVel = mean(velMag, 'omitnan');
            peakVel = max(velMag);
            if isempty(peakVel), peakVel = NaN; end

            % 8. Scan path length
            dxf = diff(X); dyf = diff(Y);
            spl = sum(sqrt(dxf.^2 + dyf.^2), 'omitnan');

            % Store features
            subj_features(tr, :) = [alphaPow, preAlpha, ...
                gazeDev, gazeStdX, gazeStdY, pupilSize, ...
                msRate, meanVel, peakVel, spl];

            % --- Peri-microsaccadic alpha envelope ---
            alpha_env_full = alpha_env; % full epoch
            % Get gaze data for full epoch for MS detection
            X_full = et_tr(1, :);
            Y_full = et_tr(2, :);
            A_full = et_tr(3, :);
            bad_full = ~isfinite(A_full) | A_full <= 0 | ...
                       X_full < 0 | X_full > screenSize(1) | ...
                       Y_full < 0 | Y_full > screenSize(2);
            X_full(bad_full) = NaN; Y_full(bad_full) = NaN;
            Y_full = screenSize(2) - Y_full; % invert Y for Cartesian convention

            % Detect MS in analysis window of full epoch
            ana_start_samp = find(ana_samp, 1, 'first');
            X_ana_full = double(et_tr(1, ana_samp));
            Y_ana_full = double(et_tr(2, ana_samp));
            bad_af = ~isfinite(et_tr(3, ana_samp)) | et_tr(3, ana_samp) <= 0;
            X_ana_full(bad_af) = NaN; Y_ana_full(bad_af) = NaN;
            Y_ana_full = screenSize(2) - Y_ana_full;
            X_ana_full = fillmissing(X_ana_full, 'linear', 'EndValues', 'nearest');
            Y_ana_full = fillmissing(Y_ana_full, 'linear', 'EndValues', 'nearest');
            ms_ana = ms_detect_engbert_local(X_ana_full, Y_ana_full, fs, ms_lambda, ms_minDurSec, ms_minISISec);

            for mi = 1:numel(ms_ana.onsets)
                % Convert analysis-window index to full-epoch index
                eeg_idx_ms = ana_start_samp + ms_ana.onsets(mi) - 1;
                lo = eeg_idx_ms - periSamp;
                hi = eeg_idx_ms + periSamp;
                if lo < 1 || hi > numel(alpha_env_full), continue; end
                epoch_ms = alpha_env_full(lo:hi);
                if any(isnan(epoch_ms)), continue; end
                subj_periMS = [subj_periMS; epoch_ms]; %#ok<AGROW>

                % Surrogate: random time in analysis window
                lo_valid = max(ana_start_samp, 1 + periSamp);
                hi_valid = min(ana_start_samp + nAnaSamp - 1, numel(alpha_env_full) - periSamp);
                if hi_valid > lo_valid
                    surr_c = randi([lo_valid, hi_valid]);
                    epoch_surr = alpha_env_full(surr_c - periSamp : surr_c + periSamp);
                    if ~any(isnan(epoch_surr))
                        subj_periSur = [subj_periSur; epoch_surr]; %#ok<AGROW>
                    end
                end
            end

            % --- Cross-correlation: alpha envelope × velocity in analysis window ---
            ae_ana = alpha_env(ana_samp);
            vel_ana = velMag;
            if numel(ae_ana) == numel(vel_ana) && all(isfinite(ae_ana)) && all(isfinite(vel_ana))
                [xc, ~] = xcorr(zscore(ae_ana), zscore(vel_ana), maxLagSamp, 'coeff');
                subj_xcorr = [subj_xcorr; xc(:)']; %#ok<AGROW>
            end

            % --- Store alpha envelope & velocity for time-resolved analysis ---
            if all(isfinite(ae_ana)) && all(isfinite(vel_ana))
                subj_alpha_ts = [subj_alpha_ts; ae_ana]; %#ok<AGROW>
                subj_vel_ts   = [subj_vel_ts;   vel_ana]; %#ok<AGROW>
            end

        end % trial loop

        % ---- 1g. Aggregate subject-level time-series results ----------

        % Peri-microsaccadic alpha (subject mean)
        nMS_subj(si) = size(subj_periMS, 1);
        if ~isempty(subj_periMS)
            periMS_subj(si, :) = mean(subj_periMS, 1, 'omitnan');
        end
        if ~isempty(subj_periSur)
            periSurr_subj(si, :) = mean(subj_periSur, 1, 'omitnan');
        end

        % Cross-correlation (subject mean)
        if ~isempty(subj_xcorr)
            xcorr_subj(si, :) = mean(subj_xcorr, 1, 'omitnan');
        end

        % Coherence (subject mean, computed from concatenated trials)
        if size(subj_alpha_ts, 1) >= 5
            ae_cat  = reshape(subj_alpha_ts', 1, []);
            vel_cat = reshape(subj_vel_ts', 1, []);
            [cxy, f_coh] = mscohere(ae_cat, vel_cat, ...
                hamming(coh_nfft), coh_novlap, coh_nfft, fs);
            if isempty(coh_f)
                coh_f = f_coh;
                coh_subj = nan(nSubj, numel(f_coh));
            end
            coh_subj(si, 1:numel(cxy)) = cxy;
        end

        % Time-resolved correlation (sliding window over analysis window)
        if size(subj_alpha_ts, 1) >= 10
            winSamp  = round(slidWin_s * fs);
            stepSamp = round(slidStep_s * fs);
            nWin = floor((nAnaSamp - winSamp) / stepSamp) + 1;
            if isempty(timeCorr_t)
                anaTimes = linspace(anaWin(1), anaWin(2), nAnaSamp);
                winCenters = nan(nWin, 1);
                for wi = 1:nWin
                    s1 = (wi-1)*stepSamp + 1;
                    s2 = s1 + winSamp - 1;
                    winCenters(wi) = mean(anaTimes([s1, min(s2, nAnaSamp)]));
                end
                timeCorr_t = winCenters;
                timeCorr_subj = nan(nSubj, nWin);
            end
            for wi = 1:nWin
                s1 = (wi-1)*stepSamp + 1;
                s2 = min(s1 + winSamp - 1, nAnaSamp);
                a_w = subj_alpha_ts(:, s1:s2);
                v_w = subj_vel_ts(:, s1:s2);
                % Correlation across trials at each window
                a_mean = mean(a_w, 2);
                v_mean = mean(v_w, 2);
                if std(a_mean) > 0 && std(v_mean) > 0
                    timeCorr_subj(si, wi) = corr(a_mean, v_mean, ...
                        'Type', 'Spearman', 'Rows', 'complete');
                end
            end
        end

        % ---- 1h. Append to group trial table --------------------------
        valid = ~all(isnan(subj_features), 2);
        nValid = sum(valid);
        subjCol  = repmat(str2double(subjectID), nValid, 1);
        condCol  = trial_cond(valid)' - condOffset;
        featMat  = subj_features(valid, :);

        % Columns: SubjID, Condition, AlphaPower, PreStimAlpha, Gaze×8
        allTrialData = [allTrialData; subjCol, condCol, featMat]; %#ok<AGROW>

        fprintf('%d trials extracted.\n', nValid);
        clear trial_eeg trial_et trial_cond subj_features
        clear subj_periMS subj_periSur subj_xcorr subj_alpha_ts subj_vel_ts
    end % subject loop

    % ---- 1i. Build MATLAB table from allTrialData --------------------
    colNames = ['SubjID', 'Condition', 'AlphaPower', 'PreStimAlpha', gazeMetricNames];
    T = array2table(allTrialData, 'VariableNames', colNames);

    % Save intermediate features
    save(fullfile(dataDir, ['trialdata_' taskName '.mat']), 'T', '-v7.3');
    save(fullfile(dataDir, ['timeseries_' taskName '.mat']), ...
        'periMS_subj', 'periSurr_subj', 'nMS_subj', ...
        'xcorr_subj', 'xcorr_lags', ...
        'coh_subj', 'coh_f', ...
        'timeCorr_subj', 'timeCorr_t', ...
        'subjects', '-v7.3');

    fprintf('\n  Features saved for %s (%d total trials).\n', taskName, height(T));

    % ==================================================================
    %% 2  ANALYSIS 1: Trial-level within-subject correlations
    % ==================================================================
    fprintf('\n--- ANALYSIS 1: Trial-level correlations (%s) ---\n', taskName);

    uniqueSubj = unique(T.SubjID);
    nS = numel(uniqueSubj);

    % For each gaze metric: within-subject Spearman r with AlphaPower
    r_ws = nan(nS, nGazeMetrics);
    for si = 1:nS
        idx = T.SubjID == uniqueSubj(si);
        a = T.AlphaPower(idx);
        for gi = 1:nGazeMetrics
            g = T.(gazeMetricNames{gi})(idx);
            vld = isfinite(a) & isfinite(g);
            if sum(vld) >= 5
                r_ws(si, gi) = corr(a(vld), g(vld), 'Type', 'Spearman');
            end
        end
    end

    % Fisher-z transform → one-sample t-test
    z_ws = atanh(r_ws);
    p_corr   = nan(1, nGazeMetrics);
    t_corr   = nan(1, nGazeMetrics);
    df_corr  = nan(1, nGazeMetrics);
    mean_r   = nan(1, nGazeMetrics);
    ci_lo    = nan(1, nGazeMetrics);
    ci_hi    = nan(1, nGazeMetrics);
    for gi = 1:nGazeMetrics
        vld = isfinite(z_ws(:, gi));
        if sum(vld) >= 3
            zv = z_ws(vld, gi);
            [~, p_corr(gi), ~, stats] = ttest(zv);
            t_corr(gi) = stats.tstat;
            df_corr(gi) = stats.df;
            mz = mean(zv);
            se = std(zv) / sqrt(numel(zv));
            mean_r(gi) = tanh(mz);
            ci_lo(gi)  = tanh(mz - 1.96*se);
            ci_hi(gi)  = tanh(mz + 1.96*se);
            fprintf('  Alpha × %-18s: r = %+.3f [%+.3f, %+.3f], t(%d) = %+.2f, p = %.4f\n', ...
                gazeMetricNames{gi}, mean_r(gi), ci_lo(gi), ci_hi(gi), ...
                df_corr(gi), t_corr(gi), p_corr(gi));
        end
    end

    % FDR correction (Benjamini-Hochberg)
    p_fdr_corr = fdr_bh_local(p_corr);

    % --- Figure: Forest plot of all within-subject correlations ---
    close all
    figure('Position', [50 50 1000 700]); hold on
    for gi = 1:nGazeMetrics
        col = cmap_corr;
        if p_fdr_corr(gi) < 0.05, col = [0.85 0.2 0.2]; end
        plot([ci_lo(gi) ci_hi(gi)], [gi gi], '-', 'Color', col, 'LineWidth', 2.5);
        plot(mean_r(gi), gi, 'o', 'MarkerSize', 10, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    end
    xline(0, 'k:', 'LineWidth', 1.5);
    set(gca, 'YTick', 1:nGazeMetrics, 'YTickLabel', gazeMetricNames, ...
        'FontSize', fontSize-4, 'YDir', 'reverse', 'TickLabelInterpreter', 'none');
    xlabel('Within-subject Spearman r (Alpha × Gaze)', 'FontSize', fontSize-2);
    title(sprintf('%s: Trial-Level Alpha-Gaze Correlations', upper(taskName)), ...
        'FontSize', fontSize);
    xlim([-0.4 0.4]);
    saveas(gcf, fullfile(figDir, ['AOC_eeglab_corr_forest_' taskName '.png']));

    % --- Figure: Scatter plots for top 4 metrics (by |mean_r|) ---
    [~, sortIdx] = sort(abs(mean_r), 'descend');
    topN = min(4, nGazeMetrics);
    figure('Position', [50 50 1800 500]);
    for pp = 1:topN
        gi = sortIdx(pp);
        subplot(1, topN, pp); hold on
        % Remove between-subject means for display
        alpha_ws  = T.AlphaPower;
        gaze_ws   = T.(gazeMetricNames{gi});
        for ssi = 1:nS
            idx = T.SubjID == uniqueSubj(ssi);
            alpha_ws(idx) = alpha_ws(idx) - mean(alpha_ws(idx), 'omitnan');
            gaze_ws(idx)  = gaze_ws(idx)  - mean(gaze_ws(idx),  'omitnan');
        end
        % Plot by condition
        uConds = sort(unique(T.Condition));
        for ci = 1:numel(uConds)
            idx = T.Condition == uConds(ci);
            scatter(alpha_ws(idx), gaze_ws(idx), 12, colors3(ci,:), ...
                'filled', 'MarkerFaceAlpha', 0.25);
        end
        % Regression line
        vld = isfinite(alpha_ws) & isfinite(gaze_ws);
        if sum(vld) > 2
            pfit = polyfit(alpha_ws(vld), gaze_ws(vld), 1);
            xl = xlim; xf = linspace(xl(1), xl(2), 100);
            plot(xf, polyval(pfit, xf), 'k-', 'LineWidth', 2.5);
        end
        xlabel('Alpha (ws-centered)'); ylabel(gazeMetricNames{gi});
        title(sprintf('r=%.3f, p=%.3f', mean_r(gi), p_corr(gi)), 'FontSize', fontSize-6);
        set(gca, 'FontSize', fontSize-6);
    end
    sgtitle(sprintf('%s: Top Alpha-Gaze Scatter Plots', upper(taskName)), 'FontSize', fontSize);
    saveas(gcf, fullfile(figDir, ['AOC_eeglab_corr_scatter_' taskName '.png']));

    % Save correlation results
    corrResults.(taskName).r_ws     = r_ws;
    corrResults.(taskName).mean_r   = mean_r;
    corrResults.(taskName).ci       = [ci_lo; ci_hi];
    corrResults.(taskName).p        = p_corr;
    corrResults.(taskName).p_fdr    = p_fdr_corr;
    corrResults.(taskName).t        = t_corr;
    corrResults.(taskName).df       = df_corr;

    % ==================================================================
    %% 3  ANALYSIS 2: Median split on alpha → gaze comparison
    % ==================================================================
    fprintf('\n--- ANALYSIS 2: Median split (%s) ---\n', taskName);

    % Within-subject median split
    hiAlpha_mean = nan(nS, nGazeMetrics);
    loAlpha_mean = nan(nS, nGazeMetrics);

    for si = 1:nS
        idx = T.SubjID == uniqueSubj(si);
        a = T.AlphaPower(idx);
        med = median(a, 'omitnan');
        hiIdx = idx & T.AlphaPower >= med;
        loIdx = idx & T.AlphaPower <  med;
        for gi = 1:nGazeMetrics
            hiAlpha_mean(si, gi) = mean(T.(gazeMetricNames{gi})(hiIdx), 'omitnan');
            loAlpha_mean(si, gi) = mean(T.(gazeMetricNames{gi})(loIdx), 'omitnan');
        end
    end

    % Paired t-tests
    p_median   = nan(1, nGazeMetrics);
    t_median   = nan(1, nGazeMetrics);
    d_median   = nan(1, nGazeMetrics); % Cohen's d
    for gi = 1:nGazeMetrics
        hi = hiAlpha_mean(:, gi);
        lo = loAlpha_mean(:, gi);
        vld = isfinite(hi) & isfinite(lo);
        if sum(vld) >= 3
            [~, p_median(gi), ~, stats] = ttest(hi(vld), lo(vld));
            t_median(gi) = stats.tstat;
            d_median(gi) = mean(hi(vld) - lo(vld)) / std(hi(vld) - lo(vld));
        end
        fprintf('  %-18s: High-Low diff t = %+.2f, p = %.4f, d = %+.3f\n', ...
            gazeMetricNames{gi}, t_median(gi), p_median(gi), d_median(gi));
    end
    p_fdr_median = fdr_bh_local(p_median);

    % --- Figure: Grouped bar plot (high vs low alpha) ---
    close all
    figure('Position', [50 50 1400 600]);
    for gi = 1:nGazeMetrics
        subplot(2, 4, gi); hold on
        hi_vals = hiAlpha_mean(:, gi);
        lo_vals = loAlpha_mean(:, gi);
        vld = isfinite(hi_vals) & isfinite(lo_vals);
        mH = mean(hi_vals(vld)); seH = std(hi_vals(vld)) / sqrt(sum(vld));
        mL = mean(lo_vals(vld)); seL = std(lo_vals(vld)) / sqrt(sum(vld));
        bar(1, mL, 0.6, 'FaceColor', colors3(1,:), 'EdgeColor', 'none');
        bar(2, mH, 0.6, 'FaceColor', colors3(3,:), 'EdgeColor', 'none');
        errorbar([1 2], [mL mH], [seL seH], 'k.', 'LineWidth', 1.5, 'CapSize', 8);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \alpha', 'High \alpha'}, ...
            'FontSize', fontSize-8);
        title(sprintf('%s\np=%.3f', gazeMetricNames{gi}, p_median(gi)), 'FontSize', fontSize-8);
        if p_fdr_median(gi) < 0.05
            text(1.5, max(mL, mH)+seH*1.5, '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
        end
    end
    sgtitle(sprintf('%s: Median Split by Alpha Power', upper(taskName)), 'FontSize', fontSize);
    saveas(gcf, fullfile(figDir, ['AOC_eeglab_mediansplit_' taskName '.png']));

    medianResults.(taskName).p     = p_median;
    medianResults.(taskName).p_fdr = p_fdr_median;
    medianResults.(taskName).t     = t_median;
    medianResults.(taskName).d     = d_median;
    medianResults.(taskName).hiMean = hiAlpha_mean;
    medianResults.(taskName).loMean = loAlpha_mean;

    % ==================================================================
    %% 4  ANALYSIS 3: GLMMs
    % ==================================================================
    fprintf('\n--- ANALYSIS 3: GLMMs (%s) ---\n', taskName);

    % Z-score alpha power within subject
    T.AlphaZ = nan(height(T), 1);
    for si = 1:nS
        idx = T.SubjID == uniqueSubj(si);
        vals = T.AlphaPower(idx);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            T.AlphaZ(idx) = (vals - mu) / sd;
        else
            T.AlphaZ(idx) = 0;
        end
    end

    T.CondCat = categorical(T.Condition);
    T.SubjCat = categorical(T.SubjID);

    glmmResults.(taskName) = struct();
    glmm_pvals = nan(1, nGazeMetrics);
    glmm_beta  = nan(1, nGazeMetrics);

    for gi = 1:nGazeMetrics
        metricName = gazeMetricNames{gi};
        try
            % Full model
            formula_full = sprintf('%s ~ AlphaZ * CondCat + (1 | SubjCat)', metricName);
            formula_red  = sprintf('%s ~ AlphaZ + CondCat + (1 | SubjCat)', metricName);
            lme_full = fitlme(T, formula_full);
            lme_red  = fitlme(T, formula_red);

            % LRT for interaction
            lrt = compare(lme_red, lme_full);
            if lrt.pValue(2) < 0.05
                lme_use = lme_full;
                intFlag = ' [interaction]';
            else
                lme_use = lme_red;
                intFlag = '';
            end

            % Extract AlphaZ coefficient
            coefNames = lme_use.CoefficientNames;
            aIdx = find(strcmp(coefNames, 'AlphaZ'));
            if ~isempty(aIdx)
                beta = lme_use.Coefficients.Estimate(aIdx);
                pval = lme_use.Coefficients.pValue(aIdx);
                glmm_pvals(gi) = pval;
                glmm_beta(gi)  = beta;
                fprintf('  GLMM %-18s: beta(AlphaZ) = %+.4f, p = %.4f%s\n', ...
                    metricName, beta, pval, intFlag);
            end

            glmmResults.(taskName).(metricName).lme  = lme_use;
            glmmResults.(taskName).(metricName).lrt  = lrt;
        catch ME
            fprintf('  GLMM %-18s: FAILED — %s\n', metricName, ME.message);
        end
    end

    p_fdr_glmm = fdr_bh_local(glmm_pvals);

    % --- Figure: GLMM coefficient plot ---
    close all
    figure('Position', [50 50 900 600]); hold on
    for gi = 1:nGazeMetrics
        col = cmap_corr;
        if isfinite(p_fdr_glmm(gi)) && p_fdr_glmm(gi) < 0.05, col = [0.85 0.2 0.2]; end
        if isfinite(glmm_beta(gi))
            bar(gi, glmm_beta(gi), 0.6, 'FaceColor', col, 'EdgeColor', 'none');
        end
    end
    set(gca, 'XTick', 1:nGazeMetrics, 'XTickLabel', gazeMetricNames, ...
        'FontSize', fontSize-6, 'XTickLabelRotation', 35, 'TickLabelInterpreter', 'none');
    ylabel('\beta (AlphaZ)', 'FontSize', fontSize-2);
    title(sprintf('%s: GLMM Coefficients (Alpha → Gaze)', upper(taskName)), 'FontSize', fontSize);
    yline(0, 'k:');
    saveas(gcf, fullfile(figDir, ['AOC_eeglab_glmm_coefs_' taskName '.png']));

    glmmResults.(taskName).p_fdr = p_fdr_glmm;

    % ==================================================================
    %% 5  ANALYSIS 4: Peri-microsaccadic alpha envelope
    % ==================================================================
    fprintf('\n--- ANALYSIS 4: Peri-microsaccadic alpha (%s) ---\n', taskName);
    fprintf('  MS per subject: %s\n', mat2str(nMS_subj'));

    % Exclude subjects without data
    validS = all(isfinite(periMS_subj), 2) & all(isfinite(periSurr_subj), 2) & nMS_subj > 0;
    PM = periMS_subj(validS, :);
    PS = periSurr_subj(validS, :);
    nS_peri = sum(validS);
    fprintf('  Valid subjects for peri-MS: %d / %d\n', nS_peri, nSubj);

    if nS_peri >= 5
        % Baseline correction (% change, baseline = -1000 to -750 ms)
        bl_idx = periTime_ms >= -1000 & periTime_ms <= -750;
        PM_bl = nan(size(PM)); PS_bl = nan(size(PS));
        for si = 1:nS_peri
            bl_r = mean(PM(si, bl_idx), 'omitnan');
            bl_s = mean(PS(si, bl_idx), 'omitnan');
            if bl_r > 0, PM_bl(si, :) = (PM(si, :) - bl_r) / bl_r * 100; end
            if bl_s > 0, PS_bl(si, :) = (PS(si, :) - bl_s) / bl_s * 100; end
        end
        validBL = all(isfinite(PM_bl), 2) & all(isfinite(PS_bl), 2);
        PM_bl = PM_bl(validBL, :); PS_bl = PS_bl(validBL, :);
        nS_bl = sum(validBL);

        ga_real = mean(PM_bl, 1); sem_real = std(PM_bl,[],1) / sqrt(nS_bl);
        ga_surr = mean(PS_bl, 1); sem_surr = std(PS_bl,[],1) / sqrt(nS_bl);
        diff_mat = PM_bl - PS_bl;
        ga_diff  = mean(diff_mat, 1); sem_diff = std(diff_mat,[],1) / sqrt(nS_bl);

        % Cluster permutation on difference
        [clusters, tvals, thr] = cluster_perm_1d_local(diff_mat, nPerm, alphaCBP);

        nSig = 0;
        for k = 1:numel(clusters)
            if clusters(k).mass >= thr.mass
                nSig = nSig + 1;
                fprintf('  Cluster %d: [%.0f, %.0f] ms, p = %.4f\n', ...
                    nSig, periTime_ms(clusters(k).idx(1)), ...
                    periTime_ms(clusters(k).idx(end)), clusters(k).p);
            end
        end
        if nSig == 0, fprintf('  No significant clusters.\n'); end

        % --- Figure: Peri-MS alpha (real vs surrogate) ---
        close all
        figure('Position', [50 50 1200 700]); hold on
        fill([periTime_ms, fliplr(periTime_ms)], ...
            [ga_surr-sem_surr, fliplr(ga_surr+sem_surr)], ...
            [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h_s = plot(periTime_ms, ga_surr, '--', 'LineWidth', 2.5, 'Color', [0.5 0.5 0.5]);
        fill([periTime_ms, fliplr(periTime_ms)], ...
            [ga_real-sem_real, fliplr(ga_real+sem_real)], ...
            colors3(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        h_r = plot(periTime_ms, ga_real, '-', 'LineWidth', 3, 'Color', colors3(1,:));
        xline(0, 'k--', 'LineWidth', 2); yline(0, 'k:');
        xlabel('Time relative to MS onset [ms]', 'FontSize', fontSize-2);
        ylabel('Alpha envelope (% change)', 'FontSize', fontSize-2);
        title(sprintf('%s: Peri-Microsaccadic Alpha', upper(taskName)), 'FontSize', fontSize);
        legend([h_r h_s], {'Microsaccade', 'Surrogate'}, 'FontSize', fontSize-6, 'Location', 'best');
        set(gca, 'FontSize', fontSize-2);
        saveas(gcf, fullfile(figDir, ['AOC_eeglab_periMS_alpha_' taskName '.png']));

        % --- Figure: Difference wave with significant clusters ---
        figure('Position', [50 50 1200 500]); hold on
        fill([periTime_ms, fliplr(periTime_ms)], ...
            [ga_diff-sem_diff, fliplr(ga_diff+sem_diff)], ...
            colors3(2,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        plot(periTime_ms, ga_diff, '-', 'LineWidth', 3, 'Color', colors3(2,:));
        for k = 1:numel(clusters)
            if clusters(k).mass >= thr.mass
                ci = clusters(k).idx;
                patch([periTime_ms(ci), fliplr(periTime_ms(ci))], ...
                    [zeros(size(ci)), fliplr(ga_diff(ci))], ...
                    colors3(2,:), 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            end
        end
        xline(0, 'k--', 'LineWidth', 2); yline(0, 'k-');
        xlabel('Time relative to MS onset [ms]', 'FontSize', fontSize-2);
        ylabel('Difference (Real - Surrogate) [%]', 'FontSize', fontSize-2);
        title(sprintf('%s: Peri-MS Alpha Difference', upper(taskName)), 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize-2);
        saveas(gcf, fullfile(figDir, ['AOC_eeglab_periMS_diff_' taskName '.png']));

        periMSResults.(taskName).ga_real   = ga_real;
        periMSResults.(taskName).ga_surr   = ga_surr;
        periMSResults.(taskName).ga_diff   = ga_diff;
        periMSResults.(taskName).clusters  = clusters;
        periMSResults.(taskName).tvals     = tvals;
        periMSResults.(taskName).thr       = thr;
    else
        fprintf('  Too few subjects for peri-MS analysis.\n');
    end

    % ==================================================================
    %% 6  ANALYSIS 5: Alpha × velocity cross-correlation
    % ==================================================================
    fprintf('\n--- ANALYSIS 5: Cross-correlation (%s) ---\n', taskName);

    validXC = all(isfinite(xcorr_subj), 2);
    if sum(validXC) >= 5
        XC = xcorr_subj(validXC, :);
        nS_xc = sum(validXC);
        ga_xc = mean(XC, 1);
        sem_xc = std(XC, [], 1) / sqrt(nS_xc);

        % One-sample t-test at each lag (vs 0), FDR corrected
        p_xc = nan(1, nLags);
        t_xc = nan(1, nLags);
        for li = 1:nLags
            [~, p_xc(li), ~, st] = ttest(XC(:, li));
            t_xc(li) = st.tstat;
        end
        p_xc_fdr = fdr_bh_local(p_xc);
        sig_lags = p_xc_fdr < 0.05;

        [~, peakIdx] = max(abs(ga_xc));
        fprintf('  Peak xcorr: r = %.3f at lag = %.0f ms (p_fdr = %.4f)\n', ...
            ga_xc(peakIdx), xcorr_lags(peakIdx), p_xc_fdr(peakIdx));

        % --- Figure: Cross-correlation ---
        close all
        figure('Position', [50 50 1100 500]); hold on
        fill([xcorr_lags, fliplr(xcorr_lags)], ...
            [ga_xc-sem_xc, fliplr(ga_xc+sem_xc)], ...
            colors3(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(xcorr_lags, ga_xc, '-', 'LineWidth', 3, 'Color', colors3(1,:));
        if any(sig_lags)
            plot(xcorr_lags(sig_lags), ga_xc(sig_lags), '.', 'Color', [0.85 0.2 0.2], 'MarkerSize', 15);
        end
        xline(0, 'k:', 'LineWidth', 1.5); yline(0, 'k:');
        xlabel('Lag [ms] (positive = alpha leads)', 'FontSize', fontSize-2);
        ylabel('Cross-correlation (normalized)', 'FontSize', fontSize-2);
        title(sprintf('%s: Alpha Envelope × Eye Velocity', upper(taskName)), 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize-2);
        saveas(gcf, fullfile(figDir, ['AOC_eeglab_xcorr_' taskName '.png']));

        xcorrResults.(taskName).ga   = ga_xc;
        xcorrResults.(taskName).sem  = sem_xc;
        xcorrResults.(taskName).p    = p_xc;
        xcorrResults.(taskName).pfdr = p_xc_fdr;
        xcorrResults.(taskName).lags = xcorr_lags;
    else
        fprintf('  Too few subjects for cross-correlation.\n');
    end

    % ==================================================================
    %% 7  ANALYSIS 6: Alpha × velocity coherence
    % ==================================================================
    fprintf('\n--- ANALYSIS 6: Coherence (%s) ---\n', taskName);

    if ~isempty(coh_f)
        validCoh = all(isfinite(coh_subj), 2);
        if sum(validCoh) >= 5
            COH = coh_subj(validCoh, :);
            nS_coh = sum(validCoh);
            ga_coh = mean(COH, 1);
            sem_coh = std(COH, [], 1) / sqrt(nS_coh);

            % Focus on 0-30 Hz
            fMask = coh_f <= 30;
            f_plot = coh_f(fMask);
            ga_coh_plot = ga_coh(fMask);
            sem_coh_plot = sem_coh(fMask);

            % Highlight alpha band
            aband = f_plot >= bandAlpha(1) & f_plot <= bandAlpha(2);
            alpha_coh = mean(ga_coh_plot(aband));
            fprintf('  Mean coherence in alpha band: %.4f\n', alpha_coh);

            % --- Figure: Coherence spectrum ---
            close all
            figure('Position', [50 50 1000 500]); hold on
            fill([f_plot', fliplr(f_plot')], ...
                [ga_coh_plot(1:numel(f_plot))-sem_coh_plot(1:numel(f_plot)), ...
                 fliplr(ga_coh_plot(1:numel(f_plot))+sem_coh_plot(1:numel(f_plot)))], ...
                colors3(2,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            plot(f_plot, ga_coh_plot(1:numel(f_plot)), '-', 'LineWidth', 3, 'Color', colors3(2,:));
            % Shade alpha band
            patch([bandAlpha(1) bandAlpha(2) bandAlpha(2) bandAlpha(1)], ...
                [0 0 max(ga_coh_plot)*1.1 max(ga_coh_plot)*1.1], ...
                [0.9 0.9 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            xlabel('Frequency [Hz]', 'FontSize', fontSize-2);
            ylabel('Magnitude-Squared Coherence', 'FontSize', fontSize-2);
            title(sprintf('%s: Alpha-Velocity Coherence', upper(taskName)), 'FontSize', fontSize);
            set(gca, 'FontSize', fontSize-2);
            saveas(gcf, fullfile(figDir, ['AOC_eeglab_coherence_' taskName '.png']));

            cohResults.(taskName).ga  = ga_coh;
            cohResults.(taskName).f   = coh_f;
            cohResults.(taskName).sem = sem_coh;
        else
            fprintf('  Too few subjects for coherence.\n');
        end
    else
        fprintf('  No coherence data available.\n');
    end

    % ==================================================================
    %% 8  ANALYSIS 7: Pre-stimulus alpha → post-stimulus gaze
    % ==================================================================
    fprintf('\n--- ANALYSIS 7: Pre-stim alpha → post-stim gaze (%s) ---\n', taskName);

    r_pre = nan(nS, nGazeMetrics);
    for si = 1:nS
        idx = T.SubjID == uniqueSubj(si);
        a = T.PreStimAlpha(idx);
        for gi = 1:nGazeMetrics
            g = T.(gazeMetricNames{gi})(idx);
            vld = isfinite(a) & isfinite(g);
            if sum(vld) >= 5
                r_pre(si, gi) = corr(a(vld), g(vld), 'Type', 'Spearman');
            end
        end
    end

    z_pre = atanh(r_pre);
    p_pre   = nan(1, nGazeMetrics);
    t_pre   = nan(1, nGazeMetrics);
    mean_r_pre = nan(1, nGazeMetrics);
    ci_lo_pre  = nan(1, nGazeMetrics);
    ci_hi_pre  = nan(1, nGazeMetrics);
    for gi = 1:nGazeMetrics
        vld = isfinite(z_pre(:, gi));
        if sum(vld) >= 3
            zv = z_pre(vld, gi);
            [~, p_pre(gi), ~, stats] = ttest(zv);
            t_pre(gi) = stats.tstat;
            mz = mean(zv); se = std(zv)/sqrt(numel(zv));
            mean_r_pre(gi) = tanh(mz);
            ci_lo_pre(gi)  = tanh(mz - 1.96*se);
            ci_hi_pre(gi)  = tanh(mz + 1.96*se);
            fprintf('  PreAlpha × %-18s: r = %+.3f, t = %+.2f, p = %.4f\n', ...
                gazeMetricNames{gi}, mean_r_pre(gi), t_pre(gi), p_pre(gi));
        end
    end
    p_fdr_pre = fdr_bh_local(p_pre);

    % --- Figure: Forest plot for pre-stim alpha ---
    close all
    figure('Position', [50 50 1000 700]); hold on
    for gi = 1:nGazeMetrics
        col = [0.3 0.5 0.3];
        if isfinite(p_fdr_pre(gi)) && p_fdr_pre(gi) < 0.05, col = [0.85 0.2 0.2]; end
        if isfinite(ci_lo_pre(gi))
            plot([ci_lo_pre(gi) ci_hi_pre(gi)], [gi gi], '-', 'Color', col, 'LineWidth', 2.5);
            plot(mean_r_pre(gi), gi, 'o', 'MarkerSize', 10, 'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        end
    end
    xline(0, 'k:', 'LineWidth', 1.5);
    set(gca, 'YTick', 1:nGazeMetrics, 'YTickLabel', gazeMetricNames, ...
        'FontSize', fontSize-4, 'YDir', 'reverse', 'TickLabelInterpreter', 'none');
    xlabel('Within-subject Spearman r (PreStim Alpha × PostStim Gaze)', 'FontSize', fontSize-2);
    title(sprintf('%s: Pre-Stimulus Alpha → Post-Stimulus Gaze', upper(taskName)), 'FontSize', fontSize);
    xlim([-0.4 0.4]);
    saveas(gcf, fullfile(figDir, ['AOC_eeglab_prestim_forest_' taskName '.png']));

    preStimResults.(taskName).mean_r = mean_r_pre;
    preStimResults.(taskName).ci     = [ci_lo_pre; ci_hi_pre];
    preStimResults.(taskName).p      = p_pre;
    preStimResults.(taskName).p_fdr  = p_fdr_pre;

    % ==================================================================
    %% 9  ANALYSIS 8: Condition covariation (subject-level)
    % ==================================================================
    fprintf('\n--- ANALYSIS 8: Condition covariation (%s) ---\n', taskName);

    uConds = sort(unique(T.Condition));
    nCond  = numel(uConds);

    % Subject-level condition means
    subj_alpha_cond = nan(nS, nCond);
    subj_gaze_cond  = nan(nS, nGazeMetrics, nCond);
    for si = 1:nS
        for ci = 1:nCond
            idx = T.SubjID == uniqueSubj(si) & T.Condition == uConds(ci);
            if sum(idx) >= 2
                subj_alpha_cond(si, ci) = mean(T.AlphaPower(idx), 'omitnan');
                for gi = 1:nGazeMetrics
                    subj_gaze_cond(si, gi, ci) = mean(T.(gazeMetricNames{gi})(idx), 'omitnan');
                end
            end
        end
    end

    % Compute change: high load - low load (for 3 conditions)
    if nCond >= 2
        delta_alpha = subj_alpha_cond(:, end) - subj_alpha_cond(:, 1);
        delta_gaze  = squeeze(subj_gaze_cond(:, :, end) - subj_gaze_cond(:, :, 1));

        r_cond = nan(1, nGazeMetrics);
        p_cond = nan(1, nGazeMetrics);
        for gi = 1:nGazeMetrics
            vld = isfinite(delta_alpha) & isfinite(delta_gaze(:, gi));
            if sum(vld) >= 5
                [r_cond(gi), p_cond(gi)] = corr(delta_alpha(vld), delta_gaze(vld, gi), ...
                    'Type', 'Spearman');
                fprintf('  Δalpha × Δ%-18s: r = %+.3f, p = %.4f\n', ...
                    gazeMetricNames{gi}, r_cond(gi), p_cond(gi));
            end
        end

        % --- Figure: Scatter plots for condition covariation (top 4) ---
        [~, sIdx] = sort(abs(r_cond), 'descend');
        topN_cov = min(4, sum(isfinite(r_cond)));
        if topN_cov > 0
            close all
            figure('Position', [50 50 1600 450]);
            for pp = 1:topN_cov
                gi = sIdx(pp);
                subplot(1, topN_cov, pp); hold on
                vld = isfinite(delta_alpha) & isfinite(delta_gaze(:, gi));
                scatter(delta_alpha(vld), delta_gaze(vld, gi), 50, cmap_corr, ...
                    'filled', 'MarkerFaceAlpha', 0.6);
                if sum(vld) > 2
                    pf = polyfit(delta_alpha(vld), delta_gaze(vld, gi), 1);
                    xl = xlim; xf = linspace(xl(1), xl(2), 100);
                    plot(xf, polyval(pf, xf), 'k-', 'LineWidth', 2);
                end
                xlabel('\Delta Alpha Power');
                ylabel(['\Delta ' gazeMetricNames{gi}]);
                title(sprintf('r=%.3f, p=%.3f', r_cond(gi), p_cond(gi)), 'FontSize', fontSize-6);
                set(gca, 'FontSize', fontSize-6);
            end
            sgtitle(sprintf('%s: Condition Covariation (High-Low Load)', upper(taskName)), ...
                'FontSize', fontSize);
            saveas(gcf, fullfile(figDir, ['AOC_eeglab_condcovar_' taskName '.png']));
        end

        condCovarResults.(taskName).r = r_cond;
        condCovarResults.(taskName).p = p_cond;
    end

    % ==================================================================
    %% 10 ANALYSIS 9: Time-resolved alpha–velocity correlation
    % ==================================================================
    fprintf('\n--- ANALYSIS 9: Time-resolved correlation (%s) ---\n', taskName);

    if ~isempty(timeCorr_t)
        validTC = all(isfinite(timeCorr_subj), 2);
        if sum(validTC) >= 5
            TC = timeCorr_subj(validTC, :);
            nS_tc = sum(validTC);
            ga_tc = mean(TC, 1, 'omitnan');
            sem_tc = std(TC, [], 1, 'omitnan') / sqrt(nS_tc);

            % One-sample t-test at each window (vs 0), FDR
            p_tc = nan(1, numel(timeCorr_t));
            for wi = 1:numel(timeCorr_t)
                [~, p_tc(wi)] = ttest(TC(:, wi));
            end
            p_tc_fdr = fdr_bh_local(p_tc);
            sig_tc = p_tc_fdr < 0.05;

            % --- Figure: Time-resolved correlation ---
            close all
            figure('Position', [50 50 1200 500]); hold on
            fill([timeCorr_t', fliplr(timeCorr_t')], ...
                [ga_tc-sem_tc, fliplr(ga_tc+sem_tc)], ...
                colors3(3,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            plot(timeCorr_t, ga_tc, '-', 'LineWidth', 3, 'Color', colors3(3,:));
            if any(sig_tc)
                plot(timeCorr_t(sig_tc), ga_tc(sig_tc), '.', ...
                    'Color', [0.85 0.2 0.2], 'MarkerSize', 15);
            end
            yline(0, 'k:'); xlabel('Time [s]', 'FontSize', fontSize-2);
            ylabel('Spearman r (Alpha × Velocity)', 'FontSize', fontSize-2);
            title(sprintf('%s: Time-Resolved Alpha-Velocity Correlation', upper(taskName)), ...
                'FontSize', fontSize);
            set(gca, 'FontSize', fontSize-2);
            saveas(gcf, fullfile(figDir, ['AOC_eeglab_timeresolved_' taskName '.png']));

            timeResults.(taskName).ga    = ga_tc;
            timeResults.(taskName).sem   = sem_tc;
            timeResults.(taskName).t     = timeCorr_t;
            timeResults.(taskName).p     = p_tc;
            timeResults.(taskName).p_fdr = p_tc_fdr;
        else
            fprintf('  Too few subjects for time-resolved analysis.\n');
        end
    else
        fprintf('  No time-resolved data.\n');
    end

end % task loop

%% ====================================================================
%% 11  SUMMARY TABLE & SAVE
%% ====================================================================
fprintf('\n%s\n  SAVING ALL RESULTS\n%s\n', repmat('=',1,60), repmat('=',1,60));

% Build summary table
summaryRows = {};
for ti = 1:nTasks
    tn = taskDefs(ti).name;
    for gi = 1:nGazeMetrics
        gn = gazeMetricNames{gi};
        row = {upper(tn), gn, ...
            sprintf('%.3f', corrResults.(tn).mean_r(gi)), ...
            sprintf('%.4f', corrResults.(tn).p(gi)), ...
            sprintf('%.4f', corrResults.(tn).p_fdr(gi)), ...
            sprintf('%.2f', medianResults.(tn).t(gi)), ...
            sprintf('%.4f', medianResults.(tn).p(gi)), ...
            sprintf('%.3f', medianResults.(tn).d(gi))};
        summaryRows = [summaryRows; row]; %#ok<AGROW>
    end
end
summaryTable = cell2table(summaryRows, 'VariableNames', ...
    {'Task', 'GazeMetric', 'Corr_r', 'Corr_p', 'Corr_pFDR', ...
     'Median_t', 'Median_p', 'Median_d'});
disp(summaryTable);

% Save all results to data directory
save(fullfile(dataDir, 'AOC_eeglab_redo_results.mat'), ...
    'corrResults', 'medianResults', 'glmmResults', ...
    'preStimResults', 'summaryTable', '-v7.3');
try
    save(fullfile(dataDir, 'AOC_eeglab_redo_results.mat'), ...
        'periMSResults', 'xcorrResults', 'cohResults', ...
        'condCovarResults', 'timeResults', '-append');
catch
end

% Save summary as CSV
writetable(summaryTable, fullfile(dataDir, 'AOC_eeglab_redo_summary.csv'));

fprintf('\n  All results and figures saved.\n');
fprintf('  Figures: %s\n', figDir);
fprintf('  Data:    %s\n', dataDir);
fprintf('\nDONE.\n');

%% ====================================================================
%% LOCAL FUNCTIONS
%% ====================================================================

function ms = ms_detect_engbert_local(X, Y, fs, lambda, minDurSec, minISISec)
%MS_DETECT_ENGBERT_LOCAL  Detect microsaccades (Engbert & Kliegl, 2003).
%   Returns struct with fields 'onsets' and 'offsets' (sample indices).
if nargin < 4 || isempty(lambda),    lambda = 6;       end
if nargin < 5 || isempty(minDurSec), minDurSec = 0.006; end
if nargin < 6 || isempty(minISISec), minISISec = 0.02;  end

X = double(X(:)');
Y = double(Y(:)');

[vx, vy] = compute_velocity_sg_local(X, Y, fs, 3);

% Robust median-based SD (MAD)
sx = 1.4826 * median(abs(vx - median(vx, 'omitnan')), 'omitnan');
sy = 1.4826 * median(abs(vy - median(vy, 'omitnan')), 'omitnan');
tx = lambda * max(sx, eps);
ty = lambda * max(sy, eps);

rad2 = (vx ./ tx).^2 + (vy ./ ty).^2;
cand = isfinite(rad2) & (rad2 > 1);

d = diff([0, cand(:)', 0]);
on  = find(d == 1);
off = find(d == -1) - 1;

minDur = max(1, round(minDurSec * fs));
keep   = (off - on + 1) >= minDur;
on  = on(keep);
off = off(keep);

% Merge close events
minISI = max(1, round(minISISec * fs));
if ~isempty(on)
    O = on(1); F = off(1);
    onM = []; offM = [];
    for i = 2:numel(on)
        if on(i) - F <= minISI
            F = off(i);
        else
            onM  = [onM, O]; offM = [offM, F]; %#ok<AGROW>
            O = on(i); F = off(i);
        end
    end
    onM = [onM, O]; offM = [offM, F];
    on = onM; off = offM;
end

ms.onsets  = on(:);
ms.offsets = off(:);
end

function [vx, vy] = compute_velocity_sg_local(X, Y, fs, polyOrd)
%COMPUTE_VELOCITY_SG_LOCAL  Savitzky-Golay velocity estimation.
Ts = 1 / fs;
X = double(X(:)');
Y = double(Y(:)');
L = numel(X);

framelen = min(21, L);
if mod(framelen, 2) == 0, framelen = framelen - 1; end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0, minLegal = minLegal + 1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L
    framelen = L;
    if mod(framelen, 2) == 0, framelen = framelen - 1; end
end
useFallback = framelen < 5 || framelen <= polyOrd;

if ~useFallback
    try
        Xs = sgolayfilt(X, polyOrd, framelen);
        Ys = sgolayfilt(Y, polyOrd, framelen);
        [~, G] = sgolay(polyOrd, framelen);
        d1 = (factorial(1) / (Ts^1)) * G(:, 2)';
        vx = conv(Xs, d1, 'same');
        vy = conv(Ys, d1, 'same');
    catch
        vx = gradient(X) * fs;
        vy = gradient(Y) * fs;
    end
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end

function [clusters, tvals, thr] = cluster_perm_1d_local(S, nPerm, alpha)
%CLUSTER_PERM_1D_LOCAL  One-sample cluster-permutation test (vs 0).
%   S: [subjects × time]. Sign-flip permutation for the null distribution.
[ns, nT] = size(S);
se    = std(S, [], 1, 'omitnan') ./ sqrt(sum(isfinite(S), 1));
m     = mean(S, 1, 'omitnan');
tvals = m ./ max(se, eps);

% Critical t (two-sided)
tcrit     = tinv(1 - 0.5 * alpha, ns - 1);
thr.tcrit = tcrit;

% Observed clusters
clusters = find_clusters_local(tvals, tcrit);

% Permutation null (max cluster mass)
maxMass = zeros(1, nPerm);
for p = 1:nPerm
    flips = (rand(ns, 1) > 0.5) * 2 - 1;
    Sp = S .* flips;
    seP = std(Sp, [], 1, 'omitnan') ./ sqrt(sum(isfinite(Sp), 1));
    mP  = mean(Sp, 1, 'omitnan');
    tP  = mP ./ max(seP, eps);
    clP = find_clusters_local(tP, tcrit);
    if isempty(clP)
        maxMass(p) = 0;
    else
        maxMass(p) = max([clP.mass]);
    end
end
maxMass = sort(maxMass);
thr.mass = maxMass(max(1, round((1 - alpha) * nPerm)));

% Cluster p-values
for k = 1:numel(clusters)
    clusters(k).p = 1 - mean(maxMass <= clusters(k).mass);
end
end

function clusters = find_clusters_local(tvals, tcrit)
%FIND_CLUSTERS_LOCAL  Find supra-threshold clusters and their mass.
above = abs(tvals) > tcrit;
clusters = struct('idx', {}, 'mass', {}, 'sign', {}, 'p', {});
if ~any(above), return; end
d   = diff([0, above, 0]);
on  = find(d == 1);
off = find(d == -1) - 1;
for c = 1:numel(on)
    idx = on(c):off(c);
    sgn = sign(mean(tvals(idx)));
    clusters(end+1).idx  = idx; %#ok<AGROW>
    clusters(end).mass   = sum(abs(tvals(idx)));
    clusters(end).sign   = sgn;
    clusters(end).p      = NaN;
end
end

function p_fdr = fdr_bh_local(p)
%FDR_BH_LOCAL  Benjamini-Hochberg FDR correction.
%   Returns adjusted p-values.
p_fdr = nan(size(p));
valid = isfinite(p);
if ~any(valid), return; end

pv = p(valid);
m  = numel(pv);
[ps, si] = sort(pv);
ranks = (1:m)';
adj = ps .* m ./ ranks;

% Enforce monotonicity (step-up)
adj = min(adj, 1);
for k = m-1:-1:1
    adj(k) = min(adj(k), adj(k+1));
end

% Unsort
p_adj = nan(size(pv));
p_adj(si) = adj;
p_fdr(valid) = p_adj;
end
