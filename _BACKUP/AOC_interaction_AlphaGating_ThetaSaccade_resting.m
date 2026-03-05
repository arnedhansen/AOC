%% AOC Interaction — Alpha Gating of Theta-Saccade Coupling: RESTING STATE
% Tests the hierarchical oscillatory control model using resting-state
% (eyes-open) continuous EEG + eye-tracking data.
%
% Theoretical framework (same as task version):
%   Buschermöhle et al. 2025 — Theta-band coherence links saccades and
%                               oculomotor brain regions
%   Wang et al. (&Hanslmayr) 2026 — Pre-stimulus alpha modulates theta
%                                     entrainment strength and memory
%
% Advantage of resting state:
%   - No task constraints: spontaneous alpha-theta-saccade dynamics
%   - Longer continuous recordings (2.5 min per condition)
%   - Within-recording comparison: fixation-cross vs blank periods
%   - Temporal precedence can be tested (alpha leads theta?)
%
% Analyses:
%   1. Sliding-window alpha-theta power comodulation
%   2. Time-lagged cross-correlation (alpha temporal precedence)
%   3. Saccade-locked theta phase coupling (ITC, Rayleigh test)
%   4. Alpha gating: median split on theta, saccade rate, coherence
%   5. Theta-velocity coherence (mscohere)
%
% Outputs (4 figures):
%   Figure 1: Alpha-Theta Temporal Coupling (3 panels)
%   Figure 2: Theta-Saccade Phase Coupling (4 panels)
%   Figure 3: Alpha Gating of Theta-Saccade System (4 panels)
%   Figure 4: Fixation-Cross vs Blank Comparison (3 panels)
%   results .mat file
%
% Prerequisites:
%   Run AOC_preprocessing_resting.m first to create dataEEG_resting.mat
%   and dataET_resting.mat with segInfo.
%   Run AOC_eeg_fex_resting_FOOOF.m to create fooof_resting.mat
%   (aperiodic-corrected oscillatory power via sliding-window FOOOF).
%
% FOOOF correction:
%   Alpha and theta power now come from the FOOOF oscillatory component
%   (model fit - aperiodic) rather than raw Hilbert envelopes.
%   This removes shared 1/f fluctuations that confound alpha-theta
%   comodulation in broadband power.

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 22;
fs       = 500;  % EEG & ET sampling rate (Hz)

% Frequency bands
thetaBand     = [4 8];       % Hz
alphaFallback = [8 14];      % Hz (when IAF unavailable)

% Sliding window parameters
winLen  = 2;     % window length in seconds
winStep = 1;     % step size in seconds (50% overlap)
winSamp = winLen * fs;
stepSamp = winStep * fs;

% Saccade-locked epoch
epochWin     = [-0.5 0.5];   % ±500 ms around saccade onset
epochSamp    = round(epochWin * fs);
nEpochSamp   = diff(epochSamp) + 1;
epochTimeVec = linspace(epochWin(1), epochWin(2), nEpochSamp);

% Segment labels
segLabels = {'fixcross', 'blank'};
segColors = [colors(1,:); colors(3,:)];  % blue-ish vs red-ish

% Lag range for cross-correlation
maxLagSec = 5;  % ±5 seconds
maxLagWin = round(maxLagSec / winStep);
lagVec    = (-maxLagWin:maxLagWin) * winStep;  % in seconds

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/_tests/AlphaGating_ThetaSaccade_resting';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Define ROI channels from first subject's channel layout
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat', 'powload2_late');
allLabels = powload2_late.label;

% Posterior ROI (alpha): O, I, P channels
roi_posterior = {};
for i = 1:numel(allLabels)
    L = allLabels{i};
    if contains(L, 'O') || contains(L, 'I') || startsWith(L, 'P')
        roi_posterior{end+1} = L; %#ok<SAGROW>
    end
end

% Frontal ROI (theta / oculomotor): F and AF channels
roi_frontal = {};
for i = 1:numel(allLabels)
    L = allLabels{i};
    if startsWith(L, 'F') || startsWith(L, 'AF')
        roi_frontal{end+1} = L; %#ok<SAGROW>
    end
end

fprintf('Posterior ROI (%d ch): %s\n', numel(roi_posterior), strjoin(roi_posterior, ', '));
fprintf('Frontal ROI (%d ch): %s\n', numel(roi_frontal), strjoin(roi_frontal, ', '));

%% ====================================================================
%%                      PRE-ALLOCATE GROUP STORAGE
%% ====================================================================
nSubj = numel(subjects);

% Per-segment, per-subject
r_alpha_theta   = nan(nSubj, 2);   % alpha-theta window-by-window Spearman r [subj x seg]
peak_lag        = nan(nSubj, 2);   % lag (s) of max xcorr [subj x seg]
xcorr_all       = nan(nSubj, 2*maxLagWin+1, 2); % full xcorr curves [subj x lag x seg]

rayleigh_z      = nan(nSubj, 2);   % Rayleigh z for theta phase at saccade onset
rayleigh_p      = nan(nSubj, 2);
mean_ITC        = nan(nSubj, 2);   % mean ITC at saccade onset (±50 ms)

theta_highA     = nan(nSubj, 2);   % theta power in high-alpha windows
theta_lowA      = nan(nSubj, 2);   % theta power in low-alpha windows
ms_highA        = nan(nSubj, 2);   % saccade rate in high-alpha windows
ms_lowA         = nan(nSubj, 2);
coh_highA       = nan(nSubj, 2);   % theta-vel coherence in high-alpha windows
coh_lowA        = nan(nSubj, 2);

% Grand-average saccade-locked traces (time-resolved)
ITC_grand       = nan(nSubj, nEpochSamp, 2);
thetaPow_grand  = nan(nSubj, nEpochSamp, 2);

% Phase at saccade onset (pooled for polar histogram)
phase_pool      = cell(2, 1);  % {seg} → concatenated phases across subjects

% Mediation storage (across subjects, pooled windows)
allWin_alpha    = cell(2,1);
allWin_theta    = cell(2,1);
allWin_msrate   = cell(2,1);
allWin_subj     = cell(2,1);

%% ====================================================================
%%                         SUBJECT LOOP
%% ====================================================================
for s = 1:nSubj
    fprintf('\n  Subject %s (%d/%d)... ', subjects{s}, s, nSubj);

    dpeeg  = fullfile(path, subjects{s}, 'eeg');
    dpgaze = fullfile(path, subjects{s}, 'gaze');

    %% Load preprocessed resting-state data + FOOOF
    try
        eegFile   = fullfile(dpeeg, 'dataEEG_resting.mat');
        etFile    = fullfile(dpgaze, 'dataET_resting.mat');
        fooofFile = fullfile(dpeeg, 'fooof_resting.mat');
        if ~exist(eegFile, 'file') || ~exist(etFile, 'file')
            fprintf('missing resting data, skipping.\n'); continue;
        end
        if ~exist(fooofFile, 'file')
            fprintf('missing FOOOF data, skipping.\n'); continue;
        end
        load(eegFile,   'dataEEG', 'segInfo');
        load(etFile,    'dataET');
        load(fooofFile, 'fooof_resting');
    catch
        fprintf('load error, skipping.\n'); continue;
    end

    %% Load IAF (try Sternberg, then N-back, then fallback)
    IAF = NaN;
    iafCandidates = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'};
    for ic = 1:numel(iafCandidates)
        iafPath = fullfile(dpeeg, iafCandidates{ic});
        if exist(iafPath, 'file')
            tmp = load(iafPath, 'IAF_subj');
            if isfinite(tmp.IAF_subj)
                IAF = tmp.IAF_subj;
                break;
            end
        end
    end
    if isfinite(IAF)
        bandAlpha = [IAF-4, IAF+2];
    else
        bandAlpha = alphaFallback;
    end

    %% Extract continuous time series
    % EEG: average across ROI channels (frontal needed for coherence + theta phase)
    front_idx = ismember(dataEEG.label, roi_frontal);

    if sum(front_idx) == 0
        fprintf('frontal ROI channels not found, skipping.\n');
        clear dataEEG dataET segInfo fooof_resting; continue;
    end

    eeg_front_raw = mean(dataEEG.trial{1}(front_idx, :), 1);  % frontal average
    t_eeg = dataEEG.time{1};
    nSamp = numel(t_eeg);

    % ET: left-eye gaze
    etLabels = dataET.label;
    xIdx = find(strcmp(etLabels, 'L-GAZE-X'), 1);
    yIdx = find(strcmp(etLabels, 'L-GAZE-Y'), 1);
    aIdx = find(strcmp(etLabels, 'L-AREA'), 1);

    if isempty(xIdx) || isempty(yIdx)
        fprintf('gaze channels not found, skipping.\n');
        clear dataEEG dataET segInfo; continue;
    end

    gazeX = double(dataET.trial{1}(xIdx, :));
    gazeY = double(dataET.trial{1}(yIdx, :));
    if ~isempty(aIdx)
        gazeA = double(dataET.trial{1}(aIdx, :));
    else
        gazeA = ones(1, numel(gazeX));
    end

    clear dataEEG dataET  % free memory

    %% Blink removal + interpolation (full recording)
    blinkMask = ~isfinite(gazeA) | (gazeA <= 0) | ~isfinite(gazeX) | ~isfinite(gazeY);
    gazeX(blinkMask) = NaN;
    gazeY(blinkMask) = NaN;
    gazeX = fillmissing(gazeX, 'linear', 'EndValues', 'nearest');
    gazeY = fillmissing(gazeY, 'linear', 'EndValues', 'nearest');
    gazeY = -gazeY;  % invert for Cartesian convention

    %% FOOOF-based oscillatory alpha and theta power (per sliding window)
    fooof_post_idx   = ismember(fooof_resting.label, roi_posterior);
    fooof_front_idx  = ismember(fooof_resting.label, roi_frontal);
    fooof_alpha_fidx = fooof_resting.freq >= bandAlpha(1) & fooof_resting.freq <= bandAlpha(2);
    fooof_theta_fidx = fooof_resting.freq >= thetaBand(1) & fooof_resting.freq <= thetaBand(2);

    if sum(fooof_post_idx) == 0
        fprintf('posterior FOOOF channels not found, skipping.\n');
        clear dataEEG dataET segInfo fooof_resting; continue;
    end

    % Mean oscillatory power across ROI channels and frequency bins per window
    fooof_alpha_ts = squeeze(mean(mean( ...
        fooof_resting.powspctrm(fooof_post_idx, fooof_alpha_fidx, :), 1, 'omitnan'), 2, 'omitnan'));
    fooof_theta_ts = squeeze(mean(mean( ...
        fooof_resting.powspctrm(fooof_front_idx, fooof_theta_fidx, :), 1, 'omitnan'), 2, 'omitnan'));
    fooof_time = fooof_resting.time;
    clear fooof_resting  % free memory

    %% Theta envelope + phase (frontal) — needed for saccade-locked analysis
    [b_t, a_t] = butter(4, thetaBand / (fs/2), 'bandpass');
    front_theta_filt = filtfilt(b_t, a_t, eeg_front_raw);
    theta_analytic   = hilbert(front_theta_filt);
    theta_env        = abs(theta_analytic);
    theta_phase      = angle(theta_analytic);

    %% Saccade detection (full recording)
    [vx_full, vy_full] = compute_velocity_sg(gazeX, gazeY, fs, 3);
    vel_full = hypot(vx_full, vy_full);

    % Outlier removal on velocity
    zv = (vel_full - nanmean(vel_full)) / (nanstd(vel_full) + eps);
    vel_full(abs(zv) > 4) = NaN;
    vel_full = fillmissing(vel_full, 'linear', 'EndValues', 'nearest');

    ms_all = ms_detect_engbert(gazeX, gazeY, fs, 6, 6e-3, 0.02);
    sacc_onsets_all = ms_all.onsets;

    %% ================================================================
    %%                      SEGMENT LOOP (fixcross, blank)
    %% ================================================================
    for seg = 1:2
        segName = segLabels{seg};
        if seg == 1
            segRange = segInfo.fixcross;
        else
            segRange = segInfo.blank;
        end

        % Add 2-second buffer from segment boundaries to avoid transients
        bufferSamp = 2 * fs;
        s1 = segRange(1) + bufferSamp;
        s2 = segRange(2) - bufferSamp;
        if s1 >= s2 || s2 > nSamp
            fprintf('  %s: segment too short, skipping.\n', segName);
            continue;
        end

        %% Sliding window analysis
        winStarts = s1:stepSamp:(s2 - winSamp + 1);
        nWin = numel(winStarts);

        if nWin < 10
            fprintf('  %s: too few windows (%d), skipping.\n', segName, nWin);
            continue;
        end

        win_alpha  = nan(nWin, 1);
        win_theta  = nan(nWin, 1);
        win_msrate = nan(nWin, 1);
        win_coh    = nan(nWin, 1);

        for w = 1:nWin
            idx = winStarts(w):(winStarts(w) + winSamp - 1);

            % FOOOF oscillatory alpha and theta power (aperiodic-corrected)
            winCenterTime = t_eeg(min(winStarts(w) + round(winSamp/2), nSamp));
            [~, fooof_widx] = min(abs(fooof_time - winCenterTime));
            win_alpha(w) = fooof_alpha_ts(fooof_widx);
            win_theta(w) = fooof_theta_ts(fooof_widx);

            % Saccade rate in this window
            sacc_in_win = sacc_onsets_all(sacc_onsets_all >= idx(1) & sacc_onsets_all <= idx(end));
            win_msrate(w) = numel(sacc_in_win) / winLen;

            % Theta-velocity coherence (mscohere)
            eeg_seg = eeg_front_raw(idx);
            vel_seg = vel_full(idx);

            if any(~isfinite(eeg_seg)) || any(~isfinite(vel_seg))
                continue;
            end

            coh_winLen = min(256, floor(numel(eeg_seg) / 2));
            if mod(coh_winLen, 2) == 0, coh_winLen = coh_winLen - 1; end
            coh_winLen = max(coh_winLen, 65);
            noverlap   = round(coh_winLen * 0.5);
            nfft       = max(512, 2^nextpow2(coh_winLen));
            win_h      = hamming(coh_winLen);

            try
                [Cxy, f] = mscohere(eeg_seg, vel_seg, win_h, noverlap, nfft, fs);
                theta_fidx = f >= thetaBand(1) & f <= thetaBand(2);
                if any(theta_fidx)
                    win_coh(w) = mean(Cxy(theta_fidx));
                end
            catch
                % skip
            end
        end

        %% Within-subject alpha-theta correlation
        vld = isfinite(win_alpha) & isfinite(win_theta);
        if sum(vld) >= 10
            r_alpha_theta(s, seg) = corr(win_alpha(vld), win_theta(vld), ...
                'Type', 'Spearman');
        end

        %% Time-lagged cross-correlation
        if sum(vld) >= 20
            a_z = (win_alpha(vld) - mean(win_alpha(vld))) / std(win_alpha(vld));
            t_z = (win_theta(vld) - mean(win_theta(vld))) / std(win_theta(vld));
            [xc, lags_raw] = xcorr(a_z, t_z, maxLagWin, 'coeff');

            % Store: positive lag = alpha leads theta
            xcorr_all(s, :, seg) = xc;
            [~, peakIdx] = max(xc);
            peak_lag(s, seg) = lagVec(peakIdx);
        end

        %% Alpha gating: median split on windows
        vld_gate = isfinite(win_alpha) & isfinite(win_theta) & isfinite(win_msrate);
        if sum(vld_gate) >= 10
            med_a = median(win_alpha(vld_gate));
            hi = vld_gate & (win_alpha >= med_a);
            lo = vld_gate & (win_alpha < med_a);

            theta_highA(s, seg) = mean(win_theta(hi), 'omitnan');
            theta_lowA(s, seg)  = mean(win_theta(lo), 'omitnan');
            ms_highA(s, seg)    = mean(win_msrate(hi), 'omitnan');
            ms_lowA(s, seg)     = mean(win_msrate(lo), 'omitnan');

            % Coherence split
            vld_coh = vld_gate & isfinite(win_coh);
            if sum(vld_coh) >= 10
                hi_c = vld_coh & (win_alpha >= med_a);
                lo_c = vld_coh & (win_alpha < med_a);
                coh_highA(s, seg) = mean(win_coh(hi_c), 'omitnan');
                coh_lowA(s, seg)  = mean(win_coh(lo_c), 'omitnan');
            end
        end

        %% Pool windows for group-level mediation
        allWin_alpha{seg}  = [allWin_alpha{seg};  win_alpha(vld_gate)];
        allWin_theta{seg}  = [allWin_theta{seg};  win_theta(vld_gate)];
        allWin_msrate{seg} = [allWin_msrate{seg}; win_msrate(vld_gate)];
        allWin_subj{seg}   = [allWin_subj{seg};   repmat(str2double(subjects{s}), sum(vld_gate), 1)];

        %% Saccade-locked theta analysis
        sacc_in_seg = sacc_onsets_all(sacc_onsets_all >= s1 & sacc_onsets_all <= s2);

        % Exclude saccades too close to segment edges for epoching
        edgeBuf = abs(epochSamp(1)) + 10;
        sacc_in_seg = sacc_in_seg(sacc_in_seg > (s1 + edgeBuf) & sacc_in_seg < (s2 - edgeBuf));

        if numel(sacc_in_seg) < 10
            continue;
        end

        % Extract theta phase and envelope epochs around each saccade
        nSacc = numel(sacc_in_seg);
        phase_epochs = nan(nEpochSamp, nSacc);
        env_epochs   = nan(nEpochSamp, nSacc);

        for si = 1:nSacc
            eidx = (sacc_in_seg(si) + epochSamp(1)):(sacc_in_seg(si) + epochSamp(2));
            if eidx(1) < 1 || eidx(end) > nSamp, continue; end
            phase_epochs(:, si) = theta_phase(eidx);
            env_epochs(:, si)   = theta_env(eidx);
        end

        % Remove invalid epochs
        validEpochs = ~any(isnan(phase_epochs), 1);
        phase_epochs = phase_epochs(:, validEpochs);
        env_epochs   = env_epochs(:, validEpochs);

        if size(phase_epochs, 2) < 10
            continue;
        end

        % Inter-trial coherence (ITC) — time-resolved
        ITC_subj = abs(mean(exp(1i * phase_epochs), 2));
        ITC_grand(s, :, seg) = ITC_subj';

        % Mean theta power around saccade
        thetaPow_grand(s, :, seg) = mean(env_epochs, 2)';

        % Mean ITC at saccade onset (±50 ms window)
        onsetIdx = epochTimeVec >= -0.05 & epochTimeVec <= 0.05;
        mean_ITC(s, seg) = mean(ITC_subj(onsetIdx), 'omitnan');

        % Phase at saccade onset
        onset_sample = find(epochTimeVec >= 0, 1, 'first');
        phases_at_onset = phase_epochs(onset_sample, :);
        phase_pool{seg} = [phase_pool{seg}; phases_at_onset(:)];

        % Rayleigh test for non-uniformity
        [rayleigh_p(s, seg), rayleigh_z(s, seg)] = rayleigh_test(phases_at_onset);

    end % segment loop

    clear eeg_front_raw theta_env theta_phase ...
          fooof_alpha_ts fooof_theta_ts fooof_time ...
          gazeX gazeY gazeA vel_full segInfo
    fprintf('done.\n');

end % subject loop

fprintf('\n  Finished processing all subjects.\n');

%% ====================================================================
%%           GROUP STATISTICS
%% ====================================================================
fprintf('\n%s\n  GROUP STATISTICS\n%s\n', repmat('=',1,60), repmat('=',1,60));

for seg = 1:2
    fprintf('\n--- %s ---\n', upper(segLabels{seg}));

    % Alpha-Theta correlation
    vr = isfinite(r_alpha_theta(:, seg));
    if sum(vr) >= 3
        z_r = atanh(r_alpha_theta(vr, seg));
        [~, p_r, ~, stats_r] = ttest(z_r);
        mean_r = tanh(mean(z_r));
        fprintf('  Alpha-Theta r: mean = %.3f, t(%d) = %.2f, p = %.4f\n', ...
            mean_r, stats_r.df, stats_r.tstat, p_r);
    end

    % Peak lag
    vl = isfinite(peak_lag(:, seg));
    if sum(vl) >= 3
        [~, p_lag, ~, stats_lag] = ttest(peak_lag(vl, seg));
        fprintf('  Peak lag: mean = %.2f s, t(%d) = %.2f, p = %.4f\n', ...
            mean(peak_lag(vl, seg)), stats_lag.df, stats_lag.tstat, p_lag);
    end

    % Rayleigh test
    vz = isfinite(rayleigh_z(:, seg));
    fprintf('  Rayleigh z: mean = %.2f (n = %d subj with z > 0)\n', ...
        nanmean(rayleigh_z(vz, seg)), sum(rayleigh_p(vz, seg) < 0.05));

    % Alpha gating: theta
    vg = isfinite(theta_highA(:, seg)) & isfinite(theta_lowA(:, seg));
    if sum(vg) >= 3
        [~, p_g, ~, stats_g] = ttest(theta_highA(vg, seg), theta_lowA(vg, seg));
        d_g = mean(theta_highA(vg, seg) - theta_lowA(vg, seg)) / std(theta_highA(vg, seg) - theta_lowA(vg, seg));
        fprintf('  Alpha gating theta: t(%d) = %.2f, p = %.4f, d = %.2f\n', ...
            stats_g.df, stats_g.tstat, p_g, d_g);
    end

    % Alpha gating: saccade rate
    vm = isfinite(ms_highA(:, seg)) & isfinite(ms_lowA(:, seg));
    if sum(vm) >= 3
        [~, p_m, ~, stats_m] = ttest(ms_highA(vm, seg), ms_lowA(vm, seg));
        d_m = mean(ms_highA(vm, seg) - ms_lowA(vm, seg)) / std(ms_highA(vm, seg) - ms_lowA(vm, seg));
        fprintf('  Alpha gating MS rate: t(%d) = %.2f, p = %.4f, d = %.2f\n', ...
            stats_m.df, stats_m.tstat, p_m, d_m);
    end

    % Alpha gating: coherence
    vc = isfinite(coh_highA(:, seg)) & isfinite(coh_lowA(:, seg));
    if sum(vc) >= 3
        [~, p_c, ~, stats_c] = ttest(coh_highA(vc, seg), coh_lowA(vc, seg));
        d_c = mean(coh_highA(vc, seg) - coh_lowA(vc, seg)) / std(coh_highA(vc, seg) - coh_lowA(vc, seg));
        fprintf('  Alpha gating coh: t(%d) = %.2f, p = %.4f, d = %.2f\n', ...
            stats_c.df, stats_c.tstat, p_c, d_c);
    end
end

%% Fixcross vs Blank comparison
fprintf('\n--- FIXCROSS vs BLANK ---\n');
for metric_name = {'r_alpha_theta', 'peak_lag', 'rayleigh_z', 'mean_ITC'}
    mn = metric_name{1};
    vals = eval(mn);
    vb = all(isfinite(vals), 2);
    if sum(vb) >= 3
        [~, p_fb, ~, stats_fb] = ttest(vals(vb,1), vals(vb,2));
        fprintf('  %s: fixcross = %.3f, blank = %.3f, t(%d) = %.2f, p = %.4f\n', ...
            mn, mean(vals(vb,1)), mean(vals(vb,2)), stats_fb.df, stats_fb.tstat, p_fb);
    end
end

%% ====================================================================
%%       MEDIATION: Alpha → Theta → Saccade Rate (pooled windows)
%% ====================================================================
fprintf('\n--- Mediation (pooled across segments) ---\n');

med_alpha_all  = [allWin_alpha{1};  allWin_alpha{2}];
med_theta_all  = [allWin_theta{1};  allWin_theta{2}];
med_msrate_all = [allWin_msrate{1}; allWin_msrate{2}];
med_subj_all   = [allWin_subj{1};   allWin_subj{2}];

beta_a = NaN; beta_b = NaN; beta_c = NaN; beta_cprime = NaN;
SE_a = NaN; SE_b = NaN;
p_a = NaN; p_b = NaN; p_c = NaN; p_cprime = NaN;
indirect = NaN; z_sobel = NaN; p_sobel = NaN;
mediation_label = 'Not computed';

try
    T = table(med_alpha_all, med_theta_all, med_msrate_all, ...
        categorical(med_subj_all), ...
        'VariableNames', {'Alpha', 'Theta', 'MSRate', 'ID'});

    % Z-score within subject
    unique_ids = unique(med_subj_all);
    T.AlphaZ = nan(height(T), 1);
    T.ThetaZ = nan(height(T), 1);
    for si = 1:numel(unique_ids)
        idx = T.ID == categorical(unique_ids(si));
        v = T.Alpha(idx); mu = mean(v,'omitnan'); sd = std(v,'omitnan');
        if sd > 0, T.AlphaZ(idx) = (v - mu)/sd; else, T.AlphaZ(idx) = 0; end
        v = T.Theta(idx); mu = mean(v,'omitnan'); sd = std(v,'omitnan');
        if sd > 0, T.ThetaZ(idx) = (v - mu)/sd; else, T.ThetaZ(idx) = 0; end
    end

    % a-path: Alpha → Theta
    lme_a = fitlme(T, 'ThetaZ ~ AlphaZ + (1 | ID)');
    beta_a = lme_a.Coefficients.Estimate(2);
    SE_a   = lme_a.Coefficients.SE(2);
    p_a    = lme_a.Coefficients.pValue(2);

    % b + c' path: MSRate ~ Alpha + Theta
    lme_bc = fitlme(T, 'MSRate ~ AlphaZ + ThetaZ + (1 | ID)');
    beta_cprime = lme_bc.Coefficients.Estimate(2);
    p_cprime    = lme_bc.Coefficients.pValue(2);
    beta_b      = lme_bc.Coefficients.Estimate(3);
    SE_b        = lme_bc.Coefficients.SE(3);
    p_b         = lme_bc.Coefficients.pValue(3);

    % c-path: total Alpha → MSRate
    lme_c = fitlme(T, 'MSRate ~ AlphaZ + (1 | ID)');
    beta_c = lme_c.Coefficients.Estimate(2);
    p_c    = lme_c.Coefficients.pValue(2);

    % Sobel test
    indirect = beta_a * beta_b;
    SE_indirect = sqrt(beta_a^2 * SE_b^2 + beta_b^2 * SE_a^2);
    z_sobel  = indirect / SE_indirect;
    p_sobel  = 2 * (1 - normcdf(abs(z_sobel)));

    if p_sobel < 0.05 && p_a < 0.05 && p_b < 0.05
        if p_cprime >= 0.05
            mediation_label = 'Full mediation';
        else
            mediation_label = 'Partial mediation';
        end
    else
        mediation_label = 'No significant mediation';
    end

    fprintf('  a-path (alpha->theta):    beta = %.4f, p = %.4f\n', beta_a, p_a);
    fprintf('  b-path (theta->MS|alpha): beta = %.4f, p = %.4f\n', beta_b, p_b);
    fprintf('  c-path (alpha->MS total): beta = %.4f, p = %.4f\n', beta_c, p_c);
    fprintf('  c''-path (alpha->MS|theta): beta = %.4f, p = %.4f\n', beta_cprime, p_cprime);
    fprintf('  Indirect: %.4f, Sobel z = %.2f, p = %.4f\n', indirect, z_sobel, p_sobel);
    fprintf('  %s\n', mediation_label);

catch ME
    fprintf('  Mediation failed: %s\n', ME.message);
end

%% ====================================================================
%%  FIGURE 1: Alpha-Theta Temporal Coupling (3 panels)
%% ====================================================================
close all
figure('Position', [50 50 1800 600]);

%% Panel A: Time-lagged cross-correlation (grand average)
subplot(1,3,1); hold on
for seg = 1:2
    xc_seg = squeeze(xcorr_all(:, :, seg));
    valid_xc = all(isfinite(xc_seg), 2);
    if sum(valid_xc) < 2, continue; end
    m_xc = mean(xc_seg(valid_xc, :), 1);
    se_xc = std(xc_seg(valid_xc, :), [], 1) / sqrt(sum(valid_xc));
    fill([lagVec, fliplr(lagVec)], [m_xc + se_xc, fliplr(m_xc - se_xc)], ...
        segColors(seg,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(lagVec, m_xc, '-', 'Color', segColors(seg,:), 'LineWidth', 2.5);
end
xline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Lag (s)  [\alpha leads \rightarrow]')
ylabel('Cross-correlation (r)')
legend(segLabels, 'Location', 'best', 'FontSize', fontSize-10)
title('Alpha-Theta Lag', 'FontSize', fontSize-4)
set(gca, 'FontSize', fontSize-6)

%% Panel B: Peak lag distribution
subplot(1,3,2); hold on
for seg = 1:2
    vl = isfinite(peak_lag(:, seg));
    if sum(vl) < 2, continue; end
    xj = seg + 0.15 * randn(sum(vl), 1);
    scatter(xj, peak_lag(vl, seg), 60, segColors(seg,:), 'filled', 'MarkerFaceAlpha', 0.6);
    m_lag = mean(peak_lag(vl, seg));
    se_lag = std(peak_lag(vl, seg)) / sqrt(sum(vl));
    plot([seg-0.2 seg+0.2], [m_lag m_lag], '-', 'Color', segColors(seg,:), 'LineWidth', 3);
    plot([seg seg], [m_lag-se_lag m_lag+se_lag], '-', 'Color', segColors(seg,:), 'LineWidth', 2);
end
yline(0, 'k:', 'LineWidth', 1.5);
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('Peak Lag (s)')
title('Temporal Precedence', 'FontSize', fontSize-4)

%% Panel C: Alpha-Theta correlation (Spearman r distribution)
subplot(1,3,3); hold on
for seg = 1:2
    vr = isfinite(r_alpha_theta(:, seg));
    if sum(vr) < 2, continue; end
    xj = seg + 0.15 * randn(sum(vr), 1);
    scatter(xj, r_alpha_theta(vr, seg), 60, segColors(seg,:), 'filled', 'MarkerFaceAlpha', 0.6);
    z_tmp = atanh(r_alpha_theta(vr, seg));
    m_r = tanh(mean(z_tmp));
    se_z = std(z_tmp) / sqrt(numel(z_tmp));
    ci_r = tanh([mean(z_tmp) - 1.96*se_z, mean(z_tmp) + 1.96*se_z]);
    plot([seg seg], ci_r, '-', 'Color', segColors(seg,:), 'LineWidth', 3);
    plot(seg, m_r, 's', 'MarkerSize', 14, 'MarkerFaceColor', segColors(seg,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);
end
yline(0, 'k:', 'LineWidth', 1.5);
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('Spearman r (\alpha_{osc} - \theta_{osc})')
title('\alpha-\theta Oscillatory Comodulation', 'FontSize', fontSize-4)

sgtitle('Resting State: Alpha-Theta Temporal Coupling (FOOOF-corrected)', 'FontSize', fontSize)
saveas(gcf, fullfile(outdir, 'AOC_Resting_1_AlphaThetaCoupling.png'));

%% ====================================================================
%%  FIGURE 2: Theta-Saccade Phase Coupling (4 panels)
%% ====================================================================
close all
figure('Position', [50 50 1600 1200]);

%% Panel A: Peri-saccadic theta power (grand average)
subplot(2,2,1); hold on
for seg = 1:2
    vt = squeeze(~any(isnan(thetaPow_grand(:,:,seg)), 2));
    if sum(vt) < 2, continue; end
    m_pow = mean(thetaPow_grand(vt, :, seg), 1);
    se_pow = std(thetaPow_grand(vt, :, seg), [], 1) / sqrt(sum(vt));
    fill([epochTimeVec, fliplr(epochTimeVec)], [m_pow + se_pow, fliplr(m_pow - se_pow)], ...
        segColors(seg,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(epochTimeVec, m_pow, '-', 'Color', segColors(seg,:), 'LineWidth', 2.5);
end
xline(0, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time from saccade onset (s)')
ylabel('\theta Power (a.u.)')
legend(segLabels, 'Location', 'best', 'FontSize', fontSize-10)
title('Peri-saccadic \theta Power', 'FontSize', fontSize-4)
set(gca, 'FontSize', fontSize-6)

%% Panel B: Theta ITC (time-resolved)
subplot(2,2,2); hold on
for seg = 1:2
    vi = squeeze(~any(isnan(ITC_grand(:,:,seg)), 2));
    if sum(vi) < 2, continue; end
    m_itc = mean(ITC_grand(vi, :, seg), 1);
    se_itc = std(ITC_grand(vi, :, seg), [], 1) / sqrt(sum(vi));
    fill([epochTimeVec, fliplr(epochTimeVec)], [m_itc + se_itc, fliplr(m_itc - se_itc)], ...
        segColors(seg,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(epochTimeVec, m_itc, '-', 'Color', segColors(seg,:), 'LineWidth', 2.5);
end
xline(0, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time from saccade onset (s)')
ylabel('Inter-Trial Coherence')
legend(segLabels, 'Location', 'best', 'FontSize', fontSize-10)
title('Saccade-Locked \theta ITC', 'FontSize', fontSize-4)
set(gca, 'FontSize', fontSize-6)

%% Panel C: Polar histogram of theta phase at saccade onset
ax_temp = subplot(2,2,3);
pos_polar = ax_temp.Position;
delete(ax_temp);
h_pol = polaraxes('Position', pos_polar);
hold(h_pol, 'on');
for seg = 1:2
    if isempty(phase_pool{seg}), continue; end
    polarhistogram(h_pol, phase_pool{seg}, 36, ...
        'FaceColor', segColors(seg,:), 'FaceAlpha', 0.4, ...
        'EdgeColor', segColors(seg,:), 'Normalization', 'probability');
end
title(h_pol, '\theta Phase at Saccade Onset', 'FontSize', fontSize-4)
set(h_pol, 'FontSize', fontSize-8)

%% Panel D: Rayleigh z distribution per subject
subplot(2,2,4); hold on
for seg = 1:2
    vz = isfinite(rayleigh_z(:, seg));
    if sum(vz) < 2, continue; end
    xj = seg + 0.15 * randn(sum(vz), 1);
    scatter(xj, rayleigh_z(vz, seg), 60, segColors(seg,:), 'filled', 'MarkerFaceAlpha', 0.6);
    m_z = mean(rayleigh_z(vz, seg));
    se_z = std(rayleigh_z(vz, seg)) / sqrt(sum(vz));
    plot([seg-0.2 seg+0.2], [m_z m_z], '-', 'Color', segColors(seg,:), 'LineWidth', 3);
    plot([seg seg], [m_z-se_z m_z+se_z], '-', 'Color', segColors(seg,:), 'LineWidth', 2);
end
yline(3, 'k--', 'p \approx .05', 'LineWidth', 1.5, 'FontSize', fontSize-10);
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('Rayleigh z')
title('\theta Phase Clustering', 'FontSize', fontSize-4)

sgtitle('Resting State: Theta-Saccade Phase Coupling', 'FontSize', fontSize)
saveas(gcf, fullfile(outdir, 'AOC_Resting_2_ThetaSaccadeCoupling.png'));

%% ====================================================================
%%  FIGURE 3: Alpha Gating of Theta-Saccade System (4 panels)
%% ====================================================================
close all
figure('Position', [50 50 1600 1200]);

%% Panel A: Alpha median split → Theta power
subplot(2,2,1); hold on
x_pos = [1 2; 4 5];  % [fixcross_lo fixcross_hi; blank_lo blank_hi]
for seg = 1:2
    vg = isfinite(theta_highA(:, seg)) & isfinite(theta_lowA(:, seg));
    if sum(vg) < 2, continue; end
    data_bar = [mean(theta_lowA(vg, seg),'omitnan'), mean(theta_highA(vg, seg),'omitnan')];
    sem_bar  = [std(theta_lowA(vg, seg),'omitnan'), std(theta_highA(vg, seg),'omitnan')] / sqrt(sum(vg));
    b = bar(x_pos(seg,:), data_bar, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(3,:);
    errorbar(x_pos(seg,:), data_bar, sem_bar, 'k.', 'LineWidth', 2, 'CapSize', 10)
    for si = find(vg)'
        plot(x_pos(seg,:), [theta_lowA(si,seg), theta_highA(si,seg)], '-', 'Color', [.6 .6 .6 .3]);
    end
end
set(gca, 'XTick', [1.5 4.5], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('Oscillatory \theta Power (FOOOF)')
set(get(gca, 'Children'), 'HandleVisibility', 'off')
h_lo = plot(nan, nan, 's', 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
h_hi = plot(nan, nan, 's', 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
legend([h_lo, h_hi], {'Low \alpha', 'High \alpha'}, 'Location', 'best', 'FontSize', fontSize-10)
title('\alpha Gating \rightarrow \theta Power', 'FontSize', fontSize-4)

%% Panel B: Alpha median split → Saccade rate
subplot(2,2,2); hold on
for seg = 1:2
    vm = isfinite(ms_highA(:, seg)) & isfinite(ms_lowA(:, seg));
    if sum(vm) < 2, continue; end
    data_bar = [mean(ms_lowA(vm, seg),'omitnan'), mean(ms_highA(vm, seg),'omitnan')];
    sem_bar  = [std(ms_lowA(vm, seg),'omitnan'), std(ms_highA(vm, seg),'omitnan')] / sqrt(sum(vm));
    b = bar(x_pos(seg,:), data_bar, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(3,:);
    errorbar(x_pos(seg,:), data_bar, sem_bar, 'k.', 'LineWidth', 2, 'CapSize', 10)
    for si = find(vm)'
        plot(x_pos(seg,:), [ms_lowA(si,seg), ms_highA(si,seg)], '-', 'Color', [.6 .6 .6 .3]);
    end
end
set(gca, 'XTick', [1.5 4.5], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('Saccade Rate [/s]')
title('\alpha Gating \rightarrow Saccade Rate', 'FontSize', fontSize-4)

%% Panel C: Alpha median split → Theta-velocity coherence
subplot(2,2,3); hold on
for seg = 1:2
    vc = isfinite(coh_highA(:, seg)) & isfinite(coh_lowA(:, seg));
    if sum(vc) < 2, continue; end
    data_bar = [mean(coh_lowA(vc, seg),'omitnan'), mean(coh_highA(vc, seg),'omitnan')];
    sem_bar  = [std(coh_lowA(vc, seg),'omitnan'), std(coh_highA(vc, seg),'omitnan')] / sqrt(sum(vc));
    b = bar(x_pos(seg,:), data_bar, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(3,:);
    errorbar(x_pos(seg,:), data_bar, sem_bar, 'k.', 'LineWidth', 2, 'CapSize', 10)
    for si = find(vc)'
        plot(x_pos(seg,:), [coh_lowA(si,seg), coh_highA(si,seg)], '-', 'Color', [.6 .6 .6 .3]);
    end
end
set(gca, 'XTick', [1.5 4.5], 'XTickLabel', segLabels, 'FontSize', fontSize-6)
ylabel('\theta-Velocity Coherence')
title('\alpha Gating \rightarrow \theta-Saccade Coh.', 'FontSize', fontSize-4)

%% Panel D: Mediation path model summary
subplot(2,2,4); axis off
xlim([0 1]); ylim([0 1]);

path_str = {
    '\bf Hierarchical Mediation Model (Resting, FOOOF) \rm'
    ''
    sprintf('a-path (\\alpha \\rightarrow \\theta):      \\beta = %.3f, p = %.4f', beta_a, p_a)
    sprintf('b-path (\\theta \\rightarrow MS|\\alpha):  \\beta = %.3f, p = %.4f', beta_b, p_b)
    ''
    sprintf('c-path (\\alpha \\rightarrow MS total):   \\beta = %.3f, p = %.4f', beta_c, p_c)
    sprintf('c''-path (\\alpha \\rightarrow MS direct):  \\beta = %.3f, p = %.4f', beta_cprime, p_cprime)
    ''
    sprintf('Indirect effect: %.4f', indirect)
    sprintf('Sobel z = %.2f, p = %.4f', z_sobel, p_sobel)
    ''
    sprintf('\\bf %s \\rm', mediation_label)
};
text(0.05, 0.95, path_str, 'FontSize', fontSize-8, 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex', 'FontName', 'FixedWidth');

sgtitle('Resting State: Alpha Gating of Theta-Saccade System (FOOOF-corrected)', 'FontSize', fontSize)
saveas(gcf, fullfile(outdir, 'AOC_Resting_3_AlphaGating.png'));

%% ====================================================================
%%  FIGURE 4: Fixation-Cross vs Blank Period Comparison (3 panels)
%% ====================================================================
close all
figure('Position', [50 50 1800 600]);

%% Panel A: Paired comparison of alpha-theta correlations
subplot(1,3,1); hold on
vb = all(isfinite(r_alpha_theta), 2);
if sum(vb) >= 2
    for si = find(vb)'
        plot([1 2], r_alpha_theta(si, :), '-', 'Color', [.6 .6 .6 .4], 'LineWidth', 1);
    end
    for seg = 1:2
        vals = r_alpha_theta(vb, seg);
        m_val = mean(vals);
        se_val = std(vals) / sqrt(numel(vals));
        errorbar(seg, m_val, se_val, 'o', 'Color', segColors(seg,:), ...
            'MarkerFaceColor', segColors(seg,:), 'MarkerSize', 12, ...
            'LineWidth', 2.5, 'CapSize', 12);
    end
    [~, p_comp] = ttest(r_alpha_theta(vb,1), r_alpha_theta(vb,2));
    title(sprintf('\\alpha-\\theta Correlation\np = %.4f', p_comp), 'FontSize', fontSize-4)
end
yline(0, 'k:', 'LineWidth', 1.5);
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'XLim', [0.5 2.5], 'FontSize', fontSize-6)
ylabel('Spearman r')

%% Panel B: Paired comparison of theta ITC at saccade onset
subplot(1,3,2); hold on
vi = all(isfinite(mean_ITC), 2);
if sum(vi) >= 2
    for si = find(vi)'
        plot([1 2], mean_ITC(si, :), '-', 'Color', [.6 .6 .6 .4], 'LineWidth', 1);
    end
    for seg = 1:2
        vals = mean_ITC(vi, seg);
        m_val = mean(vals);
        se_val = std(vals) / sqrt(numel(vals));
        errorbar(seg, m_val, se_val, 'o', 'Color', segColors(seg,:), ...
            'MarkerFaceColor', segColors(seg,:), 'MarkerSize', 12, ...
            'LineWidth', 2.5, 'CapSize', 12);
    end
    [~, p_itc] = ttest(mean_ITC(vi,1), mean_ITC(vi,2));
    title(sprintf('\\theta ITC at Saccade\np = %.4f', p_itc), 'FontSize', fontSize-4)
end
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'XLim', [0.5 2.5], 'FontSize', fontSize-6)
ylabel('ITC (at onset)')

%% Panel C: Paired comparison of Rayleigh z
subplot(1,3,3); hold on
vz = all(isfinite(rayleigh_z), 2);
if sum(vz) >= 2
    for si = find(vz)'
        plot([1 2], rayleigh_z(si, :), '-', 'Color', [.6 .6 .6 .4], 'LineWidth', 1);
    end
    for seg = 1:2
        vals = rayleigh_z(vz, seg);
        m_val = mean(vals);
        se_val = std(vals) / sqrt(numel(vals));
        errorbar(seg, m_val, se_val, 'o', 'Color', segColors(seg,:), ...
            'MarkerFaceColor', segColors(seg,:), 'MarkerSize', 12, ...
            'LineWidth', 2.5, 'CapSize', 12);
    end
    [~, p_rz] = ttest(rayleigh_z(vz,1), rayleigh_z(vz,2));
    title(sprintf('\\theta Phase Clustering\np = %.4f', p_rz), 'FontSize', fontSize-4)
end
yline(3, 'k--', 'p \approx .05', 'LineWidth', 1.5, 'FontSize', fontSize-10);
set(gca, 'XTick', [1 2], 'XTickLabel', segLabels, 'XLim', [0.5 2.5], 'FontSize', fontSize-6)
ylabel('Rayleigh z')

sgtitle('Resting State: Fixation-Cross vs Blank Period', 'FontSize', fontSize)
saveas(gcf, fullfile(outdir, 'AOC_Resting_4_FixcrossVsBlank.png'));

%% ====================================================================
%%                          SAVE RESULTS
%% ====================================================================
results = struct();
results.subjects        = subjects;
results.segLabels       = segLabels;
results.r_alpha_theta   = r_alpha_theta;
results.peak_lag        = peak_lag;
results.xcorr_all       = xcorr_all;
results.lagVec          = lagVec;
results.rayleigh_z      = rayleigh_z;
results.rayleigh_p      = rayleigh_p;
results.mean_ITC        = mean_ITC;
results.ITC_grand       = ITC_grand;
results.thetaPow_grand  = thetaPow_grand;
results.epochTimeVec    = epochTimeVec;
results.phase_pool      = phase_pool;
results.theta_highA     = theta_highA;
results.theta_lowA      = theta_lowA;
results.ms_highA        = ms_highA;
results.ms_lowA         = ms_lowA;
results.coh_highA       = coh_highA;
results.coh_lowA        = coh_lowA;
results.mediation.beta_a      = beta_a;
results.mediation.beta_b      = beta_b;
results.mediation.beta_c      = beta_c;
results.mediation.beta_cprime = beta_cprime;
results.mediation.indirect    = indirect;
results.mediation.z_sobel     = z_sobel;
results.mediation.p_sobel     = p_sobel;
results.mediation.label       = mediation_label;

save(fullfile(outdir, 'AOC_Resting_AlphaGating_ThetaSaccade_results.mat'), 'results');
fprintf('\nDone! Results and figures saved to:\n  %s\n', outdir);

%% ====================================================================
%%                         LOCAL FUNCTIONS
%% ====================================================================

function ms = ms_detect_engbert(X, Y, fs, lambda, minDurSec, minISISec)
%MS_DETECT_ENGBERT  Detect microsaccades using Engbert & Kliegl (2003).
if nargin < 4 || isempty(lambda),    lambda = 6;        end
if nargin < 5 || isempty(minDurSec), minDurSec = 0.006; end
if nargin < 6 || isempty(minISISec), minISISec = 0.02;  end

[vx, vy] = compute_velocity_sg(X, Y, fs, 3);

sx = 1.4826 * median(abs(vx - median(vx, 'omitnan')), 'omitnan');
sy = 1.4826 * median(abs(vy - median(vy, 'omitnan')), 'omitnan');
tx = lambda * max(sx, eps);
ty = lambda * max(sy, eps);

rad2 = (vx ./ tx).^2 + (vy ./ ty).^2;
cand = isfinite(rad2) & (rad2 > 1);

d = diff([0; cand(:); 0]);
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
            onM  = [onM; O]; offM = [offM; F]; %#ok<AGROW>
            O = on(i); F = off(i);
        end
    end
    onM  = [onM; O]; offM = [offM; F];
    on = onM; off = offM;
end

ms.onsets  = on;
ms.offsets = off;
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
%COMPUTE_VELOCITY_SG  Savitzky-Golay velocity estimation.
X = double(X); Y = double(Y);
Ts = 1/fs;
L  = numel(X);
framelen = min(21, L);
if mod(framelen, 2) == 0, framelen = framelen - 1; end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0, minLegal = minLegal + 1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L, framelen = L - mod(L, 2) + 1; end
useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / (Ts^1)) * G(:, 2)';
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end

function [p, z] = rayleigh_test(theta)
%RAYLEIGH_TEST  Test for non-uniformity of circular data.
%   Uses Rayleigh's R statistic. Small p → non-uniform distribution.
n = numel(theta);
if n < 2
    p = 1; z = 0; return;
end
R = abs(mean(exp(1i * theta(:))));
z = n * R^2;
% Approximation from Mardia & Jupp (2000)
p = exp(-z) * (1 + (2*z - z^2) / (4*n) - (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n^2));
p = max(0, min(1, p));
end
