%% AOC Interaction — Alpha Gating of Theta-Saccade Coupling
% Tests the hierarchical oscillatory control model:
%   Alpha (state/gating) → Theta (sampling/coordination) → Saccades + Memory
%
% Theoretical framework:
%   Buschermöhle et al. 2025 — Theta-band coherence links saccades and
%                               oculomotor brain regions
%   Wang et al. (&Hanslmayr) 2026 — Pre-stimulus alpha modulates theta
%                                     entrainment strength and memory
%
% Predictions tested:
%   1. Low pre-stimulus alpha → stronger retention-period theta power
%      (alpha release enables theta entrainment)
%   2. Higher retention theta → more (micro)saccades
%      (theta organises rhythmic visuo-motor sampling)
%   3. Pre-stimulus alpha gates theta-saccade coupling
%      (low alpha → stronger theta-velocity coherence)
%   4. Mediation: Alpha → Theta → Saccade rate
%      (theta mediates the alpha-saccade relationship)
%
% Data: Sternberg WM (loads 2/4/6), N-back (1/2/3-back)
%       Combined EEG + eye-tracking, trial-level analysis
%
% Methods:
%   Pre-stimulus alpha: posterior ROI, bandpass [IAF-4, IAF+2] Hz, Hilbert
%   Retention theta: frontal ROI, bandpass [4 8] Hz, Hilbert
%   Microsaccade rate: Engbert & Kliegl (2003) algorithm
%   Theta-velocity coherence: mscohere, theta-band (4-8 Hz)
%   Statistics: within-subject Spearman, median split, GLMMs, Sobel mediation
%
% Outputs per task (Sternberg + N-back):
%   Figure 1: Pre-stimulus Alpha → Retention Theta (3 panels)
%   Figure 2: Retention Theta → Saccade Rate (3 panels)
%   Figure 3: Hierarchical Model Summary (4 panels)
%   Figure 4: Condition-Resolved Effects (3 panels)
%   results .mat file

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 22;
fs       = 500;  % EEG & ET sampling rate (Hz)

% Frequency bands
thetaBand     = [4 8];       % Hz
alphaFallback = [8 14];      % Hz (when IAF unavailable)

% Time windows
preStimWin = [-0.5 -0.25];   % pre-stimulus alpha window

% Task definitions
tasks      = {'sternberg', 'nback'};
eegFiles   = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback'};
etFiles    = {'dataET_sternberg', 'dataET_nback'};
iafFiles   = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'};
retWins    = {[1 2], [0 2]};   % Sternberg: late retention; N-back: full stimulus
condLabels = {{'Load 2', 'Load 4', 'Load 6'}, {'1-back', '2-back', '3-back'}};

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AlphaGating_ThetaSaccade';
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

%% Results container
results = struct();

%% ====================================================================
%%                         TASK LOOP
%% ====================================================================
for taskIdx = 1:2
    taskName = tasks{taskIdx};
    retWin   = retWins{taskIdx};

    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName), repmat('=',1,60));

    % Trial-level storage (across all subjects)
    all_preAlpha = [];    % pre-stimulus alpha power (posterior)
    all_retTheta = [];    % retention theta power (frontal)
    all_retAlpha = [];    % retention alpha power (posterior, control)
    all_msrate   = [];    % microsaccade rate in retention window
    all_thetaCoh = [];    % theta-velocity coherence per trial
    all_subj     = [];
    all_cond     = [];

    % Subject-level storage for median-split coherence
    coh_highAlpha_subj = nan(numel(subjects), 1);
    coh_lowAlpha_subj  = nan(numel(subjects), 1);

    %% ================================================================
    %%                       SUBJECT LOOP
    %% ================================================================
    for s = 1:numel(subjects)
        fprintf('  %s | Subject %s (%d/%d)... ', upper(taskName), subjects{s}, s, numel(subjects));

        dpeeg  = fullfile(path, subjects{s}, 'eeg');
        dpgaze = fullfile(path, subjects{s}, 'gaze');

        %% Load raw epoched EEG
        try
            cd(dpeeg)
            load(eegFiles{taskIdx})  % → dataTFR
        catch
            fprintf('missing EEG, skipping.\n'); continue;
        end

        %% Load IAF
        IAF_file = fullfile(dpeeg, iafFiles{taskIdx});
        if exist(IAF_file, 'file')
            tmp = load(IAF_file, 'IAF_subj');
            IAF = tmp.IAF_subj;
        else
            IAF = NaN;
        end
        if isfinite(IAF)
            bandAlpha = [IAF-4, IAF+2];
        else
            bandAlpha = alphaFallback;
        end

        %% Posterior alpha envelope (Hilbert)
        cfg = []; cfg.channel = roi_posterior;
        postData = ft_selectdata(cfg, dataTFR);
        cfg = []; cfg.avgoverchan = 'yes';
        postAvg = ft_selectdata(cfg, postData); clear postData

        cfg            = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = bandAlpha;
        cfg.demean     = 'yes';
        cfg.hilbert    = 'abs';
        alphaEnv = ft_preprocessing(cfg, postAvg); clear postAvg

        %% Frontal data: theta envelope + broadband (for coherence)
        cfg = []; cfg.channel = roi_frontal;
        frontData = ft_selectdata(cfg, dataTFR);
        cfg = []; cfg.avgoverchan = 'yes';
        frontAvg = ft_selectdata(cfg, frontData); clear frontData

        cfg            = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = thetaBand;
        cfg.demean     = 'yes';
        cfg.hilbert    = 'abs';
        thetaEnv = ft_preprocessing(cfg, frontAvg);
        % Keep frontAvg for theta-velocity coherence

        %% Load gaze
        try
            cd(dpgaze)
            load(etFiles{taskIdx}, 'dataETlong')
        catch
            fprintf('missing gaze, skipping.\n');
            clear dataTFR frontAvg alphaEnv thetaEnv; continue;
        end

        %% Match trial counts
        nT = min([numel(alphaEnv.trial), numel(thetaEnv.trial), ...
                  numel(frontAvg.trial),  size(dataETlong.trialinfo, 1)]);

        if nT == 0
            fprintf('no trials, skipping.\n');
            clear dataTFR frontAvg alphaEnv thetaEnv dataETlong; continue;
        end

        %% Pre-allocate trial-level vectors
        preAlpha_tr = nan(nT, 1);
        retTheta_tr = nan(nT, 1);
        retAlpha_tr = nan(nT, 1);
        msrate_tr   = nan(nT, 1);
        thetaCoh_tr = nan(nT, 1);
        condCode_tr = dataTFR.trialinfo(1:nT, 1);

        %% ============================================================
        %%                       TRIAL LOOP
        %% ============================================================
        for tr = 1:nT
            t_eeg  = alphaEnv.time{tr};
            t_gaze = dataETlong.time{tr};

            %% Pre-stimulus alpha (posterior envelope)
            preIdx = t_eeg >= preStimWin(1) & t_eeg <= preStimWin(2);
            ae = double(alphaEnv.trial{tr}(1, :));
            ae = fillmissing(ae, 'linear', 'EndValues', 'nearest');
            preAlpha_tr(tr) = mean(ae(preIdx), 'omitnan');

            %% Retention theta (frontal envelope)
            retIdx_eeg = t_eeg >= retWin(1) & t_eeg <= retWin(2);
            te = double(thetaEnv.trial{tr}(1, :));
            te = fillmissing(te, 'linear', 'EndValues', 'nearest');
            retTheta_tr(tr) = mean(te(retIdx_eeg), 'omitnan');

            %% Retention alpha (posterior envelope, control measure)
            retAlpha_tr(tr) = mean(ae(retIdx_eeg), 'omitnan');

            %% Gaze: retention window — microsaccade detection
            retIdx_gaze = t_gaze >= retWin(1) & t_gaze <= retWin(2);
            X = double(dataETlong.trial{tr}(1, retIdx_gaze));
            Y = double(dataETlong.trial{tr}(2, retIdx_gaze));

            if numel(X) < 50, continue; end

            % Blink removal + interpolation
            if size(dataETlong.trial{tr}, 1) >= 3
                A = double(dataETlong.trial{tr}(3, retIdx_gaze));
                blink = ~isfinite(A) | (A <= 0);
            else
                blink = ~isfinite(X) | ~isfinite(Y);
            end
            if any(blink)
                X(blink) = NaN; Y(blink) = NaN;
                if isnan(X(1)),   i1 = find(isfinite(X), 1, 'first'); if ~isempty(i1), X(1:i1-1) = X(i1); end, end
                if isnan(X(end)), i2 = find(isfinite(X), 1, 'last');  if ~isempty(i2), X(i2+1:end) = X(i2); end, end
                if isnan(Y(1)),   i1 = find(isfinite(Y), 1, 'first'); if ~isempty(i1), Y(1:i1-1) = Y(i1); end, end
                if isnan(Y(end)), i2 = find(isfinite(Y), 1, 'last');  if ~isempty(i2), Y(i2+1:end) = Y(i2); end, end
                X = fillmissing(X, 'linear');
                Y = fillmissing(Y, 'linear');
            end
            Y = -Y;  % invert for Cartesian convention

            % Microsaccade detection (Engbert & Kliegl)
            ms = ms_detect_engbert(X, Y, fs, 6, 6e-3, 0.02);
            T_dur = numel(X) / fs;
            msrate_tr(tr) = numel(ms.onsets) / T_dur;

            %% Theta-velocity coherence (mscohere, theta band)
            % Use broadband frontal EEG and eye velocity
            eeg_broad = double(frontAvg.trial{tr}(1, :));
            eeg_broad = fillmissing(eeg_broad, 'linear', 'EndValues', 'nearest');

            % Full-trial gaze for velocity computation
            X_full = double(dataETlong.trial{tr}(1, :));
            Y_full = double(dataETlong.trial{tr}(2, :));
            if size(dataETlong.trial{tr}, 1) >= 3
                A_full = double(dataETlong.trial{tr}(3, :));
                blink_full = ~isfinite(A_full) | (A_full <= 0);
            else
                blink_full = ~isfinite(X_full) | ~isfinite(Y_full);
            end
            if any(blink_full)
                X_full(blink_full) = NaN; Y_full(blink_full) = NaN;
                X_full = fillmissing(X_full, 'linear', 'EndValues', 'nearest');
                Y_full = fillmissing(Y_full, 'linear', 'EndValues', 'nearest');
            end

            [vx, vy] = compute_velocity_sg(X_full, -Y_full, fs, 3);
            vel_full = hypot(vx, vy);

            % Handle velocity outliers
            zv = (vel_full - nanmean(vel_full)) / (nanstd(vel_full) + eps);
            vel_full(abs(zv) > 4) = NaN;
            vel_full = fillmissing(vel_full, 'linear', 'EndValues', 'nearest');

            % Align gaze velocity to EEG time grid
            if numel(vel_full) ~= numel(eeg_broad)
                gaze_time_full = linspace(t_gaze(1), t_gaze(end), numel(vel_full));
                vel_full = interp1(gaze_time_full, vel_full, t_eeg, 'linear', 'extrap');
            end

            eeg_ret = eeg_broad(retIdx_eeg);
            vel_ret = vel_full(retIdx_eeg);

            if numel(eeg_ret) < 200 || any(~isfinite(eeg_ret)) || any(~isfinite(vel_ret))
                continue;
            end

            % mscohere: theta-band coherence
            coh_winLen = min(256, floor(numel(eeg_ret) / 2));
            if mod(coh_winLen, 2) == 0, coh_winLen = coh_winLen - 1; end
            coh_winLen = max(coh_winLen, 65);
            noverlap   = round(coh_winLen * 0.5);
            nfft       = max(512, 2^nextpow2(coh_winLen));
            win        = hamming(coh_winLen);

            try
                [Cxy, f] = mscohere(eeg_ret, vel_ret, win, noverlap, nfft, fs);
                theta_fidx = f >= thetaBand(1) & f <= thetaBand(2);
                if any(theta_fidx)
                    coh_val = mean(Cxy(theta_fidx));
                    thetaCoh_tr(tr) = max(0, min(0.999, coh_val));
                end
            catch
                % skip
            end
        end % trial loop

        %% Subject-level median split (coherence × alpha)
        valid_coh = isfinite(preAlpha_tr) & isfinite(thetaCoh_tr);
        if sum(valid_coh) >= 10
            med_alpha = nanmedian(preAlpha_tr(valid_coh));
            hi = valid_coh & (preAlpha_tr >= med_alpha);
            lo = valid_coh & (preAlpha_tr < med_alpha);
            coh_highAlpha_subj(s) = mean(thetaCoh_tr(hi), 'omitnan');
            coh_lowAlpha_subj(s)  = mean(thetaCoh_tr(lo), 'omitnan');
        end

        %% Append to grand storage
        valid = isfinite(preAlpha_tr) & isfinite(retTheta_tr) & isfinite(msrate_tr);
        all_preAlpha = [all_preAlpha; preAlpha_tr(valid)];   %#ok<AGROW>
        all_retTheta = [all_retTheta; retTheta_tr(valid)];   %#ok<AGROW>
        all_retAlpha = [all_retAlpha; retAlpha_tr(valid)];   %#ok<AGROW>
        all_msrate   = [all_msrate;   msrate_tr(valid)];     %#ok<AGROW>
        all_thetaCoh = [all_thetaCoh; thetaCoh_tr(valid)];   %#ok<AGROW>
        all_subj     = [all_subj;     repmat(str2double(subjects{s}), sum(valid), 1)]; %#ok<AGROW>
        all_cond     = [all_cond;     condCode_tr(valid)];   %#ok<AGROW>

        %% Clean up
        clear dataTFR frontAvg alphaEnv thetaEnv dataETlong
        fprintf('done (%d valid trials).\n', sum(valid));

    end % subject loop

    %% ================================================================
    %%  Convert condition codes
    %% ================================================================
    all_cond = all_cond - 20;  % Sternberg: 22/24/26 → 2/4/6; N-back: 21/22/23 → 1/2/3

    unique_ids   = unique(all_subj);
    nS           = numel(unique_ids);
    unique_conds = sort(unique(all_cond));

    fprintf('\n  Total: %d subjects, %d trials\n', nS, numel(all_preAlpha));

    %% ================================================================
    %%       ANALYSIS 1: Pre-stimulus Alpha → Retention Theta
    %% ================================================================
    fprintf('\n--- Analysis 1: Pre-stim Alpha → Retention Theta (%s) ---\n', taskName);

    r_alpha_theta = nan(nS, 1);
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        a = all_preAlpha(idx);
        t = all_retTheta(idx);
        if sum(isfinite(a) & isfinite(t)) >= 5
            r_alpha_theta(si) = corr(a, t, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    valid_r1 = isfinite(r_alpha_theta);
    z_at = atanh(r_alpha_theta(valid_r1));
    [~, p_at, ~, stats_at] = ttest(z_at);
    mean_r_at = tanh(mean(z_at));
    se_z_at   = std(z_at) / sqrt(numel(z_at));
    ci_r_at   = tanh([mean(z_at) - 1.96*se_z_at, mean(z_at) + 1.96*se_z_at]);

    fprintf('  Alpha→Theta r: mean = %.3f [%.3f, %.3f], t(%d) = %.2f, p = %.4f\n', ...
        mean_r_at, ci_r_at(1), ci_r_at(2), stats_at.df, stats_at.tstat, p_at);

    %% ================================================================
    %%       ANALYSIS 2: Retention Theta → Saccade Rate
    %% ================================================================
    fprintf('\n--- Analysis 2: Retention Theta → Saccade Rate (%s) ---\n', taskName);

    r_theta_ms = nan(nS, 1);
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        t = all_retTheta(idx);
        m = all_msrate(idx);
        if sum(isfinite(t) & isfinite(m)) >= 5
            r_theta_ms(si) = corr(t, m, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    valid_r2 = isfinite(r_theta_ms);
    z_tm = atanh(r_theta_ms(valid_r2));
    [~, p_tm, ~, stats_tm] = ttest(z_tm);
    mean_r_tm = tanh(mean(z_tm));
    se_z_tm   = std(z_tm) / sqrt(numel(z_tm));
    ci_r_tm   = tanh([mean(z_tm) - 1.96*se_z_tm, mean(z_tm) + 1.96*se_z_tm]);

    fprintf('  Theta→MS r: mean = %.3f [%.3f, %.3f], t(%d) = %.2f, p = %.4f\n', ...
        mean_r_tm, ci_r_tm(1), ci_r_tm(2), stats_tm.df, stats_tm.tstat, p_tm);

    %% Direct: Alpha → MS (for comparison / mediation)
    r_alpha_ms = nan(nS, 1);
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        a = all_preAlpha(idx);
        m = all_msrate(idx);
        if sum(isfinite(a) & isfinite(m)) >= 5
            r_alpha_ms(si) = corr(a, m, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end
    valid_r3 = isfinite(r_alpha_ms);
    mean_r_am = tanh(mean(atanh(r_alpha_ms(valid_r3))));

    %% ================================================================
    %%       ANALYSIS 3: Alpha Gating — Median Split
    %% ================================================================
    fprintf('\n--- Analysis 3: Alpha Gating of Theta & Saccades (%s) ---\n', taskName);

    theta_high = nan(nS, 1); theta_low = nan(nS, 1);
    ms_high    = nan(nS, 1); ms_low    = nan(nS, 1);

    for si = 1:nS
        idx = all_subj == unique_ids(si);
        a = all_preAlpha(idx);
        t = all_retTheta(idx);
        m = all_msrate(idx);

        vld = isfinite(a) & isfinite(t) & isfinite(m);
        a = a(vld); t = t(vld); m = m(vld);
        if numel(a) < 10, continue; end

        med_a = median(a);
        hi = a >= med_a; lo = a < med_a;

        theta_high(si) = mean(t(hi), 'omitnan');
        theta_low(si)  = mean(t(lo), 'omitnan');
        ms_high(si)    = mean(m(hi), 'omitnan');
        ms_low(si)     = mean(m(lo), 'omitnan');
    end

    % Paired t-tests
    vs = isfinite(theta_high) & isfinite(theta_low);
    [~, p_theta_split, ~, stats_theta_split] = ttest(theta_high(vs), theta_low(vs));
    d_theta = mean(theta_high(vs) - theta_low(vs)) / std(theta_high(vs) - theta_low(vs));

    vs_ms = isfinite(ms_high) & isfinite(ms_low);
    [~, p_ms_split, ~, stats_ms_split] = ttest(ms_high(vs_ms), ms_low(vs_ms));
    d_ms = mean(ms_high(vs_ms) - ms_low(vs_ms)) / std(ms_high(vs_ms) - ms_low(vs_ms));

    fprintf('  Theta (high vs low alpha): t(%d) = %.2f, p = %.4f, d = %.3f\n', ...
        stats_theta_split.df, stats_theta_split.tstat, p_theta_split, d_theta);
    fprintf('  MS rate (high vs low alpha): t(%d) = %.2f, p = %.4f, d = %.3f\n', ...
        stats_ms_split.df, stats_ms_split.tstat, p_ms_split, d_ms);

    %% ================================================================
    %%       ANALYSIS 4: Theta-Velocity Coherence × Alpha Split
    %% ================================================================
    fprintf('\n--- Analysis 4: Theta-Velocity Coherence × Alpha (%s) ---\n', taskName);

    coh_hi = coh_highAlpha_subj(isfinite(coh_highAlpha_subj) & isfinite(coh_lowAlpha_subj));
    coh_lo = coh_lowAlpha_subj(isfinite(coh_highAlpha_subj) & isfinite(coh_lowAlpha_subj));

    if numel(coh_hi) > 2
        [~, p_coh_split, ~, stats_coh_split] = ttest(coh_hi, coh_lo);
        d_coh = mean(coh_hi - coh_lo) / std(coh_hi - coh_lo);
        fprintf('  Theta-Vel Coh (high vs low alpha): t(%d) = %.2f, p = %.4f, d = %.3f\n', ...
            stats_coh_split.df, stats_coh_split.tstat, p_coh_split, d_coh);
        fprintf('  High alpha: M = %.4f; Low alpha: M = %.4f\n', ...
            mean(coh_hi, 'omitnan'), mean(coh_lo, 'omitnan'));
    else
        p_coh_split = NaN; d_coh = NaN;
        fprintf('  Insufficient data for coherence split.\n');
    end

    %% ================================================================
    %%       ANALYSIS 5: Mediation — Alpha → Theta → Saccade Rate
    %% ================================================================
    fprintf('\n--- Analysis 5: Mediation & GLMMs (%s) ---\n', taskName);

    T = table(all_preAlpha, all_retTheta, all_msrate, ...
        categorical(all_cond), categorical(all_subj), ...
        'VariableNames', {'PreAlpha', 'RetTheta', 'MSRate', 'Condition', 'ID'});

    % Z-score within subject
    T.PreAlphaZ = nan(height(T), 1);
    T.RetThetaZ = nan(height(T), 1);
    for si = 1:nS
        idx = T.ID == categorical(unique_ids(si));
        vals = T.PreAlpha(idx);
        mu = mean(vals, 'omitnan'); sd = std(vals, 'omitnan');
        if sd > 0, T.PreAlphaZ(idx) = (vals - mu) / sd; else, T.PreAlphaZ(idx) = 0; end

        vals = T.RetTheta(idx);
        mu = mean(vals, 'omitnan'); sd = std(vals, 'omitnan');
        if sd > 0, T.RetThetaZ(idx) = (vals - mu) / sd; else, T.RetThetaZ(idx) = 0; end
    end

    % Mediation analysis (Baron & Kenny + Sobel test)
    beta_a = NaN; beta_b = NaN; beta_c = NaN; beta_cprime = NaN;
    SE_a = NaN; SE_b = NaN;
    p_a = NaN; p_b = NaN; p_c = NaN; p_cprime = NaN;
    indirect = NaN; z_sobel = NaN; p_sobel = NaN;
    mediation_label = 'Not computed';

    try
        % a-path: PreAlpha → RetTheta
        lme_a = fitlme(T, 'RetThetaZ ~ PreAlphaZ + (1 | ID)');
        beta_a = lme_a.Coefficients.Estimate(2);
        SE_a   = lme_a.Coefficients.SE(2);
        p_a    = lme_a.Coefficients.pValue(2);

        % b-path + c'-path: MSRate ~ PreAlpha + RetTheta
        lme_bc = fitlme(T, 'MSRate ~ PreAlphaZ + RetThetaZ + (1 | ID)');
        beta_cprime = lme_bc.Coefficients.Estimate(2);  % direct alpha→MS
        p_cprime    = lme_bc.Coefficients.pValue(2);
        beta_b      = lme_bc.Coefficients.Estimate(3);  % theta→MS (controlling alpha)
        SE_b        = lme_bc.Coefficients.SE(3);
        p_b         = lme_bc.Coefficients.pValue(3);

        % c-path: total effect PreAlpha → MSRate
        lme_c = fitlme(T, 'MSRate ~ PreAlphaZ + (1 | ID)');
        beta_c = lme_c.Coefficients.Estimate(2);
        p_c    = lme_c.Coefficients.pValue(2);

        % Indirect effect & Sobel test
        indirect = beta_a * beta_b;
        SE_indirect = sqrt(beta_a^2 * SE_b^2 + beta_b^2 * SE_a^2);
        z_sobel  = indirect / SE_indirect;
        p_sobel  = 2 * (1 - normcdf(abs(z_sobel)));

        % Interpret mediation
        if p_sobel < 0.05 && p_a < 0.05 && p_b < 0.05
            if p_cprime >= 0.05
                mediation_label = 'Full mediation';
            else
                mediation_label = 'Partial mediation';
            end
        else
            mediation_label = 'No significant mediation';
        end

        fprintf('\n  Mediation results:\n');
        fprintf('    a-path (alpha→theta):    beta = %.4f, p = %.4f\n', beta_a, p_a);
        fprintf('    b-path (theta→MS|alpha): beta = %.4f, p = %.4f\n', beta_b, p_b);
        fprintf('    c-path (alpha→MS total): beta = %.4f, p = %.4f\n', beta_c, p_c);
        fprintf('    c''-path (alpha→MS|theta): beta = %.4f, p = %.4f\n', beta_cprime, p_cprime);
        fprintf('    Indirect effect: %.4f, Sobel z = %.2f, p = %.4f\n', indirect, z_sobel, p_sobel);
        fprintf('    Interpretation: %s\n', mediation_label);

        % Condition-modulated models
        fprintf('\n  Condition-modulated GLMMs:\n');
        lme_at_cond = fitlme(T, 'RetThetaZ ~ PreAlphaZ * Condition + (1 | ID)');
        fprintf('    RetTheta ~ PreAlpha * Condition:\n');
        disp(lme_at_cond);

        lme_tm_cond = fitlme(T, 'MSRate ~ RetThetaZ * Condition + (1 | ID)');
        fprintf('    MSRate ~ RetTheta * Condition:\n');
        disp(lme_tm_cond);

    catch ME
        fprintf('  GLMM/mediation failed: %s\n', ME.message);
    end

    %% ================================================================
    %%                     FIGURE 1: Pre-stim Alpha → Retention Theta
    %% ================================================================
    close all
    figure('Position', [50 50 1800 600]);

    %% Panel A: Scatter (within-subject centered)
    subplot(1,3,1); hold on

    alpha_ws = all_preAlpha;
    theta_ws = all_retTheta;
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        alpha_ws(idx) = all_preAlpha(idx) - mean(all_preAlpha(idx), 'omitnan');
        theta_ws(idx) = all_retTheta(idx) - mean(all_retTheta(idx), 'omitnan');
    end

    h_sc = gobjects(numel(unique_conds), 1);
    for ci = 1:numel(unique_conds)
        idx = all_cond == unique_conds(ci);
        h_sc(ci) = scatter(alpha_ws(idx), theta_ws(idx), 15, colors(ci,:), ...
            'filled', 'MarkerFaceAlpha', 0.2);
    end

    vld = isfinite(alpha_ws) & isfinite(theta_ws);
    p_fit = polyfit(alpha_ws(vld), theta_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg = plot(xfit, polyval(p_fit, xfit), 'k-', 'LineWidth', 3);

    xlabel('Pre-stim \alpha Power (w-s centered)')
    ylabel('Retention \theta Power (w-s centered)')
    legend([h_sc; h_reg], [condLabels{taskIdx}; {'Regression'}], ...
        'Location', 'best', 'FontSize', fontSize-10)
    title(sprintf('r = %.3f [%.3f, %.3f], p = %.4f', ...
        mean_r_at, ci_r_at(1), ci_r_at(2), p_at), 'FontSize', fontSize-6)
    set(gca, 'FontSize', fontSize-6)

    %% Panel B: Median split bar plot
    subplot(1,3,2); hold on
    valid_s = isfinite(theta_high) & isfinite(theta_low);
    data_bar = [mean(theta_low(valid_s), 'omitnan'), mean(theta_high(valid_s), 'omitnan')];
    sem_bar  = [std(theta_low(valid_s), 'omitnan'), std(theta_high(valid_s), 'omitnan')] / sqrt(sum(valid_s));

    b = bar(1:2, data_bar, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(3,:);
    errorbar(1:2, data_bar, sem_bar, 'k.', 'LineWidth', 2, 'CapSize', 12)

    for si = find(valid_s)'
        plot([1 2], [theta_low(si), theta_high(si)], '-', 'Color', [.6 .6 .6 .4], 'LineWidth', 1);
    end

    set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \alpha', 'High \alpha'}, 'FontSize', fontSize-6)
    ylabel('Retention \theta Power')
    title(sprintf('Median Split: t(%d) = %.2f, p = %.4f, d = %.2f', ...
        stats_theta_split.df, stats_theta_split.tstat, p_theta_split, d_theta), ...
        'FontSize', fontSize-6)

    %% Panel C: Within-subject correlation raincloud
    subplot(1,3,3); hold on
    vals = r_alpha_theta(valid_r1);
    xj = 0.85 + 0.3 * rand(size(vals));
    scatter(xj, vals, 60, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot([1 1], ci_r_at, 'k-', 'LineWidth', 3)
    plot(1, mean_r_at, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(1,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2)
    yline(0, 'k:', 'LineWidth', 1.5)
    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('\\alpha \\rightarrow \\theta Correlation\nt(%d) = %.2f, p = %.4f', ...
        stats_at.df, stats_at.tstat, p_at), 'FontSize', fontSize-6)
    set(gca, 'XTick', [], 'FontSize', fontSize-6)

    sgtitle(sprintf('%s: Pre-stimulus Alpha \\rightarrow Retention Theta', upper(taskName)), ...
        'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_AlphaGating_1_AlphaTheta_' taskName '.png']));

    %% ================================================================
    %%                     FIGURE 2: Retention Theta → Saccade Rate
    %% ================================================================
    close all
    figure('Position', [50 50 1800 600]);

    %% Panel A: Scatter
    subplot(1,3,1); hold on
    theta_ws2 = all_retTheta;
    ms_ws     = all_msrate;
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        theta_ws2(idx) = all_retTheta(idx) - mean(all_retTheta(idx), 'omitnan');
        ms_ws(idx)     = all_msrate(idx) - mean(all_msrate(idx), 'omitnan');
    end

    h_sc2 = gobjects(numel(unique_conds), 1);
    for ci = 1:numel(unique_conds)
        idx = all_cond == unique_conds(ci);
        h_sc2(ci) = scatter(theta_ws2(idx), ms_ws(idx), 15, colors(ci,:), ...
            'filled', 'MarkerFaceAlpha', 0.2);
    end

    vld = isfinite(theta_ws2) & isfinite(ms_ws);
    p_fit = polyfit(theta_ws2(vld), ms_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg = plot(xfit, polyval(p_fit, xfit), 'k-', 'LineWidth', 3);

    xlabel('Retention \theta Power (w-s centered)')
    ylabel('MS Rate [/s] (w-s centered)')
    legend([h_sc2; h_reg], [condLabels{taskIdx}; {'Regression'}], ...
        'Location', 'best', 'FontSize', fontSize-10)
    title(sprintf('r = %.3f [%.3f, %.3f], p = %.4f', ...
        mean_r_tm, ci_r_tm(1), ci_r_tm(2), p_tm), 'FontSize', fontSize-6)
    set(gca, 'FontSize', fontSize-6)

    %% Panel B: Median split (theta on MS rate)
    subplot(1,3,2); hold on
    ms_hiTheta = nan(nS, 1); ms_loTheta = nan(nS, 1);
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        t = all_retTheta(idx); m = all_msrate(idx);
        vld_t = isfinite(t) & isfinite(m);
        t = t(vld_t); m = m(vld_t);
        if numel(t) < 10, continue; end
        med_t = median(t);
        ms_hiTheta(si) = mean(m(t >= med_t), 'omitnan');
        ms_loTheta(si) = mean(m(t < med_t), 'omitnan');
    end

    vs_tm = isfinite(ms_hiTheta) & isfinite(ms_loTheta);
    data_bar2 = [mean(ms_loTheta(vs_tm), 'omitnan'), mean(ms_hiTheta(vs_tm), 'omitnan')];
    sem_bar2  = [std(ms_loTheta(vs_tm), 'omitnan'), std(ms_hiTheta(vs_tm), 'omitnan')] / sqrt(sum(vs_tm));

    b = bar(1:2, data_bar2, 0.6);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(1,:);
    b.CData(2,:) = colors(3,:);
    errorbar(1:2, data_bar2, sem_bar2, 'k.', 'LineWidth', 2, 'CapSize', 12)
    for si = find(vs_tm)'
        plot([1 2], [ms_loTheta(si), ms_hiTheta(si)], '-', 'Color', [.6 .6 .6 .4], 'LineWidth', 1);
    end

    [~, p_tm_split] = ttest(ms_hiTheta(vs_tm), ms_loTheta(vs_tm));
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \theta', 'High \theta'}, 'FontSize', fontSize-6)
    ylabel('MS Rate [/s]')
    title(sprintf('Median Split: p = %.4f', p_tm_split), 'FontSize', fontSize-6)

    %% Panel C: Correlation distribution
    subplot(1,3,3); hold on
    vals2 = r_theta_ms(valid_r2);
    xj = 0.85 + 0.3 * rand(size(vals2));
    scatter(xj, vals2, 60, colors(2,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot([1 1], ci_r_tm, 'k-', 'LineWidth', 3)
    plot(1, mean_r_tm, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(2,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2)
    yline(0, 'k:', 'LineWidth', 1.5)
    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('\\theta \\rightarrow Saccade Correlation\nt(%d) = %.2f, p = %.4f', ...
        stats_tm.df, stats_tm.tstat, p_tm), 'FontSize', fontSize-6)
    set(gca, 'XTick', [], 'FontSize', fontSize-6)

    sgtitle(sprintf('%s: Retention Theta \\rightarrow Saccade Rate', upper(taskName)), ...
        'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_AlphaGating_2_ThetaSaccade_' taskName '.png']));

    %% ================================================================
    %%                     FIGURE 3: Full Hierarchical Model
    %% ================================================================
    close all
    figure('Position', [50 50 1600 1200]);

    %% Panel A: Alpha → Theta (bars)
    subplot(2,2,1); hold on
    b = bar(1:2, [mean(theta_low(vs), 'omitnan'), mean(theta_high(vs), 'omitnan')], 0.6);
    b.FaceColor = 'flat'; b.CData(1,:) = colors(1,:); b.CData(2,:) = colors(3,:);
    errorbar(1:2, ...
        [mean(theta_low(vs), 'omitnan'), mean(theta_high(vs), 'omitnan')], ...
        [std(theta_low(vs), 'omitnan'), std(theta_high(vs), 'omitnan')] / sqrt(sum(vs)), ...
        'k.', 'LineWidth', 2, 'CapSize', 12)
    for si = find(vs)', plot([1 2], [theta_low(si), theta_high(si)], '-', 'Color', [.6 .6 .6 .4]); end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \alpha', 'High \alpha'}, 'FontSize', fontSize-6)
    ylabel('Retention \theta Power')
    title(sprintf('\\alpha \\rightarrow \\theta\np = %.4f, d = %.2f', p_theta_split, d_theta), ...
        'FontSize', fontSize-4)

    %% Panel B: Alpha → Saccade Rate (bars)
    subplot(2,2,2); hold on
    b = bar(1:2, [mean(ms_low(vs_ms), 'omitnan'), mean(ms_high(vs_ms), 'omitnan')], 0.6);
    b.FaceColor = 'flat'; b.CData(1,:) = colors(1,:); b.CData(2,:) = colors(3,:);
    errorbar(1:2, ...
        [mean(ms_low(vs_ms), 'omitnan'), mean(ms_high(vs_ms), 'omitnan')], ...
        [std(ms_low(vs_ms), 'omitnan'), std(ms_high(vs_ms), 'omitnan')] / sqrt(sum(vs_ms)), ...
        'k.', 'LineWidth', 2, 'CapSize', 12)
    for si = find(vs_ms)', plot([1 2], [ms_low(si), ms_high(si)], '-', 'Color', [.6 .6 .6 .4]); end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \alpha', 'High \alpha'}, 'FontSize', fontSize-6)
    ylabel('MS Rate [/s]')
    title(sprintf('\\alpha \\rightarrow Saccades\np = %.4f, d = %.2f', p_ms_split, d_ms), ...
        'FontSize', fontSize-4)

    %% Panel C: Theta-Velocity Coherence × Alpha
    subplot(2,2,3); hold on
    if numel(coh_hi) > 2
        b = bar(1:2, [mean(coh_lo, 'omitnan'), mean(coh_hi, 'omitnan')], 0.6);
        b.FaceColor = 'flat'; b.CData(1,:) = colors(1,:); b.CData(2,:) = colors(3,:);
        errorbar(1:2, ...
            [mean(coh_lo, 'omitnan'), mean(coh_hi, 'omitnan')], ...
            [std(coh_lo, 'omitnan'), std(coh_hi, 'omitnan')] / sqrt(numel(coh_lo)), ...
            'k.', 'LineWidth', 2, 'CapSize', 12)
        for si = 1:numel(coh_lo)
            plot([1 2], [coh_lo(si), coh_hi(si)], '-', 'Color', [.6 .6 .6 .4]);
        end
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Low \alpha', 'High \alpha'}, 'FontSize', fontSize-6)
        ylabel('\theta-Velocity Coherence')
        title(sprintf('\\alpha \\rightarrow \\theta-Saccade Coh.\np = %.4f, d = %.2f', ...
            p_coh_split, d_coh), 'FontSize', fontSize-4)
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'FontSize', fontSize-4)
        axis off
    end

    %% Panel D: Mediation Path Model (text summary)
    subplot(2,2,4); axis off
    xlim([0 1]); ylim([0 1]);

    path_str = {
        '\bf Hierarchical Mediation Model \rm'
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

    sgtitle(sprintf('%s: Hierarchical Oscillatory Control Model', upper(taskName)), ...
        'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_AlphaGating_3_HierarchicalModel_' taskName '.png']));

    %% ================================================================
    %%                     FIGURE 4: Condition-Resolved Effects
    %% ================================================================
    close all
    figure('Position', [50 50 1800 600]);

    %% Panel A: Alpha → Theta correlation by condition
    subplot(1,3,1); hold on
    r_by_cond_at = nan(nS, numel(unique_conds));
    for ci = 1:numel(unique_conds)
        for si = 1:nS
            idx = all_subj == unique_ids(si) & all_cond == unique_conds(ci);
            a = all_preAlpha(idx); t = all_retTheta(idx);
            if sum(isfinite(a) & isfinite(t)) >= 5
                r_by_cond_at(si, ci) = corr(a, t, 'Type', 'Spearman', 'Rows', 'complete');
            end
        end
    end

    for ci = 1:numel(unique_conds)
        vals_c = r_by_cond_at(isfinite(r_by_cond_at(:,ci)), ci);
        xj = ci + 0.15 * randn(size(vals_c));
        scatter(xj, vals_c, 40, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.5);
        m_c = mean(vals_c); se_c = std(vals_c) / sqrt(numel(vals_c));
        plot([ci-0.2 ci+0.2], [m_c m_c], '-', 'Color', colors(ci,:), 'LineWidth', 3);
        plot([ci ci], [m_c-se_c m_c+se_c], '-', 'Color', colors(ci,:), 'LineWidth', 2);
    end
    yline(0, 'k:', 'LineWidth', 1.5)
    set(gca, 'XTick', 1:numel(unique_conds), 'XTickLabel', condLabels{taskIdx}, 'FontSize', fontSize-6)
    ylabel('Within-subject r (\alpha \rightarrow \theta)')
    title('\alpha \rightarrow \theta by Condition', 'FontSize', fontSize-4)

    %% Panel B: Theta → MS correlation by condition
    subplot(1,3,2); hold on
    r_by_cond_tm = nan(nS, numel(unique_conds));
    for ci = 1:numel(unique_conds)
        for si = 1:nS
            idx = all_subj == unique_ids(si) & all_cond == unique_conds(ci);
            t = all_retTheta(idx); m = all_msrate(idx);
            if sum(isfinite(t) & isfinite(m)) >= 5
                r_by_cond_tm(si, ci) = corr(t, m, 'Type', 'Spearman', 'Rows', 'complete');
            end
        end
    end

    for ci = 1:numel(unique_conds)
        vals_c = r_by_cond_tm(isfinite(r_by_cond_tm(:,ci)), ci);
        xj = ci + 0.15 * randn(size(vals_c));
        scatter(xj, vals_c, 40, colors(ci,:), 'filled', 'MarkerFaceAlpha', 0.5);
        m_c = mean(vals_c); se_c = std(vals_c) / sqrt(numel(vals_c));
        plot([ci-0.2 ci+0.2], [m_c m_c], '-', 'Color', colors(ci,:), 'LineWidth', 3);
        plot([ci ci], [m_c-se_c m_c+se_c], '-', 'Color', colors(ci,:), 'LineWidth', 2);
    end
    yline(0, 'k:', 'LineWidth', 1.5)
    set(gca, 'XTick', 1:numel(unique_conds), 'XTickLabel', condLabels{taskIdx}, 'FontSize', fontSize-6)
    ylabel('Within-subject r (\theta \rightarrow MS)')
    title('\theta \rightarrow Saccade by Condition', 'FontSize', fontSize-4)

    %% Panel C: Condition × Alpha split interaction on Theta
    subplot(1,3,3); hold on
    theta_cond_split = nan(nS, numel(unique_conds), 2);  % [subj × cond × (low/high)]
    for ci = 1:numel(unique_conds)
        for si = 1:nS
            idx = all_subj == unique_ids(si) & all_cond == unique_conds(ci);
            a = all_preAlpha(idx); t = all_retTheta(idx);
            vld_t = isfinite(a) & isfinite(t);
            a = a(vld_t); t = t(vld_t);
            if numel(a) < 6, continue; end
            med_a = median(a);
            theta_cond_split(si, ci, 1) = mean(t(a < med_a), 'omitnan');
            theta_cond_split(si, ci, 2) = mean(t(a >= med_a), 'omitnan');
        end
    end

    x_offset = [-0.15 0.15];
    split_labels = {'Low \alpha', 'High \alpha'};
    split_colors = [colors(1,:); colors(3,:)];
    h_split = gobjects(2, 1);
    for split_i = 1:2
        vals_mean = squeeze(nanmean(theta_cond_split(:, :, split_i), 1));
        vals_sem  = squeeze(nanstd(theta_cond_split(:, :, split_i), [], 1)) / sqrt(nS);
        h_split(split_i) = errorbar((1:numel(unique_conds)) + x_offset(split_i), ...
            vals_mean, vals_sem, 'o-', 'Color', split_colors(split_i,:), ...
            'MarkerFaceColor', split_colors(split_i,:), 'LineWidth', 2, ...
            'MarkerSize', 10, 'CapSize', 10);
    end

    set(gca, 'XTick', 1:numel(unique_conds), 'XTickLabel', condLabels{taskIdx}, 'FontSize', fontSize-6)
    ylabel('Retention \theta Power')
    legend(h_split, split_labels, 'Location', 'best', 'FontSize', fontSize-8)
    title('Condition \times \alpha Split', 'FontSize', fontSize-4)

    sgtitle(sprintf('%s: Condition-Resolved Effects', upper(taskName)), 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_AlphaGating_4_ConditionEffects_' taskName '.png']));

    %% Store results for this task
    results.(taskName).r_alpha_theta  = r_alpha_theta;
    results.(taskName).mean_r_at      = mean_r_at;
    results.(taskName).ci_r_at        = ci_r_at;
    results.(taskName).p_at           = p_at;
    results.(taskName).stats_at       = stats_at;
    results.(taskName).r_theta_ms     = r_theta_ms;
    results.(taskName).mean_r_tm      = mean_r_tm;
    results.(taskName).ci_r_tm        = ci_r_tm;
    results.(taskName).p_tm           = p_tm;
    results.(taskName).stats_tm       = stats_tm;
    results.(taskName).r_alpha_ms     = r_alpha_ms;
    results.(taskName).mean_r_am      = mean_r_am;
    results.(taskName).theta_high     = theta_high;
    results.(taskName).theta_low      = theta_low;
    results.(taskName).ms_high        = ms_high;
    results.(taskName).ms_low         = ms_low;
    results.(taskName).p_theta_split  = p_theta_split;
    results.(taskName).d_theta        = d_theta;
    results.(taskName).p_ms_split     = p_ms_split;
    results.(taskName).d_ms           = d_ms;
    results.(taskName).coh_highAlpha  = coh_hi;
    results.(taskName).coh_lowAlpha   = coh_lo;
    results.(taskName).p_coh_split    = p_coh_split;
    results.(taskName).d_coh          = d_coh;
    results.(taskName).mediation.beta_a      = beta_a;
    results.(taskName).mediation.beta_b      = beta_b;
    results.(taskName).mediation.beta_c      = beta_c;
    results.(taskName).mediation.beta_cprime = beta_cprime;
    results.(taskName).mediation.indirect    = indirect;
    results.(taskName).mediation.z_sobel     = z_sobel;
    results.(taskName).mediation.p_sobel     = p_sobel;
    results.(taskName).mediation.label       = mediation_label;
    results.(taskName).r_by_cond_at   = r_by_cond_at;
    results.(taskName).r_by_cond_tm   = r_by_cond_tm;

end % task loop

%% Save all results
save(fullfile(outdir, 'AOC_AlphaGating_ThetaSaccade_results.mat'), 'results');
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
