%% AOC Interaction — Alpha × Microsaccade Time-Course Coupling
% Correlates trial-averaged alpha power and microsaccade rate time courses
% within each participant (collapsed across conditions).
%
% Analyses:
%   A) Per-subject Pearson r between trial-averaged alpha and MS rate time
%      courses. Group-level: Fisher-z one-sample t-test. Per-subject
%      significance via circular-shift null distribution.
%   B) Sliding-window correlation (500 ms window, 50 ms step).
%      Grand average +/- SEM with cluster-based permutation test.
%
% Data: Sternberg (loads 2/4/6) and N-back (1/2/3-back), all conditions
%       collapsed.
%
% Artifact control:
%   1) Convolution padding — MS detection + Gaussian smoothing are performed
%      on a wider extraction window (±150 ms) to eliminate boundary artifacts;
%      the result is then trimmed to the analysis epoch.
%   2) Retention-only analysis — Analyses A & B are restricted to a clean
%      retention window [0.5, 1.8] s, excluding the stimulus-onset transient
%      (ERD, MS inhibition/rebound) and epoch-end artifacts.
%   3) Group-mean residualisation — The grand-average time course is
%      subtracted from each subject before computing correlations, so that
%      only individual deviations (not shared task-evoked structure) drive
%      the coupling measure.
%
% Methods:
%   Alpha: Posterior ROI average, bandpass [IAF-4, IAF+2] Hz, Hilbert
%   MS rate: Engbert & Kliegl (2003) detection, Gaussian-smoothed
%            (sigma = 50 ms) to obtain a continuous rate estimate
%   Time course: Trial-averaged per subject at 500 Hz, epoch [-0.5, 2] s
%
% Outputs per task:
%   Figure 1 — Dual-axis overlay: grand-average alpha + MS rate (full epoch)
%   Figure 2 — Per-subject correlation raincloud (retention window, residualised)
%   Figure 3 — Sliding-window r with cluster permutation (retention, residualised)
%   Figure 4 — Individual subject time-course overlays (residualised)

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 22;
fs       = 500;  % Hz

% Task definitions
tasks        = {'sternberg', 'nback'};
taskEEGFiles = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback'};
taskETFiles  = {'dataET_sternberg', 'dataET_nback'};
taskIAFFiles = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'};

% Analysis epoch (for extraction & plotting)
anaWin = [-0.5 2];
tVec   = anaWin(1) : 1/fs : anaWin(2);
nSamp  = numel(tVec);

% Gaussian smoothing kernel for MS rate (sigma = 50 ms)
sigma_ms   = 0.050;                       % seconds
sigma_samp = sigma_ms * fs;               % 25 samples
halfKern   = round(3 * sigma_samp);       % 3-sigma half-width
t_kern     = -halfKern:halfKern;
kernel     = exp(-t_kern.^2 / (2 * sigma_samp^2));
kernel     = kernel / sum(kernel);         % normalise to sum = 1

% FIX 1: Convolution padding — extract gaze from wider window
padTime    = 3 * sigma_ms;                % 150 ms padding per side
extractWin = [anaWin(1) - padTime, anaWin(2) + padTime];

% FIX 2: Retention-only correlation window (avoids onset transient + edges)
corrWin = [0.5 1.8];
corrIdx = tVec >= corrWin(1) & tVec <= corrWin(2);  % logical into tVec

% Sliding-window parameters
winWidth = 0.500;  % seconds
winStep  = 0.050;  % seconds
winSamp  = round(winWidth * fs);

% Permutation parameters
nPerm    = 2000;   % cluster permutations (Analysis B)
nPermA   = 1000;   % circular-shift permutations (Analysis A)
alphaThr = 0.05;

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/MSxAlpha_TimeCourse';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Posterior ROI channels (O, I, P)
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat', 'powload2_late');
allLabs = powload2_late.label;
roi_labels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I') || startsWith(L, 'P')
        roi_labels{end+1} = L; %#ok<SAGROW>
    end
end
fprintf('Posterior ROI (%d ch): %s\n', numel(roi_labels), strjoin(roi_labels, ', '));

%% Results container
results = struct();

%% ====================================================================
%%                         TASK LOOP
%% ====================================================================
for taskIdx = 1:2
    taskName = tasks{taskIdx};
    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName), repmat('=',1,60));

    %% Storage: subject-level time courses [nSubjects x nSamp]
    alpha_tc_all = nan(numel(subjects), nSamp);
    ms_tc_all    = nan(numel(subjects), nSamp);
    nTrials_subj = zeros(numel(subjects), 1);
    valid_subj   = false(numel(subjects), 1);

    %% ================================================================
    %%                       SUBJECT LOOP
    %% ================================================================
    for s = 1:numel(subjects)
        fprintf('  %s | Subject %s (%d/%d)... ', ...
            upper(taskName), subjects{s}, s, numel(subjects));

        dpeeg  = fullfile(path, subjects{s}, 'eeg');
        dpgaze = fullfile(path, subjects{s}, 'gaze');

        %% Load EEG
        try
            cd(dpeeg)
            load(taskEEGFiles{taskIdx})  % -> dataTFR
        catch
            fprintf('missing EEG, skipping.\n'); continue;
        end

        %% Subject-specific IAF
        IAF_file = fullfile(dpeeg, taskIAFFiles{taskIdx});
        if exist(IAF_file, 'file')
            tmp = load(IAF_file, 'IAF_subj');
            IAF = tmp.IAF_subj;
        else
            IAF = NaN;
        end
        if isfinite(IAF)
            bandAlpha = [IAF-4, IAF+2];
        else
            bandAlpha = [8 14];
        end

        %% Posterior alpha envelope (ROI average -> bandpass -> Hilbert)
        cfg = []; cfg.channel = roi_labels;
        roi = ft_selectdata(cfg, dataTFR);
        cfg = []; cfg.avgoverchan = 'yes';
        roiAvg = ft_selectdata(cfg, roi);
        clear roi

        cfg            = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = bandAlpha;
        cfg.demean     = 'yes';
        cfg.hilbert    = 'abs';
        alphaEnv = ft_preprocessing(cfg, roiAvg);
        clear roiAvg

        %% Load gaze
        try
            cd(dpgaze)
            load(taskETFiles{taskIdx}, 'dataETlong')
        catch
            fprintf('missing gaze, skipping.\n');
            clear dataTFR alphaEnv; continue;
        end

        %% Match trial counts
        nT = min(numel(alphaEnv.trial), size(dataETlong.trialinfo, 1));
        if nT == 0
            fprintf('no trials, skipping.\n');
            clear dataTFR alphaEnv dataETlong; continue;
        end

        %% Trial loop: accumulate time courses
        alpha_trials = nan(nT, nSamp);
        ms_trials    = nan(nT, nSamp);

        for tr = 1:nT
            t_eeg  = alphaEnv.time{tr};
            t_gaze = dataETlong.time{tr};

            %% ---- Alpha envelope ----
            ae = double(alphaEnv.trial{tr}(1, :));
            ae = fillmissing(ae, 'linear', 'EndValues', 'nearest');

            idx_eeg = t_eeg >= anaWin(1) & t_eeg <= anaWin(2);
            ae_seg  = ae(idx_eeg);

            if numel(ae_seg) >= nSamp
                alpha_trials(tr, :) = ae_seg(1:nSamp);
            else
                alpha_trials(tr, 1:numel(ae_seg)) = ae_seg;
            end

            %% ---- Microsaccade rate (with padded convolution) ----
            % FIX 1: Extract gaze from wider window to avoid edge artifacts
            idx_gaze_pad = t_gaze >= extractWin(1) & t_gaze <= extractWin(2);
            X_pad = double(dataETlong.trial{tr}(1, idx_gaze_pad));
            Y_pad = double(dataETlong.trial{tr}(2, idx_gaze_pad));
            t_gaze_pad = t_gaze(idx_gaze_pad);

            if numel(X_pad) < 50, continue; end

            % Blink removal + linear interpolation (on padded data)
            if size(dataETlong.trial{tr}, 1) >= 3
                A_pad = double(dataETlong.trial{tr}(3, idx_gaze_pad));
                blink = ~isfinite(A_pad) | (A_pad <= 0);
            else
                blink = ~isfinite(X_pad) | ~isfinite(Y_pad);
            end

            if any(blink)
                X_pad(blink) = NaN; Y_pad(blink) = NaN;
                if isnan(X_pad(1)),   i1 = find(isfinite(X_pad), 1, 'first'); if ~isempty(i1), X_pad(1:i1-1) = X_pad(i1); end, end
                if isnan(X_pad(end)), i2 = find(isfinite(X_pad), 1, 'last');  if ~isempty(i2), X_pad(i2+1:end) = X_pad(i2); end, end
                if isnan(Y_pad(1)),   i1 = find(isfinite(Y_pad), 1, 'first'); if ~isempty(i1), Y_pad(1:i1-1) = Y_pad(i1); end, end
                if isnan(Y_pad(end)), i2 = find(isfinite(Y_pad), 1, 'last');  if ~isempty(i2), Y_pad(i2+1:end) = Y_pad(i2); end, end
                X_pad = fillmissing(X_pad, 'linear');
                Y_pad = fillmissing(Y_pad, 'linear');
            end
            Y_pad = -Y_pad;  % Cartesian convention

            % Detect microsaccades on padded window
            ms = ms_detect_engbert(X_pad, Y_pad, fs, 6, 6e-3, 0.02);

            % Binary onset vector -> Gaussian-smoothed rate (Hz)
            ms_bin_pad = zeros(1, numel(X_pad));
            if ~isempty(ms.onsets)
                on = ms.onsets;
                on = on(on >= 1 & on <= numel(ms_bin_pad));
                ms_bin_pad(on) = 1;
            end
            ms_rate_pad = conv(ms_bin_pad, kernel, 'same') * fs;

            % Trim convolved signal back to analysis window
            trim_idx = t_gaze_pad >= anaWin(1) & t_gaze_pad <= anaWin(2);
            ms_rate  = ms_rate_pad(trim_idx);

            % Match to nSamp
            if numel(ms_rate) >= nSamp
                ms_trials(tr, :) = ms_rate(1:nSamp);
            else
                ms_trials(tr, 1:numel(ms_rate)) = ms_rate;
            end
        end % trial loop

        %% Average across trials -> subject-level time courses
        alpha_tc_all(s, :) = mean(alpha_trials, 1, 'omitnan');
        ms_tc_all(s, :)    = mean(ms_trials, 1, 'omitnan');
        nTrials_subj(s)    = nT;

        if all(isfinite(alpha_tc_all(s, :))) && all(isfinite(ms_tc_all(s, :)))
            valid_subj(s) = true;
        end

        clear dataTFR alphaEnv dataETlong
        fprintf('done (%d trials).\n', nT);
    end % subject loop

    %% Select valid subjects
    alpha_tc = alpha_tc_all(valid_subj, :);
    ms_tc    = ms_tc_all(valid_subj, :);
    subj_ids = subjects(valid_subj);
    nS       = sum(valid_subj);
    fprintf('\n  Valid subjects: %d / %d\n', nS, numel(subjects));

    %% FIX 3: Residualise — subtract group-mean time course
    ga_alpha_raw = mean(alpha_tc, 1, 'omitnan');
    ga_ms_raw    = mean(ms_tc, 1, 'omitnan');

    alpha_tc_resid = alpha_tc - ga_alpha_raw;
    ms_tc_resid    = ms_tc    - ga_ms_raw;

    fprintf('  Residualised: subtracted group-mean time course from each subject.\n');
    fprintf('  Correlation window: [%.1f, %.1f] s (%d samples)\n', ...
        corrWin(1), corrWin(2), sum(corrIdx));

    %% ================================================================
    %%  ANALYSIS A: Per-subject time-course correlation
    %%              (residualised, retention window only)
    %% ================================================================
    fprintf('\n--- Analysis A: Per-subject correlation (%s) ---\n', taskName);

    r_subj = nan(nS, 1);
    for si = 1:nS
        a = alpha_tc_resid(si, corrIdx);
        m = ms_tc_resid(si, corrIdx);
        vld = isfinite(a) & isfinite(m);
        if sum(vld) >= 50
            r_subj(si) = corr(a(vld)', m(vld)', 'Type', 'Pearson');
        end
    end

    % Group-level: Fisher-z -> one-sample t-test
    valid_r = isfinite(r_subj);
    z_subj  = atanh(r_subj(valid_r));
    [~, p_group, ~, stats_group] = ttest(z_subj);
    mean_r = tanh(mean(z_subj));
    se_z   = std(z_subj) / sqrt(numel(z_subj));
    ci_r   = tanh([mean(z_subj) - 1.96*se_z, mean(z_subj) + 1.96*se_z]);

    fprintf('  Group r: %.3f [%.3f, %.3f], t(%d) = %.2f, p = %.4f\n', ...
        mean_r, ci_r(1), ci_r(2), stats_group.df, stats_group.tstat, p_group);

    %% Circular-shift null distribution (per subject)
    p_subj = nan(nS, 1);
    for si = 1:nS
        a = alpha_tc_resid(si, corrIdx);
        m = ms_tc_resid(si, corrIdx);
        vld = isfinite(a) & isfinite(m);
        a_v = a(vld); m_v = m(vld);
        nV = numel(a_v);
        if nV < 50 || ~isfinite(r_subj(si)), continue; end

        r_null = nan(1, nPermA);
        for p = 1:nPermA
            shift = randi([round(nV * 0.1), round(nV * 0.9)]);
            r_null(p) = corr(a_v', circshift(m_v, shift)', 'Type', 'Pearson');
        end
        p_subj(si) = mean(abs(r_null) >= abs(r_subj(si)));
    end

    nSig = sum(p_subj < 0.05, 'omitnan');
    fprintf('  Significant (circ-shift, p < .05): %d / %d\n', ...
        nSig, sum(isfinite(p_subj)));

    %% ================================================================
    %%  ANALYSIS B: Sliding-window correlation
    %%              (residualised, windows restricted to corrWin)
    %% ================================================================
    fprintf('\n--- Analysis B: Sliding-window correlation (%s) ---\n', taskName);

    halfWin      = floor(winSamp / 2);
    stepSamp     = round(winStep * fs);
    corrSampStart = find(tVec >= corrWin(1), 1, 'first');
    corrSampEnd   = find(tVec <= corrWin(2), 1, 'last');

    % Window centres: entire window must fall within corrWin
    winCenters = (corrSampStart + halfWin) : stepSamp : (corrSampEnd - halfWin);
    nWin       = numel(winCenters);
    tWin       = tVec(winCenters);

    fprintf('  Sliding windows: %d (centres %.2f – %.2f s)\n', ...
        nWin, tWin(1), tWin(end));

    r_sliding = nan(nS, nWin);
    for si = 1:nS
        a = alpha_tc_resid(si, :);
        m = ms_tc_resid(si, :);
        for wi = 1:nWin
            idx = (winCenters(wi) - halfWin) : (winCenters(wi) + halfWin);
            a_w = a(idx);
            m_w = m(idx);
            vld = isfinite(a_w) & isfinite(m_w);
            if sum(vld) >= 50
                r_sliding(si, wi) = corr(a_w(vld)', m_w(vld)', 'Type', 'Pearson');
            end
        end
    end

    % Grand average
    ga_r  = mean(r_sliding, 1, 'omitnan');
    sem_r = std(r_sliding, [], 1, 'omitnan') / sqrt(nS);

    % Cluster-permutation on Fisher-z r values (one-sample vs 0)
    z_sliding = atanh(r_sliding);
    z_sliding(~isfinite(z_sliding)) = 0;
    [clusters, tvals_sw, thr] = cluster_permutation_1d(z_sliding, nPerm, alphaThr);

    nSigClus = 0;
    for k = 1:numel(clusters)
        if clusters(k).mass >= thr.mass
            nSigClus = nSigClus + 1;
            fprintf('  Cluster %d: [%.2f, %.2f] s, mass = %.1f, p = %.4f\n', ...
                nSigClus, tWin(clusters(k).idx(1)), ...
                tWin(clusters(k).idx(end)), clusters(k).mass, clusters(k).p);
        end
    end
    if nSigClus == 0, fprintf('  No significant clusters.\n'); end

    %% ================================================================
    %%  FIGURE 1: Dual-axis time-course overlay (raw, full epoch)
    %% ================================================================
    close all
    figure('Position', [50 50 1400 700]);

    sem_alpha = std(alpha_tc, [], 1, 'omitnan') / sqrt(nS);
    sem_ms_raw = std(ms_tc, [], 1, 'omitnan') / sqrt(nS);

    yyaxis left
    fill([tVec, fliplr(tVec)], ...
        [ga_alpha_raw - sem_alpha, fliplr(ga_alpha_raw + sem_alpha)], ...
        colors(1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
    plot(tVec, ga_alpha_raw, '-', 'LineWidth', 2.5, 'Color', colors(1,:))
    ylabel('Alpha Power [\muV]')
    set(gca, 'YColor', colors(1,:))

    yyaxis right
    fill([tVec, fliplr(tVec)], ...
        [ga_ms_raw - sem_ms_raw, fliplr(ga_ms_raw + sem_ms_raw)], ...
        colors(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(tVec, ga_ms_raw, '-', 'LineWidth', 2.5, 'Color', colors(2,:))
    ylabel('Microsaccade Rate [Hz]')
    set(gca, 'YColor', colors(2,:))

    % Mark correlation window
    xline(corrWin(1), 'k:', 'LineWidth', 1.5)
    xline(corrWin(2), 'k:', 'LineWidth', 1.5)
    xline(0, 'k--', 'LineWidth', 2)

    % Shade excluded regions
    yl = ylim;
    patch([anaWin(1) corrWin(1) corrWin(1) anaWin(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor', 'none');
    patch([corrWin(2) anaWin(2) anaWin(2) corrWin(2)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0 0 0], ...
        'FaceAlpha', 0.06, 'EdgeColor', 'none');

    xlabel('Time [s]')
    title(sprintf('%s: Alpha Power & Microsaccade Rate (shaded = excluded from correlation)', ...
        upper(taskName)), 'FontSize', fontSize - 2)
    set(gca, 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ...
        ['AOC_MSxAlpha_TC_overlay_' taskName '.png']));

    %% ================================================================
    %%  FIGURE 2: Per-subject correlations (raincloud, residualised)
    %% ================================================================
    figure('Position', [50 50 800 700]); hold on

    vals = r_subj(valid_r);
    xj = 0.85 + 0.3 * rand(size(vals));

    % Colour by circular-shift significance
    sig_mask = false(size(vals));
    p_valid = p_subj(valid_r);
    sig_mask(isfinite(p_valid)) = p_valid(isfinite(p_valid)) < 0.05;

    scatter(xj(~sig_mask), vals(~sig_mask), 60, [0.6 0.6 0.6], ...
        'filled', 'MarkerFaceAlpha', 0.6);
    scatter(xj(sig_mask), vals(sig_mask), 80, colors(1,:), ...
        'filled', 'MarkerFaceAlpha', 0.8);

    plot([1 1], ci_r, 'k-', 'LineWidth', 3)
    plot(1, mean_r, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(1,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2)
    yline(0, 'k:', 'LineWidth', 1.5)

    xlim([0.5 1.5])
    ylabel('Pearson r (\alpha \times MS rate)')
    title(sprintf(['%s: Per-Subject Correlation (residualised, [%.1f–%.1f] s)\n' ...
        'r = %.3f [%.3f, %.3f], t(%d) = %.2f, p = %.4f\n' ...
        'Significant: %d/%d (circ-shift, p < .05)'], ...
        upper(taskName), corrWin(1), corrWin(2), ...
        mean_r, ci_r(1), ci_r(2), ...
        stats_group.df, stats_group.tstat, p_group, ...
        nSig, sum(isfinite(p_subj))), 'FontSize', fontSize - 6)
    set(gca, 'XTick', [], 'FontSize', fontSize)
    legend({'n.s.', 'p < .05 (circ-shift)'}, 'Location', 'best', ...
        'FontSize', fontSize - 6)
    saveas(gcf, fullfile(outdir, ...
        ['AOC_MSxAlpha_TC_persubj_r_' taskName '.png']));

    %% ================================================================
    %%  FIGURE 3: Sliding-window r with cluster permutation (residualised)
    %% ================================================================
    figure('Position', [50 50 1400 700]); hold on

    fill([tWin, fliplr(tWin)], [ga_r - sem_r, fliplr(ga_r + sem_r)], ...
        colors(2,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(tWin, ga_r, '-', 'LineWidth', 3, 'Color', colors(2,:))

    % Shade significant clusters
    for k = 1:numel(clusters)
        if clusters(k).mass >= thr.mass
            cidx = clusters(k).idx;
            xx = [tWin(cidx), fliplr(tWin(cidx))];
            yy = [zeros(size(cidx)), fliplr(ga_r(cidx))];
            patch(xx, yy, colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end

    yline(0, 'k-'); xline(0, 'k--', 'LineWidth', 2)
    xlabel('Time [s] (window centre)')
    ylabel('Pearson r (\alpha \times MS rate)')
    title(sprintf(['%s: Sliding-Window Correlation ' ...
        '(%d ms window, %d ms step, residualised)'], ...
        upper(taskName), winWidth*1000, winStep*1000), 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ...
        ['AOC_MSxAlpha_TC_slidingwin_' taskName '.png']));

    %% ================================================================
    %%  FIGURE 4: Individual subject time-course overlays (residualised)
    %% ================================================================
    nPlot = min(nS, 16);
    nRows = ceil(sqrt(nPlot));
    nCols = ceil(nPlot / nRows);

    figure('Position', [50 50 1800 1200]);
    for si = 1:nPlot
        subplot(nRows, nCols, si); hold on

        % Z-score residualised time courses for overlay
        a_z = zscore(alpha_tc_resid(si, :));
        m_z = zscore(ms_tc_resid(si, :));

        plot(tVec, a_z, '-', 'LineWidth', 1.5, 'Color', colors(1,:))
        plot(tVec, m_z, '-', 'LineWidth', 1.5, 'Color', colors(2,:))
        xline(0, 'k--', 'LineWidth', 1)

        % Shade correlation window
        yl_sub = ylim;
        patch([corrWin(1) corrWin(2) corrWin(2) corrWin(1)], ...
            [yl_sub(1) yl_sub(1) yl_sub(2) yl_sub(2)], ...
            [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        % Re-plot lines on top
        plot(tVec, a_z, '-', 'LineWidth', 1.5, 'Color', colors(1,:))
        plot(tVec, m_z, '-', 'LineWidth', 1.5, 'Color', colors(2,:))

        if isfinite(r_subj(si))
            starStr = '';
            if isfinite(p_subj(si)) && p_subj(si) < 0.05
                starStr = '*';
            end
            title(sprintf('S%s: r=%.2f%s', subj_ids{si}, ...
                r_subj(si), starStr), 'FontSize', 10)
        else
            title(sprintf('S%s', subj_ids{si}), 'FontSize', 10)
        end

        if si == 1
            legend({'\alpha (z)', 'MS (z)'}, 'FontSize', 8, ...
                'Location', 'best')
        end
        xlim(anaWin)
        set(gca, 'FontSize', 9)
    end
    sgtitle(sprintf('%s: Individual Residualised Time Courses (z-scored, green = analysis window)', ...
        upper(taskName)), 'FontSize', fontSize - 2)
    saveas(gcf, fullfile(outdir, ...
        ['AOC_MSxAlpha_TC_individual_' taskName '.png']));

    %% ================================================================
    %%  Store results for this task
    %% ================================================================
    results.(taskName).alpha_tc       = alpha_tc;
    results.(taskName).ms_tc          = ms_tc;
    results.(taskName).alpha_tc_resid = alpha_tc_resid;
    results.(taskName).ms_tc_resid    = ms_tc_resid;
    results.(taskName).tVec           = tVec;
    results.(taskName).corrWin        = corrWin;
    results.(taskName).subj_ids       = subj_ids;
    results.(taskName).nTrials        = nTrials_subj(valid_subj);
    results.(taskName).r_subj         = r_subj;
    results.(taskName).p_subj         = p_subj;
    results.(taskName).mean_r         = mean_r;
    results.(taskName).ci_r           = ci_r;
    results.(taskName).p_group        = p_group;
    results.(taskName).stats_group    = stats_group;
    results.(taskName).r_sliding      = r_sliding;
    results.(taskName).tWin           = tWin;
    results.(taskName).ga_r           = ga_r;
    results.(taskName).sem_r          = sem_r;
    results.(taskName).clusters       = clusters;
    results.(taskName).tvals          = tvals_sw;
    results.(taskName).thr            = thr;

end % task loop

%% Save all results
save(fullfile(outdir, 'AOC_MSxAlpha_TimeCourse_results.mat'), 'results');
fprintf('\nDone! All results and figures saved to:\n  %s\n', outdir);

%% ====================================================================
%%                         LOCAL FUNCTIONS
%% ====================================================================

function ms = ms_detect_engbert(X, Y, fs, lambda, minDurSec, minISISec)
%MS_DETECT_ENGBERT  Detect microsaccades using Engbert & Kliegl (2003).
%   Returns struct with fields 'onsets' and 'offsets' (sample indices).
if nargin < 4 || isempty(lambda),    lambda = 6;        end
if nargin < 5 || isempty(minDurSec), minDurSec = 0.006; end
if nargin < 6 || isempty(minISISec), minISISec = 0.02;  end

[vx, vy] = compute_velocity_sg(X, Y, fs, 3);

% Robust median-based SD
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
Ts = 1 / fs;
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

function [clusters, tvals, thr] = cluster_permutation_1d(S, nPerm, alpha)
%CLUSTER_PERMUTATION_1D  One-sample cluster-permutation test (vs 0).
%   S: subjects x time points. Sign-flip permutation for null.
[ns, ~] = size(S);
se    = std(S, [], 1, 'omitnan') ./ sqrt(sum(isfinite(S), 1));
m     = mean(S, 1, 'omitnan');
tvals = m ./ se;

% Critical t (two-sided)
tcrit     = tinv(1 - 0.5 * alpha, ns - 1);
thr.tcrit = tcrit;

% Observed clusters
clusters = compute_clusters(tvals, tcrit);

% Permutation null distribution (max cluster mass)
maxMass = zeros(1, nPerm);
for p = 1:nPerm
    flips = (rand(ns, 1) > 0.5) * 2 - 1;
    Sprm  = S .* flips;
    seP   = std(Sprm, [], 1, 'omitnan') ./ sqrt(sum(isfinite(Sprm), 1));
    mP    = mean(Sprm, 1, 'omitnan');
    tP    = mP ./ seP;
    clP   = compute_clusters(tP, tcrit);
    if isempty(clP)
        maxMass(p) = 0;
    else
        maxMass(p) = max([clP.mass]);
    end
end
maxMass = sort(maxMass);
thr.mass = maxMass(round((1 - alpha) * nPerm));

% Assign cluster p-values
for k = 1:numel(clusters)
    clusters(k).p = 1 - mean(maxMass <= clusters(k).mass);
end
end

function clusters = compute_clusters(tvals, tcrit)
%COMPUTE_CLUSTERS  Find supra-threshold clusters and their mass.
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
