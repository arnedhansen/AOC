%% AOC Interaction — Microsaccade Rate × Alpha Power
% Two analyses for each task (Sternberg & N-back):
%   A) Trial-by-trial covariation of posterior alpha power and microsaccade rate
%      — Within-subject Spearman correlations (Fisher-z, group t-test)
%      — GLMM: MSRate ~ AlphaZ * Condition + (1|ID) via fitlme
%   B) Peri-microsaccadic alpha envelope (±1000 ms, event-locked)
%      — Grand average real vs surrogate (random events matched per trial)
%      — Cluster-permutation test on baseline-corrected difference
%
% Hypothesis: Lower alpha → fewer microsaccades (less fixation stability needed).
% Exploratory analysis.
%
% Key outputs (per task):
%   Scatter plot (within-subject), correlation raincloud,
%   peri-MS alpha curve, difference wave with significant clusters

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 26;

fs       = 500;           % EEG & ET sampling rate (Hz)
anaWin   = [0 2];         % retention window (s)

% Peri-microsaccadic parameters
periWin_s  = 1.0;                             % ±1000 ms
periSamp   = round(periWin_s * fs);           % 500 samples
nPeriSamp  = 2 * periSamp + 1;               % 1001 samples
periTime   = linspace(-periWin_s, periWin_s, nPeriSamp) * 1000; % ms

% Cluster-permutation parameters
nPerm    = 2000;
alphaThr = 0.05;

% Task definitions
tasks        = {'sternberg', 'nback'};
taskEEGFiles = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback'};
taskETFiles  = {'dataET_sternberg', 'dataET_nback'};
taskIAFFiles = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'};

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/MSxAlpha';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Posterior ROI labels (channels containing I, O, or starting with P)
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat');
allLabs = powload2_late.label;
roi_labels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I') || startsWith(L, 'P')
        roi_labels{end+1} = L; %#ok<SAGROW>
    end
end
fprintf('ROI channels (%d): %s\n', numel(roi_labels), strjoin(roi_labels, ', '));

%% Results container
results = struct();

%% ====================================================================
%%                         TASK LOOP
%% ====================================================================
for taskIdx = 1:2
    taskName = tasks{taskIdx};
    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName), repmat('=',1,60));

    %% Holders — Analysis A (trial-by-trial)
    all_alpha  = [];
    all_msrate = [];
    all_subj   = [];
    all_cond   = [];

    %% Holders — Analysis B (peri-microsaccadic)
    periMS_allSubs   = cell(numel(subjects), 1);   % subject-mean peri-MS envelope
    periSurr_allSubs = cell(numel(subjects), 1);   % subject-mean surrogate envelope
    nMS_perSubj      = zeros(numel(subjects), 1);

    %% ================================================================
    %%                       SUBJECT LOOP
    %% ================================================================
    for s = 1:numel(subjects)
        
        fprintf('Task: %s | Subject %s (%d/%d)\n', taskName, subjects{s}, s, numel(subjects));

        dpeeg  = fullfile(path, subjects{s}, 'eeg');
        dpgaze = fullfile(path, subjects{s}, 'gaze');

        %% Load EEG
        cd(dpeeg)
        load(taskEEGFiles{taskIdx})  % → dataTFR

        % Subject-specific IAF
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

        % Select posterior ROI channels, average across channels
        cfg = []; cfg.channel = roi_labels;
        roi = ft_selectdata(cfg, dataTFR);
        cfg = []; cfg.avgoverchan = 'yes';
        roiAvg = ft_selectdata(cfg, roi);

        % Bandpass filter in alpha band → Hilbert envelope (amplitude)
        cfg            = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = bandAlpha;
        cfg.demean     = 'yes';
        cfg.hilbert    = 'abs';
        alphaEnv = ft_preprocessing(cfg, roiAvg);

        %% Load gaze
        cd(dpgaze)
        load(taskETFiles{taskIdx}, 'dataETlong')

        % Match trial counts (EEG vs gaze)
        nT = min(numel(alphaEnv.trial), size(dataETlong.trialinfo, 1));

        %% Per-trial storage
        subj_alpha_trial = nan(nT, 1);
        subj_ms_trial    = nan(nT, 1);
        subj_cond_trial  = nan(nT, 1);

        periMS_epochs   = [];   % [nEpochs × nPeriSamp]
        periSurr_epochs = [];

        %% ============================================================
        %%                       TRIAL LOOP
        %% ============================================================
        for tr = 1:nT
            t_eeg  = alphaEnv.time{tr};
            t_gaze = dataETlong.time{tr};

            % --- Alpha envelope (full trial for epoch extraction) ---
            ae_full = double(alphaEnv.trial{tr}(1, :));
            ae_full = fillmissing(ae_full, 'linear', 'EndValues', 'nearest');

            % Retention window indices (EEG)
            idx_ana_eeg  = t_eeg >= anaWin(1) & t_eeg <= anaWin(2);
            ae_ana       = ae_full(idx_ana_eeg);

            % --- Gaze: retention window ---
            idx_ana_gaze = t_gaze >= anaWin(1) & t_gaze <= anaWin(2);
            X = double(dataETlong.trial{tr}(1, idx_ana_gaze));
            Y = double(dataETlong.trial{tr}(2, idx_ana_gaze));
            A = double(dataETlong.trial{tr}(3, idx_ana_gaze));
            Y = -Y;  % invert Y for Cartesian convention

            % Blink removal + linear interpolation
            blink = ~isfinite(A) | (A <= 0);
            if any(blink)
                X(blink) = NaN; Y(blink) = NaN;
                % Pad edges before interpolation
                if isnan(X(1)),   i1 = find(isfinite(X), 1, 'first'); if ~isempty(i1), X(1:i1-1) = X(i1); end, end
                if isnan(X(end)), i2 = find(isfinite(X), 1, 'last');  if ~isempty(i2), X(i2+1:end) = X(i2); end, end
                if isnan(Y(1)),   i1 = find(isfinite(Y), 1, 'first'); if ~isempty(i1), Y(1:i1-1) = Y(i1); end, end
                if isnan(Y(end)), i2 = find(isfinite(Y), 1, 'last');  if ~isempty(i2), Y(i2+1:end) = Y(i2); end, end
                X = fillmissing(X, 'linear');
                Y = fillmissing(Y, 'linear');
            end

            % Detect microsaccades (Engbert & Kliegl)
            ms = ms_detect_engbert(X, Y, fs, 6, 6e-3, 0.02);

            %% Analysis A: trial-level alpha power and microsaccade rate
            subj_alpha_trial(tr) = mean(ae_ana, 'omitnan');
            T_dur = numel(X) / fs;  % 2 seconds
            subj_ms_trial(tr)    = numel(ms.onsets) / T_dur;  % events/s
            subj_cond_trial(tr)  = dataETlong.trialinfo(tr, 1);

            %% Analysis B: peri-microsaccadic alpha epochs
            t_gaze_ana = t_gaze(idx_ana_gaze);

            for mi = 1:numel(ms.onsets)
                % Time of this microsaccade onset
                ms_time = t_gaze_ana(ms.onsets(mi));

                % Find closest EEG sample in full trial
                [~, eeg_idx] = min(abs(t_eeg - ms_time));

                % Check that ±500 ms fits within trial
                lo = eeg_idx - periSamp;
                hi = eeg_idx + periSamp;
                if lo < 1 || hi > numel(ae_full), continue; end

                epoch = ae_full(lo:hi);
                if any(isnan(epoch)), continue; end
                periMS_epochs = [periMS_epochs; epoch]; %#ok<AGROW>
            end

            % Surrogate: same number of random time points in retention window
            nMS_tr = numel(ms.onsets);
            if nMS_tr > 0
                % Valid EEG sample range (within retention, with margins)
                ana_start = find(t_eeg >= anaWin(1), 1, 'first');
                ana_end   = find(t_eeg <= anaWin(2), 1, 'last');
                lo_valid  = max(ana_start, 1 + periSamp);
                hi_valid  = min(ana_end, numel(ae_full) - periSamp);

                if hi_valid > lo_valid
                    surr_idx = randi([lo_valid, hi_valid], nMS_tr, 1);
                    for mi = 1:nMS_tr
                        epoch = ae_full(surr_idx(mi) - periSamp : surr_idx(mi) + periSamp);
                        if any(isnan(epoch)), continue; end
                        periSurr_epochs = [periSurr_epochs; epoch]; %#ok<AGROW>
                    end
                end
            end
        end % trial loop

        %% Append subject data (Analysis A)
        valid = isfinite(subj_alpha_trial) & isfinite(subj_ms_trial);
        all_alpha  = [all_alpha;  subj_alpha_trial(valid)]; %#ok<AGROW>
        all_msrate = [all_msrate; subj_ms_trial(valid)];    %#ok<AGROW>
        all_subj   = [all_subj;   repmat(str2double(subjects{s}), sum(valid), 1)]; %#ok<AGROW>
        all_cond   = [all_cond;   subj_cond_trial(valid)];  %#ok<AGROW>

        %% Store subject-level means (Analysis B)
        if ~isempty(periMS_epochs)
            periMS_allSubs{s} = mean(periMS_epochs, 1, 'omitnan');
        else
            periMS_allSubs{s} = nan(1, nPeriSamp);
        end
        if ~isempty(periSurr_epochs)
            periSurr_allSubs{s} = mean(periSurr_epochs, 1, 'omitnan');
        else
            periSurr_allSubs{s} = nan(1, nPeriSamp);
        end
        nMS_perSubj(s) = size(periMS_epochs, 1);
    end % subject loop

    %% ================================================================
    %%       ANALYSIS A: Trial-by-trial covariation
    %% ================================================================
    fprintf('\n--- Analysis A: Trial-by-trial covariation (%s) ---\n', taskName);

    % Convert condition codes to readable levels
    all_cond = all_cond - 20;  % Sternberg: 22/24/26 → 2/4/6; N-back: 21/22/23 → 1/2/3

    %% Within-subject Spearman correlations
    unique_ids = unique(all_subj);
    nS = numel(unique_ids);
    r_ws = nan(nS, 1);

    for si = 1:nS
        idx = all_subj == unique_ids(si);
        a = all_alpha(idx);
        m = all_msrate(idx);
        if sum(isfinite(a) & isfinite(m)) >= 5
            r_ws(si) = corr(a, m, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    % Group-level: Fisher-z transform → one-sample t-test vs 0
    valid_r = isfinite(r_ws);
    z_ws = atanh(r_ws(valid_r));
    [~, p_group, ~, stats_group] = ttest(z_ws);
    mean_r = tanh(mean(z_ws));
    se_z   = std(z_ws) / sqrt(numel(z_ws));
    ci_r   = tanh([mean(z_ws) - 1.96*se_z, mean(z_ws) + 1.96*se_z]);

    fprintf('  Within-subject r: mean = %.3f, 95%% CI = [%.3f, %.3f]\n', mean_r, ci_r(1), ci_r(2));
    fprintf('  t(%d) = %.2f, p = %.4f\n', stats_group.df, stats_group.tstat, p_group);

    %% GLMM: MSRate ~ AlphaZ * Condition + (1|ID)
    T = table(all_alpha, all_msrate, categorical(all_cond), categorical(all_subj), ...
        'VariableNames', {'AlphaPower', 'MSRate', 'Condition', 'ID'});

    % Z-score alpha power within subject (remove between-subject variance)
    T.AlphaZ = nan(height(T), 1);
    for si = 1:nS
        idx = T.ID == categorical(unique_ids(si));
        vals = T.AlphaPower(idx);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            T.AlphaZ(idx) = (vals - mu) / sd;
        else
            T.AlphaZ(idx) = 0;
        end
    end

    try
        lme_full    = fitlme(T, 'MSRate ~ AlphaZ * Condition + (1 | ID)');
        lme_reduced = fitlme(T, 'MSRate ~ AlphaZ + Condition + (1 | ID)');
        lrt = compare(lme_reduced, lme_full);

        fprintf('\n  GLMM: MSRate ~ AlphaZ * Condition + (1|ID)\n');
        fprintf('  LRT for interaction: Chi2 = %.2f, p = %.4f\n', lrt.LRStat(2), lrt.pValue(2));

        if lrt.pValue(2) < 0.05
            fprintf('  → Significant interaction; reporting full model\n');
            lme_final = lme_full;
        else
            fprintf('  → No significant interaction; reporting reduced model\n');
            lme_final = lme_reduced;
        end
        disp(lme_final)
    catch ME
        fprintf('  GLMM failed: %s\n', ME.message);
        lme_final = [];
    end

    %% Figure A1: Within-subject scatter (alpha vs MS rate)
    close all
    figure('Position', [0 0 1000 800]); hold on

    % Remove between-subject means for display
    alpha_ws  = all_alpha;
    msrate_ws = all_msrate;
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        alpha_ws(idx)  = all_alpha(idx)  - mean(all_alpha(idx),  'omitnan');
        msrate_ws(idx) = all_msrate(idx) - mean(all_msrate(idx), 'omitnan');
    end

    % Scatter by condition
    unique_conds = sort(unique(all_cond));
    nCond = numel(unique_conds);
    h_sc = gobjects(nCond, 1);
    for ci = 1:nCond
        idx = all_cond == unique_conds(ci);
        h_sc(ci) = scatter(alpha_ws(idx), msrate_ws(idx), 20, colors(ci,:), ...
            'filled', 'MarkerFaceAlpha', 0.3);
    end

    % Overall regression line
    vld = isfinite(alpha_ws) & isfinite(msrate_ws);
    p_fit = polyfit(alpha_ws(vld), msrate_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg = plot(xfit, polyval(p_fit, xfit), 'k-', 'LineWidth', 3);

    xlabel('Alpha Power (within-subject centered)')
    ylabel('Microsaccade Rate [MS/s] (within-subject centered)')

    if strcmp(taskName, 'sternberg')
        condLabels = arrayfun(@(x) sprintf('Load %d', x), unique_conds, 'UniformOutput', false);
    else
        condLabels = arrayfun(@(x) sprintf('%d-back', x), unique_conds, 'UniformOutput', false);
    end
    legend([h_sc; h_reg], [condLabels; {'Regression'}], 'Location', 'best', 'FontSize', fontSize - 8)

    title(sprintf('%s: Alpha Power vs Microsaccade Rate\nr = %.3f [%.3f, %.3f], p = %.4f', ...
        upper(taskName), mean_r, ci_r(1), ci_r(2), p_group), 'FontSize', fontSize - 4)
    set(gca, 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_MSxAlpha_scatter_' taskName '.png']));

    %% Figure A2: Raincloud of within-subject correlations
    figure('Position', [0 0 700 700]); hold on
    vals = r_ws(valid_r);

    % Jittered points
    xj = 0.85 + 0.3 * rand(size(vals));
    scatter(xj, vals, 60, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.6);

    % Mean + 95% CI
    plot([1 1], ci_r, 'k-', 'LineWidth', 3)
    plot(1, mean_r, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(1,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2)
    yline(0, 'k:', 'LineWidth', 1.5)

    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('%s: Within-Subject Correlations (\\alpha vs MS)\nr = %.3f, t(%d) = %.2f, p = %.4f', ...
        upper(taskName), mean_r, stats_group.df, stats_group.tstat, p_group), 'FontSize', fontSize - 6)
    set(gca, 'XTick', [], 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_MSxAlpha_corr_raincloud_' taskName '.png']));

    %% ================================================================
    %%       ANALYSIS B: Peri-microsaccadic alpha envelope
    %% ================================================================
    fprintf('\n--- Analysis B: Peri-microsaccadic alpha (%s) ---\n', taskName);
    fprintf('  MS per subject: %s\n', mat2str(nMS_perSubj'));

    % Assemble subject-level matrices [nSubjects × nPeriSamp]
    periMS_mat   = cell2mat(periMS_allSubs);
    periSurr_mat = cell2mat(periSurr_allSubs);

    % Exclude subjects with insufficient data
    validS = all(isfinite(periMS_mat), 2) & all(isfinite(periSurr_mat), 2);
    periMS_mat   = periMS_mat(validS, :);
    periSurr_mat = periSurr_mat(validS, :);
    nS_valid = sum(validS);
    fprintf('  Valid subjects for peri-MS: %d / %d\n', nS_valid, numel(subjects));

    % Baseline correction (% change from -1000 to -750 ms pre-event baseline)
    bl_idx = periTime >= -1000 & periTime <= -750;
    periMS_bl   = nan(size(periMS_mat));
    periSurr_bl = nan(size(periSurr_mat));
    for si = 1:nS_valid
        bl_real = mean(periMS_mat(si, bl_idx), 'omitnan');
        bl_surr = mean(periSurr_mat(si, bl_idx), 'omitnan');
        if bl_real > 0
            periMS_bl(si, :)   = (periMS_mat(si, :)   - bl_real) / bl_real * 100;
        end
        if bl_surr > 0
            periSurr_bl(si, :) = (periSurr_mat(si, :) - bl_surr) / bl_surr * 100;
        end
    end

    % Remove subjects with invalid baseline
    validBL = all(isfinite(periMS_bl), 2) & all(isfinite(periSurr_bl), 2);
    periMS_bl   = periMS_bl(validBL, :);
    periSurr_bl = periSurr_bl(validBL, :);
    nS_bl = sum(validBL);

    % Grand averages
    ga_real = mean(periMS_bl, 1, 'omitnan');
    sem_real = std(periMS_bl, [], 1, 'omitnan') / sqrt(nS_bl);
    ga_surr = mean(periSurr_bl, 1, 'omitnan');
    sem_surr = std(periSurr_bl, [], 1, 'omitnan') / sqrt(nS_bl);

    % Difference (real - surrogate) per subject
    diff_mat = periMS_bl - periSurr_bl;
    ga_diff  = mean(diff_mat, 1, 'omitnan');
    sem_diff = std(diff_mat, [], 1, 'omitnan') / sqrt(nS_bl);

    % Cluster-permutation test (paired, one-sample on difference vs 0)
    [clusters, tvals, thr] = cluster_permutation_1d(diff_mat, nPerm, alphaThr);

    % Print cluster results
    nSigClus = 0;
    for k = 1:numel(clusters)
        if clusters(k).mass >= thr.mass
            nSigClus = nSigClus + 1;
            fprintf('  Significant cluster %d: [%.0f, %.0f] ms, p = %.4f\n', ...
                nSigClus, periTime(clusters(k).idx(1)), periTime(clusters(k).idx(end)), clusters(k).p);
        end
    end
    if nSigClus == 0
        fprintf('  No significant clusters found.\n');
    end

    %% Figure B1: Grand average peri-microsaccadic alpha envelope
    close all
    figure('Position', [0 0 1200 800]); hold on

    % Surrogate (gray, dashed)
    fill([periTime, fliplr(periTime)], [ga_surr - sem_surr, fliplr(ga_surr + sem_surr)], ...
        [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h_surr = plot(periTime, ga_surr, '--', 'LineWidth', 2.5, 'Color', [0.5 0.5 0.5]);

    % Real (colored, solid)
    fill([periTime, fliplr(periTime)], [ga_real - sem_real, fliplr(ga_real + sem_real)], ...
        colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h_real = plot(periTime, ga_real, '-', 'LineWidth', 3, 'Color', colors(1,:));

    xline(0, 'k--', 'LineWidth', 2)
    yline(0, 'k:')
    xlabel('Time relative to microsaccade onset [ms]')
    ylabel('Alpha envelope (% change from baseline)')
    title(sprintf('%s: Peri-Microsaccadic Alpha Power', upper(taskName)), 'FontSize', fontSize)
    legend([h_real, h_surr], {'Microsaccade', 'Surrogate'}, 'Location', 'best', 'FontSize', fontSize - 6)
    set(gca, 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_periMS_alpha_' taskName '.png']));

    %% Figure B2: Difference wave (real - surrogate) with significant clusters
    figure('Position', [0 0 1200 600]); hold on

    fill([periTime, fliplr(periTime)], [ga_diff - sem_diff, fliplr(ga_diff + sem_diff)], ...
        colors(2,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(periTime, ga_diff, '-', 'LineWidth', 3, 'Color', colors(2,:))

    % Shade significant clusters
    for k = 1:numel(clusters)
        if clusters(k).mass >= thr.mass
            idx = clusters(k).idx;
            xx = [periTime(idx), fliplr(periTime(idx))];
            yy = [zeros(size(idx)), fliplr(ga_diff(idx))];
            patch(xx, yy, colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end

    xline(0, 'k--', 'LineWidth', 2)
    yline(0, 'k-')
    xlabel('Time relative to microsaccade onset [ms]')
    ylabel('Difference (Real - Surrogate) [%]')
    title(sprintf('%s: Peri-MS Alpha Difference (Real - Surrogate)', upper(taskName)), 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_periMS_alpha_diff_' taskName '.png']));

    %% Figure B3: Per-lag t-curve with significant clusters
    figure('Position', [0 0 1200 600]); hold on
    plot(periTime, tvals, 'LineWidth', 2, 'Color', colors(2,:))
    yline(thr.tcrit, 'k--'); yline(-thr.tcrit, 'k--');
    xline(0, 'k:', 'LineWidth', 1.5)
    xlabel('Time relative to microsaccade onset [ms]')
    ylabel('t-statistic (Real - Surrogate vs 0)')

    for k = 1:numel(clusters)
        if clusters(k).mass >= thr.mass
            idx = clusters(k).idx;
            xx = [periTime(idx), fliplr(periTime(idx))];
            yy = [zeros(size(idx)), fliplr(tvals(idx))];
            patch(xx, yy, colors(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    title(sprintf('%s: Peri-MS Alpha t-curve (cluster permutation)', upper(taskName)), 'FontSize', fontSize)
    set(gca, 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_periMS_alpha_tcurve_' taskName '.png']));

    %% Store results for this task
    results.(taskName).r_within      = r_ws;
    results.(taskName).mean_r        = mean_r;
    results.(taskName).ci_r          = ci_r;
    results.(taskName).p_group       = p_group;
    results.(taskName).t_group       = stats_group.tstat;
    results.(taskName).df_group      = stats_group.df;
    results.(taskName).nMS_perSubj   = nMS_perSubj;
    results.(taskName).periMS_subj   = periMS_mat;
    results.(taskName).periSurr_subj = periSurr_mat;
    results.(taskName).periMS_bl     = periMS_bl;
    results.(taskName).periSurr_bl   = periSurr_bl;
    results.(taskName).periTime      = periTime;
    results.(taskName).clusters      = clusters;
    results.(taskName).tvals         = tvals;
    results.(taskName).thr           = thr;
    if ~isempty(lme_final)
        results.(taskName).lme = lme_final;
    end

end % task loop

%% Save all results
save(fullfile(outdir, 'AOC_MSxAlpha_results.mat'), 'results');
fprintf('\nDone! Results and figures saved to:\n  %s\n', outdir);

%% ====================================================================
%%                         LOCAL FUNCTIONS
%% ====================================================================

function ms = ms_detect_engbert(X, Y, fs, lambda, minDurSec, minISISec)
%MS_DETECT_ENGBERT  Detect microsaccades using Engbert & Kliegl (2003).
%   Returns struct with fields 'onsets' and 'offsets' (sample indices).
if nargin < 4 || isempty(lambda),    lambda = 6;       end
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
%COMPUTE_VELOCITY_SG  Savitzky–Golay velocity estimation.
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
%   S: subjects × time points. Sign-flip permutation for null.
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
