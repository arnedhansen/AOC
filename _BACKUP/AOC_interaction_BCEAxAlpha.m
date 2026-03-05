%% AOC Interaction — BCEA × Alpha Power
% Exploratory analysis testing the relationship between fixation stability
% (BCEA, Bivariate Contour Ellipse Area) and posterior alpha power at the
% trial level, for both Sternberg and N-back tasks.
%
% BCEA (95% confidence):
%   BCEA = 2 * k95 * pi * sigma_x * sigma_y * sqrt(1 - rho^2)
%   where k95 = -log(1 - 0.95). Larger BCEA = less stable fixation.
%
% BCEA Lateralization:
%   BCEALat = (mean(x) - screenCentreX) / screenCentreX
%   Range: -1 (all left) to +1 (all right); 0 = centred.
%
% Posterior Alpha Power:
%   Mean power in IAF-band over occipital channels (O, I).
%
% Alpha Lateralization (Stroganova et al., 2007):
%   ALI = (R_alpha - L_alpha) / (R_alpha + L_alpha)
%
% Analyses (per task: Sternberg, N-back):
%   1. BCEA area × Alpha power (trial-level)
%      — Within-subject Spearman correlations (Fisher-z, group t-test)
%      — GLMM: AlphaPower ~ BCEAZ * Condition + (1 | ID)
%   2. BCEA lateralization × Alpha lateralization (trial-level)
%      — Within-subject Spearman correlations (Fisher-z, group t-test)
%      — GLMM: AlphaLat ~ BCEALatZ * Condition + (1 | ID)
%   3. Topographic split: high vs low BCEA trials (alpha topo median split)
%   4. Topographic split: left vs right gaze (BCEA lateralization direction)
%   5. Condition comparison panels (BCEA and alpha by load level)
%
% Visualizations (per task):
%   Fig 1: Scatter — BCEA vs alpha power (within-subject centred, by condition)
%   Fig 2: Raincloud — within-subject Spearman r (BCEA × alpha)
%   Fig 3: Scatter — BCEA lateralization vs alpha lateralization
%   Fig 4: Raincloud — within-subject Spearman r (BCEALat × alphaLat)
%   Fig 5: Topography — alpha power by BCEA median split (low / high / diff)
%   Fig 6: Topography — alpha power by gaze direction (left / right / diff)
%   Fig 7: Condition comparison (BCEA, alpha, lateralizations by load)
%
% Prerequisites:
%   - merged_data_sternberg_trials.mat / merged_data_nback_trials.mat
%   - dataEEG_TFR_sternberg.mat / dataEEG_TFR_nback.mat (for topos)
%   - IAF_*_subject.mat (for subject-specific alpha band)

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 26;

fs       = 500;           % EEG & ET sampling rate (Hz)
anaWin   = [0 2];         % retention / analysis window (s)
screenCX = 400;           % screen centre X (pixels, 800x600)

% Task definitions
tasks        = {'sternberg', 'nback'};
mergeFiles   = {'merged_data_sternberg_trials.mat', 'merged_data_nback_trials.mat'};
mergeVars    = {'merged_data_sternberg_trials', 'merged_data_nback_trials'};
taskEEGFiles = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback'};
taskIAFFiles = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'};

% Alpha band (fallback; overridden by IAF when available)
alphaRange = [8 14];

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/BCEAxAlpha';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Define posterior ROI (occipital: O, I channels)
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat', 'powload2_late');
allLabs = powload2_late.label;

occ_channels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I')
        occ_channels{end+1} = L; %#ok<SAGROW>
    end
end

fprintf('Occipital ROI (%d channels): %s\n', numel(occ_channels), strjoin(occ_channels, ', '));

%% Results container
results = struct();

%% ====================================================================
%%                         TASK LOOP
%% ====================================================================
for taskIdx = 1:numel(tasks)
    taskName = tasks{taskIdx};
    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName), repmat('=',1,60));

    %% Load merged trial-level data
    mergeFile = fullfile(path, mergeFiles{taskIdx});
    tmp = load(mergeFile);
    T = tmp.(mergeVars{taskIdx});
    fprintf('  Loaded %d trials from %s\n', height(T), mergeFiles{taskIdx});

    % Remove rows with missing key variables
    valid = isfinite(T.BCEAFull) & isfinite(T.AlphaPowerFull) & ...
            isfinite(T.BCEALatFull) & isfinite(T.Lateralization);
    T = T(valid, :);
    fprintf('  Valid trials (all key vars finite): %d\n', height(T));

    unique_ids   = unique(T.ID);
    nS           = numel(unique_ids);
    unique_conds = sort(unique(T.Condition));
    nCond        = numel(unique_conds);

    if strcmp(taskName, 'sternberg')
        condLabels = arrayfun(@(x) sprintf('Load %d', x), unique_conds, 'UniformOutput', false);
    else
        condLabels = arrayfun(@(x) sprintf('%d-back', x), unique_conds, 'UniformOutput', false);
    end

    %% ================================================================
    %%    ANALYSIS 1: BCEA × Alpha Power
    %% ================================================================
    fprintf('\n--- Analysis 1: BCEA x Alpha Power ---\n');

    %% Within-subject Spearman correlations
    r_bcea_alpha = nan(nS, 1);
    for si = 1:nS
        idx = T.ID == unique_ids(si);
        b = T.BCEAFull(idx);
        a = T.AlphaPowerFull(idx);
        if sum(isfinite(b) & isfinite(a)) >= 5
            r_bcea_alpha(si) = corr(b, a, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    % Group-level: Fisher-z -> one-sample t-test vs 0
    valid_r1 = isfinite(r_bcea_alpha);
    z_ws1 = atanh(r_bcea_alpha(valid_r1));
    if numel(z_ws1) >= 3
        [~, p_group1, ~, stats_group1] = ttest(z_ws1);
        mean_r1 = tanh(mean(z_ws1));
        se_z1   = std(z_ws1) / sqrt(numel(z_ws1));
        ci_r1   = tanh([mean(z_ws1) - 1.96*se_z1, mean(z_ws1) + 1.96*se_z1]);
    else
        p_group1 = NaN; mean_r1 = NaN; ci_r1 = [NaN NaN];
        stats_group1 = struct('df', NaN, 'tstat', NaN);
    end

    fprintf('  Within-subject r: mean = %.3f, 95%% CI = [%.3f, %.3f]\n', mean_r1, ci_r1(1), ci_r1(2));
    fprintf('  t(%d) = %.2f, p = %.4f\n', stats_group1.df, stats_group1.tstat, p_group1);

    %% GLMM: AlphaPower ~ BCEAZ * Condition + (1 | ID)
    Tglm1 = table(T.AlphaPowerFull, T.BCEAFull, categorical(T.Condition), categorical(T.ID), ...
        'VariableNames', {'AlphaPower', 'BCEA', 'Condition', 'ID'});

    % Z-score BCEA within subject
    Tglm1.BCEAZ = nan(height(Tglm1), 1);
    for si = 1:nS
        idx = Tglm1.ID == categorical(unique_ids(si));
        vals = Tglm1.BCEA(idx);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            Tglm1.BCEAZ(idx) = (vals - mu) / sd;
        else
            Tglm1.BCEAZ(idx) = 0;
        end
    end

    try
        lme1_full    = fitlme(Tglm1, 'AlphaPower ~ BCEAZ * Condition + (1 | ID)');
        lme1_reduced = fitlme(Tglm1, 'AlphaPower ~ BCEAZ + Condition + (1 | ID)');
        lrt1 = compare(lme1_reduced, lme1_full);

        fprintf('\n  GLMM: AlphaPower ~ BCEAZ * Condition + (1|ID)\n');
        fprintf('  LRT for interaction: Chi2 = %.2f, p = %.4f\n', lrt1.LRStat(2), lrt1.pValue(2));

        if lrt1.pValue(2) < 0.05
            fprintf('  -> Significant interaction; reporting full model\n');
            lme1_final = lme1_full;
        else
            fprintf('  -> No significant interaction; reporting reduced model\n');
            lme1_final = lme1_reduced;
        end
        disp(lme1_final)
    catch ME
        fprintf('  GLMM failed: %s\n', ME.message);
        lme1_final = [];
    end

    %% ================================================================
    %%    ANALYSIS 2: BCEA Lateralization × Alpha Lateralization
    %% ================================================================
    fprintf('\n--- Analysis 2: BCEA Lat x Alpha Lat ---\n');

    %% Within-subject Spearman correlations
    r_lat = nan(nS, 1);
    for si = 1:nS
        idx = T.ID == unique_ids(si);
        bl = T.BCEALatFull(idx);
        al = T.Lateralization(idx);
        if sum(isfinite(bl) & isfinite(al)) >= 5
            r_lat(si) = corr(bl, al, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    % Group-level: Fisher-z -> one-sample t-test vs 0
    valid_r2 = isfinite(r_lat);
    z_ws2 = atanh(r_lat(valid_r2));
    if numel(z_ws2) >= 3
        [~, p_group2, ~, stats_group2] = ttest(z_ws2);
        mean_r2 = tanh(mean(z_ws2));
        se_z2   = std(z_ws2) / sqrt(numel(z_ws2));
        ci_r2   = tanh([mean(z_ws2) - 1.96*se_z2, mean(z_ws2) + 1.96*se_z2]);
    else
        p_group2 = NaN; mean_r2 = NaN; ci_r2 = [NaN NaN];
        stats_group2 = struct('df', NaN, 'tstat', NaN);
    end

    fprintf('  Within-subject r: mean = %.3f, 95%% CI = [%.3f, %.3f]\n', mean_r2, ci_r2(1), ci_r2(2));
    fprintf('  t(%d) = %.2f, p = %.4f\n', stats_group2.df, stats_group2.tstat, p_group2);

    %% GLMM: AlphaLat ~ BCEALatZ * Condition + (1 | ID)
    Tglm2 = table(T.Lateralization, T.BCEALatFull, categorical(T.Condition), categorical(T.ID), ...
        'VariableNames', {'AlphaLat', 'BCEALat', 'Condition', 'ID'});

    % Z-score BCEA lateralization within subject
    Tglm2.BCEALatZ = nan(height(Tglm2), 1);
    for si = 1:nS
        idx = Tglm2.ID == categorical(unique_ids(si));
        vals = Tglm2.BCEALat(idx);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            Tglm2.BCEALatZ(idx) = (vals - mu) / sd;
        else
            Tglm2.BCEALatZ(idx) = 0;
        end
    end

    try
        lme2_full    = fitlme(Tglm2, 'AlphaLat ~ BCEALatZ * Condition + (1 | ID)');
        lme2_reduced = fitlme(Tglm2, 'AlphaLat ~ BCEALatZ + Condition + (1 | ID)');
        lrt2 = compare(lme2_reduced, lme2_full);

        fprintf('\n  GLMM: AlphaLat ~ BCEALatZ * Condition + (1|ID)\n');
        fprintf('  LRT for interaction: Chi2 = %.2f, p = %.4f\n', lrt2.LRStat(2), lrt2.pValue(2));

        if lrt2.pValue(2) < 0.05
            fprintf('  -> Significant interaction; reporting full model\n');
            lme2_final = lme2_full;
        else
            fprintf('  -> No significant interaction; reporting reduced model\n');
            lme2_final = lme2_reduced;
        end
        disp(lme2_final)
    catch ME
        fprintf('  GLMM failed: %s\n', ME.message);
        lme2_final = [];
    end

    %% ================================================================
    %%    TOPOGRAPHIC SPLIT: Load EEG per subject
    %% ================================================================
    fprintf('\n--- Topographic split (loading per-subject EEG) ---\n');

    topo_highBCEA_allSubj  = [];
    topo_lowBCEA_allSubj   = [];
    topo_gazeLeft_allSubj  = [];
    topo_gazeRight_allSubj = [];
    topo_template          = [];

    for s = 1:numel(subjects)
        fprintf('  Topo: Subject %s (%d/%d)\n', subjects{s}, s, numel(subjects));
        subjID = str2double(subjects{s});
        dpeeg  = fullfile(path, subjects{s}, 'eeg');

        % Load TFR data
        eegFile = fullfile(dpeeg, [taskEEGFiles{taskIdx} '.mat']);
        if ~exist(eegFile, 'file')
            fprintf('    Missing EEG data, skipping.\n');
            continue
        end
        cd(dpeeg);
        try
            load(taskEEGFiles{taskIdx}, 'dataTFR');
        catch
            fprintf('    Load error (EEG), skipping.\n');
            continue
        end

        % Subject-specific IAF
        IAF = NaN;
        iafFile = fullfile(dpeeg, taskIAFFiles{taskIdx});
        if exist(iafFile, 'file')
            tmp = load(iafFile, 'IAF_subj');
            if isfinite(tmp.IAF_subj)
                IAF = tmp.IAF_subj;
            end
        end
        if isfinite(IAF)
            bandAlpha = [IAF-4, IAF+2];
        else
            bandAlpha = alphaRange;
        end

        % Frequency analysis on retention window
        cfg            = [];
        cfg.method     = 'mtmfft';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 1:0.5:30;
        cfg.keeptrials = 'yes';

        cfgsel         = [];
        cfgsel.latency = anaWin;
        dataSel        = ft_selectdata(cfgsel, dataTFR);
        freq_trials    = ft_freqanalysis(cfg, dataSel);

        alpha_fidx = freq_trials.freq >= bandAlpha(1) & freq_trials.freq <= bandAlpha(2);
        nT_eeg = size(freq_trials.powspctrm, 1);

        % Match trials to merged table by global trial ID
        eeg_trialIDs = dataTFR.trialinfo(:, 2);
        subj_rows = T(T.ID == subjID, :);

        trialBCEA    = nan(nT_eeg, 1);
        trialBCEALat = nan(nT_eeg, 1);
        for tr = 1:min(nT_eeg, numel(eeg_trialIDs))
            match = subj_rows.Trial == eeg_trialIDs(tr);
            if any(match)
                trialBCEA(tr)    = subj_rows.BCEAFull(find(match, 1));
                trialBCEALat(tr) = subj_rows.BCEALatFull(find(match, 1));
            end
        end

        % --- Median split by BCEA area ---
        validBCEA = isfinite(trialBCEA);
        if sum(validBCEA) < 10
            clear dataTFR dataSel freq_trials
            continue
        end
        medBCEA = median(trialBCEA(validBCEA));

        highIdx = validBCEA & (trialBCEA >= medBCEA);
        lowIdx  = validBCEA & (trialBCEA < medBCEA);

        if sum(highIdx) >= 5 && sum(lowIdx) >= 5
            highPow = squeeze(mean(mean(freq_trials.powspctrm(highIdx, :, alpha_fidx), 3), 1));
            lowPow  = squeeze(mean(mean(freq_trials.powspctrm(lowIdx,  :, alpha_fidx), 3), 1));
            topo_highBCEA_allSubj = cat(2, topo_highBCEA_allSubj, highPow(:));
            topo_lowBCEA_allSubj  = cat(2, topo_lowBCEA_allSubj,  lowPow(:));
        end

        % --- Split by gaze direction (BCEA lateralization) ---
        gazeLeftIdx  = isfinite(trialBCEALat) & (trialBCEALat < 0);
        gazeRightIdx = isfinite(trialBCEALat) & (trialBCEALat > 0);

        if sum(gazeLeftIdx) >= 5
            leftPow = squeeze(mean(mean(freq_trials.powspctrm(gazeLeftIdx, :, alpha_fidx), 3), 1));
            topo_gazeLeft_allSubj = cat(2, topo_gazeLeft_allSubj, leftPow(:));
        end
        if sum(gazeRightIdx) >= 5
            rightPow = squeeze(mean(mean(freq_trials.powspctrm(gazeRightIdx, :, alpha_fidx), 3), 1));
            topo_gazeRight_allSubj = cat(2, topo_gazeRight_allSubj, rightPow(:));
        end

        % Save FieldTrip template for topoplot (first valid subject)
        if isempty(topo_template) && exist('freq_trials', 'var')
            topo_template = freq_trials;
            topo_template = rmfield(topo_template, 'powspctrm');
            topo_template.dimord = 'chan_freq';
            topo_template.freq   = mean(bandAlpha);
        end

        clear dataTFR dataSel freq_trials
    end % subject loop (topo)

    %% ================================================================
    %%    FIGURE 1: Scatter — BCEA vs Alpha Power
    %% ================================================================
    close all
    figure('Position', [0 0 1000 800]); hold on

    % Within-subject centering for display
    bcea_ws  = T.BCEAFull;
    alpha_ws = T.AlphaPowerFull;
    for si = 1:nS
        idx = T.ID == unique_ids(si);
        bcea_ws(idx)  = T.BCEAFull(idx)       - mean(T.BCEAFull(idx),       'omitnan');
        alpha_ws(idx) = T.AlphaPowerFull(idx)  - mean(T.AlphaPowerFull(idx), 'omitnan');
    end

    h_sc = gobjects(nCond, 1);
    for ci = 1:nCond
        idx = T.Condition == unique_conds(ci);
        cIdx = min(ci, size(colors, 1));
        h_sc(ci) = scatter(bcea_ws(idx), alpha_ws(idx), 20, colors(cIdx,:), ...
            'filled', 'MarkerFaceAlpha', 0.3);
    end

    % Regression line
    vld = isfinite(bcea_ws) & isfinite(alpha_ws);
    p_fit = polyfit(bcea_ws(vld), alpha_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg = plot(xfit, polyval(p_fit, xfit), 'k-', 'LineWidth', 3);

    xlabel('BCEA [px^2] (within-subject centred)')
    ylabel('Alpha Power [\muV^2/Hz] (within-subject centred)')
    legend([h_sc; h_reg], [condLabels; {'Regression'}], 'Location', 'best', 'FontSize', fontSize - 8)
    title(sprintf('%s: BCEA vs Alpha Power\nr = %.3f [%.3f, %.3f], p = %.4f', ...
        upper(taskName), mean_r1, ci_r1(1), ci_r1(2), p_group1), 'FontSize', fontSize - 4)
    set(gca, 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_BCEAxAlpha_scatter_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 2: Raincloud — Within-Subject Correlations (BCEA × Alpha)
    %% ================================================================
    figure('Position', [0 0 700 700]); hold on
    vals1 = r_bcea_alpha(valid_r1);

    % Jittered points
    xj = 0.85 + 0.3 * rand(size(vals1));
    scatter(xj, vals1, 60, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.6);

    % Mean + 95% CI
    if ~isnan(mean_r1)
        plot([1 1], ci_r1, 'k-', 'LineWidth', 3)
        plot(1, mean_r1, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(1,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2)
    end
    yline(0, 'k:', 'LineWidth', 1.5)

    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('%s: Within-Subject Correlations\nBCEA \\times \\alpha Power\nr = %.3f, t(%d) = %.2f, p = %.4f', ...
        upper(taskName), mean_r1, stats_group1.df, stats_group1.tstat, p_group1), 'FontSize', fontSize - 6)
    set(gca, 'XTick', [], 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_BCEAxAlpha_corr_raincloud_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 3: Scatter — BCEA Lateralization vs Alpha Lateralization
    %% ================================================================
    figure('Position', [0 0 1000 800]); hold on

    % Within-subject centering
    bcealat_ws  = T.BCEALatFull;
    alphalat_ws = T.Lateralization;
    for si = 1:nS
        idx = T.ID == unique_ids(si);
        bcealat_ws(idx)  = T.BCEALatFull(idx)   - mean(T.BCEALatFull(idx),   'omitnan');
        alphalat_ws(idx) = T.Lateralization(idx) - mean(T.Lateralization(idx), 'omitnan');
    end

    h_sc2 = gobjects(nCond, 1);
    for ci = 1:nCond
        idx = T.Condition == unique_conds(ci);
        cIdx = min(ci, size(colors, 1));
        h_sc2(ci) = scatter(bcealat_ws(idx), alphalat_ws(idx), 20, colors(cIdx,:), ...
            'filled', 'MarkerFaceAlpha', 0.3);
    end

    % Regression line
    vld = isfinite(bcealat_ws) & isfinite(alphalat_ws);
    p_fit2 = polyfit(bcealat_ws(vld), alphalat_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg2 = plot(xfit, polyval(p_fit2, xfit), 'k-', 'LineWidth', 3);

    xlabel('BCEA Lateralization [L-R] (within-subject centred)')
    ylabel('Alpha Lateralization (within-subject centred)')
    legend([h_sc2; h_reg2], [condLabels; {'Regression'}], 'Location', 'best', 'FontSize', fontSize - 8)
    title(sprintf('%s: BCEA Lat vs Alpha Lat\nr = %.3f [%.3f, %.3f], p = %.4f', ...
        upper(taskName), mean_r2, ci_r2(1), ci_r2(2), p_group2), 'FontSize', fontSize - 4)
    set(gca, 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_BCEALatxAlphaLat_scatter_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 4: Raincloud — Within-Subject Correlations (BCEALat × AlphaLat)
    %% ================================================================
    figure('Position', [0 0 700 700]); hold on
    vals2 = r_lat(valid_r2);

    % Jittered points
    xj = 0.85 + 0.3 * rand(size(vals2));
    scatter(xj, vals2, 60, colors(2,:), 'filled', 'MarkerFaceAlpha', 0.6);

    % Mean + 95% CI
    if ~isnan(mean_r2)
        plot([1 1], ci_r2, 'k-', 'LineWidth', 3)
        plot(1, mean_r2, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(2,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2)
    end
    yline(0, 'k:', 'LineWidth', 1.5)

    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('%s: Within-Subject Correlations\nBCEA Lat \\times \\alpha Lat\nr = %.3f, t(%d) = %.2f, p = %.4f', ...
        upper(taskName), mean_r2, stats_group2.df, stats_group2.tstat, p_group2), 'FontSize', fontSize - 6)
    set(gca, 'XTick', [], 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_BCEALatxAlphaLat_corr_raincloud_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 5: Topography — Alpha by BCEA (High vs Low Median Split)
    %% ================================================================
    if ~isempty(topo_highBCEA_allSubj) && ~isempty(topo_lowBCEA_allSubj) && ~isempty(topo_template)
        topo_H    = mean(topo_highBCEA_allSubj, 2);
        topo_L    = mean(topo_lowBCEA_allSubj, 2);
        topo_diff = topo_H - topo_L;

        topo_struct = topo_template;

        figure('Position', [0 0 1500 500]);

        cfg_topo = [];
        cfg_topo.xlim     = topo_template.freq * [1 1];
        cfg_topo.zlim     = 'maxabs';
        cfg_topo.layout   = 'EEG1005.lay';
        cfg_topo.colorbar = 'yes';
        cfg_topo.comment  = 'no';
        cfg_topo.style    = 'straight';

        % Panel 1: Low BCEA (stable fixation)
        subplot(1,3,1)
        topo_struct.powspctrm = topo_L;
        ft_topoplotER(cfg_topo, topo_struct);
        title('Low BCEA (stable fixation)', 'FontSize', fontSize - 4)

        % Panel 2: High BCEA (unstable fixation)
        subplot(1,3,2)
        topo_struct.powspctrm = topo_H;
        ft_topoplotER(cfg_topo, topo_struct);
        title('High BCEA (unstable fixation)', 'FontSize', fontSize - 4)

        % Panel 3: Difference (High - Low)
        subplot(1,3,3)
        topo_struct.powspctrm = topo_diff;
        cfg_topo.zlim = 'maxabs';
        ft_topoplotER(cfg_topo, topo_struct);
        title('Difference (High - Low)', 'FontSize', fontSize - 4)

        sgtitle(sprintf('%s: Alpha Topography by BCEA Median Split', upper(taskName)), 'FontSize', fontSize)
        saveas(gcf, fullfile(outdir, ['AOC_BCEAxAlpha_topo_' taskName '.png']));
    end

    %% ================================================================
    %%    FIGURE 6: Topography — Alpha by Gaze Direction (BCEA Lat)
    %% ================================================================
    if ~isempty(topo_gazeLeft_allSubj) && ~isempty(topo_gazeRight_allSubj) && ~isempty(topo_template)
        topo_gL   = mean(topo_gazeLeft_allSubj, 2);
        topo_gR   = mean(topo_gazeRight_allSubj, 2);
        topo_gDiff = topo_gR - topo_gL;

        topo_struct = topo_template;

        figure('Position', [0 0 1500 500]);

        cfg_topo = [];
        cfg_topo.xlim     = topo_template.freq * [1 1];
        cfg_topo.zlim     = 'maxabs';
        cfg_topo.layout   = 'EEG1005.lay';
        cfg_topo.colorbar = 'yes';
        cfg_topo.comment  = 'no';
        cfg_topo.style    = 'straight';

        % Panel 1: Gaze Left trials
        subplot(1,3,1)
        topo_struct.powspctrm = topo_gL;
        ft_topoplotER(cfg_topo, topo_struct);
        title('Gaze Left Trials', 'FontSize', fontSize - 4)

        % Panel 2: Gaze Right trials
        subplot(1,3,2)
        topo_struct.powspctrm = topo_gR;
        ft_topoplotER(cfg_topo, topo_struct);
        title('Gaze Right Trials', 'FontSize', fontSize - 4)

        % Panel 3: Difference (Right - Left)
        subplot(1,3,3)
        topo_struct.powspctrm = topo_gDiff;
        cfg_topo.zlim = 'maxabs';
        ft_topoplotER(cfg_topo, topo_struct);
        title('Difference (R - L gaze)', 'FontSize', fontSize - 4)

        sgtitle(sprintf('%s: Alpha Topography by BCEA Gaze Direction', upper(taskName)), 'FontSize', fontSize)
        saveas(gcf, fullfile(outdir, ['AOC_BCEALatxAlpha_topo_gaze_' taskName '.png']));
    end

    %% ================================================================
    %%    FIGURE 7: Condition Comparison (BCEA, Alpha, Lateralizations)
    %% ================================================================
    figure('Position', [0 0 1800 600]);

    % Panel A: BCEA by condition
    subplot(1,3,1); hold on
    for ci = 1:nCond
        cond_vals = nan(nS, 1);
        for si = 1:nS
            idx = T.ID == unique_ids(si) & T.Condition == unique_conds(ci);
            if sum(idx) >= 3
                cond_vals(si) = mean(T.BCEAFull(idx), 'omitnan');
            end
        end
        vld_c = isfinite(cond_vals);
        xj = ci + 0.15 * randn(sum(vld_c), 1);
        cIdx = min(ci, size(colors, 1));
        scatter(xj, cond_vals(vld_c), 60, colors(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.6);
        m_val = mean(cond_vals(vld_c));
        se_val = std(cond_vals(vld_c)) / sqrt(sum(vld_c));
        plot([ci-0.2, ci+0.2], [m_val m_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 3);
        plot([ci ci], [m_val-se_val m_val+se_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', fontSize - 4)
    ylabel('BCEA [px^2]')
    title('BCEA by Condition', 'FontSize', fontSize - 4)

    % Panel B: Alpha Power by condition
    subplot(1,3,2); hold on
    for ci = 1:nCond
        cond_vals = nan(nS, 1);
        for si = 1:nS
            idx = T.ID == unique_ids(si) & T.Condition == unique_conds(ci);
            if sum(idx) >= 3
                cond_vals(si) = mean(T.AlphaPowerFull(idx), 'omitnan');
            end
        end
        vld_c = isfinite(cond_vals);
        xj = ci + 0.15 * randn(sum(vld_c), 1);
        cIdx = min(ci, size(colors, 1));
        scatter(xj, cond_vals(vld_c), 60, colors(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.6);
        m_val = mean(cond_vals(vld_c));
        se_val = std(cond_vals(vld_c)) / sqrt(sum(vld_c));
        plot([ci-0.2, ci+0.2], [m_val m_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 3);
        plot([ci ci], [m_val-se_val m_val+se_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', fontSize - 4)
    ylabel('Alpha Power [\muV^2/Hz]')
    title('\alpha Power by Condition', 'FontSize', fontSize - 4)

    % Panel C: Lateralization comparison (subject-level scatter)
    subplot(1,3,3); hold on
    for ci = 1:nCond
        bcea_lat_cond  = nan(nS, 1);
        alpha_lat_cond = nan(nS, 1);
        for si = 1:nS
            idx = T.ID == unique_ids(si) & T.Condition == unique_conds(ci);
            if sum(idx) >= 3
                bcea_lat_cond(si)  = mean(T.BCEALatFull(idx), 'omitnan');
                alpha_lat_cond(si) = mean(T.Lateralization(idx), 'omitnan');
            end
        end
        vld_c = isfinite(bcea_lat_cond) & isfinite(alpha_lat_cond);
        cIdx = min(ci, size(colors, 1));
        scatter(bcea_lat_cond(vld_c), alpha_lat_cond(vld_c), 60, colors(cIdx,:), ...
            'filled', 'MarkerFaceAlpha', 0.6);
    end
    xlabel('BCEA Lateralization [L-R]')
    ylabel('\alpha Lateralization')
    legend(condLabels, 'Location', 'best', 'FontSize', fontSize - 8)
    title('Lateralization by Condition', 'FontSize', fontSize - 4)
    set(gca, 'FontSize', fontSize - 4)

    sgtitle(sprintf('%s: Condition Comparison', upper(taskName)), 'FontSize', fontSize)
    saveas(gcf, fullfile(outdir, ['AOC_BCEAxAlpha_conditions_' taskName '.png']));

    %% ================================================================
    %%    Store Results
    %% ================================================================
    results.(taskName).r_bcea_alpha    = r_bcea_alpha;
    results.(taskName).mean_r1         = mean_r1;
    results.(taskName).ci_r1           = ci_r1;
    results.(taskName).p_group1        = p_group1;
    results.(taskName).t_group1        = stats_group1.tstat;
    results.(taskName).df_group1       = stats_group1.df;
    results.(taskName).r_lat           = r_lat;
    results.(taskName).mean_r2         = mean_r2;
    results.(taskName).ci_r2           = ci_r2;
    results.(taskName).p_group2        = p_group2;
    results.(taskName).t_group2        = stats_group2.tstat;
    results.(taskName).df_group2       = stats_group2.df;
    results.(taskName).nTrials         = height(T);
    results.(taskName).nSubjects       = nS;
    if exist('lme1_final', 'var') && ~isempty(lme1_final)
        results.(taskName).lme_bcea_alpha = lme1_final;
    end
    if exist('lme2_final', 'var') && ~isempty(lme2_final)
        results.(taskName).lme_lat = lme2_final;
    end

end % task loop

%% Save all results
save(fullfile(outdir, 'AOC_BCEAxAlpha_results.mat'), 'results');
fprintf('\nDone! Results and figures saved to:\n  %s\n', outdir);
