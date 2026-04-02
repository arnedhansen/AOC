%% Reciprocal subject×load LMEs: Gaze% ~ Load*Alpha and Alpha ~ Load*Gaze%
%
% Long-format rows: one per subject × WM load (2/4/6). Alpha = occipital 8–14 Hz
% TFR power [dB] (FOOOF-baselined TFR), averaged 0.5–1.5 s. Gaze% = mean trial-level
% gaze deviation percentage change [%]: 100*(task/baseline - 1).
% Baseline [-0.5 -0.25] s, task [0.5 1.5] s.
%
% Load is coded as (Load-4)/2 so 0 = load 4, +1 = load 6 (same as AOC_split_AlphaLoads_500_1500ms).
% In each model the *moderator* is grand-mean centered across all rows (Alpha_c or GazePct_c) so
% the interaction is interpretable; DVs stay on the raw scale (GazePct, Alpha).
%
% Models (random slope for Load by subject, with intercept-only fallback):
%   (1) GazePct ~ Load * Alpha_c + (Load|Subject)
%   (2) Alpha ~ Load * GazePct_c + (Load|Subject)
%
% These are associational within-subject×load data; they do not identify causal direction.

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;

if ispc
    base_data = 'W:\Students\Arne\AOC';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
end
feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'splits', 'ReciprocalAlphaGazeLME_Pct');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

cond_vals = [2 4 6];
t_analysis_eeg = [0.5 1.5];
t_base_gaze = [-0.5 -0.25];
t_task_gaze = [0.5 1.5];

winsor_cfg = struct('enable', true, 'prctile', [2 98]);

log_dir = fullfile(base_data, 'data', 'controls', 'logs');
if ~isfolder(log_dir)
    mkdir(log_dir);
end
cmdlog_file = fullfile(log_dir, sprintf('AOC_alpha_gaze_reciprocal_LME_Pct_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
diary('off');
diary(cmdlog_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Command window log file: %s\n', cmdlog_file);

fprintf('\n=== Reciprocal LMEs: Gaze%%(GazePct) ~ Load*Alpha_c, Alpha ~ Load*GazePct_c ===\n');
fprintf('Figure/log directory: %s\n', fig_dir);

%% Load TFR per subject
cfg_bl = [];
cfg_bl.baseline = [-.5 -.25];
cfg_bl.baselinetype = 'absolute';
load2 = cell(1, numel(subjects));
load4 = cell(1, numel(subjects));
load6 = cell(1, numel(subjects));

for subj = 1:numel(subjects)
    fprintf('LOADING TFR Subject %d / %d\n', subj, numel(subjects));
    eeg_dir = fullfile(path, subjects{subj}, 'eeg');
    f = fullfile(eeg_dir, 'tfr_stern.mat');
    if ~isfile(f)
        error('Missing: %s', f);
    end
    datTFR = load(f);
    if isfield(datTFR, 'tfr2_fooof_bl')
        load2{subj} = datTFR.tfr2_fooof_bl;
        load4{subj} = datTFR.tfr4_fooof_bl;
        load6{subj} = datTFR.tfr6_fooof_bl;
    else
        load2{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr2_fooof);
        load4{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr4_fooof);
        load6{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr6_fooof);
    end
end
disp('TFR LOADING FINISHED');

%% Occipital channels
occ_channels = {};
labels = load2{1}.label;
for i = 1:numel(labels)
    label = labels{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label; %#ok<AGROW>
    end
end
channels = occ_channels;

%% Grand average + winsorize
cfg = [];
cfg.keepindividual = 'yes';
ga2 = ft_freqgrandaverage(cfg, load2{:});
ga4 = ft_freqgrandaverage(cfg, load4{:});
ga6 = ft_freqgrandaverage(cfg, load6{:});
if winsor_cfg.enable
    ga2 = winsorize_freq_subjects(ga2, winsor_cfg.prctile);
    ga4 = winsorize_freq_subjects(ga4, winsor_cfg.prctile);
    ga6 = winsorize_freq_subjects(ga6, winsor_cfg.prctile);
end

cfg = [];
cfg.frequency = [8 14];
cfg.latency = t_analysis_eeg;
cfg.channel = channels;
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
val2 = ft_selectdata(cfg, ga2);
val4 = ft_selectdata(cfg, ga4);
val6 = ft_selectdata(cfg, ga6);
alpha2 = val2.powspctrm;
alpha4 = val4.powspctrm;
alpha6 = val6.powspctrm;

nSubj = numel(subjects);

% Canonical figure size requested
fig_pos = [0 0 1512 982];

%% Gaze percentage change per load (ET)
fprintf('\nComputing gaze deviation percentage change [%%] per load from dataET_sternberg...\n');
[DEV2, DEV4, DEV6] = sternberg_gaze_deviation_percent_change_by_load(subjects, path, t_base_gaze, t_task_gaze);

%% Long table: subject × load
Subject_col = [];
Load_raw = [];
GazePct_col = [];
Alpha_col = [];

idx_ok = find(all(isfinite([alpha2(:), alpha4(:), alpha6(:), DEV2(:), DEV4(:), DEV6(:)]), 2));
fprintf('\nSubjects with finite alpha (2/4/6) and gaze %% change (2/4/6): %d / %d\n', numel(idx_ok), nSubj);

for k = 1:numel(idx_ok)
    ii = idx_ok(k);
    av = [alpha2(ii), alpha4(ii), alpha6(ii)];
    gv = [DEV2(ii), DEV4(ii), DEV6(ii)];
    for c = 1:3
        Subject_col(end+1, 1) = ii; %#ok<AGROW>
        Load_raw(end+1, 1) = cond_vals(c); %#ok<AGROW>
        Alpha_col(end+1, 1) = av(c); %#ok<AGROW>
        GazePct_col(end+1, 1) = gv(c); %#ok<AGROW>
    end
end

if numel(Subject_col) < 15
    error('Too few observations for LME (need at least ~15 rows).');
end

mu_alpha_grand = mean(Alpha_col);
mu_gaze_grand = mean(GazePct_col);
Alpha_c = Alpha_col - mu_alpha_grand;
GazePct_c = GazePct_col - mu_gaze_grand;

Load_ctr = (Load_raw - 4) / 2;

tbl = table(Subject_col, Load_ctr, GazePct_col, Alpha_col, Alpha_c, GazePct_c, ...
    'VariableNames', {'Subject', 'Load', 'GazePct', 'Alpha', 'Alpha_c', 'GazePct_c'});
tbl.Subject = categorical(tbl.Subject);

fprintf('\n--- Data summary ---\n');
fprintf('Rows (subject×load): %d; subjects: %d\n', height(tbl), numel(unique(Subject_col)));
fprintf('Grand means: Alpha=%.6g [dB], Gaze%%=%.6g [%%] (subtracted from moderator in each model).\n', ...
    mu_alpha_grand, mu_gaze_grand);
fprintf('Load coding: Load = (WM_load - 4) / 2  (0 = load 4; +1 = load 6; -1 = load 2).\n');
fprintf('EEG window [%.2f %.2f] s; gaze task [%.2f %.2f] s; gaze baseline [%.2f %.2f] s.\n', ...
    t_analysis_eeg(1), t_analysis_eeg(2), t_task_gaze(1), t_task_gaze(2), t_base_gaze(1), t_base_gaze(2));

%% Model 1: GazePct ~ Load * Alpha_c + (Load|Subject)
fprintf('\n========== MODEL 1: Gaze%% ~ Load * Alpha_c + (Load|Subject) ==========\n');
fprintf('DV: Gaze%% [%] at each load. Moderator: Alpha_c (Alpha grand-mean centered).\n');
try
    [lme1, f1] = fitlme_with_random_slope_fallback(tbl, 'GazePct ~ Load * Alpha_c');
    fprintf('Fitted: %s\n', f1);
    disp(anova(lme1));
    fprintf('Fixed effects:\n');
    disp(lme1.Coefficients);
    try
        fprintf('\n95%% CI (fixed):\n');
        disp(coefCI(lme1));
    catch %#ok<CTCH>
    end
catch ME
    fprintf('Model 1 failed: %s\n', ME.message);
end

%% Model 2: Alpha ~ Load * GazePct_c + (Load|Subject)
fprintf('\n========== MODEL 2: Alpha ~ Load * GazePct_c + (Load|Subject) ==========\n');
fprintf('DV: Alpha [dB] at each load. Moderator: GazePct_c (Gaze%% grand-mean centered).\n');
try
    [lme2, f2] = fitlme_with_random_slope_fallback(tbl, 'Alpha ~ Load * GazePct_c');
    fprintf('Fitted: %s\n', f2);
    disp(anova(lme2));
    fprintf('Fixed effects:\n');
    disp(lme2.Coefficients);
    try
        fprintf('\n95%% CI (fixed):\n');
        disp(coefCI(lme2));
    catch %#ok<CTCH>
    end
catch ME
    fprintf('Model 2 failed: %s\n', ME.message);
end

%% ========================= Visualizations =========================
fprintf('\n=== Visualizations: predicted interactions + subject-level slope coupling ===\n');

% Subject-level slopes across loads (using load coding -1,0,+1)
load_step = [-1; 0; 1];
alphaSlope_subj = nan(nSubj, 1);
gazePctSlope_subj = nan(nSubj, 1);
for ii = 1:nSubj
    if all(isfinite([alpha2(ii), alpha4(ii), alpha6(ii)])) && all(isfinite([DEV2(ii), DEV4(ii), DEV6(ii)]))
        yA = [alpha2(ii); alpha4(ii); alpha6(ii)];
        yG = [DEV2(ii); DEV4(ii); DEV6(ii)];
        Xmat = [ones(3,1), load_step];
        bA = Xmat \ yA;
        bG = Xmat \ yG;
        alphaSlope_subj(ii) = bA(2);      % change per +2 items (because load_step is in units of +2 items)
        gazePctSlope_subj(ii) = bG(2);   % change per +2 items
    end
end
maskSlope = isfinite(alphaSlope_subj) & isfinite(gazePctSlope_subj);
nSlope = sum(maskSlope);
fprintf('Subject-level slope coupling usable pairs: %d\n', nSlope);

% Descriptives for the coupling model sample
x_slope = alphaSlope_subj(maskSlope);
y_slope = gazePctSlope_subj(maskSlope);
if nSlope > 0
    fprintf('  Alpha slope (dB per +2 items): min=%.6g, max=%.6g, mean=%.6g, SD=%.6g, median=%.6g, n=%d\n', ...
        min(x_slope), max(x_slope), mean(x_slope), std(x_slope), median(x_slope), nSlope);
    fprintf('  Gaze%% slope (%% per +2 items): min=%.6g, max=%.6g, mean=%.6g, SD=%.6g, median=%.6g, n=%d\n', ...
        min(y_slope), max(y_slope), mean(y_slope), std(y_slope), median(y_slope), nSlope);
end

% Linear regression (for CW output)
p_fit_slopeCouple = [];
if nSlope >= 2
    p_fit_slopeCouple = polyfit(x_slope, y_slope, 1); % [slope intercept]
    yhat = polyval(p_fit_slopeCouple, x_slope);
    SSE = sum((y_slope - yhat).^2, 'omitnan');
    SST = sum((y_slope - mean(y_slope)).^2, 'omitnan');
    R2 = 1 - SSE / max(SST, eps);
    fprintf('  OLS (unstandardized): Gaze%%_slope = %.6g + %.6g * Alpha_slope; R^2=%.5f\n', ...
        p_fit_slopeCouple(2), p_fit_slopeCouple(1), R2);

    if nSlope >= 3
        df = nSlope - 2;
        x_bar = mean(x_slope);
        Sxx = sum((x_slope - x_bar).^2, 'omitnan');
        s2 = SSE / max(df, 1);
        SE_slope = sqrt(s2 / max(Sxx, eps));
        t_slope = p_fit_slopeCouple(1) / max(SE_slope, eps);
        p_two = 2 * tcdf(-abs(t_slope), df);
        tcrit = tinv(1 - 0.025, df);
        ci_slope = p_fit_slopeCouple(1) + [-1 1] * tcrit * SE_slope;
        fprintf('  Slope test: t(%d)=%.4g, p(two-sided)=%.6g, slope SE=%.6g, 95%% CI [%.6g, %.6g]\n', ...
            df, t_slope, p_two, SE_slope, ci_slope(1), ci_slope(2));
    else
        fprintf('  Slope inference skipped (need n>=3).\n');
    end
end

if nSlope >= 3
    [r_p, p_p] = corr(alphaSlope_subj(maskSlope), gazePctSlope_subj(maskSlope), 'Type', 'Pearson', 'rows', 'complete');
    [r_s, p_s] = corr(alphaSlope_subj(maskSlope), gazePctSlope_subj(maskSlope), 'Type', 'Spearman', 'rows', 'complete');
    fprintf('Slope coupling: Pearson r=%.5f (p=%.6g); Spearman rho=%.5f (p=%.6g)\n', r_p, p_p, r_s, p_s);
else
    fprintf('Slope coupling correlation skipped (need >= 3 pairs).\n');
end

% Plot: alpha slope vs gaze% slope
figure('Position', fig_pos, 'Color', 'w');
hold on;
if nSlope > 0
    scatter(alphaSlope_subj(maskSlope), gazePctSlope_subj(maskSlope), 35, [0.2 0.4 0.8], 'filled', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.2 0.4 0.8] * 0.5, 'LineWidth', 0.6);
    if nSlope >= 2 && ~isempty(p_fit_slopeCouple)
        xl = linspace(min(alphaSlope_subj(maskSlope)), max(alphaSlope_subj(maskSlope)), 100);
        plot(xl, polyval(p_fit_slopeCouple, xl), '-', 'Color', [0.2 0.4 0.8] * 0.7, 'LineWidth', 2);
    end
end
xline(0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5);
yline(0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5);
grid on; box off;
xlabel('Alpha slope across load (dB per +2 items)');
ylabel('Gaze% slope across load (% per +2 items)');
title('Subject-level coupling: alpha slope vs gaze% slope');
saveas(gcf, fullfile(fig_dir, 'AOC_reciprocal_SubjectSlopeCoupling_AlphaSlope_vs_GazePctSlope.png'));
close(gcf);

% Plot: predicted interaction curves for Model 1 and Model 2
% Model 1: GazePct ~ Load * Alpha_c
if exist('lme1', 'var')
    loadLevels = unique(tbl.Load);
    % Map coded load to labels
    % Load coded: (2-4)/2=-1, (4-4)/2=0, (6-4)/2=+1
    loadLabels = cell(size(loadLevels));
    for iL = 1:numel(loadLevels)
        if abs(loadLevels(iL) + 1) < 1e-9
            loadLabels{iL} = 'WM load 2';
        elseif abs(loadLevels(iL) - 0) < 1e-9
            loadLabels{iL} = 'WM load 4';
        elseif abs(loadLevels(iL) - 1) < 1e-9
            loadLabels{iL} = 'WM load 6';
        else
            loadLabels{iL} = sprintf('Load code %.3g', loadLevels(iL));
        end
    end

    alphaGrid = linspace(min(tbl.Alpha_c), max(tbl.Alpha_c), 120)';
    nGrid = numel(alphaGrid);
    nLevels = numel(loadLevels);
    subjDummy = tbl.Subject(1);

    figure('Position', fig_pos, 'Color', 'w');
    hold on;
    for iL = 1:nLevels
        tblPred = table( ...
            repmat(subjDummy, nGrid, 1), ...
            repmat(loadLevels(iL), nGrid, 1), ...
            alphaGrid, ...
            'VariableNames', {'Subject', 'Load', 'Alpha_c'});
        try
            yhat = predict(lme1, tblPred, 'Conditional', false);
        catch
            yhat = predict(lme1, tblPred);
        end
        plot(alphaGrid, yhat, 'LineWidth', 2, 'DisplayName', loadLabels{iL});
    end
    yline(0, 'k:', 'LineWidth', 1.5);
    grid on; box off;
    xlabel('Alpha_c (grand-mean centered alpha [dB])');
    ylabel('Predicted Gaze%');
    title('Model 1 prediction: Gaze% ~ Load * Alpha_c');
    legend('Location', 'best', 'Box', 'off');
    saveas(gcf, fullfile(fig_dir, 'AOC_reciprocal_Model1_PredictedInteraction_Load_by_Alpha_c_GazePct.png'));
    close(gcf);
end

% Model 2: Alpha ~ Load * GazePct_c
if exist('lme2', 'var')
    loadLevels = unique(tbl.Load);
    loadLabels = cell(size(loadLevels));
    for iL = 1:numel(loadLevels)
        if abs(loadLevels(iL) + 1) < 1e-9
            loadLabels{iL} = 'WM load 2';
        elseif abs(loadLevels(iL) - 0) < 1e-9
            loadLabels{iL} = 'WM load 4';
        elseif abs(loadLevels(iL) - 1) < 1e-9
            loadLabels{iL} = 'WM load 6';
        else
            loadLabels{iL} = sprintf('Load code %.3g', loadLevels(iL));
        end
    end

    gazeGrid = linspace(min(tbl.GazePct_c), max(tbl.GazePct_c), 120)';
    nGrid = numel(gazeGrid);
    nLevels = numel(loadLevels);
    subjDummy = tbl.Subject(1);

    figure('Position', fig_pos, 'Color', 'w');
    hold on;
    for iL = 1:nLevels
        tblPred = table( ...
            repmat(subjDummy, nGrid, 1), ...
            repmat(loadLevels(iL), nGrid, 1), ...
            gazeGrid, ...
            'VariableNames', {'Subject', 'Load', 'GazePct_c'});
        try
            yhat = predict(lme2, tblPred, 'Conditional', false);
        catch
            yhat = predict(lme2, tblPred);
        end
        plot(gazeGrid, yhat, 'LineWidth', 2, 'DisplayName', loadLabels{iL});
    end
    yline(0, 'k:', 'LineWidth', 1.5);
    grid on; box off;
    xlabel('GazePct_c (grand-mean centered gaze% [%%])');
    ylabel('Predicted Alpha [dB]');
    title('Model 2 prediction: Alpha ~ Load * GazePct_c');
    legend('Location', 'best', 'Box', 'off');
    saveas(gcf, fullfile(fig_dir, 'AOC_reciprocal_Model2_PredictedInteraction_Load_by_GazePct_c_Alpha.png'));
    close(gcf);
end

fprintf('\n=== Done. Associational models; interpret Load:moderator as moderation at fixed grand mean of the other variable. ===\n');
fprintf('Log: %s\n', cmdlog_file);

%% --- Local functions ---

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm)
    return
end
P = freq_in.powspctrm;
if size(P, 1) < 3
    return
end
P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function [lme, used_formula] = fitlme_with_random_slope_fallback(tbl_in, fixed_formula)
formula_rs = sprintf('%s + (Load|Subject)', fixed_formula);
formula_ri = sprintf('%s + (1|Subject)', fixed_formula);
try
    lme = fitlme(tbl_in, formula_rs);
    used_formula = formula_rs;
catch
    lme = fitlme(tbl_in, formula_ri);
    used_formula = formula_ri;
end
end

function [DEV2, DEV4, DEV6] = sternberg_gaze_deviation_percent_change_by_load(subjects, base_path, t_base, t_task)
% Per subject: mean across trials of percentage change 100*(gd_task/gd_baseline - 1)
% for loads 2/4/6.
screenW = 800; screenH = 600;
centreX = 400; centreY = 300;
blink_win = 50;
min_valid_samples = 100;
bounds_x = [0 screenW];
bounds_y = [0 screenH];

nSubj = numel(subjects);
DEV2 = nan(nSubj, 1);
DEV4 = nan(nSubj, 1);
DEV6 = nan(nSubj, 1);

for si = 1:nSubj
    fprintf('Gaze ET Subject %d / %d\n', si, nSubj);
    et_file = fullfile(base_path, subjects{si}, 'gaze', 'dataET_sternberg.mat');
    if ~isfile(et_file)
        et_file = fullfile(base_path, subjects{si}, 'gaze', 'dataET_sternberg');
    end
    if ~isfile(et_file)
        warning('Missing ET file for subject %s', subjects{si});
        continue
    end
    tmp = load(et_file);
    dataETlong = select_dataetlong(tmp, subjects{si});
    nTrials = size(dataETlong.trialinfo, 1);
    cond_code = dataETlong.trialinfo(:, 1) - 20;

    gd_by_cond = cell(3, 1);
    for c = 1:3
        gd_by_cond{c} = [];
    end

    for trl = 1:nTrials
        raw_dat = dataETlong.trial{trl};
        t = dataETlong.time{trl};
        if isempty(raw_dat) || isempty(t)
            continue
        end
        raw_dat = raw_dat(1:3, :);
        raw_dat(2, :) = screenH - raw_dat(2, :);

        inb = raw_dat(1, :) >= bounds_x(1) & raw_dat(1, :) <= bounds_x(2) & ...
            raw_dat(2, :) >= bounds_y(1) & raw_dat(2, :) <= bounds_y(2);
        raw_dat(:, ~inb) = NaN;

        raw_dat = remove_blinks(raw_dat, blink_win);

        idx_b = t >= t_base(1) & t <= t_base(2);
        idx_t = t >= t_task(1) & t <= t_task(2);
        dat_b = raw_dat(:, idx_b);
        dat_t = raw_dat(:, idx_t);

        ok_b = sum(all(isfinite(dat_b(1:2, :)), 1)) >= min_valid_samples;
        ok_t = sum(all(isfinite(dat_t(1:2, :)), 1)) >= min_valid_samples;
        if ~ok_b || ~ok_t
            continue
        end

        dx_b = dat_b(1, :) - centreX;
        dy_b = dat_b(2, :) - centreY;
        gaze_dev_b = nanmean(sqrt(dx_b.^2 + dy_b.^2));

        dx_t = dat_t(1, :) - centreX;
        dy_t = dat_t(2, :) - centreY;
        gaze_dev_t = nanmean(sqrt(dx_t.^2 + dy_t.^2));

        if ~(isfinite(gaze_dev_b) && gaze_dev_b > 0 && isfinite(gaze_dev_t) && gaze_dev_t > 0)
            continue
        end
        gd_pct = 100 * (gaze_dev_t / gaze_dev_b - 1);
        if ~isfinite(gd_pct)
            continue
        end

        cc = cond_code(trl);
        if cc == 2
            gd_by_cond{1}(end+1, 1) = gd_pct;
        elseif cc == 4
            gd_by_cond{2}(end+1, 1) = gd_pct;
        elseif cc == 6
            gd_by_cond{3}(end+1, 1) = gd_pct;
        end
    end

    if ~isempty(gd_by_cond{1})
        DEV2(si) = mean(gd_by_cond{1}, 'omitnan');
    end
    if ~isempty(gd_by_cond{2})
        DEV4(si) = mean(gd_by_cond{2}, 'omitnan');
    end
    if ~isempty(gd_by_cond{3})
        DEV6(si) = mean(gd_by_cond{3}, 'omitnan');
    end
end
end

function dataETlong = select_dataetlong(tmp, subj_label)
if isfield(tmp, 'dataETlong')
    dataETlong = tmp.dataETlong;
elseif isfield(tmp, 'dataet')
    dataETlong = tmp.dataet;
elseif isfield(tmp, 'dataET')
    dataETlong = tmp.dataET;
else
    error('No ET data structure found in dataET_sternberg for subject %s.', subj_label);
end
if ~isfield(dataETlong, 'trial') || ~isfield(dataETlong, 'trialinfo') || ~isfield(dataETlong, 'time')
    error('ET structure for subject %s is missing trial/trialinfo/time fields.', subj_label);
end
end
