%% AOC Split Sternberg Alpha (Subject-Level)
% Subject-level split (fixed across conditions) using merged_data_sternberg:
%   mean AlphaPower_FOOOF_bl across WM2/4/6 (baselined, FOOOFed alpha)
%   < -cutoff  -> reduction group
%   > cutoff   -> amplification group
%   |alpha| <= cutoff -> excluded (percentage-based band around zero)
%
% Cutoff = alpha_zero_pct of the alpha_ref_percentile of |alpha| (default: 5% of 95th percentile).
%
% Uses baselined+FOOOFed alpha (power spectra, TFR, topoplots) and baselined gaze metrics
% (SPL, Dev, BCEA, Vel; all % change from baseline).
%
% Outliers excluded via Tukey 1.5*IQR (per metric per condition) before visualization/analyses.
%
% Generates:
% - Alpha split inclusion figure (participants by alpha, thresholds, group assignment)
% - Power spectra (3 conditions) for both groups
% - TFRs per condition for both groups + group differences
% - Topoplots per condition for both groups + group-difference topoplots
% - Rainclouds for Alpha, SPL, Vel, Dev, BCEA
% - Time-course panels (SPL, velocity, gaze deviation) with effect-size strips

%% Setup
startup
[subjects, pathAOC, colors, headmodel] = setup('AOC');

if ispc
    base_data = 'W:\Students\Arne\AOC';
    seb_path = 'W:\Students\Arne\toolboxes\shadedErrorBar';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
    seb_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar';
end
addpath(seb_path);

feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'interactions', 'SplitAlphaAmpRed');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
fprintf('\n=== AOC Split Sternberg Alpha ===\n');
fprintf('Figure directory: %s\n', fig_dir);

% Keep canonical figure size requested.
fig_pos = [0 0 1512 982];

% Conditions and labels
cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

fontSize = 20;
rng(42);
alpha_zero_pct = 5; % Exclude near-zero alpha within this % of robust |alpha| reference
alpha_ref_percentile = 95; % Use this percentile of |alpha| for reference (outlier-resistant)

%% Load subject-level merged data and define alpha split
fprintf('\n=== Loading merged data ===\n');
merged_file = fullfile(feat_dir, 'merged_data_sternberg.mat');
if ~isfile(merged_file)
    error('Missing file: %s', merged_file);
end
S = load(merged_file, 'merged_data_sternberg');
if ~isfield(S, 'merged_data_sternberg')
    error('Variable merged_data_sternberg not found in %s', merged_file);
end
T = struct2table(S.merged_data_sternberg);

% Compute subject-level split value from full-window FOOOF-baselined alpha.
uIDs = unique(T.ID);
nSubj = numel(uIDs);
alpha_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    alpha_mean(i) = mean(T.AlphaPower_FOOOF_bl(mask), 'omitnan');
end

% Percentage-based exclusion band around zero (outlier-resistant reference).
valid_alpha_mean = alpha_mean(isfinite(alpha_mean));
if isempty(valid_alpha_mean)
    error('No finite subject-level alpha values found for split.');
end
alpha_abs_ref = prctile(abs(valid_alpha_mean), alpha_ref_percentile);
alpha_zero_margin = (alpha_zero_pct / 100) * alpha_abs_ref;

reduction_ids = uIDs(alpha_mean < -alpha_zero_margin);
amplification_ids = uIDs(alpha_mean > alpha_zero_margin);
zero_ids = uIDs(abs(alpha_mean) <= alpha_zero_margin);

fprintf('\n=== Split Summary (AlphaPower_FOOOF_bl, full window) ===\n');
fprintf('Subjects total: %d\n', nSubj);
fprintf('Zero-band setting: %.2f%% of %dth-percentile |alpha| (ref=%.4f, cutoff=%.4f)\n', ...
    alpha_zero_pct, alpha_ref_percentile, alpha_abs_ref, alpha_zero_margin);
fprintf('Reduction (< -%.4f): %d\n', alpha_zero_margin, numel(reduction_ids));
fprintf('Amplification (> %.4f): %d\n', alpha_zero_margin, numel(amplification_ids));
fprintf('Excluded (|alpha| <= %.4f): %d\n', alpha_zero_margin, numel(zero_ids));

if numel(reduction_ids) < 2 || numel(amplification_ids) < 2
    error('Insufficient subjects per split group.');
end

%% Alpha split inclusion figure
fprintf('\n=== Plotting alpha split inclusion figure ===\n');
figure('Position', fig_pos, 'Color', 'w');
hold on
% Grey yline at 0
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
% Red ylines at thresholds
yline(alpha_zero_margin, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
yline(-alpha_zero_margin, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
% Plot subjects: reduction blue, amplification red, excluded grey
x_vals = (1:nSubj)';
idx_red = ismember(uIDs, reduction_ids);
idx_amp = ismember(uIDs, amplification_ids);
idx_excl = ismember(uIDs, zero_ids);
scatter(x_vals(idx_excl), alpha_mean(idx_excl), 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
scatter(x_vals(idx_red), alpha_mean(idx_red), 80, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.8);
scatter(x_vals(idx_amp), alpha_mean(idx_amp), 80, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
xlabel('Participant (index)');
ylabel('Alpha Power [dB]');
title('Alpha Split: Reduction (blue), Amplification (red), Excluded (grey)');
legend({'Excluded', 'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize - 2);
set(gca, 'FontSize', fontSize);
box on
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlpha_inclusion.png'));
close(gcf);

%% Preallocate containers
% EEG power spectra/topography
pow_red = cell(1, 3);
pow_amp = cell(1, 3);
for c = 1:3
    pow_red{c} = {};
    pow_amp{c} = {};
end

% EEG TFR
tfr_red = cell(1, 3);
tfr_amp = cell(1, 3);
for c = 1:3
    tfr_red{c} = {};
    tfr_amp{c} = {};
end

% Gaze summary metrics from merged subject table (pre-outlier exclusion)
metrics_raw = struct();
metrics_raw.Alpha = nan(nSubj, 3);
metrics_raw.SPL = nan(nSubj, 3);           % ScanPathLengthFullBL [% change]
metrics_raw.Dev = nan(nSubj, 3);          % GazeDeviationFullBL [% change]
metrics_raw.BCEA = nan(nSubj, 3);         % BCEAFullBL [% change]
metrics_raw.Vel = nan(nSubj, 3);          % velocity [% change], computed below

% Time courses per subject x condition (full resolution, pre-outlier exclusion)
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
spl_tc_raw = nan(nSubj, 3, Tf);
dev_tc_raw = nan(nSubj, 3, Tf);
vel_tc_raw = nan(nSubj, 3, Tf);
bcea_tc_raw = nan(nSubj, 3, Tf);

missing_eeg = {};
missing_tfr = {};
missing_gaze = {};

%% Aggregate per-subject data
fprintf('\n=== Aggregating per-subject data (%d subjects) ===\n', nSubj);
for s = 1:nSubj
    disp(s)
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end
    subj_rows = T(T.ID == sid, :);

    % Subject-level condition metrics from merged_data_sternberg
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics_raw.Alpha(s, c) = mean(subj_rows.AlphaPower_FOOOF_bl(cmask), 'omitnan');
            metrics_raw.SPL(s, c) = mean(subj_rows.ScanPathLengthFullBL(cmask), 'omitnan');
            metrics_raw.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
            metrics_raw.BCEA(s, c) = mean(subj_rows.BCEAFullBL(cmask), 'omitnan');
        end
    end

    % EEG power spectra and topography sources (baselined, FOOOFed)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    try
        P = load(fullfile(eeg_dir, 'power_stern_fooof.mat'), 'pow2_fooof_bl', 'pow4_fooof_bl', 'pow6_fooof_bl');
        pow_conds = {P.pow2_fooof_bl, P.pow4_fooof_bl, P.pow6_fooof_bl};
        for c = 1:3
            if ismember(sid, reduction_ids)
                pow_red{c}{end+1} = pow_conds{c}; %#ok<AGROW>
            elseif ismember(sid, amplification_ids)
                pow_amp{c}{end+1} = pow_conds{c}; %#ok<AGROW>
            end
        end
    catch
        missing_eeg{end+1} = sid_str; %#ok<AGROW>
    end

    % EEG TFR sources
    try
        % Use FOOOFed + baselined TFR to avoid raw-power scale skew.
        R = load(fullfile(eeg_dir, 'tfr_stern.mat'), 'tfr2_fooof_bl', 'tfr4_fooof_bl', 'tfr6_fooof_bl');
        tfr_conds = {R.tfr2_fooof_bl, R.tfr4_fooof_bl, R.tfr6_fooof_bl};
        for c = 1:3
            if ismember(sid, reduction_ids)
                tfr_red{c}{end+1} = tfr_conds{c}; %#ok<AGROW>
            elseif ismember(sid, amplification_ids)
                tfr_amp{c}{end+1} = tfr_conds{c}; %#ok<AGROW>
            end
        end
    catch
        missing_tfr{end+1} = sid_str; %#ok<AGROW>
    end

    % Gaze series for time courses and velocity summary
    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', 'gaze_series_sternberg_trials.mat');
    if ~isfile(gaze_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    G = load(gaze_file);
    if ~isfield(G, 'trialinfo')
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    conds = parse_trialinfo_conds(G.trialinfo);
    if isempty(conds)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    % SPL time course from ScanPathSeries
    if isfield(G, 'ScanPathSeries') && isfield(G, 'ScanPathSeriesT')
        for c = 1:3
            tr_mask = conds == cond_codes(c);
            tr_idx = find(tr_mask);
            if isempty(tr_idx)
                continue
            end
            mat = nan(numel(tr_idx), Tf);
            for k = 1:numel(tr_idx)
                tr = tr_idx(k);
                srl = G.ScanPathSeries{tr};
                tt = G.ScanPathSeriesT{tr};
                if isempty(srl) || isempty(tt)
                    continue
                end
                try
                    mat(k, :) = interp1(tt, srl, t_plot, 'linear', NaN);
                catch
                end
            end
            spl_tc_raw(s, c, :) = nanmean(mat, 1);
        end
    end

    % Deviation and velocity time course from gaze x/y if available
    has_xy = isfield(G, 'gaze_x') && isfield(G, 'gaze_y');
    if has_xy
        for c = 1:3
            tr_mask = conds == cond_codes(c);
            tr_idx = find(tr_mask);
            if isempty(tr_idx)
                continue
            end

            dev_mat = nan(numel(tr_idx), Tf);
            vel_mat = nan(numel(tr_idx), Tf);
            bcea_mat = nan(numel(tr_idx), Tf);
            vel_full_trials = nan(numel(tr_idx), 1);
            vel_bl_trials = nan(numel(tr_idx), 1);

            for k = 1:numel(tr_idx)
                tr = tr_idx(k);
                x = double(G.gaze_x{tr});
                y = double(G.gaze_y{tr});
                if isempty(x) || isempty(y) || numel(x) ~= numel(y)
                    continue
                end

                tt = linspace(-0.5, 2, numel(x));
                dev = sqrt((x - 400).^2 + (y - 300).^2);
                [vx, vy] = compute_velocity_sg(x, y, fs, 3);
                [vx, vy] = clean_velocity_components(vx, vy);
                vel = hypot(vx, vy);

                try
                    dev_mat(k, :) = interp1(tt, dev, t_plot, 'linear', NaN);
                    vel_mat(k, :) = interp1(tt, vel, t_plot, 'linear', NaN);
                    bcea_mat(k, :) = compute_bcea_timecourse(x, y, tt, t_plot, 0.1);
                catch
                end

                idx_full = tt >= 0 & tt <= 2;
                idx_bl = tt >= -0.5 & tt <= -0.25;
                v_full = mean(vel(idx_full), 'omitnan');
                v_bl = mean(vel(idx_bl), 'omitnan');
                vel_full_trials(k) = v_full;
                if isfinite(v_full) && isfinite(v_bl) && v_bl > 0
                    vel_bl_trials(k) = 100 * (v_full - v_bl) / v_bl;
                end
            end

            dev_tc_raw(s, c, :) = nanmean(dev_mat, 1);
            vel_tc_raw(s, c, :) = nanmean(vel_mat, 1);
            bcea_tc_raw(s, c, :) = nanmean(bcea_mat, 1);
            metrics_raw.Vel(s, c) = mean(vel_bl_trials, 'omitnan');
        end
    end
end

%% Define groups in row index space
is_red = ismember(uIDs, reduction_ids);
is_amp = ismember(uIDs, amplification_ids);

%% Exclude outliers (Tukey 1.5*IQR) before visualization and analyses
fprintf('\n=== Outlier exclusion (Tukey 1.5*IQR, per metric per condition) ===\n');
[metrics, spl_tc, dev_tc, vel_tc, bcea_tc] = exclude_outliers_tukey(metrics_raw, spl_tc_raw, dev_tc_raw, vel_tc_raw, bcea_tc_raw);

%% Determine occipital channels (from first available power file)
channels = {};
for c = 1:3
    if ~isempty(pow_red{c})
        p0 = pow_red{c}{1};
        channels = occ_channels_from_labels(p0.label);
        break
    elseif ~isempty(pow_amp{c})
        p0 = pow_amp{c}{1};
        channels = occ_channels_from_labels(p0.label);
        break
    end
end
if isempty(channels)
    error('No power data available for channel definition.');
end

%% Power spectra (both groups, single figure)
close all
fprintf('\n=== Plotting power spectra ===\n');
plot_group_power_spectrum_combined(pow_red, pow_amp, channels, colors, cond_labels, ...
    fullfile(fig_dir, 'AOC_splitAlpha_powspctrm.png'), fig_pos, fontSize);

%% TFRs per condition (both groups + diff, 3x3)
fprintf('\n=== Plotting TFRs ===\n');
color_map_tfr = interp1(linspace(0,1,5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0,1,64));

plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, ...
    color_map_tfr, fig_dir, fig_pos, fontSize);

%% Topoplots per condition (both groups + differences)
fprintf('\n=== Plotting topoplots ===\n');
plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir, fig_pos, fontSize);

%% Rainclouds
fprintf('\n=== Plotting rainclouds ===\n');
plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% EEG and Gaze time courses with effect-size
close all
Tf = size(spl_tc, 3);
t_full_eeg = linspace(-0.5, 2, Tf);
addEEG_TC = true;
eeg_tc = [];
if addEEG_TC
    fprintf('\n=== Extracting EEG alpha time course from TFR ===\n');
    eeg_tc = extract_alpha_timecourse_tfr(tfr_red, tfr_amp, is_red, is_amp, channels, Tf, t_full_eeg);
end

fprintf('\n=== Plotting time courses ===\n');
plot_timecourse_with_effect_CBPT(dev_tc, is_red, is_amp, cond_labels, colors, ...
    'Gaze Deviation [px]', 'gaze_deviation', fig_dir, fig_pos, fontSize, fs, eeg_tc, addEEG_TC);
% plot_timecourse_with_effect_CBPT(spl_tc, is_red, is_amp, cond_labels, colors, ...
%     'Scan Path Length [px]', 'spl', fig_dir, fig_pos, fontSize, fs, eeg_tc, addEEG_TC);
% plot_timecourse_with_effect_CBPT(vel_tc, is_red, is_amp, cond_labels, colors, ...
%     'Eye Velocity [px/s]', 'velocity', fig_dir, fig_pos, fontSize, fs, eeg_tc, addEEG_TC);
% plot_timecourse_with_effect_CBPT(bcea_tc, is_red, is_amp, cond_labels, colors, ...
%     'BCEA [px^2]', 'bcea', fig_dir, fig_pos, fontSize, fs, eeg_tc, addEEG_TC);

%% Sanity checks
fprintf('\n=== Data Diagnostics ===\n');
fprintf('Missing EEG power (power_stern_fooof.mat): %d\n', numel(unique(missing_eeg)));
fprintf('Missing TFR files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

%% Group comparison: Reduction vs Amplification (gaze metrics)
fprintf('\n=== Parametric group comparison (reduction > amp) ===\n');
metric_names = {'SPL', 'Dev', 'BCEA', 'Vel'};
metric_labels = {'Scan path length [%%]', 'Gaze deviation [%%]', 'BCEA [%%]', 'Velocity [%%]'};

for m = 1:numel(metric_names)
    key = metric_names{m};
    X = metrics.(key);
    subj_mean = mean(X, 2, 'omitnan');

    red_vals = subj_mean(is_red);
    amp_vals = subj_mean(is_amp);
    red_vals = red_vals(isfinite(red_vals));
    amp_vals = amp_vals(isfinite(amp_vals));

    if numel(red_vals) < 2 || numel(amp_vals) < 2
        fprintf('  %s: insufficient data (n_red=%d, n_amp=%d)\n', key, numel(red_vals), numel(amp_vals));
        continue
    end

    [~, p_one, ~, stats] = ttest2(red_vals, amp_vals, 'Tail', 'right');  % H1: reduction > amp
    [~, p_two] = ttest2(red_vals, amp_vals);
    d = (mean(red_vals) - mean(amp_vals)) / sqrt(((numel(red_vals)-1)*var(red_vals) + (numel(amp_vals)-1)*var(amp_vals)) / (numel(red_vals)+numel(amp_vals)-2));
    fprintf('  %s: Red mean=%.2f, Amp mean=%.2f; diff=%.2f; one-tailed p(reduction>amp)=%.4f; two-tailed p=%.4f; d=%.3f; n_red=%d, n_amp=%d\n', ...
        key, mean(red_vals), mean(amp_vals), mean(red_vals)-mean(amp_vals), p_one, p_two, d, numel(red_vals), numel(amp_vals));
end

%% Mixed ANOVA: Group × Condition
fprintf('\n=== Mixed ANOVA (Group × Condition) on gaze metrics ===\n');
% Repeated measures: Condition (within-Ss), Group (between-Ss). Use fitrm + ranova.
incl = is_red | is_amp;
subj_idx = find(incl);
n_incl = numel(subj_idx);

for m = 1:numel(metric_names)
    key = metric_names{m};
    X = metrics.(key);
    X_incl = X(incl, :);

    % Exclude subjects with < 2 valid conditions (fitrm needs enough repeated measures)
    n_valid_per_subj = sum(isfinite(X_incl), 2);
    keep = n_valid_per_subj >= 2;
    if sum(keep) < 4
        fprintf('  %s: insufficient subjects for repeated-measures ANOVA (n=%d)\n', key, sum(keep));
        continue
    end

    X_use = X_incl(keep, :);
    grp_use = 1 + is_amp(subj_idx(keep));  % 1=reduction, 2=amplification

    % fitrm expects wide format: rows = subjects, cols = condition values
    t = array2table(X_use, 'VariableNames', {'Cond_2', 'Cond_4', 'Cond_6'});
    t.Group = categorical(grp_use, [1 2], {'Reduction', 'Amplification'});
    within = table([2; 4; 6], 'VariableNames', {'Load'});

    try
        rm = fitrm(t, 'Cond_2-Cond_6~Group', 'WithinDesign', within);
        % Between-Ss (Group): from anova(rm), which uses Within/Between columns
        anovatbl = anova(rm);
        idx_grp = strcmp(string(anovatbl.Between), 'Group');
        p_grp = anovatbl.pValue(idx_grp);
        if isempty(p_grp), p_grp = NaN; elseif numel(p_grp) > 1, p_grp = p_grp(1); end
        % Within-Ss (Load, Load:Group): from ranova(rm), row names like (Intercept):Load, Group:Load
        ranovatbl = ranova(rm);
        rn = ranovatbl.Properties.RowNames;
        p_cond = ranovatbl.pValue(1);  % (Intercept):Load = main effect of Load
        idx_int = cellfun(@(x) contains(x, 'Group') && contains(x, ':'), rn);
        p_int = ranovatbl.pValue(idx_int);
        if isempty(p_int), p_int = NaN; elseif numel(p_int) > 1, p_int = p_int(1); end
        fprintf('  %s: Condition (Load) p=%.4f; Group p=%.4f; Condition×Group p=%.4f; n=%d\n', ...
            key, p_cond, p_grp, p_int, sum(keep));
    catch ME
        fprintf('  %s: Mixed ANOVA failed (%s)\n', key, ME.message);
    end
end

%% Linear mixed models (LMM): Group × Condition + (1|Subject)
fprintf('\n=== Linear mixed models (Value ~ Group * Condition + (1|Subject)) ===\n');
% Long-format table: Subject, Condition, Group, Value. Random intercept by subject.
incl = is_red | is_amp;
subj_idx = find(incl);

for m = 1:numel(metric_names)
    key = metric_names{m};
    X = metrics.(key);
    X_incl = X(incl, :);

    % Stack to long format
    rows = [];
    for s = 1:numel(subj_idx)
        for c = 1:3
            val = X_incl(s, c);
            if isfinite(val)
                grp = 1 + is_amp(subj_idx(s));
                rows(end+1, :) = [subj_idx(s), c, grp, val]; %#ok<AGROW>
            end
        end
    end

    if size(rows, 1) < 15
        fprintf('  %s: insufficient data for LMM (n_obs=%d)\n', key, size(rows, 1));
        continue
    end

    tbl = array2table(rows, 'VariableNames', {'Subject', 'Condition', 'Group', 'Value'});
    tbl.Subject = categorical(tbl.Subject);
    tbl.Condition = categorical(tbl.Condition, [1 2 3], cond_labels);
    tbl.Group = categorical(tbl.Group, [1 2], {'Reduction', 'Amplification'});

    try
        lme = fitlme(tbl, 'Value ~ Group * Condition + (1|Subject)');
        aov = anova(lme);
        terms = aov.Term;
        get_p = @(t) aov.pValue(strcmp(terms, t));
        p_grp = get_p('Group');
        p_cond = get_p('Condition');
        p_int = get_p('Group:Condition');
        if isempty(p_grp), p_grp = NaN; end
        if isempty(p_cond), p_cond = NaN; end
        if isempty(p_int), p_int = NaN; end
        fprintf('  %s: Group p=%.4f; Condition p=%.4f; Group×Condition p=%.4f; n_obs=%d, n_subj=%d\n', ...
            key, p_grp, p_cond, p_int, size(rows, 1), numel(unique(tbl.Subject)));
    catch ME
        fprintf('  %s: LMM failed (%s)\n', key, ME.message);
    end
end

%% ========================= Local Functions =========================
function [metrics_out, spl_tc_out, dev_tc_out, vel_tc_out, bcea_tc_out] = exclude_outliers_tukey(metrics_in, spl_tc_in, dev_tc_in, vel_tc_in, bcea_tc_in)
% Apply Tukey 1.5*IQR rule per metric per condition. Set outliers to NaN.
% Also excludes corresponding time-course data for SPL, Dev, Vel, BCEA.
metrics_out = metrics_in;
spl_tc_out = spl_tc_in;
dev_tc_out = dev_tc_in;
vel_tc_out = vel_tc_in;
bcea_tc_out = bcea_tc_in;

metric_fields = {'Alpha', 'SPL', 'Vel', 'Dev', 'BCEA'};
tc_map = containers.Map({'SPL', 'Dev', 'Vel', 'BCEA'}, {1, 2, 3, 4});  % index for spl/dev/vel/bcea
tc_cells = {spl_tc_out, dev_tc_out, vel_tc_out, bcea_tc_out};

for m = 1:numel(metric_fields)
    key = metric_fields{m};
    X = metrics_out.(key);
    nOut = 0;
    for c = 1:3
        vals = X(:, c);
        vals_f = vals(isfinite(vals));
        if numel(vals_f) < 4
            continue
        end
        q1 = prctile(vals_f, 25);
        q3 = prctile(vals_f, 75);
        iqr_val = q3 - q1;
        if iqr_val <= 0
            continue
        end
        lo = q1 - 1.5 * iqr_val;
        hi = q3 + 1.5 * iqr_val;
        out_mask = (vals < lo) | (vals > hi);
        nOut = nOut + sum(out_mask);
        X(out_mask, c) = NaN;
        % Exclude corresponding time course for gaze metrics
        if isKey(tc_map, key)
            tc_idx = tc_map(key);
            tc = tc_cells{tc_idx};
            tc(out_mask, c, :) = NaN;
            tc_cells{tc_idx} = tc;
        end
    end
    metrics_out.(key) = X;
    if nOut > 0
        fprintf('  %s: %d outlier(s) excluded\n', key, nOut);
    end
end
spl_tc_out = tc_cells{1};
dev_tc_out = tc_cells{2};
vel_tc_out = tc_cells{3};
bcea_tc_out = tc_cells{4};
end

function conds = parse_trialinfo_conds(trialinfo)
conds = [];
if isempty(trialinfo)
    return
end
if isvector(trialinfo)
    conds = trialinfo(:);
elseif size(trialinfo, 2) >= 1
    if size(trialinfo, 2) == 2
        conds = trialinfo(:, 1);
    elseif size(trialinfo, 1) == 2
        conds = trialinfo(1, :)';
    else
        conds = trialinfo(:, 1);
    end
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'}) || contains(lab, {'PO'})
        ch{end+1} = lab; %#ok<AGROW>
    end
end
if isempty(ch)
    ch = labels;
end
end

function subj_folder = resolve_subject_folder(subjects, sid)
subj_folder = '';
for i = 1:numel(subjects)
    sval = str2double(subjects{i});
    if isfinite(sval) && sval == sid
        subj_folder = subjects{i};
        return
    end
end
end

function [vx, vy] = clean_velocity_components(vx, vy)
halfwin = 11; % SG edge region for framelen=21
if numel(vx) > 2*halfwin
    vx(1:halfwin) = NaN;
    vx(end-halfwin+1:end) = NaN;
    vy(1:halfwin) = NaN;
    vy(end-halfwin+1:end) = NaN;
end
zvx = (vx - nanmean(vx)) ./ (nanstd(vx) + eps);
zvy = (vy - nanmean(vy)) ./ (nanstd(vy) + eps);
bad = abs(zvx) > 4 | abs(zvy) > 4;
vx(bad) = NaN;
vy(bad) = NaN;
end

function bcea_t = compute_bcea_timecourse(x, y, tt, t_plot, win_half)
% BCEA time course: 95% bivariate contour ellipse area in sliding windows.
% x, y: gaze coordinates; tt: time vector; t_plot: output time points;
% win_half: half-window in s (e.g. 0.1 for 200 ms total).
% Returns row vector of BCEA [px²] at each t_plot (NaN if insufficient samples).
k95 = -log(1 - 0.95);
nT = numel(t_plot);
bcea_t = nan(1, nT);
for j = 1:nT
    t0 = t_plot(j);
    idx = tt >= t0 - win_half & tt <= t0 + win_half;
    xw = x(idx);
    yw = y(idx);
    valid = isfinite(xw) & isfinite(yw);
    xw = xw(valid);
    yw = yw(valid);
    if numel(xw) < 10
        continue
    end
    sx = std(xw);
    sy = std(yw);
    if sx <= 0 || sy <= 0
        continue
    end
    r = corrcoef(xw, yw);
    rho = r(1, 2);
    if ~isfinite(rho)
        continue
    end
    rho2 = max(0, min(1, rho^2));  % guard for numerical issues
    bcea_t(j) = 2 * k95 * pi * sx * sy * sqrt(1 - rho2);
end
end

function plot_group_power_spectrum_combined(pow_red, pow_amp, channels, colors, cond_labels, out_file, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping power spectrum (incomplete condition data).');
    return
end
ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], pow_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], pow_amp{c}{:});
end

elecs = ismember(ga_red{1}.label, channels);
freqs = ga_red{1}.freq;

% Compute spectra and determine shared ylim (FOOOFed data in dB)
% Restrict to 5-30 Hz for ylim (matches xlim)
freq_idx_plot = freqs >= 5 & freqs <= 30;
all_m = [];
all_se = [];
pw = {pow_red, pow_amp};
for grp = 1:2
    pow_cells = pw{grp};
    for c = 1:3
        subj_n = numel(pow_cells{c});
        subj_spec = nan(subj_n, numel(freqs));
        for s = 1:subj_n
            p = pow_cells{c}{s};
            subj_spec(s, :) = mean(p.powspctrm(elecs, :), 1, 'omitnan');
        end
        m = mean(subj_spec, 1, 'omitnan');
        se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
        all_m = [all_m; m]; all_se = [all_se; se]; %#ok<AGROW>
    end
end
% Absolute max centered around 0: ylim = [-ymax, ymax] (5-30 Hz only)
m_plot = all_m(:, freq_idx_plot);
se_plot = all_se(:, freq_idx_plot);
all_vals = [reshape(m_plot, [], 1); reshape(m_plot - se_plot, [], 1); reshape(m_plot + se_plot, [], 1)];
all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals)
    ymax_abs = 1;
else
    ymax_abs = max(abs(all_vals));
    ymax_abs = max(ymax_abs, 0.1);
end
ylim_shared = [-ymax_abs, ymax_abs];

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

pw2 = {pow_red, pow_amp};
ttls = {'Reduction', 'Amplification'};
for grp = 1:2
    nexttile; hold on
    pow_cells = pw2{grp};
    ttl = ttls{grp};
    for c = 1:3
        subj_n = numel(pow_cells{c});
        subj_spec = nan(subj_n, numel(freqs));
        for s = 1:subj_n
            p = pow_cells{c}{s};
            subj_spec(s, :) = mean(p.powspctrm(elecs, :), 1, 'omitnan');
        end
        m = mean(subj_spec, 1, 'omitnan');
        se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
        eb(c) = shadedErrorBar(freqs, m, se, 'lineProps', {'-'}, 'transparent', true);
        set(eb(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
        set(eb(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb(c).edge(1), 'Color', colors(c, :));
        set(eb(c).edge(2), 'Color', colors(c, :));
    end
    xlim([5 30]);
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ylim(ylim_shared);
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    title(ttl, 'FontSize', fsz + 4);
    legend([eb.mainLine], cond_labels, 'FontSize', fsz - 2, 'Location', 'best');
    set(gca, 'FontSize', fsz);
    box on
end
saveas(gcf, out_file);
close(gcf);
end

function plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, color_map, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping TFR plot (incomplete condition data).');
    return
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], tfr_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], tfr_amp{c}{:});
end

% Shared clim for Reduction and Amplification
[~, ch_idx] = ismember(channels, ga_red{1}.label);
freq_idx = ga_red{1}.freq >= 5 & ga_red{1}.freq <= 30;
time_idx = ga_red{1}.time >= -0.5 & ga_red{1}.time <= 2;
mx = 0;
for c = 1:3
    Ared = squeeze(mean(ga_red{c}.powspctrm(ch_idx, :, :), 1, 'omitnan'));
    Aamp = squeeze(mean(ga_amp{c}.powspctrm(ch_idx, :, :), 1, 'omitnan'));
    mx = max(mx, max(abs(Ared(freq_idx, time_idx)), [], 'all'));
    mx = max(mx, max(abs(Aamp(freq_idx, time_idx)), [], 'all'));
end
clim_abs = [-mx mx];

cfg = [];
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;

% Single 3x2 figure: column 1 = Reduction, column 2 = Amplification
fsz_tfr = round(fsz * 0.8);  % ~20% reduction
figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    % Column 1: Reduction (rows 1–3)
    ax = subplot(3, 2, (c-1)*2 + 1);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_red{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    colorbar(ax);
    xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
    set(ax, 'FontSize', fsz_tfr);
    title(ax, sprintf('Reduction - %s', cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');

    % Column 2: Amplification (rows 1–3)
    ax = subplot(3, 2, (c-1)*2 + 2);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_amp{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    colorbar(ax);
    xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
    set(ax, 'FontSize', fsz_tfr);
    title(ax, sprintf('Amplification - %s', cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');
end
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlpha_tfr_all.png'));
close(gcf);
end

function plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping topoplots (missing power data).');
    return
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], pow_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], pow_amp{c}{:});
end

cfg = [];
cfg.layout = headmodel.layANThead;
allch = cfg.layout.label;
cfg.channel = allch(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.figure = 'gcf';
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];

cfg.colormap = rdbu_cmap(64);  % Blue-white-red diverging
freq_idx = ga_red{1}.freq >= 8 & ga_red{1}.freq <= 14;
all_alpha = [];
for c = 1:3
    Ared = mean(ga_red{c}.powspctrm(:, freq_idx), 2, 'omitnan');
    Aamp = mean(ga_amp{c}.powspctrm(:, freq_idx), 2, 'omitnan');
    all_alpha = [all_alpha; Ared(:); Aamp(:)]; %#ok<AGROW>
end
all_alpha = all_alpha(isfinite(all_alpha));
if isempty(all_alpha)
    cfg.zlim = 'maxabs';
else
    mx = prctile(abs(all_alpha), 95);
    if mx <= 0
        mx = max(abs(all_alpha));
    end
    cfg.zlim = [-mx mx];  % Symmetric for blue-white-red diverging
end

% Single 2x3 figure: row 1 = Reduction, row 2 = Amplification
figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    % Row 1: Reduction (cols 1–3)
    ax = subplot(2, 3, c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_red{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Reduction - %s', cond_labels{c}), 'Interpreter', 'none');

    % Row 2: Amplification (cols 1–3)
    ax = subplot(2, 3, 3 + c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_amp{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Amplification - %s', cond_labels{c}), 'Interpreter', 'none');
end
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlpha_topo_all.png'));
close(gcf);
end

function plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fsz)
metric_defs = { ...
    'Alpha', 'AlphaPower_FOOOF_bl', false; ...
    'SPL', 'Scan path length', true; ...
    'Vel', 'Velocity', true; ...
    'Dev', 'Gaze deviation', true; ...
    'BCEA', 'BCEA', true; ...
    };

for m = 1:size(metric_defs, 1)
    key = metric_defs{m, 1};
    varname = metric_defs{m, 2};
    is_pct = metric_defs{m, 3};
    X = metrics.(key);

    % Shared ymax for both groups, centered around 0; add margin for density curve
    all_vals = X(isfinite(X));
    if isempty(all_vals)
        ymax = 1;
    else
        ymax = max(abs(all_vals));
        if ymax <= 0
            ymax = 1;
        end
        ymax = ymax * 1.15;  % accommodate full density curve
    end
    ylim_shared = [-ymax ymax];

    if is_pct
        ylab = [varname ' [%]'];
    else
        ylab = varname;
    end

    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    nexttile; hold on
    for c = 1:3
        draw_one_cloud(X(is_red, c), c, colors(c,:), 0.3, 96, 0.45);
    end
    yline(0, '--', 'Color', [0.4 0.4 0.4]);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title('Reduction', 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    box off

    nexttile; hold on
    for c = 1:3
        draw_one_cloud(X(is_amp, c), c, colors(c,:), 0.3, 96, 0.45);
    end
    yline(0, '--', 'Color', [0.4 0.4 0.4]);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title('Amplification', 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    box off

    sgtitle(varname, 'FontSize', fsz+2, 'Interpreter', 'none');
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_raincloud_%s.png', lower(key))));
    close(gcf);
end
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha)
y = yvals(isfinite(yvals));
if numel(y) < 3
    return
end
[f, xi] = ksdensity(y, 'NumPoints', 120);
f = f / max(f) * 0.35;
fill([xpos - f, fliplr(repmat(xpos, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25);
q3 = prctile(y, 75);
med = median(y);
p5 = prctile(y, 5);
p95 = prctile(y, 95);
plot([xpos xpos], [p5 q1], '-k', 'LineWidth', 1.2);
plot([xpos xpos], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [xpos-box_w/2, q1, box_w, q3-q1], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(xpos + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = 0.10*(rand(numel(y),1)-0.5);
scatter(xpos + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

% Cluster-based permutation test
function plot_timecourse_with_effect_CBPT(tc, is_red, is_amp, cond_labels, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, eeg_tc, addEEG_TC)
if nargin < 12, eeg_tc = []; end
if nargin < 13, addEEG_TC = ~isempty(eeg_tc); end
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 2, nT);
tc(~isfinite(tc)) = NaN;
Rall = squeeze(mean(tc(is_red, :, :), 2, 'omitnan'));
Aall = squeeze(mean(tc(is_amp, :, :), 2, 'omitnan'));
win_sm = max(1, round(0.05 * fs));
if win_sm > 1
    Rall = movmean(Rall, win_sm, 2, 'omitnan');
    Aall = movmean(Aall, win_sm, 2, 'omitnan');
end
figure('Position', fig_pos, 'Color', 'w');
nR = size(Rall, 1);
nA = size(Aall, 1);
has_eeg = addEEG_TC && ~isempty(eeg_tc) && size(eeg_tc, 1) >= (nR + nA) && size(eeg_tc, 2) >= 1 && size(eeg_tc, 3) >= nT;
if has_eeg
    tiledlayout(4, 1, 'TileSpacing', 'compact');
else
    tiledlayout(3, 1, 'TileSpacing', 'compact');
end
nexttile([2 1]); hold on
mR = mean(Rall, 1, 'omitnan');
mA = mean(Aall, 1, 'omitnan');
nR_fin = sum(isfinite(Rall), 1);
nA_fin = sum(isfinite(Aall), 1);
sR = std(Rall, 0, 1, 'omitnan') ./ max(sqrt(nR_fin), 1);
sA = std(Aall, 0, 1, 'omitnan') ./ max(sqrt(nA_fin), 1);
sR(~isfinite(sR)) = NaN;
sA(~isfinite(sA)) = NaN;
e1 = shadedErrorBar(t_plot, mR, sR, 'lineProps', {'-'}, 'transparent', true);
e2 = shadedErrorBar(t_plot, mA, sA, 'lineProps', {'-'}, 'transparent', true);
set(e1.mainLine, 'Color', colors(1,:), 'LineWidth', 3.5);
set(e2.mainLine, 'Color', colors(3,:), 'LineWidth', 3.5);
set(e1.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.25);
set(e2.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.25);
xline(0, '--k');
ylabel(ylab);
xlim([-0.5 2]);
box on
set(gca, 'FontSize', fsz-4);
% Legend with colored patch boxes (clearer than thin lines)
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1 leg_p2], {'Reduction', 'Amplification'}, 'Location', 'northeast', 'FontSize', fsz-2);
if strcmp(save_tag, 'gaze_deviation')
    yl = ylim;
    ylim([max(0, yl(1)) yl(2)]);
end
nexttile; hold on
n_perm = 1000;
min_per_group = 3;
d = nan(1, nT);
for t = 1:nT
    x = Rall(:, t); y = Aall(:, t);
    x = x(isfinite(x)); y = y(isfinite(y));
    if numel(x) < min_per_group || numel(y) < min_per_group, continue, end
    sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / max(numel(x)+numel(y)-2, 1));
    d(t) = (mean(y) - mean(x)) / max(sp, eps);
end
ds_factor = 10;
Rall_ds = Rall(:, 1:ds_factor:end);
Aall_ds = Aall(:, 1:ds_factor:end);
nT_ds = size(Rall_ds, 2);
t_plot_ds = t_plot(1:ds_factor:end);
dt_ds = ds_factor * dt;
[clusters, tvals_cl, thr, maxMassNull, maxExtentNull] = ft_cluster_permutation_1d(Rall_ds, Aall_ds, n_perm, 0.05, 'onetail_pos', t_plot_ds);
nAbove = sum(tvals_cl > thr.tcrit & isfinite(tvals_cl));
maxClMass = max([0, arrayfun(@(k) clusters(k).mass, 1:numel(clusters))]);
maxClExtent = max([0, arrayfun(@(k) clusters(k).extent, 1:numel(clusters))]);
fprintf('  [%s] n_red=%d n_amp=%d tcrit=%.2f (one-tailed Red>Amp) t>tcrit at %d timepts; max cluster mass=%.1f; max cluster extent=%d\n', ...
    save_tag, nR, nA, thr.tcrit, nAbove, maxClMass, maxClExtent);
if ~isempty(maxMassNull)
    fprintf('    CBPT gaze: nT_ds=%d (%.0f ms bins); null mass: median=%.1f, 90th=%.1f, 95th=%.1f; null extent: median=%d, 90th=%d, 95th=%d\n', ...
        nT_ds, ds_factor*1000/fs, median(maxMassNull), prctile(maxMassNull, 90), prctile(maxMassNull, 95), ...
        median(maxExtentNull), prctile(maxExtentNull, 90), prctile(maxExtentNull, 95));
else
    fprintf('    CBPT gaze: nT_ds=%d (%.0f ms bins); FieldTrip ft_timelockstatistics (cluster mass)\n', nT_ds, ds_factor*1000/fs);
end
if ~isempty(clusters)
    for k = 1:numel(clusters)
        t_lo = t_plot_ds(clusters(k).idx(1));
        t_hi = t_plot_ds(clusters(k).idx(end));
        pass_p = clusters(k).p < 0.05;
        status = 'n.s.';
        if pass_p
            status = 'SIGNIFICANT';
        end
        fprintf('    Cluster %d: t=[%.2f, %.2f]s, extent=%d pts, mass=%.1f; p=%.4f; %s\n', ...
            k, t_lo, t_hi, clusters(k).extent, clusters(k).mass, clusters(k).p, status);
    end
    if isfinite(thr.mass)
        fprintf('    Gap to significance: mass need +%.1f (have %.1f, need %.1f); extent need +%d (have %d, need %d)\n', ...
            thr.mass - maxClMass, maxClMass, thr.mass, thr.extent - maxClExtent, maxClExtent, thr.extent);
    end
else
    t_fin = tvals_cl(isfinite(tvals_cl));
    if ~isempty(t_fin)
        [t_max, idx_max] = max(tvals_cl);
        fprintf('    No clusters formed (no contiguous t>tcrit); largest t=%.2f at t=%.2fs\n', ...
            t_max, t_plot_ds(idx_max));
    else
        fprintf('    No clusters; no valid t-values\n');
    end
end
sig_cluster = false(1, nT_ds);
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        sig_cluster(clusters(k).idx) = true;
    end
end
sig_uncorr = (tvals_cl > thr.tcrit) & isfinite(tvals_cl);
sig = sig_cluster;
if ~any(sig) && any(sig_uncorr)
    sig = sig_uncorr;
    fprintf('  [%s] (cluster n.s.; shading uncorrected t>tcrit)\n', save_tag);
end
d_ds = d(1:ds_factor:end);
if any(sig_cluster)
    sig = sig & isfinite(d_ds);
end
run_start = [false, diff(sig) == 1];
run_end = [diff(sig) == -1, false];
if sig(1), run_start(1) = true; end
if sig(end), run_end(end) = true; end
starts = find(run_start);
ends = find(run_end);
d_fin = d(isfinite(d));
mx = max(abs(d_fin), [], 'omitnan');
if isempty(d_fin) || ~isfinite(mx) || mx == 0
    ylims = [-0.6 0.6];
else
    ylims = [-max(mx + 0.1, 0.6), max(mx + 0.1, 0.6)];
end
patch_alpha = 0.4 * ~any(sig_cluster) + 0.25 * any(sig_cluster);
for k = 1:numel(starts)
    t1 = max(0, t_plot_ds(starts(k)) - dt_ds/2);
    t2 = t_plot_ds(ends(k)) + dt_ds/2;
    patch([t1 t2 t2 t1], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.5 0.5 0.5], 'FaceAlpha', patch_alpha, 'EdgeColor', 'none');
end
ylim(ylims);
plot(t_plot, d, 'k-', 'LineWidth', 3.5);
yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Cohen''s d');
xlim([-0.5 2]);
box on
set(gca, 'FontSize', fsz-4);
if ~any(sig_cluster) && any(sig_uncorr)
    text(0.75, ylims(2) - 0.08*diff(ylims), 'WARNING: No significant clusters; shading shows uncorrected t>t_{crit}', ...
        'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'HorizontalAlignment', 'center');
end

% Alpha power time course panel (4th row when has_eeg) - two lines with shaded error bars
if has_eeg
    EegR = squeeze(mean(eeg_tc(is_red, :, :), 2, 'omitnan'));
    EegA = squeeze(mean(eeg_tc(is_amp, :, :), 2, 'omitnan'));
    if win_sm > 1
        EegR = movmean(EegR, win_sm, 2, 'omitnan');
        EegA = movmean(EegA, win_sm, 2, 'omitnan');
    end
    mR_eeg = mean(EegR, 1, 'omitnan');
    mA_eeg = mean(EegA, 1, 'omitnan');
    nR_fin_eeg = sum(isfinite(EegR), 1);
    nA_fin_eeg = sum(isfinite(EegA), 1);
    sR_eeg = std(EegR, 0, 1, 'omitnan') ./ max(sqrt(nR_fin_eeg), 1);
    sA_eeg = std(EegA, 0, 1, 'omitnan') ./ max(sqrt(nA_fin_eeg), 1);
    sR_eeg(~isfinite(sR_eeg)) = NaN;
    sA_eeg(~isfinite(sA_eeg)) = NaN;
    nexttile; hold on
    ebR = shadedErrorBar(t_plot, mR_eeg, sR_eeg, 'lineProps', {'-'}, 'transparent', true);
    ebA = shadedErrorBar(t_plot, mA_eeg, sA_eeg, 'lineProps', {'-'}, 'transparent', true);
    set(ebR.mainLine, 'Color', colors(1,:), 'LineWidth', 3.5);
    set(ebA.mainLine, 'Color', colors(3,:), 'LineWidth', 3.5);
    set(ebR.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.25);
    set(ebA.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.25);
    yline(0, '--');
    xline(0, '--k');
    all_eeg_vals = [mR_eeg(:); mA_eeg(:); mR_eeg(:)-sR_eeg(:); mR_eeg(:)+sR_eeg(:); mA_eeg(:)-sA_eeg(:); mA_eeg(:)+sA_eeg(:)];
    all_eeg_vals = all_eeg_vals(isfinite(all_eeg_vals));
    ymax_eeg = max(abs(all_eeg_vals), [], 'omitnan') * 1.1;
    if isempty(all_eeg_vals) || ~isfinite(ymax_eeg) || ymax_eeg == 0
        ymax_eeg = 0.1;
    end
    ylim([-ymax_eeg ymax_eeg]);
    ylabel('Alpha Power [dB]');
    xlabel('Time [s]');
    xlim([-0.5 2]);
    box on
    set(gca, 'FontSize', fsz-4);
    leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
    leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
end

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_timecourse_%s_CBPT.png', save_tag)));
end

function [clusters, tvals, thr, maxMassNull, maxExtentNull] = ft_cluster_permutation_1d(Rall, Aall, nPerm, alpha, tail, t_plot_ds)
% FT_CLUSTER_PERMUTATION_1D  Cluster-based permutation test using FieldTrip.
%
% Uses ft_timelockstatistics (Maris & Oostenveld, 2007) for two-sample
% independent-samples cluster permutation on 1D time-series data.
%
% Inputs:
%   Rall      - nR x nT matrix (group 1, e.g. Reduction)
%   Aall      - nA x nT matrix (group 2, e.g. Amplification)
%   nPerm     - number of permutations (e.g. 2000)
%   alpha     - significance level (e.g. 0.05)
%   tail      - 'onetail_pos' (H1: group1 > group2, i.e. t > 0; use for Red > Amp)
%               'onetail_neg' (H1: group1 < group2, i.e. t < 0)
%               or 'twotail' (H1: group1 ~= group2)
%   t_plot_ds - (optional) 1 x nT actual time vector. If provided, uses cfg.latency=[0 2]
%               to restrict clustering to post-stimulus only (avoids baseline clusters).
%
% Outputs (compatible with previous cluster_permutation_2sample_1d):
%   clusters     - struct array with .idx, .mass, .extent, .p, .p_extent
%   tvals        - 1 x nT t-values (group1 - group2; FT indepsamplesT)
%   thr          - struct with .tcrit, .mass, .extent
%   maxMassNull  - []; FieldTrip does not expose null distribution
%   maxExtentNull - []
%
% Requires: FieldTrip (startup must add it to path).
%
% Reference: Maris & Oostenveld (2007). J Neurosci Methods 164:177-190.

if nargin < 5
    tail = 'twotail';
end
if nargin < 6
    t_plot_ds = [];
end

nR = size(Rall, 1);
nA = size(Aall, 1);
nT = size(Rall, 2);
df = nR + nA - 2;

% Replace NaN with condition mean at that time (conservative imputation for FieldTrip)
for t = 1:nT
    r = Rall(:, t);
    a = Aall(:, t);
    mr = mean(r(isfinite(r)), 'omitnan');
    ma = mean(a(isfinite(a)), 'omitnan');
    if ~isfinite(mr), mr = 0; end
    if ~isfinite(ma), ma = 0; end
    r(~isfinite(r)) = mr;
    a(~isfinite(a)) = ma;
    Rall(:, t) = r;
    Aall(:, t) = a;
end

% Build FieldTrip timelock structures (rpt_chan_time)
chan_label = 'metric';
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    time_vec = t_plot_ds;
else
    time_vec = (0:nT-1) / 500;
end

tl1 = struct();
tl1.label = {chan_label};
tl1.time = time_vec;
tl1.dimord = 'rpt_chan_time';
tl1.trial = reshape(Rall, [nR, 1, nT]);  % nR x 1 x nT

tl2 = struct();
tl2.label = {chan_label};
tl2.time = time_vec;
tl2.dimord = 'rpt_chan_time';
tl2.trial = reshape(Aall, [nA, 1, nT]);

% Neighbours: single channel, no spatial neighbours (temporal clustering only)
cfg_neigh = struct();
cfg_neigh(1).label = chan_label;
cfg_neigh(1).neighblabel = {};

% Statistics config
cfg = struct();
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = cfg_neigh;
cfg.numrandomization = nPerm;
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    cfg.latency = [0 2];  % restrict to post-stimulus only
else
    cfg.latency = 'all';
end
cfg.channel = chan_label;

if strcmpi(tail, 'onetail_pos')
    % H1: group1 > group2  =>  t = (mean1 - mean2)/se > 0
    cfg.tail = 1;
    cfg.clustertail = 1;
    cfg.alpha = alpha;  % one-tailed
elseif strcmpi(tail, 'onetail_neg')
    % H1: group1 < group2  =>  t = (mean1 - mean2)/se < 0
    cfg.tail = -1;
    cfg.clustertail = -1;
    cfg.alpha = alpha;  % one-tailed
else
    cfg.tail = 0;
    cfg.clustertail = 0;
    cfg.alpha = alpha / 2;  % two-tailed: split alpha
end

cfg.design = [ones(1, nR), 2*ones(1, nA)];
cfg.ivar = 1;

% Run FieldTrip cluster permutation
stat = ft_timelockstatistics(cfg, tl1, tl2);

% Extract t-values (stat is chan x time)
tvals = stat.stat(1, :);

% When we used cfg.latency=[0 2], stat has fewer time points; map back to full range
post_idx = [];
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    post_idx = find(t_plot_ds >= 0);
    if size(tvals, 2) < nT && ~isempty(post_idx)
        tvals_full = nan(1, nT);
        n_sel = min(size(tvals, 2), numel(post_idx));
        tvals_full(post_idx(1:n_sel)) = tvals(1:n_sel);
        tvals = tvals_full;
    end
end

% Map to our cluster format (FieldTrip uses chan x time; we have 1 chan)
clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {}, 'p_extent', {});

% FieldTrip may omit negclusters/posclusters when no clusters found
has_neg = isfield(stat, 'negclusters') && ~isempty(stat.negclusters);
has_pos = isfield(stat, 'posclusters') && ~isempty(stat.posclusters);

if strcmpi(tail, 'onetail_pos')
    if has_pos
        clist = stat.posclusters;
        lmat = stat.posclusterslabelmat;
    else
        clist = struct('prob', {});
        lmat = zeros(1, nT);
    end
elseif strcmpi(tail, 'onetail_neg')
    if has_neg
        clist = stat.negclusters;
        lmat = stat.negclusterslabelmat;
    else
        clist = struct('prob', {});
        lmat = zeros(1, nT);
    end
else
    % Two-tailed: collect from both pos and neg (use FT structs directly to match fields)
    clist = [];
    lmat = zeros(1, nT);
    if has_pos
        clist = stat.posclusters(:);
        lmat = stat.posclusterslabelmat;
    end
    if has_neg
        lmat_neg = stat.negclusterslabelmat;
        npos = numel(clist);
        lmat(lmat_neg > 0) = npos + lmat_neg(lmat_neg > 0);
        if isempty(clist)
            clist = stat.negclusters(:);
        else
            clist = [clist; stat.negclusters(:)];
        end
    end
end

for k = 1:numel(clist)
    idx = find(lmat(1, :) == k);
    if isempty(idx)
        continue
    end
    idx = idx(1):idx(end);  % ensure contiguous
    if ~isempty(post_idx)
        valid = idx <= numel(post_idx);
        idx = post_idx(idx(valid));  % map to full time indices
        if isempty(idx), continue; end
    end
    clusters(end+1).idx = idx; %#ok<AGROW>
    clusters(end).mass = sum(abs(tvals(idx)), 'omitnan');
    clusters(end).extent = numel(idx);
    clusters(end).p = clist(k).prob;
    clusters(end).p_extent = clist(k).prob;
end

% Threshold for uncorrected shading
if strcmpi(tail, 'onetail_pos') || strcmpi(tail, 'onetail_neg')
    thr.tcrit = tinv(1 - alpha, df);  % one-tailed upper critical value
else
    thr.tcrit = tinv(1 - alpha/2, df);  % two-tailed
end
% FT does not expose explicit mass/extent thresholds
thr.mass = NaN;
thr.extent = NaN;

maxMassNull = [];
maxExtentNull = [];
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
Ts = 1 / fs;
L = numel(X);
framelen = min(21, L);
if mod(framelen, 2) == 0
    framelen = framelen - 1;
end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0
    minLegal = minLegal + 1;
end
if framelen < minLegal
    framelen = minLegal;
end
if framelen > L
    framelen = L - mod(L, 2) + 1;
end
useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / Ts) * G(:, 2)';
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end

function cmap = rdbu_cmap(n)
% RDBU_CMAP - Red-Blue diverging colormap (no external deps, ColorBrewer RdBu-like).
% Usage: cmap = rdbu_cmap(n);  returns n x 3 RGB matrix.
rdbu_11 = [33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; ...
    253 219 199; 244 165 130; 214 96 77; 178 24 43] / 255;
x = linspace(0, 1, size(rdbu_11, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, rdbu_11, xi, 'linear');
end

function eeg_tc = extract_alpha_timecourse_tfr(tfr_red, tfr_amp, is_red, is_amp, channels, Tf, t_target)
% Extract alpha power (8-14 Hz) over occipital channels from TFR, per subject x condition.
% Returns eeg_tc: nSubj x 3 x Tf.
nSubj = numel(is_red);
eeg_tc = nan(nSubj, 3, Tf);
red_list = find(is_red);
amp_list = find(is_amp);
alpha_freq = [8 14];
for c = 1:3
    for j = 1:numel(red_list)
        if j > numel(tfr_red{c}), continue, end
        T = tfr_red{c}{j};
        if ~isfield(T, 'powspctrm') || ~isfield(T, 'label'), continue, end
        ch_idx = find(ismember(T.label, channels));
        freq_idx = T.freq >= alpha_freq(1) & T.freq <= alpha_freq(2);
        if isempty(ch_idx) || ~any(freq_idx), continue, end
        alpha_raw = squeeze(mean(mean(T.powspctrm(ch_idx, freq_idx, :), 1, 'omitnan'), 2, 'omitnan'));
        if numel(alpha_raw) ~= numel(T.time), alpha_raw = alpha_raw(:); end
        try
            eeg_tc(red_list(j), c, :) = interp1(T.time(:), alpha_raw(:), t_target(:), 'linear', NaN);
        catch
        end
    end
    for j = 1:numel(amp_list)
        if j > numel(tfr_amp{c}), continue, end
        T = tfr_amp{c}{j};
        if ~isfield(T, 'powspctrm') || ~isfield(T, 'label'), continue, end
        ch_idx = find(ismember(T.label, channels));
        freq_idx = T.freq >= alpha_freq(1) & T.freq <= alpha_freq(2);
        if isempty(ch_idx) || ~any(freq_idx), continue, end
        alpha_raw = squeeze(mean(mean(T.powspctrm(ch_idx, freq_idx, :), 1, 'omitnan'), 2, 'omitnan'));
        if numel(alpha_raw) ~= numel(T.time), alpha_raw = alpha_raw(:); end
        try
            eeg_tc(amp_list(j), c, :) = interp1(T.time(:), alpha_raw(:), t_target(:), 'linear', NaN);
        catch
        end
    end
end
end
