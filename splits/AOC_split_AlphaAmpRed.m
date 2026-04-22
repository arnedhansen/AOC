%% AOC Split Alpha Amp/Red (Subject-Level) — Sternberg
% Subject-level split (fixed across conditions) using merged_data_<task>:
%   mean AlphaPower_FOOOF_bl across load levels (baselined, FOOOFed alpha)
%   split0 mode: < 0 -> reduction, >= 0 -> amplification
%   fixed mode:  < threshold -> reduction, >= threshold -> amplification
%   median mode: < median -> reduction, >= median -> amplification
%
% Uses split-pipeline FOOOF sources:
%   Sternberg: power_stern_fooof_TFR.mat
% plus merged gaze/behavioral metrics.
%
% Generates (per task):
% - Alpha split inclusion figure (participants by alpha, thresholds, group assignment)
% - Power spectra (3 conditions) for both groups
% - TFRs per condition for both groups + group differences
% - Topoplots per condition for both groups + group-difference topoplots
% - Rainclouds for Alpha and gaze deviation
% - Time-course panels for gaze deviation (percent baseline): collapsed and per condition

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
pathAOC = paths.features;

if ispc
    base_data = 'W:\Students\Arne\AOC';
    seb_path = 'W:\Students\Arne\toolboxes\shadedErrorBar';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
    seb_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar';
end
addpath(seb_path);

feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'splits', 'SplitAlphaAmpRed');
fig_dir_nback = fullfile(base_data, 'figures', 'splits', 'SplitAlphaAmpRed', 'Nback');
stats_dir = fullfile(base_data, 'data', 'stats', 'splits');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
if ~isfolder(fig_dir_nback)
    mkdir(fig_dir_nback);
end
if ~isfolder(stats_dir)
    mkdir(stats_dir);
end
fprintf('\n=== AOC Split Alpha Amp/Red (Sternberg + N-back variants) ===\n');
fprintf('Figure directory: %s\n', fig_dir);

% Keep canonical figure size requested.
fig_pos = [0 0 1512 982];

% Task/split definitions:
% 1) Sternberg split at 0
% 2) N-back split
tasks(1).tag = 'sternberg';
tasks(1).split_mode = 'zero';
tasks(1).split_label = 'split0';
tasks(1).merged_file = 'AOC_merged_data_sternberg.mat';
tasks(1).merged_var = 'merged_data_sternberg';
tasks(1).cond_vals = [2 4 6];
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};
tasks(1).power_fname = 'power_stern_fooof_TFR.mat';
tasks(1).pow_vars = {'pow2_fooof_bl', 'pow4_fooof_bl', 'pow6_fooof_bl'};
tasks(1).tfr_fname = 'tfr_stern.mat';
tasks(1).tfr_vars = {'tfr2_fooof_bl', 'tfr4_fooof_bl', 'tfr6_fooof_bl'};
tasks(1).gaze_fname = 'gaze_series_sternberg_trials.mat';
tasks(1).power_missing_label = 'power_stern_fooof_TFR.mat';

tasks(2).tag = 'nback';
% N-back split configuration:
nback_split_mode = 'fixed';
nback_split_fixed_value = -0.05;
tasks(2).split_mode = nback_split_mode;
if strcmpi(nback_split_mode, 'fixed')
    tasks(2).split_label = sprintf('splitFixed_%s', strrep(sprintf('%.3f', nback_split_fixed_value), '.', 'p'));
    tasks(2).split_fixed_value = nback_split_fixed_value;
elseif strcmpi(nback_split_mode, 'median')
    tasks(2).split_label = 'splitMedian';
else
    error('Unsupported nback_split_mode: %s', nback_split_mode);
end
tasks(2).merged_file = 'AOC_merged_data_nback.mat';
tasks(2).merged_var = 'merged_data_nback';
tasks(2).cond_vals = [1 2 3];
tasks(2).cond_codes = [21 22 23];
tasks(2).cond_labels = {'1-back', '2-back', '3-back'};
tasks(2).power_fname = 'power_nback_fooof.mat';
tasks(2).pow_vars = {'pow1_fooof_bl', 'pow2_fooof_bl', 'pow3_fooof_bl'};
tasks(2).tfr_fname = 'tfr_nback.mat';
tasks(2).tfr_vars = {'tfr1_fooof_bl', 'tfr2_fooof_bl', 'tfr3_fooof_bl'};
tasks(2).gaze_fname = 'gaze_series_nback_trials.mat';
tasks(2).power_missing_label = 'power_nback_fooof.mat';

fontSize = 20;
tfr_winsor_cfg = struct();
tfr_winsor_cfg.enable = true;
tfr_winsor_cfg.prctile = [2 98]; % subject-level clipping per TF bin

for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
split_mode = tk.split_mode;
split_label = tk.split_label;
cond_vals = tk.cond_vals;
cond_codes = tk.cond_codes;
cond_labels = tk.cond_labels;
fig_prefix = sprintf('AOC_splitAlphaAmpRed_%s_%s', task_tag, split_label);
fig_dir_task = fig_dir;
if strcmpi(task_tag, 'nback')
    fig_dir_task = fig_dir_nback;
end

fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));

%% Load subject-level merged data and define alpha split
fprintf('\n=== Loading merged data (%s) ===\n', task_tag);
merged_file = fullfile(feat_dir, tk.merged_file);
if ~isfile(merged_file)
    warning('Skipping task %s: missing file %s', task_tag, merged_file);
    continue
end
S = load(merged_file, tk.merged_var);
if ~isfield(S, tk.merged_var)
    warning('Skipping task %s: variable %s not found in %s', task_tag, tk.merged_var, merged_file);
    continue
end
T = struct2table(S.(tk.merged_var));

% Compute subject-level split value from full-window FOOOF-baselined alpha.
uIDs = unique(T.ID);
nSubj = numel(uIDs);
alpha_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    alpha_mean(i) = mean(T.AlphaPower_FOOOF_bl(mask), 'omitnan');
end

% Robust filtering for split reference:
% remove non-finite values and pathological subject-level alpha outliers
% before deriving the zero band. This prevents single corrupt values from
% inflating the percentile-based cutoff.
split_valid = isfinite(alpha_mean);
vals = alpha_mean(split_valid);
if numel(vals) >= 4
    q1_split = prctile(vals, 25);
    q3_split = prctile(vals, 75);
    iqr_split = q3_split - q1_split;
    if isfinite(iqr_split) && iqr_split > 0
        lo_split = q1_split - 3 * iqr_split;
        hi_split = q3_split + 3 * iqr_split;
        split_valid = split_valid & (alpha_mean >= lo_split) & (alpha_mean <= hi_split);
    end
end

valid_alpha_mean = alpha_mean(split_valid);
if isempty(valid_alpha_mean)
    error('No finite subject-level alpha values found for split.');
end
alpha_split_threshold = NaN;
split_info_str = '';
switch lower(split_mode)
    case 'median'
        alpha_split_threshold = median(valid_alpha_mean, 'omitnan');
        reduction_ids = uIDs(split_valid & (alpha_mean < alpha_split_threshold));
        amplification_ids = uIDs(split_valid & (alpha_mean >= alpha_split_threshold));
        zero_ids = [];
        invalid_ids = uIDs(~split_valid);
        split_info_str = sprintf('Median split at %.4f', alpha_split_threshold);
    case 'fixed'
        if ~isfield(tk, 'split_fixed_value') || ~isfinite(tk.split_fixed_value)
            error('Task %s uses split_mode=fixed but split_fixed_value is missing/invalid.', task_tag);
        end
        alpha_split_threshold = tk.split_fixed_value;
        reduction_ids = uIDs(split_valid & (alpha_mean < alpha_split_threshold));
        amplification_ids = uIDs(split_valid & (alpha_mean >= alpha_split_threshold));
        zero_ids = [];
        invalid_ids = uIDs(~split_valid);
        split_info_str = sprintf('Fixed split at %.4f', alpha_split_threshold);
    otherwise
        alpha_split_threshold = 0;
        reduction_ids = uIDs(split_valid & (alpha_mean < alpha_split_threshold));
        amplification_ids = uIDs(split_valid & (alpha_mean >= alpha_split_threshold));
        zero_ids = [];
        invalid_ids = uIDs(~split_valid);
        split_info_str = 'Split at 0.0000 (no near-zero exclusion)';
end

fprintf('\n=== Split Summary [%s | %s] (AlphaPower_FOOOF_bl, full window) ===\n', task_tag, split_mode);
fprintf('Subjects total: %d\n', nSubj);
fprintf('%s\n', split_info_str);
if strcmpi(split_mode, 'median')
    fprintf('Reduction (< median %.4f): %d\n', alpha_split_threshold, numel(reduction_ids));
    fprintf('Amplification (>= median %.4f): %d\n', alpha_split_threshold, numel(amplification_ids));
elseif strcmpi(split_mode, 'fixed')
    fprintf('Reduction (< fixed %.4f): %d\n', alpha_split_threshold, numel(reduction_ids));
    fprintf('Amplification (>= fixed %.4f): %d\n', alpha_split_threshold, numel(amplification_ids));
else
    fprintf('Reduction (< 0.0000): %d\n', numel(reduction_ids));
    fprintf('Amplification (>= 0.0000): %d\n', numel(amplification_ids));
end
if ~isempty(invalid_ids)
    fprintf('Excluded (invalid/pathological alpha): %d\n', numel(invalid_ids));
end

if numel(reduction_ids) < 2 || numel(amplification_ids) < 2
    warning('Task %s: insufficient subjects per split group (red=%d, amp=%d). Skipping task.', ...
        task_tag, numel(reduction_ids), numel(amplification_ids));
    continue
end

%% Alpha split inclusion figure
fprintf('\n=== Plotting alpha split inclusion figure ===\n');
figure('Position', fig_pos, 'Color', 'w');
hold on
% Grey yline at 0
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
if strcmpi(split_mode, 'median') || strcmpi(split_mode, 'fixed')
    yline(alpha_split_threshold, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
else
    yline(0, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
end
% Plot subjects: reduction blue, amplification red, excluded grey
x_vals = (1:nSubj)';
idx_red = ismember(uIDs, reduction_ids);
idx_amp = ismember(uIDs, amplification_ids);
idx_excl = ismember(uIDs, zero_ids) | ismember(uIDs, invalid_ids);
h_excl = scatter(x_vals(idx_excl), alpha_mean(idx_excl), 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
h_red = scatter(x_vals(idx_red), alpha_mean(idx_red), 80, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.8);
h_amp = scatter(x_vals(idx_amp), alpha_mean(idx_amp), 80, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
xlabel('Participant (index)');
ylabel('Alpha Power [dB]');
title(sprintf('Alpha Split [%s | %s]', task_tag, split_mode), 'Interpreter', 'none');
legend([h_excl, h_red, h_amp], {'Excluded', 'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize - 2, 'Box', 'off');
set(gca, 'FontSize', fontSize);
box on
saveas(gcf, fullfile(fig_dir_task, sprintf('%s_inclusion.png', fig_prefix)));
close(gcf);

%% Preallocate containers
% EEG power spectra/topography
pow_red = cell(1, 3);
pow_amp = cell(1, 3);
for c = 1:3
    pow_red{c} = {};
    pow_amp{c} = {};
end

% EEG TFR (tfr_*_subj: subject indices corresponding to each TFR, for correct variance mapping)
tfr_red = cell(1, 3);
tfr_amp = cell(1, 3);
tfr_red_subj = [];
tfr_amp_subj = [];
for c = 1:3
    tfr_red{c} = {};
    tfr_amp{c} = {};
end

% Gaze summary metrics from merged subject table (pre-outlier exclusion)
metrics = struct();
metrics.Alpha = nan(nSubj, 3);
metrics.Dev = nan(nSubj, 3);          % GazeDeviationFullBL [% change]

% Time courses per subject x condition (full resolution, pre-outlier exclusion)
fs = 500;
t_full = -0.5:1/fs:3;
t_plot = t_full(2:end);
Tf = numel(t_plot);
dev_tc = nan(nSubj, 3, Tf);

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

    % Subject-level condition metrics from merged_data_<task>
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics.Alpha(s, c) = mean(subj_rows.AlphaPower_FOOOF_bl(cmask), 'omitnan');
            metrics.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
        end
    end

    % EEG power spectra and topography sources (baselined, FOOOFed)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    try
        P = load(fullfile(eeg_dir, tk.power_fname), tk.pow_vars{:});
        pow_conds = {P.(tk.pow_vars{1}), P.(tk.pow_vars{2}), P.(tk.pow_vars{3})};
        for c = 1:3
            if ismember(sid, reduction_ids)
                pow_red{c}{end+1} = pow_conds{c};
            elseif ismember(sid, amplification_ids)
                pow_amp{c}{end+1} = pow_conds{c};
            end
        end
    catch
        missing_eeg{end+1} = sid_str;
    end

    % EEG TFR sources (store subject index for correct time-course extraction)
    try
        % Use FOOOFed + baselined TFR to avoid raw-power scale skew.
        R = load(fullfile(eeg_dir, tk.tfr_fname), tk.tfr_vars{:});
        tfr_conds = {R.(tk.tfr_vars{1}), R.(tk.tfr_vars{2}), R.(tk.tfr_vars{3})};
        if ismember(sid, reduction_ids)
            for c = 1:3
                tfr_red{c}{end+1} = tfr_conds{c};
            end
            tfr_red_subj(end+1) = s;
        elseif ismember(sid, amplification_ids)
            for c = 1:3
                tfr_amp{c}{end+1} = tfr_conds{c};
            end
            tfr_amp_subj(end+1) = s;
        end
    catch
        missing_tfr{end+1} = sid_str;
    end

    % Gaze series for deviation time courses
    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', tk.gaze_fname);
    if ~isfile(gaze_file)
        missing_gaze{end+1} = sid_str;
        continue
    end

    G = load(gaze_file);
    if ~isfield(G, 'trialinfo')
        missing_gaze{end+1} = sid_str;
        continue
    end

    conds = parse_trialinfo_conds(G.trialinfo);
    if isempty(conds)
        missing_gaze{end+1} = sid_str;
        continue
    end

    % ET coverage diagnostic for the requested 1-3 s window (per subject).
    if isfield(G, 'time') && ~isempty(G.time)
        tmax_trials = nan(1, numel(G.time));
        for tr = 1:numel(G.time)
            tt_cov = G.time{tr};
            if isempty(tt_cov), continue, end
            tmax_trials(tr) = tt_cov(end);
        end
        n_reach_3s = sum(tmax_trials >= 3, 'omitnan');
        fprintf('ET coverage %s: %d/%d trials reach 3 s (tmax median=%.3f s, min=%.3f s)\n', ...
            sid_str, n_reach_3s, numel(tmax_trials), median(tmax_trials, 'omitnan'), min(tmax_trials, [], 'omitnan'));
    end

    % Deviation time course from gaze x/y if available
    has_xy = isfield(G, 'gaze_x') && isfield(G, 'gaze_y');
    if has_xy
        for c = 1:3
            tr_mask = conds == cond_codes(c);
            tr_idx = find(tr_mask);
            if isempty(tr_idx)
                continue
            end

            dev_mat = nan(numel(tr_idx), Tf);
            for k = 1:numel(tr_idx)
                tr = tr_idx(k);
                x = double(G.gaze_x{tr});
                y = double(G.gaze_y{tr});
                if isempty(x) || isempty(y) || numel(x) ~= numel(y)
                    continue
                end

                if isfield(G, 'time') && numel(G.time) >= tr && ~isempty(G.time{tr})
                    tt = double(G.time{tr});
                    if numel(tt) ~= numel(x)
                        tt = linspace(tt(1), tt(end), numel(x));
                    end
                else
                    tt = linspace(-0.5, 3, numel(x));
                end
                dev = sqrt((x - 400).^2 + (y - 300).^2);

                try
                    dev_mat(k, :) = interp1(tt, dev, t_plot, 'linear', NaN);
                catch
                end
            end

            dev_tc(s, c, :) = nanmean(dev_mat, 1);
        end
    end
end

%% Define groups in row index space
is_red = ismember(uIDs, reduction_ids);
is_amp = ismember(uIDs, amplification_ids);

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
    fullfile(fig_dir_task, sprintf('%s_powspctrm.png', fig_prefix)), fig_pos, fontSize);

%% TFRs per condition (both groups + diff, 3x3)
fprintf('\n=== Plotting TFRs ===\n');
color_map_tfr = customcolormap_preset('red-white-blue');

plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, ...
    color_map_tfr, fig_dir_task, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg);

%% TFRs collapsed over conditions (both groups)
fprintf('\n=== Plotting collapsed TFRs (across conditions) ===\n');
plot_group_tfrs_collapsed(tfr_red, tfr_amp, channels, headmodel, color_map_tfr, ...
    fig_dir_task, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg);

%% Topoplots per condition (both groups + differences)
fprintf('\n=== Plotting topoplots ===\n');
plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir_task, fig_prefix, fig_pos, fontSize);

%% Topoplots collapsed over conditions (both groups)
fprintf('\n=== Plotting collapsed topoplots (across conditions) ===\n');
plot_group_topos_collapsed(pow_red, pow_amp, channels, headmodel, fig_dir_task, fig_prefix, fig_pos, fontSize);

%% Rainclouds
fprintf('\n=== Plotting rainclouds ===\n');
plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir_task, fig_prefix, fig_pos, fontSize);

%% EEG and Gaze time courses with effect-size
close all
Tf = size(dev_tc, 3);
t_full_eeg = linspace(-0.5, 3, Tf);
fprintf('\n=== Extracting EEG alpha time course from TFR ===\n');
eeg_tc = extract_alpha_timecourse_tfr(tfr_red, tfr_amp, tfr_red_subj, tfr_amp_subj, nSubj, channels, Tf, t_full_eeg);

% Baselined gaze deviation time course: percent change from baseline (-0.5 to -0.25 s).
fprintf('\n=== Computing baselined time courses ===\n');
t_vec = linspace(-0.5, 3, Tf);
bl_idx = (t_vec >= -0.5) & (t_vec <= -0.25);

dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
% Percent baseline: (value/baseline - 1) * 100.
dev_tc_pct = (dev_tc ./ dev_bl_3d - 1) * 100;
dev_tc_pct(~isfinite(dev_tc_pct)) = NaN;

% Plot time courses (always save both: gaze-only and combined EEG+gaze)
close all
fontSizeTC = 25;
rng(123)
ds_factor = 10; % downsampling to ds_factor*2ms windows
tc_viz_smooth_sec = 0.05; % slight display-only smoothing for time-course plots

% Exclude subjects with incomplete gaze time courses in analysis window.
% Criterion: for each condition, signal must be finite at t=3 s and have
% >=99.5% finite samples in [0, 3] s.
tc_window_idx = (t_vec >= 0) & (t_vec <= 3);
tc_complete_min_frac = 0.995;
tc_finite_frac = squeeze(mean(isfinite(dev_tc(:, :, tc_window_idx)), 3));
tc_has_endpoint = squeeze(isfinite(dev_tc(:, :, end)));
tc_complete_by_cond = (tc_finite_frac >= tc_complete_min_frac) & tc_has_endpoint;
tc_complete_subj = all(tc_complete_by_cond, 2);
tc_excluded_subj = ~tc_complete_subj;
dev_tc_pct(tc_excluded_subj, :, :) = NaN;
is_red_tc = is_red & tc_complete_subj;
is_amp_tc = is_amp & tc_complete_subj;
fprintf('Time-course completeness filter: finite frac >= %.3f in [0,3]s + finite endpoint at 3s\n', tc_complete_min_frac);
fprintf('Excluded incomplete gaze time-course subjects: %d\n', sum(tc_excluded_subj & (is_red | is_amp)));
if any(tc_excluded_subj & (is_red | is_amp))
    excl_ids = uIDs(tc_excluded_subj & (is_red | is_amp));
    fprintf('Excluded IDs: %s\n', sprintf('%d ', excl_ids));
end

% Keep collapsed-across-conditions time course output
plot_timecourse_with_effect_CBPT(dev_tc_pct, is_red_tc, is_amp_tc, colors, ...
    'Gaze Deviation [%]', sprintf('%s_%s_gaze_deviation_pct', task_tag, split_label), fig_dir_task, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc, false, 'Alpha Power [dB]', tc_viz_smooth_sec);
plot_timecourse_individuals(dev_tc_pct, is_red_tc, is_amp_tc, colors, ...
    'Gaze Deviation [%]', 'Collapsed over conditions', ...
    sprintf('%s_%s_gaze_deviation_pct_individuals_collapsed', task_tag, split_label), ...
    fig_dir_task, fig_pos, fontSizeTC, fs, tc_viz_smooth_sec);

% Additional condition-wise time course outputs
for c = 1:numel(cond_vals)
    tc_cond = dev_tc_pct(:, c, :);
    eeg_tc_cond = eeg_tc(:, c, :);
    save_tag_cond = sprintf('%s_%s_gaze_deviation_pct_%s', task_tag, split_label, sanitize_label_for_fname(cond_labels{c}));
    plot_timecourse_with_effect_CBPT(tc_cond, is_red_tc, is_amp_tc, colors, ...
        'Gaze Deviation [%]', save_tag_cond, fig_dir_task, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc_cond, false, 'Alpha Power [dB]', tc_viz_smooth_sec);
    plot_timecourse_individuals(tc_cond, is_red_tc, is_amp_tc, colors, ...
        'Gaze Deviation [%]', cond_labels{c}, ...
        sprintf('%s_individuals', save_tag_cond), ...
        fig_dir_task, fig_pos, fontSizeTC, fs, tc_viz_smooth_sec);
end

%% Sanity checks
fprintf('\n=== Data Diagnostics [%s] ===\n', task_tag);
fprintf('Missing EEG power (%s): %d\n', tk.power_missing_label, numel(unique(missing_eeg)));
fprintf('Missing TFR files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir_task);

%% Export CSV for Python statistics script
fprintf('\n=== Exporting CSV for Python stats ===\n');
t_win_lo = 0;
t_win_hi = 2;
task_idx_gaze = (t_vec >= t_win_lo) & (t_vec <= t_win_hi);
dev_pct_by_load = squeeze(mean(dev_tc_pct(:, :, task_idx_gaze), 3, 'omitnan'));

n_rows = nSubj * numel(cond_vals);
ID_col = nan(n_rows, 1);
LoadValue_col = nan(n_rows, 1);
LoadLabel_col = strings(n_rows, 1);
Group_col = strings(n_rows, 1);
Included_col = false(n_rows, 1);
Alpha_col = nan(n_rows, 1);
GazeDev_col = nan(n_rows, 1);
GazeDevPct_col = nan(n_rows, 1);

r = 0;
for s = 1:nSubj
    for c = 1:numel(cond_vals)
        r = r + 1;
        ID_col(r) = uIDs(s);
        LoadValue_col(r) = cond_vals(c);
        LoadLabel_col(r) = string(cond_labels{c});
        Alpha_col(r) = metrics.Alpha(s, c);
        GazeDev_col(r) = metrics.Dev(s, c);
        GazeDevPct_col(r) = dev_pct_by_load(s, c);

        if is_red(s)
            Group_col(r) = "Reduction";
            Included_col(r) = is_red_tc(s);
        elseif is_amp(s)
            Group_col(r) = "Amplification";
            Included_col(r) = is_amp_tc(s);
        else
            Group_col(r) = "Excluded";
            Included_col(r) = false;
        end
    end
end

stats_tbl = table( ...
    ID_col, LoadValue_col, LoadLabel_col, Group_col, Included_col, ...
    Alpha_col, GazeDev_col, GazeDevPct_col, ...
    'VariableNames', { ...
    'ID', 'LoadValue', 'LoadLabel', 'Group', 'Included', ...
    'AlphaPower_FOOOF_bl', 'GazeDeviationFullBL', 'GazeDev_pct_0_2s'});

csv_out = fullfile(stats_dir, sprintf('AOC_splitAmpRed_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV saved to: %s\n', csv_out);

end % task loop

%% ========================= Local Functions =========================
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
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end+1} = lab;
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

function txt = sanitize_label_for_fname(label)
txt = lower(label);
txt = regexprep(txt, '[^a-z0-9]+', '_');
txt = regexprep(txt, '_+', '_');
txt = regexprep(txt, '^_|_$', '');
end


function plot_group_power_spectrum_combined(pow_red, pow_amp, channels, colors, cond_labels, out_file, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping power spectrum (incomplete condition data).');
    return
end
ga_red = cell(1, 3);
ga_amp = cell(1, 3);
guard_cfg = struct('hard_abs', 1e4, 'winsor_prc', [1 99]);
for c = 1:3
    red_clean = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_red{c}, 'UniformOutput', false);
    amp_clean = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_amp{c}, 'UniformOutput', false);
    ga_red{c} = ft_freqgrandaverage([], red_clean{:});
    ga_amp{c} = ft_freqgrandaverage([], amp_clean{:});
    pow_red{c} = red_clean;
    pow_amp{c} = amp_clean;
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
        all_m = [all_m; m]; all_se = [all_se; se];
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
    legend([eb.mainLine], cond_labels, 'FontSize', fsz - 2, 'Location', 'best', 'Box', 'off');
    set(gca, 'FontSize', fsz);
    box on
end
saveas(gcf, out_file);
close(gcf);
end

function plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping TFR plot (incomplete condition data).');
    return
end

if nargin < 11 || isempty(winsor_cfg)
    winsor_cfg = struct('enable', false, 'prctile', [2 98]);
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    cfg = [];
    cfg.keepindividual = 'yes';
    ga_red{c} = ft_freqgrandaverage(cfg, tfr_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage(cfg, tfr_amp{c}{:});
    if winsor_cfg.enable
        ga_red{c} = winsorize_freq_subjects(ga_red{c}, winsor_cfg.prctile);
        ga_amp{c} = winsorize_freq_subjects(ga_amp{c}, winsor_cfg.prctile);
    end
end

% Shared clim for Reduction and Amplification
[~, ch_idx] = ismember(channels, ga_red{1}.label);
ch_idx = ch_idx(ch_idx > 0);
if isempty(ch_idx)
    warning('Skipping TFR plot (none of the requested channels found in TFR labels).');
    return
end
freq_idx = ga_red{1}.freq >= 5 & ga_red{1}.freq <= 30;
time_idx = ga_red{1}.time >= -0.5 & ga_red{1}.time <= 3;
mx = 0;
for c = 1:3
    Ared = mean_over_channels_tfr(ga_red{c}, ch_idx);
    Aamp = mean_over_channels_tfr(ga_amp{c}, ch_idx);
    mx = max(mx, max(abs(Ared(freq_idx, time_idx)), [], 'all'));
    mx = max(mx, max(abs(Aamp(freq_idx, time_idx)), [], 'all'));
end
clim_abs = [-mx*0.9 mx*0.9]

cfg = [];
cfg.channel = channels;
cfg.colorbar = 'no';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 3];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;

% Single 3x2 figure: column 1 = Reduction, column 2 = Amplification
fsz_tfr = round(fsz * 0.8);  % ~20% reduction
% target_ticks = 9;
% tick_step = 0.05 * ceil(((2 * mx) / max(target_ticks - 1, 1)) / 0.05);
% tick_step = max(tick_step, 0.05);
% tick_max = tick_step * ceil(mx / tick_step);
% cb_ticks = -tick_max:tick_step:tick_max;
% clim_abs = [-tick_max tick_max];
figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    % Column 1: Reduction (rows 1–3)
    ax = subplot(3, 2, (c-1)*2 + 1);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_red{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    cbar = colorbar(ax);
    xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
    set(ax, 'FontSize', fsz_tfr);
    cbar.FontSize = fsz_tfr - 4;
    %cbar.Ticks = cb_ticks;
    cbar.Label.String = 'Power [dB]';
    cbar.Label.FontSize = fsz_tfr;
    title(ax, sprintf('Reduction - %s', cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');

    % Column 2: Amplification (rows 1–3)
    ax = subplot(3, 2, (c-1)*2 + 2);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_amp{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    cbar = colorbar(ax);
    xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
    set(ax, 'FontSize', fsz_tfr);
    cbar.FontSize = fsz_tfr - 4;
    %cbar.Ticks = cb_ticks;
    cbar.Label.String = 'Power [dB]';
    cbar.Label.FontSize = fsz_tfr;
    title(ax, sprintf('Amplification - %s', cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');
end
saveas(gcf, fullfile(fig_dir, sprintf('%s_tfr_all.png', fig_prefix)));
close(gcf);
end

function A = mean_over_channels_tfr(T, ch_idx)
% Return freq x time map averaged across selected channels and subjects (if present).
P = T.powspctrm;
if isempty(P)
    A = nan(numel(T.freq), numel(T.time));
    return
end

Psize = size(P);
chan_dim = find(Psize == numel(T.label), 1, 'first');
if isempty(chan_dim)
    if ndims(P) >= 2
        chan_dim = 2; % Common FieldTrip keepindividual layout: rpt x chan x freq x time
    else
        chan_dim = 1;
    end
end

switch ndims(P)
    case 4
        if chan_dim == 1
            A = squeeze(mean(P(ch_idx, :, :, :), 1, 'omitnan')); % rpt x freq x time
        elseif chan_dim == 2
            A = squeeze(mean(P(:, ch_idx, :, :), 2, 'omitnan')); % rpt x freq x time
        else
            error('Unsupported channel dimension in 4D TFR data.');
        end
        if ndims(A) == 3
            A = squeeze(mean(A, 1, 'omitnan')); % freq x time
        end
    case 3
        if chan_dim == 1
            A = squeeze(mean(P(ch_idx, :, :), 1, 'omitnan')); % freq x time
        elseif chan_dim == 2
            A = squeeze(mean(P(:, ch_idx, :), 2, 'omitnan'));
            if ndims(A) == 2 && size(A, 1) == numel(T.time)
                A = A'; % enforce freq x time
            end
        else
            error('Unsupported channel dimension in 3D TFR data.');
        end
    otherwise
        error('Unsupported TFR powspctrm dimensionality.');
end
end

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm)
    return
end

P = freq_in.powspctrm;
if size(P, 1) < 3
    % Too few subjects for stable percentile clipping.
    return
end

P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir, fig_prefix, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping topoplots (missing power data).');
    return
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
guard_cfg = struct('hard_abs', 1e4, 'winsor_prc', [1 99]);
for c = 1:3
    red_clean = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_red{c}, 'UniformOutput', false);
    amp_clean = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_amp{c}, 'UniformOutput', false);
    ga_red{c} = ft_freqgrandaverage([], red_clean{:});
    ga_amp{c} = ft_freqgrandaverage([], amp_clean{:});
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
    all_alpha = [all_alpha; Ared(:); Aamp(:)];
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
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_all.png', fig_prefix)));
close(gcf);
end

function plot_group_tfrs_collapsed(tfr_red, tfr_amp, channels, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping collapsed TFR plot (incomplete condition data).');
    return
end

if nargin < 9 || isempty(winsor_cfg)
    winsor_cfg = struct('enable', false, 'prctile', [2 98]);
end

% Pool all subject-level TFRs across conditions within each group.
tfr_red_all = [tfr_red{1}, tfr_red{2}, tfr_red{3}];
tfr_amp_all = [tfr_amp{1}, tfr_amp{2}, tfr_amp{3}];
if isempty(tfr_red_all) || isempty(tfr_amp_all)
    warning('Skipping collapsed TFR plot (no pooled TFR data).');
    return
end

cfg = [];
cfg.keepindividual = 'yes';
ga_red = ft_freqgrandaverage(cfg, tfr_red_all{:});
ga_amp = ft_freqgrandaverage(cfg, tfr_amp_all{:});
if winsor_cfg.enable
    ga_red = winsorize_freq_subjects(ga_red, winsor_cfg.prctile);
    ga_amp = winsorize_freq_subjects(ga_amp, winsor_cfg.prctile);
end

[~, ch_idx] = ismember(channels, ga_red.label);
ch_idx = ch_idx(ch_idx > 0);
if isempty(ch_idx)
    warning('Skipping collapsed TFR plot (none of requested channels found).');
    return
end
freq_idx = ga_red.freq >= 5 & ga_red.freq <= 30;
time_idx = ga_red.time >= -0.5 & ga_red.time <= 3;
Ared = mean_over_channels_tfr(ga_red, ch_idx);
Aamp = mean_over_channels_tfr(ga_amp, ch_idx);
mx = max([max(abs(Ared(freq_idx, time_idx)), [], 'all'), max(abs(Aamp(freq_idx, time_idx)), [], 'all')]);
if ~isfinite(mx) || mx <= 0
    mx = 0.1;
end
clim_abs = [-0.9 * mx, 0.9 * mx];

cfg = [];
cfg.channel = channels;
cfg.colorbar = 'no';
cfg.zlim = 'maxabs';
cfg.xlim = [-0.5 3];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;

fsz_tfr = round(fsz * 0.8);
figure('Position', [0 0 1512*0.666 982], 'Color', 'w');
ax = subplot(2, 1, 1);
cfg.figure = ax;
ft_singleplotTFR(cfg, ga_red);
colormap(ax, color_map);
set(ax, 'CLim', clim_abs);
cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr);
rectangle('Position', [0 8 3 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.FontSize = fsz_tfr - 4;
cbar.Label.String = 'Power [dB]';
cbar.Label.FontSize = fsz_tfr;
title(ax, 'Reduction', 'FontSize', fsz_tfr+4, 'Interpreter', 'none');

ax = subplot(2, 1, 2);
cfg.figure = ax;
ft_singleplotTFR(cfg, ga_amp);
colormap(ax, color_map);
set(ax, 'CLim', clim_abs);
cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr);
rectangle('Position', [0 8 3 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.FontSize = fsz_tfr - 4;
cbar.Label.String = 'Power [dB]';
cbar.Label.FontSize = fsz_tfr;
title(ax, 'Amplification', 'FontSize', fsz_tfr+4, 'Interpreter', 'none');

saveas(gcf, fullfile(fig_dir, sprintf('%s_tfr_collapsedConditions.png', fig_prefix)));
end

function plot_group_topos_collapsed(pow_red, pow_amp, channels, headmodel, fig_dir, fig_prefix, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping collapsed topoplots (missing power data).');
    return
end

% Pool all subject-level power spectra across conditions within each group.
pow_red_all = [pow_red{1}, pow_red{2}, pow_red{3}];
pow_amp_all = [pow_amp{1}, pow_amp{2}, pow_amp{3}];
if isempty(pow_red_all) || isempty(pow_amp_all)
    warning('Skipping collapsed topoplots (no pooled power data).');
    return
end

guard_cfg = struct('hard_abs', 1e4, 'winsor_prc', [1 99]);
pow_red_all = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_red_all, 'UniformOutput', false);
pow_amp_all = cellfun(@(x) sanitize_pow_struct(x, guard_cfg), pow_amp_all, 'UniformOutput', false);
ga_red = ft_freqgrandaverage([], pow_red_all{:});
ga_amp = ft_freqgrandaverage([], pow_amp_all{:});

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
cfg.colormap = rdbu_cmap(64);

freq_idx = ga_red.freq >= 8 & ga_red.freq <= 14;
Ared = mean(ga_red.powspctrm(:, freq_idx), 2, 'omitnan');
Aamp = mean(ga_amp.powspctrm(:, freq_idx), 2, 'omitnan');
all_alpha = [Ared(:); Aamp(:)];
all_alpha = all_alpha(isfinite(all_alpha));
if isempty(all_alpha)
    cfg.zlim = 'maxabs';
else
    mx = prctile(abs(all_alpha), 95);
    if mx <= 0
        mx = max(abs(all_alpha));
    end
    if ~isfinite(mx) || mx <= 0
        cfg.zlim = 'maxabs';
    else
        cfg.zlim = [-mx mx];
    end
end

figure('Position', fig_pos, 'Color', 'w');
ax = subplot(1, 2, 1);
cfg.figure = ax;
ft_topoplotER(cfg, ga_red);
colorbar(ax);
set(ax, 'FontSize', fsz);
title(ax, 'Reduction - collapsed over conditions', 'Interpreter', 'none');

ax = subplot(1, 2, 2);
cfg.figure = ax;
ft_topoplotER(cfg, ga_amp);
colorbar(ax);
set(ax, 'FontSize', fsz);
title(ax, 'Amplification - collapsed over conditions', 'Interpreter', 'none');

saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_prefix, fig_pos, fsz)
metric_defs = { ...
    'Alpha', 'AlphaPower_FOOOF_bl', false; ...
    'Dev', 'Gaze deviation', false; ...
    };

for m = 1:size(metric_defs, 1)
    key = metric_defs{m, 1};
    varname = metric_defs{m, 2};
    is_dB = metric_defs{m, 3};
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

    if is_dB
        ylab = [varname ' [dB]'];
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
    saveas(gcf, fullfile(fig_dir, sprintf('%s_raincloud_%s.png', fig_prefix, lower(key))));
    close(gcf);
end
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha, density_offset, box_offset, dot_jitter_halfwidth)
y = yvals(isfinite(yvals));
if numel(y) < 3
    return
end
[f, xi] = ksdensity(y, 'NumPoints', 120);
if nargin < 7 || isempty(density_offset), density_offset = 0.08; end
if nargin < 8 || isempty(box_offset), box_offset = 0.03; end
if nargin < 9 || isempty(dot_jitter_halfwidth), dot_jitter_halfwidth = box_w / 2; end
if max(f) > 0
    f = f / max(f) * 0.35;
else
    f = zeros(size(f));
end
x_den = xpos - density_offset;
x_box = xpos + box_offset;
fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25);
q3 = prctile(y, 75);
med = median(y);
p5 = prctile(y, 5);
p95 = prctile(y, 95);
plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [x_box-box_w/2, q1, box_w, q3-q1], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = dot_jitter_halfwidth * 2 * (rand(numel(y),1)-0.5);
scatter(x_box + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

function plot_timecourse_individuals(tc, is_red, is_amp, colors, ylab, title_tag, save_tag, fig_dir, fig_pos, fsz, fs, smooth_sec)
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 3, nT);
if nargin < 12 || isempty(smooth_sec)
    smooth_sec = 0;
end
win_sm = max(1, round(smooth_sec * fs));

% Collapse over conditions if multiple condition slices are present.
tc_subj = squeeze(mean(tc, 2, 'omitnan')); % nSubj x nT
if isvector(tc_subj)
    tc_subj = tc_subj(:)';
end

R = tc_subj(is_red, :);
A = tc_subj(is_amp, :);
R(~isfinite(R)) = NaN;
A(~isfinite(A)) = NaN;
if win_sm > 1
    R = movmean(R, win_sm, 2, 'omitnan');
    A = movmean(A, win_sm, 2, 'omitnan');
end
colR_light = colors(1, :) * 0.35 + 0.65;
colA_light = colors(3, :) * 0.35 + 0.65;

figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile; hold on
if ~isempty(R)
    plot(t_plot, R', 'Color', colR_light, 'LineWidth', 0.8);
    plot(t_plot, mean(R, 1, 'omitnan'), 'Color', colors(1, :), 'LineWidth', 3);
end
xline(0, '--k');
xlim([-0.5 3]);
ylabel(ylab);
title(sprintf('Reduction (n=%d) - %s', size(R, 1), title_tag), 'Interpreter', 'none');
set(gca, 'FontSize', fsz - 6);
box on

nexttile; hold on
if ~isempty(A)
    plot(t_plot, A', 'Color', colA_light, 'LineWidth', 0.8);
    plot(t_plot, mean(A, 1, 'omitnan'), 'Color', colors(3, :), 'LineWidth', 3);
end
xline(0, '--k');
xlim([-0.5 3]);
xlabel('Time [s]');
ylabel(ylab);
title(sprintf('Amplification (n=%d) - %s', size(A, 1), title_tag), 'Interpreter', 'none');
set(gca, 'FontSize', fsz - 6);
box on

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlphaAmpRed_timecourse_%s.png', save_tag)));
close(gcf);
end

% Cluster-based permutation test
function plot_timecourse_with_effect_CBPT(tc, is_red, is_amp, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, eeg_tc, addEEG_TC, eeg_ylab, smooth_sec)
if nargin < 11 || isempty(ds_factor), ds_factor = 50; end
if nargin < 12, eeg_tc = []; end
if nargin < 13, addEEG_TC = ~isempty(eeg_tc); end
if nargin < 14, eeg_ylab = 'Alpha Power [dB]'; end
if nargin < 15 || isempty(smooth_sec), smooth_sec = 0.05; end
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 3, nT);
tc(~isfinite(tc)) = NaN;
Rall = squeeze(mean(tc(is_red, :, :), 2, 'omitnan'));
Aall = squeeze(mean(tc(is_amp, :, :), 2, 'omitnan'));
win_sm = max(1, round(smooth_sec * fs));
if win_sm > 1
    Rall = movmean(Rall, win_sm, 2, 'omitnan');
    Aall = movmean(Aall, win_sm, 2, 'omitnan');
end
figure('Position', fig_pos, 'Color', 'w');
nR = size(Rall, 1);
nA = size(Aall, 1);
show_eeg = addEEG_TC && ~isempty(eeg_tc) && size(eeg_tc, 1) >= (nR + nA) && size(eeg_tc, 2) >= 1 && size(eeg_tc, 3) >= nT;
if show_eeg
    tiledlayout(5, 1, 'TileSpacing', 'compact');
else
    tiledlayout(3, 1, 'TileSpacing', 'compact');
end
ax_gaze = nexttile([2 1]); hold on
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
xlim([-0.5 3]);
box on
set(gca, 'FontSize', fsz-4);
% Legend with colored patch boxes (clearer than thin lines)
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1 leg_p2], {'Reduction', 'Amplification'}, 'Location', 'northeast', 'FontSize', fsz-2, 'Box', 'off');
if contains(save_tag, 'gaze_deviation_pct')
    % ylim from shaded envelopes (mean ± SEM), not just the main lines
    upper = max(mR + sR, mA + sA);
    lower = min(mR - sR, mA - sA);
    y_hi = max(upper, [], 'omitnan');
    y_lo = min(lower, [], 'omitnan');
    if isfinite(y_hi) && isfinite(y_lo) && y_hi > y_lo
        pad = max(0.15, 0.05 * (y_hi - y_lo));
        ylim([y_lo - pad, y_hi + pad]);
    end
end
nexttile; hold on
n_perm = 10000;
min_per_group = 3;
d = nan(1, nT);
for t = 1:nT
    x = Rall(:, t); y = Aall(:, t);
    x = x(isfinite(x)); y = y(isfinite(y));
    if numel(x) < min_per_group || numel(y) < min_per_group, continue, end
    sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / max(numel(x)+numel(y)-2, 1));
    d(t) = (mean(y) - mean(x)) / max(sp, eps);
end
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
    bin_width_ms = dt_ds * 1000;
    for k = 1:numel(clusters)
        t_lo = t_plot_ds(clusters(k).idx(1));
        t_hi = t_plot_ds(clusters(k).idx(end));
        t_start = t_lo - dt_ds/2;  % leading edge of first bin (post-stimulus)
        t_end = t_hi + dt_ds/2;    % trailing edge of last bin
        duration_ms = clusters(k).extent * bin_width_ms;
        pass_p = clusters(k).p < 0.05;
        status = 'n.s.';
        if pass_p
            status = 'SIGNIFICANT';
        end
        fprintf('    Cluster %d: window [%.3f, %.3f]s post-stim (duration=%d ms, %d bins @ %.0f ms); mass=%.1f; p=%.4f; %s\n', ...
            k, t_start, t_end, round(duration_ms), clusters(k).extent, bin_width_ms, clusters(k).mass, clusters(k).p, status);
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
ylabel('Effect Size [Cohen''s d]');
xlim([-0.5 3]);
box on
set(gca, 'FontSize', fsz-4);
if ~any(sig_cluster) && any(sig_uncorr)
    text(0.75, ylims(2) - 0.08*diff(ylims), 'WARNING: No significant clusters; shading shows uncorrected t>t_{crit}', ...
        'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'HorizontalAlignment', 'center');
end

% Alpha power time course panel in combined figure.
% Error bars: between-subject SEM = std(sample) / sqrt(n) at each time point.
if show_eeg
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
    ax_eeg = gca;
    ylabel(eeg_ylab);
    xlabel('Time [s]');
    xlim([-0.5 3]);
    box on
    set(gca, 'FontSize', fsz-4);
    leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
    leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
    % Align y-labels to same horizontal position
    drawnow;
    pos_g = get(ax_gaze, 'Position');
    pos_e = get(ax_eeg, 'Position');
    left_min = min(pos_g(1), pos_e(1));
    set(ax_gaze, 'Position', [left_min, pos_g(2), pos_g(3) + (pos_g(1) - left_min), pos_g(4)]);
    set(ax_eeg,  'Position', [left_min, pos_e(2), pos_e(3) + (pos_e(1) - left_min), pos_e(4)]);

    % EEG effect-size strip (same logic as gaze panel above)
    nexttile; hold on
    d_eeg = nan(1, nT);
    for t = 1:nT
        x = EegR(:, t); y = EegA(:, t);
        x = x(isfinite(x)); y = y(isfinite(y));
        if numel(x) < min_per_group || numel(y) < min_per_group, continue, end
        sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / max(numel(x)+numel(y)-2, 1));
        d_eeg(t) = (mean(y) - mean(x)) / max(sp, eps);
    end

    EegR_ds = EegR(:, 1:ds_factor:end);
    EegA_ds = EegA(:, 1:ds_factor:end);
    [clusters_eeg, tvals_cl_eeg, thr_eeg] = ft_cluster_permutation_1d(EegR_ds, EegA_ds, n_perm, 0.05, 'onetail_pos', t_plot_ds);
    sig_cluster_eeg = false(1, nT_ds);
    for k = 1:numel(clusters_eeg)
        if clusters_eeg(k).p < 0.05
            sig_cluster_eeg(clusters_eeg(k).idx) = true;
        end
    end
    sig_uncorr_eeg = (tvals_cl_eeg > thr_eeg.tcrit) & isfinite(tvals_cl_eeg);
    sig_eeg = sig_cluster_eeg;
    if ~any(sig_eeg) && any(sig_uncorr_eeg)
        sig_eeg = sig_uncorr_eeg;
        fprintf('  [%s EEG] (cluster n.s.; shading uncorrected t>tcrit)\n', save_tag);
    end

    d_eeg_ds = d_eeg(1:ds_factor:end);
    if any(sig_cluster_eeg)
        sig_eeg = sig_eeg & isfinite(d_eeg_ds);
    end
    run_start_eeg = [false, diff(sig_eeg) == 1];
    run_end_eeg = [diff(sig_eeg) == -1, false];
    if sig_eeg(1), run_start_eeg(1) = true; end
    if sig_eeg(end), run_end_eeg(end) = true; end
    starts_eeg = find(run_start_eeg);
    ends_eeg = find(run_end_eeg);

    d_eeg_fin = d_eeg(isfinite(d_eeg));
    mx_eeg = max(abs(d_eeg_fin), [], 'omitnan');
    if isempty(d_eeg_fin) || ~isfinite(mx_eeg) || mx_eeg == 0
        ylims_eeg = [-0.6 0.6];
    else
        ylims_eeg = [-max(mx_eeg + 0.1, 0.6), max(mx_eeg + 0.1, 0.6)];
    end

    patch_alpha_eeg = 0.4 * ~any(sig_cluster_eeg) + 0.25 * any(sig_cluster_eeg);
    for k = 1:numel(starts_eeg)
        t1 = max(0, t_plot_ds(starts_eeg(k)) - dt_ds/2);
        t2 = t_plot_ds(ends_eeg(k)) + dt_ds/2;
        patch([t1 t2 t2 t1], [ylims_eeg(1) ylims_eeg(1) ylims_eeg(2) ylims_eeg(2)], [0.5 0.5 0.5], 'FaceAlpha', patch_alpha_eeg, 'EdgeColor', 'none');
    end

    ylim(ylims_eeg);
    plot(t_plot, d_eeg, 'k-', 'LineWidth', 3.5);
    yline(0, '--');
    xline(0, '--k');
    xlabel('Time [s]');
    ylabel('EEG Cohen''s d');
    xlim([-0.5 3]);
    box on
    set(gca, 'FontSize', fsz-4);
end

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlphaAmpRed_timecourse_%s_CBPT.png', save_tag)));
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
%   t_plot_ds - (optional) 1 x nT actual time vector. If provided, uses cfg.latency=[0 3]
%               to restrict clustering to post-stimulus only (avoids baseline clusters).
%
% Outputs:
%   clusters     - struct array with .idx, .mass, .extent, .p, .p_extent
%   tvals        - 1 x nT t-values (group1 - group2; FT indepsamplesT)
%   thr          - struct with .tcrit, .mass, .extent
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
    cfg.latency = [0 3];  % restrict to post-stimulus only
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

% When we used cfg.latency=[0 3], stat has fewer time points; map back to full range
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
    clusters(end+1).idx = idx;
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

function cmap = rdbu_cmap(n)
% RDBU_CMAP - Red-Blue diverging colormap (no external deps, ColorBrewer RdBu-like).
% Usage: cmap = rdbu_cmap(n);  returns n x 3 RGB matrix.
rdbu_11 = [33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; ...
    253 219 199; 244 165 130; 214 96 77; 178 24 43] / 255;
x = linspace(0, 1, size(rdbu_11, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, rdbu_11, xi, 'linear');
end

function eeg_tc = extract_alpha_timecourse_tfr(tfr_red, tfr_amp, tfr_red_subj, tfr_amp_subj, nSubj, channels, Tf, t_target)
% Extract alpha power (8-14 Hz) over occipital channels from TFR, per subject x condition.
% tfr_*_subj: subject indices (1..nSubj) for each TFR in the cell arrays; required for correct
%   variance when some subjects have missing TFR (otherwise j-th TFR may not match j-th in group).
% Returns eeg_tc: nSubj x 3 x Tf.
eeg_tc = nan(nSubj, 3, Tf);
alpha_freq = [8 14];
for c = 1:3
    for j = 1:numel(tfr_red_subj)
        if j > numel(tfr_red{c}), break, end
        T = tfr_red{c}{j};
        s = tfr_red_subj(j);
        if ~isfield(T, 'powspctrm') || ~isfield(T, 'label'), continue, end
        ch_idx = find(ismember(T.label, channels));
        freq_idx = T.freq >= alpha_freq(1) & T.freq <= alpha_freq(2);
        if isempty(ch_idx) || ~any(freq_idx), continue, end
        alpha_raw = squeeze(mean(mean(T.powspctrm(ch_idx, freq_idx, :), 1, 'omitnan'), 2, 'omitnan'));
        if numel(alpha_raw) ~= numel(T.time), alpha_raw = alpha_raw(:); end
        try
            eeg_tc(s, c, :) = interp1(T.time(:), alpha_raw(:), t_target(:), 'linear', NaN);
        catch
        end
    end
    for j = 1:numel(tfr_amp_subj)
        if j > numel(tfr_amp{c}), break, end
        T = tfr_amp{c}{j};
        s = tfr_amp_subj(j);
        if ~isfield(T, 'powspctrm') || ~isfield(T, 'label'), continue, end
        ch_idx = find(ismember(T.label, channels));
        freq_idx = T.freq >= alpha_freq(1) & T.freq <= alpha_freq(2);
        if isempty(ch_idx) || ~any(freq_idx), continue, end
        alpha_raw = squeeze(mean(mean(T.powspctrm(ch_idx, freq_idx, :), 1, 'omitnan'), 2, 'omitnan'));
        if numel(alpha_raw) ~= numel(T.time), alpha_raw = alpha_raw(:); end
        try
            eeg_tc(s, c, :) = interp1(T.time(:), alpha_raw(:), t_target(:), 'linear', NaN);
        catch
        end
    end
end
end

function S = sanitize_pow_struct(S, cfg)
if ~isfield(S, 'powspctrm') || isempty(S.powspctrm)
    return
end
X = S.powspctrm;
X(~isfinite(X)) = NaN;
X(abs(X) > cfg.hard_abs) = NaN;
vals = X(isfinite(X));
if numel(vals) >= 20
    lo = prctile(vals, cfg.winsor_prc(1));
    hi = prctile(vals, cfg.winsor_prc(2));
    X = min(max(X, lo), hi);
end
S.powspctrm = X;
end
