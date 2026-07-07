%% AOC Split ERS/ERD (Subject-Level) — MS (Sternberg + N-back)
% MS-only outputs: inclusion figure, MS rainclouds, MS time courses (no EEG plots).
% MS time courses use the same ET pipeline as AOC_gaze_microsaccades_*.m
% (detect_microsaccades, 50 ms Gaussian smoothing, baseline [-1.5 -0.5] s).
%
% Generates (per task):
% - ERSD split inclusion figure
% - Rainclouds for MS rate (MSRateFullBL)
% - Time-course panels for MS rate (percent baseline): collapsed and per condition

%% Setup
startup
[subjects, paths, colors] = setup('AOC');

addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir_root = fullfile(paths.figures, 'splits', 'SplitERSERD');
stats_dir = paths.splits_stats;
if ~isfolder(fig_dir_root)
    mkdir(fig_dir_root);
end
if ~isfolder(stats_dir)
    mkdir(stats_dir);
end
fprintf('\n=== AOC Split ERS/ERD — Microsaccades (Sternberg + N-back) ===\n');
fprintf('Main figure directory: %s\n', fig_dir_root);

% Keep canonical figure size requested.
fig_pos = [0 0 1512 982];

% Sternberg task (split0 on mean ERSD_late)
tasks(1).tag = 'sternberg';
tasks(1).split_label = 'split0';
tasks(1).merged_file = 'AOC_merged_data_sternberg.mat';
tasks(1).merged_var = 'merged_data_sternberg';
tasks(1).cond_vals = [2 4 6];
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};
tasks(1).et_fname = 'dataET_sternberg';
tasks(1).ersd_var = 'ERSD_late';
tasks(1).fig_subdir = 'SternbergMS';
tasks(1).group_lbl_low = 'ERD';
tasks(1).group_lbl_high = 'ERS';

% N-back task (median split on pooled ERSD_full)
tasks(2).tag = 'nback';
tasks(2).split_label = 'splitMedian';
tasks(2).merged_file = 'AOC_merged_data_nback.mat';
tasks(2).merged_var = 'merged_data_nback';
tasks(2).cond_vals = [1 2 3];
tasks(2).cond_codes = [21 22 23];
tasks(2).cond_labels = {'1-back', '2-back', '3-back'};
tasks(2).et_fname = 'dataET_nback';
tasks(2).ersd_var = 'ERSD_full';
tasks(2).fig_subdir = 'NbackMS';
tasks(2).group_lbl_low = 'More ERD';
tasks(2).group_lbl_high = 'Less ERD';

fontSize = 40;

for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
split_label = tk.split_label;
cond_vals = tk.cond_vals;
cond_codes = tk.cond_codes;
cond_labels = tk.cond_labels;
fig_prefix = sprintf('AOC_splitERSERD_MS_%s_%s', task_tag, split_label);
fig_dir_task = fullfile(fig_dir_root, tk.fig_subdir);
if ~isfolder(fig_dir_task)
    mkdir(fig_dir_task);
end
fprintf('Supplementary figure directory: %s\n', fig_dir_task);

fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));

%% Load subject-level merged data and define ERSD split
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
if ~ismember(tk.ersd_var, T.Properties.VariableNames)
    warning('Skipping task %s: variable %s not found in merged table', task_tag, tk.ersd_var);
    continue
end

% Compute subject-level split value from ERSD (occipital 8-14 Hz, dB).
uIDs = unique(T.ID);
nSubj = numel(uIDs);
ersd_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    ersd_mean(i) = mean(T.(tk.ersd_var)(mask), 'omitnan');
end

% Robust filtering for split reference:
% remove non-finite values and pathological subject-level ERSD outliers
% before deriving the zero band. This prevents single corrupt values from
% inflating the percentile-based cutoff.
split_valid = isfinite(ersd_mean);
vals = ersd_mean(split_valid);
if numel(vals) >= 4
    q1_split = prctile(vals, 25);
    q3_split = prctile(vals, 75);
    iqr_split = q3_split - q1_split;
    if isfinite(iqr_split) && iqr_split > 0
        lo_split = q1_split - 3 * iqr_split;
        hi_split = q3_split + 3 * iqr_split;
        split_valid = split_valid & (ersd_mean >= lo_split) & (ersd_mean <= hi_split);
    end
end

valid_ersd_mean = ersd_mean(split_valid);
if isempty(valid_ersd_mean)
    error('No finite subject-level ERSD values found for split.');
end
ersd_split_threshold = NaN;
split_info_str = '';
zero_ids = [];
if strcmpi(task_tag, 'sternberg')
    ersd_split_threshold = 0;
    split_summary_suffix = sprintf('%s, [1 2]s', tk.ersd_var);
    split_info_str = 'Split at 0.0000 (no near-zero exclusion)';
else
    ersd_split_threshold = median(valid_ersd_mean, 'omitnan');
    split_summary_suffix = sprintf('%s, pooled loads', tk.ersd_var);
    split_info_str = sprintf('Median split at %.4f', ersd_split_threshold);
end
erd_ids = uIDs(split_valid & (ersd_mean < ersd_split_threshold));
ers_ids = uIDs(split_valid & (ersd_mean >= ersd_split_threshold));
invalid_ids = uIDs(~split_valid);

fprintf('\n=== Split Summary [%s | %s] (%s) ===\n', task_tag, split_label, split_summary_suffix);
fprintf('Subjects total: %d\n', nSubj);
fprintf('%s\n', split_info_str);
if strcmpi(task_tag, 'sternberg')
    fprintf('ERD (< 0.0000): %d\n', numel(erd_ids));
    fprintf('ERS (>= 0.0000): %d\n', numel(ers_ids));
else
    fprintf('ERD (< %.4f): %d\n', ersd_split_threshold, numel(erd_ids));
    fprintf('ERS (>= %.4f): %d\n', ersd_split_threshold, numel(ers_ids));
end
if ~isempty(invalid_ids)
    fprintf('Excluded (invalid/pathological ERSD): %d\n', numel(invalid_ids));
end

if numel(erd_ids) < 2 || numel(ers_ids) < 2
    warning('Task %s: insufficient subjects per split group (ERD=%d, ERS=%d). Skipping task.', ...
        task_tag, numel(erd_ids), numel(ers_ids));
    continue
end

%% ERSD split inclusion figure
fprintf('\n=== Plotting ERSD split inclusion figure ===\n');
figure('Position', fig_pos, 'Color', 'w');
hold on
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
if strcmpi(task_tag, 'sternberg')
    yline(0, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
else
    yline(ersd_split_threshold, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
end
x_vals = (1:nSubj)';
idx_erd = ismember(uIDs, erd_ids);
idx_ers = ismember(uIDs, ers_ids);
idx_excl = ismember(uIDs, zero_ids) | ismember(uIDs, invalid_ids);
h_excl = scatter(x_vals(idx_excl), ersd_mean(idx_excl), 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
h_erd = scatter(x_vals(idx_erd), ersd_mean(idx_erd), 80, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.8);
h_ers = scatter(x_vals(idx_ers), ersd_mean(idx_ers), 80, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
legend_excl = sprintf('Excluded (n=%d)', sum(idx_excl));
legend_erd = sprintf('%s (n=%d)', tk.group_lbl_low, sum(idx_erd));
legend_ers = sprintf('%s (n=%d)', tk.group_lbl_high, sum(idx_ers));
xlabel('Participant (index)');
ylabel('ERSD [dB]');
if strcmpi(task_tag, 'nback')
    inclusion_title = sprintf('More vs Less ERD Split [%s | %s]', task_tag, split_label);
else
    inclusion_title = sprintf('ERS/ERD Split [%s | %s]', task_tag, split_label);
end
title(inclusion_title, 'Interpreter', 'none');
legend([h_excl, h_erd, h_ers], {legend_excl, legend_erd, legend_ers}, 'Location', 'best', 'FontSize', fontSize - 2, 'Box', 'off');
set(gca, 'FontSize', fontSize);
box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir_task, sprintf('%s_inclusion.png', fig_prefix)));
close(gcf);

%% Preallocate containers
metrics = struct();
metrics.ERSD = nan(nSubj, 3);         % for CSV export only
metrics.MS = nan(nSubj, 3);           % MSRateFullBL [% change]

% MS time-course pipeline (matches AOC_gaze_microsaccades_sternberg.m)
ms_cfg = init_ms_tc_cfg(cond_vals);
ms_tc_pct = nan(nSubj, 3, ms_cfg.n_samp);

missing_et = {};

%% Aggregate per-subject data
fprintf('\n=== Aggregating per-subject data (%d subjects) ===\n', nSubj);
for s = 1:nSubj
    clc; fprintf('[SPLIT ERS/ERD MS - %s] Aggregating data for Subject %d / %d \n', upper(task_tag), s, nSubj);
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
            metrics.ERSD(s, c) = mean(subj_rows.(tk.ersd_var)(cmask), 'omitnan');
            metrics.MS(s, c) = mean(subj_rows.MSRateFullBL(cmask), 'omitnan');
        end
    end

    et_file = fullfile(feat_dir, subj_folder, 'gaze', [tk.et_fname, '.mat']);
    if ~isfile(et_file)
        missing_et{end+1} = sid_str;
        continue
    end

    try
        subj_tc = build_ms_tc_pct_from_et(et_file, cond_vals, ms_cfg);
        ms_tc_pct(s, :, :) = subj_tc;
    catch ME
        missing_et{end+1} = sid_str;
        fprintf('  Warning: MS time course failed for %s (%s)\n', sid_str, ME.message);
    end
end

%% Define groups in row index space
is_red = ismember(uIDs, erd_ids);
is_amp = ismember(uIDs, ers_ids);

%% Rainclouds (MS only)
fprintf('\n=== Plotting MS rainclouds ===\n');
plot_ms_rainclouds(metrics.MS, is_red, is_amp, cond_labels, colors, fig_dir_task, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);

%% MS time courses with effect-size
close all
fprintf('\n=== Preparing microsaccade time courses ===\n');
t_vec = ms_cfg.t_vec;
idx_viable = (t_vec >= 0) & (t_vec <= 2);

% Subject-level QC (matches AOC_gaze_microsaccades_*.m summarize_subject_tc)
keep_tc = true(nSubj, 1);
for c = 1:3
    Xc = reshape(ms_tc_pct(:, c, :), nSubj, ms_cfg.n_samp);
    [Xc, keep_c] = preprocess_ms_subject_tc(Xc, idx_viable, ms_cfg);
    ms_tc_pct(:, c, :) = reshape(Xc, [nSubj, 1, ms_cfg.n_samp]);
    keep_tc = keep_tc & keep_c;
end

ms_ylabel = sprintf('Microsaccade\nRate Change [%%]');
ms_tag_base = 'MS_pct';
ms_csv_varname = 'MS_pct_0_2s';

close all
fontSizeTC = 40;
rng(123)
fs_ms = ms_cfg.fsample;
ds_factor = 10; % downsampling to 20 ms bins (500 Hz / 10)
tc_viz_smooth_sec = 0.05; % slight display-only smoothing for time-course plots

% Exclude subjects with incomplete MS time courses in analysis window.
tc_window_idx = idx_viable;
tc_complete_min_frac = 0.995;
tc_finite_frac = squeeze(mean(isfinite(ms_tc_pct(:, :, tc_window_idx)), 3));
tc_has_endpoint = squeeze(isfinite(ms_tc_pct(:, :, end)));
tc_complete_by_cond = (tc_finite_frac >= tc_complete_min_frac) & tc_has_endpoint;
tc_complete_subj = all(tc_complete_by_cond, 2) & keep_tc;
tc_excluded_subj = ~tc_complete_subj;
ms_tc_pct(tc_excluded_subj, :, :) = NaN;
is_red_tc = is_red & tc_complete_subj;
is_amp_tc = is_amp & tc_complete_subj;
fprintf('Time-course completeness filter: finite frac >= %.3f in [0,2]s + finite endpoint at 2s\n', tc_complete_min_frac);
fprintf('Excluded incomplete MS time-course subjects: %d\n', sum(tc_excluded_subj & (is_red | is_amp)));
if any(tc_excluded_subj & (is_red | is_amp))
    excl_ids = uIDs(tc_excluded_subj & (is_red | is_amp));
    fprintf('Excluded IDs: %s\n', sprintf('%d ', excl_ids));
end

cbpt_report_file = fullfile(stats_dir, sprintf('AOC_splitERSERD_MS_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_split_AlphaAmpRed_MS.m', ...
    'task_tag', task_tag, ...
    'split_label', split_label, ...
    'split_info', split_info_str, ...
    'ersd_var', tk.ersd_var, ...
    'group_lbl_low', tk.group_lbl_low, ...
    'group_lbl_high', tk.group_lbl_high, ...
    'n_low_split', numel(erd_ids), ...
    'n_high_split', numel(ers_ids), ...
    'n_low_tc', sum(is_red_tc), ...
    'n_high_tc', sum(is_amp_tc), ...
    'metric', 'Microsaccade rate change [% baseline]'));

plot_timecourse_with_effect_CBPT(ms_tc_pct, is_red_tc, is_amp_tc, colors, ...
    ms_ylabel, sprintf('%s_%s_%s', task_tag, split_label, ms_tag_base), fig_dir_root, fig_pos, fontSizeTC, fs_ms, ds_factor, t_vec, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
plot_timecourse_individuals(ms_tc_pct, is_red_tc, is_amp_tc, colors, ...
    ms_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_%s_%s_individuals_collapsed', task_tag, split_label, ms_tag_base), ...
    fig_dir_task, fig_pos, fontSizeTC, fs_ms, t_vec, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high);

% Additional condition-wise time course outputs
for c = 1:numel(cond_vals)
    tc_cond = ms_tc_pct(:, c, :);
    save_tag_cond = sprintf('%s_%s_%s_%s', task_tag, split_label, ms_tag_base, sanitize_label_for_fname(cond_labels{c}));
    plot_timecourse_with_effect_CBPT(tc_cond, is_red_tc, is_amp_tc, colors, ...
        ms_ylabel, save_tag_cond, fig_dir_task, fig_pos, fontSizeTC, fs_ms, ds_factor, t_vec, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
    plot_timecourse_individuals(tc_cond, is_red_tc, is_amp_tc, colors, ...
        ms_ylabel, cond_labels{c}, ...
        sprintf('%s_individuals', save_tag_cond), ...
        fig_dir_task, fig_pos, fontSizeTC, fs_ms, t_vec, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high);
end

%% Sanity checks
fprintf('\n=== Data Diagnostics [%s] ===\n', task_tag);
fprintf('Missing ET files for MS time courses: %d\n', numel(unique(missing_et)));
fprintf('Main figures saved to: %s\n', fig_dir_root);
fprintf('Supplementary figures saved to: %s\n', fig_dir_task);

%% Export CSV for Python statistics script
fprintf('\n=== Exporting CSV for Python stats ===\n');
t_win_lo = 0;
t_win_hi = 2;
task_idx_ms = (t_vec >= t_win_lo) & (t_vec <= t_win_hi);
ms_summary_by_load = squeeze(mean(ms_tc_pct(:, :, task_idx_ms), 3, 'omitnan'));

n_rows = nSubj * numel(cond_vals);
ID_col = nan(n_rows, 1);
LoadValue_col = nan(n_rows, 1);
LoadLabel_col = strings(n_rows, 1);
Group_col = strings(n_rows, 1);
Included_col = false(n_rows, 1);
ERSD_col = nan(n_rows, 1);
MS_col = nan(n_rows, 1);
MSSummary_col = nan(n_rows, 1);

r = 0;
for s = 1:nSubj
    for c = 1:numel(cond_vals)
        r = r + 1;
        ID_col(r) = uIDs(s);
        LoadValue_col(r) = cond_vals(c);
        LoadLabel_col(r) = string(cond_labels{c});
        ERSD_col(r) = metrics.ERSD(s, c);
        MS_col(r) = metrics.MS(s, c);
        MSSummary_col(r) = ms_summary_by_load(s, c);

        if is_red(s)
            Group_col(r) = tk.group_lbl_low;
            Included_col(r) = is_red_tc(s);
        elseif is_amp(s)
            Group_col(r) = tk.group_lbl_high;
            Included_col(r) = is_amp_tc(s);
        else
            Group_col(r) = "Excluded";
            Included_col(r) = false;
        end
    end
end

stats_tbl = table( ...
    ID_col, LoadValue_col, LoadLabel_col, Group_col, Included_col, ...
    ERSD_col, MS_col, MSSummary_col, ...
    'VariableNames', { ...
    'ID', 'LoadValue', 'LoadLabel', 'Group', 'Included', ...
    tk.ersd_var, 'MSRateFullBL', ms_csv_varname});

csv_out = fullfile(stats_dir, sprintf('AOC_splitERSERD_MS_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV saved to: %s\n', csv_out);
fprintf('CBPT report saved to: %s\n', cbpt_report_file);

end % task loop


%% ========================= Local Functions =========================
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

function ms_cfg = init_ms_tc_cfg(cond_vals)
% Pipeline parameters aligned with AOC_gaze_microsaccades_*.m
ms_cfg.fsample = 500;
ms_cfg.screenW = 800;
ms_cfg.screenH = 600;
ms_cfg.blink_win = 50;
ms_cfg.sigma_ms = 50;
sigma_samp = round(ms_cfg.sigma_ms / (1000 / ms_cfg.fsample));
kHalf = 3 * sigma_samp;
x_kern = -kHalf:kHalf;
ms_cfg.gKernel = exp(-x_kern.^2 / (2 * sigma_samp^2));
ms_cfg.gKernel = ms_cfg.gKernel / sum(ms_cfg.gKernel);
ms_cfg.t_comp = [-1.5 2.5];
ms_cfg.n_comp = round(diff(ms_cfg.t_comp) * ms_cfg.fsample) + 1;
ms_cfg.t_comp_vec = linspace(ms_cfg.t_comp(1), ms_cfg.t_comp(2), ms_cfg.n_comp);
ms_cfg.t_win = [-0.5 2];
[~, crop_start] = min(abs(ms_cfg.t_comp_vec - ms_cfg.t_win(1)));
[~, crop_end] = min(abs(ms_cfg.t_comp_vec - ms_cfg.t_win(2)));
ms_cfg.crop_idx = crop_start:crop_end;
ms_cfg.n_samp = numel(ms_cfg.crop_idx);
ms_cfg.t_vec = ms_cfg.t_comp_vec(ms_cfg.crop_idx);
ms_cfg.bl_win = [-1.5 -0.5];
ms_cfg.bl_idx_comp = ms_cfg.t_comp_vec >= ms_cfg.bl_win(1) & ms_cfg.t_comp_vec <= ms_cfg.bl_win(2);
ms_cfg.cond_vals = cond_vals;
ms_cfg.min_trials_per_cond = 3;
ms_cfg.outlier_k_iqr = 1.5;
ms_cfg.max_interp_gap_sec = 0.20;
ms_cfg.min_subject_coverage = 0.85;
ms_cfg.smooth_sec = 0.05;
ms_cfg.win_sm = max(1, round(ms_cfg.smooth_sec * ms_cfg.fsample));
end

function subj_tc = build_ms_tc_pct_from_et(et_file, cond_vals, ms_cfg)
nConds = numel(cond_vals);
subj_tc = nan(1, nConds, ms_cfg.n_samp);

S = load(et_file);
if isfield(S, 'dataETlong')
    et = S.dataETlong;
elseif isfield(S, 'dataet')
    et = S.dataet;
else
    error('No ET struct in %s.', et_file);
end

fsample = ms_cfg.fsample;
n_comp = ms_cfg.n_comp;
condSpikes = cell(nConds, 1);

for trl = 1:numel(et.trial)
    raw = et.trial{trl};
    t = et.time{trl};

    raw = raw(1:min(3, size(raw, 1)), :);
    raw(2, :) = ms_cfg.screenH - raw(2, :);

    oob = raw(1,:) < 0 | raw(1,:) > ms_cfg.screenW | ...
          raw(2,:) < 0 | raw(2,:) > ms_cfg.screenH;
    raw(:, oob) = NaN;

    raw = remove_blinks(raw, ms_cfg.blink_win);

    idx_win = t >= ms_cfg.t_comp(1) & t <= ms_cfg.t_comp(2);
    gx = raw(1, idx_win);
    gy = raw(2, idx_win);

    if sum(isfinite(gx) & isfinite(gy)) < 50
        continue
    end

    velData = [gx; gy];
    [~, msDetails] = detect_microsaccades(fsample, velData, length(gx));

    spikeVec = zeros(1, length(gx));
    if ~isempty(msDetails.Onset)
        onsets = msDetails.Onset(msDetails.Onset >= 1 & msDetails.Onset <= length(gx));
        spikeVec(onsets) = 1;
    end

    if length(spikeVec) >= n_comp
        spikeVec = spikeVec(1:n_comp);
    else
        spikeVec(end+1:n_comp) = 0;
    end

    if size(et.trialinfo, 2) > 1
        cond = et.trialinfo(trl, 1) - 20;
    else
        cond = et.trialinfo(trl) - 20;
    end

    cIdx = find(cond_vals == cond, 1);
    if isempty(cIdx)
        continue
    end

    condSpikes{cIdx}(end+1, :) = spikeVec;
end

for c = 1:nConds
    if size(condSpikes{c}, 1) >= ms_cfg.min_trials_per_cond
        rate = mean(condSpikes{c}, 1) * fsample;
        smoothed = conv(rate, ms_cfg.gKernel, 'same');
        bl_mean = nanmean(smoothed(ms_cfg.bl_idx_comp));
        if isfinite(bl_mean) && bl_mean > 0
            smoothed_pct = (smoothed - bl_mean) ./ bl_mean * 100;
        else
            smoothed_pct = nan(size(smoothed));
        end
        subj_tc(1, c, :) = smoothed_pct(ms_cfg.crop_idx);
    end
end
end

function [X, keep_subj] = preprocess_ms_subject_tc(X, idx_viable, ms_cfg)
X(~isfinite(X)) = NaN;
med_metric = median(X(:, idx_viable), 2, 'omitnan');
[X, keep_subj] = exclude_outlier_trajectories(X, med_metric, ms_cfg.outlier_k_iqr);
max_interp_gap_smp = max(1, round(ms_cfg.max_interp_gap_sec * ms_cfg.fsample));
for s = 1:size(X, 1)
    X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
end
subj_cov = mean(isfinite(X(:, idx_viable)), 2);
keep_subj = keep_subj & (subj_cov >= ms_cfg.min_subject_coverage);
X(~keep_subj, :) = NaN;
if ms_cfg.win_sm > 1
    X = movmean(X, ms_cfg.win_sm, 2, 'omitnan');
end
end

function [X_keep, keep_subj] = exclude_outlier_trajectories(X, metric, k_iqr)
med_m = median(metric, 'omitnan');
iqr_m = iqr(metric);
if ~isfinite(iqr_m) || iqr_m == 0
    low_m = -inf;
    high_m = inf;
else
    low_m = med_m - k_iqr * iqr_m;
    high_m = med_m + k_iqr * iqr_m;
end
keep_subj = isfinite(metric) & (metric >= low_m) & (metric <= high_m);
X_keep = X;
X_keep(~keep_subj, :) = NaN;
end

function x = fill_short_nan_gaps(x, max_gap_smp)
if isempty(x)
    return
end
valid = isfinite(x);
if all(~valid) || all(valid)
    return
end

n = numel(x);
i = 1;
while i <= n
    if valid(i)
        i = i + 1;
        continue
    end
    j = i;
    while j <= n && ~valid(j)
        j = j + 1;
    end
    gap_start = i;
    gap_end = j - 1;
    gap_len = gap_end - gap_start + 1;

    left = gap_start - 1;
    right = gap_end + 1;
    if left >= 1 && right <= n && valid(left) && valid(right) && gap_len <= max_gap_smp
        x(gap_start:gap_end) = interp1([left right], [x(left) x(right)], gap_start:gap_end);
    end
    i = j;
end
end

function plot_ms_rainclouds(X, is_red, is_amp, cond_labels, colors, fig_dir, fig_prefix, fig_pos, fsz, group_lbl_low, group_lbl_high)
if nargin < 10 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 11 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
all_vals = X(isfinite(X));
if isempty(all_vals)
    ymax = 1;
else
    ymax = max(abs(all_vals));
    if ymax <= 0
        ymax = 1;
    end
    ymax = ymax * 1.15;
end
ylim_shared = [-ymax ymax];

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile; hold on
for c = 1:3
    draw_one_cloud(X(is_red, c), c, colors(c,:), 0.3, 96, 0.45);
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim(ylim_shared);
set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
title(group_lbl_low, 'Interpreter', 'none');
ylabel('Microsaccade rate', 'Interpreter', 'none');
box off

nexttile; hold on
for c = 1:3
    draw_one_cloud(X(is_amp, c), c, colors(c,:), 0.3, 96, 0.45);
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim(ylim_shared);
set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
title(group_lbl_high, 'Interpreter', 'none');
ylabel('Microsaccade rate', 'Interpreter', 'none');
box off

sgtitle('Microsaccade rate', 'FontSize', fsz+2, 'Interpreter', 'none');
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_raincloud_ms.png', fig_prefix)));
close(gcf);
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

function plot_timecourse_individuals(tc, is_red, is_amp, colors, ylab, title_tag, save_tag, fig_dir, fig_pos, fsz, fs, t_vec, smooth_sec, group_lbl_low, group_lbl_high)
if nargin < 14 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 15 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
nT = size(tc, 3);
if nargin < 12 || isempty(t_vec)
    dt = 1 / fs;
    t_plot = linspace(-0.5 + dt/2, 2 - dt/2, nT);
else
    t_plot = t_vec(:)';
    if numel(t_plot) ~= nT
        error('t_vec length (%d) must match time dimension (%d).', numel(t_plot), nT);
    end
end
if nargin < 13 || isempty(smooth_sec)
    smooth_sec = 0;
end
win_sm = max(1, round(smooth_sec * fs));
R = group_tc_by_time(tc(is_red, :, :));
A = group_tc_by_time(tc(is_amp, :, :));
R = as_subjects_by_time(R, nT);
A = as_subjects_by_time(A, nT);
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
    plot(t_plot, mean(R, 1, 'omitnan'), 'Color', colors(1, :), 'LineWidth', 2.5);
end
xline(0, '--k');
xlim([-0.5 2]);
ylabel(ylab);
title(sprintf('%s (n=%d) - %s', group_lbl_low, size(R, 1), title_tag), 'Interpreter', 'none');
set(gca, 'FontSize', fsz - 6);
box off
nexttile; hold on
if ~isempty(A)
    plot(t_plot, A', 'Color', colA_light, 'LineWidth', 0.8);
    plot(t_plot, mean(A, 1, 'omitnan'), 'Color', colors(3, :), 'LineWidth', 2.5);
end
xline(0, '--k');
xlim([-0.5 2]);
xlabel('Time [s]');
ylabel(ylab);
title(sprintf('%s (n=%d) - %s', group_lbl_high, size(A, 1), title_tag), 'Interpreter', 'none');
set(gca, 'FontSize', fsz - 6);
box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_MS_timecourse_%s.png', save_tag)));
close(gcf);
end

function plot_timecourse_with_effect_CBPT(tc, is_red, is_amp, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, t_vec, smooth_sec, group_lbl_low, group_lbl_high, cbpt_report_file)
if nargin < 11 || isempty(ds_factor), ds_factor = 10; end
if nargin < 13 || isempty(smooth_sec), smooth_sec = 0.05; end
if nargin < 14 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 15 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
if nargin < 16, cbpt_report_file = ''; end
nT = size(tc, 3);
if nargin < 12 || isempty(t_vec)
    dt = 1 / fs;
    t_plot = linspace(-0.5 + dt/2, 2 - dt/2, nT);
else
    t_plot = t_vec(:)';
    if numel(t_plot) ~= nT
        error('t_vec length (%d) must match time dimension (%d).', numel(t_plot), nT);
    end
end
dt = mean(diff(t_plot), 'omitnan');
tc(~isfinite(tc)) = NaN;
Rall = group_tc_by_time(tc(is_red, :, :));
Aall = group_tc_by_time(tc(is_amp, :, :));
Rall = as_subjects_by_time(Rall, nT);
Aall = as_subjects_by_time(Aall, nT);
win_sm = max(1, round(smooth_sec * fs));
if win_sm > 1
    Rall = movmean(Rall, win_sm, 2, 'omitnan');
    Aall = movmean(Aall, win_sm, 2, 'omitnan');
end
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');
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
set(e1.mainLine, 'Color', colors(1,:), 'LineWidth', 2.5);
set(e2.mainLine, 'Color', colors(3,:), 'LineWidth', 2.5);
set(e1.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.20);
set(e2.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.20);
set(e1.edge(1), 'Color', 'none');
set(e1.edge(2), 'Color', 'none');
set(e2.edge(1), 'Color', 'none');
set(e2.edge(2), 'Color', 'none');
xline(0, '--k');
ylabel(ylab);
xlim([-0.5 2]);
box off
set(gca, 'FontSize', fsz-4);
leg_p1 = patch(NaN, NaN, colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(1,:), 'LineWidth', 1.5);
leg_p2 = patch(NaN, NaN, colors(3,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(3,:), 'LineWidth', 1.5);
legend([leg_p1 leg_p2], {[' ' group_lbl_low], [' ' group_lbl_high]}, 'Location', 'best', 'FontSize', fsz*0.75, 'Box', 'off');
if contains(save_tag, 'MS_pct')
    upper = max(mR + sR, mA + sA);
    lower = min(mR - sR, mA - sA);
    y_hi = max(upper, [], 'omitnan');
    y_lo = min(lower, [], 'omitnan');
    if isfinite(y_hi) && isfinite(y_lo) && y_hi > y_lo
        pad = max(0.15, 0.05 * (y_hi - y_lo));
        ylim([y_lo - pad, y_hi + pad]);
    end
end
ax_d = nexttile; hold on
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
if numel(t_plot_ds) ~= nT_ds
    nT_ds = min(numel(t_plot_ds), nT_ds);
    t_plot_ds = t_plot_ds(1:nT_ds);
    Rall_ds = Rall_ds(:, 1:nT_ds);
    Aall_ds = Aall_ds(:, 1:nT_ds);
end
post_idx = t_plot_ds >= 0 & t_plot_ds <= 2;
Rall_cbpt = as_subjects_by_time(Rall_ds(:, post_idx), sum(post_idx));
Aall_cbpt = as_subjects_by_time(Aall_ds(:, post_idx), sum(post_idx));
t_cbpt = t_plot_ds(post_idx);
dt_ds = ds_factor * dt;
alpha_cbpt = 0.05;
tail_cbpt = 'twotail';
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d(Rall_cbpt, Aall_cbpt, n_perm, alpha_cbpt, tail_cbpt, t_cbpt);
report_cfg = struct( ...
    'tag', save_tag, ...
    'modality', 'MS', ...
    'nR', size(Rall_cbpt, 1), ...
    'nA', size(Aall_cbpt, 1), ...
    'lbl_low', group_lbl_low, ...
    'lbl_high', group_lbl_high, ...
    'n_perm', n_perm, ...
    'alpha', alpha_cbpt, ...
    'tail', tail_cbpt, ...
    'nT_ds', numel(t_cbpt), ...
    'bin_ms', ds_factor * 1000 / fs, ...
    'fs', fs, ...
    'ds_factor', ds_factor, ...
    'clusters', clusters, ...
    'tvals', tvals_cl, ...
    'thr', thr, ...
    't_plot', t_cbpt, ...
    'dt_ds', dt_ds, ...
    'maxMassNull', [], ...
    'maxExtentNull', []);
log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(report_cfg));
sig_cluster = false(1, nT_ds);
post_pos = find(post_idx);
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        sig_cluster(post_pos(clusters(k).idx)) = true;
    end
end
sig_uncorr = false(1, nT_ds);
sig_uncorr(post_pos) = (abs(tvals_cl) > thr.tcrit) & isfinite(tvals_cl);
sig = sig_cluster;
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
box off
set(gca, 'FontSize', fsz-4);
drawnow;
align_stacked_tc_panels(ax_d, ax_gaze);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_MS_timecourse_%s_CBPT.png', save_tag)));
close(gcf);
end

function align_stacked_tc_panels(ax_ref, ax_top)
% Match ylabel horizontal position to the reference axis (Cohen's d panel).
% TiledChartLayout owns axis positions, so only the label is moved.
drawnow;
align_ylabel_to_reference(ax_ref, ax_top);
end

function align_ylabel_to_reference(ax_ref, ax_targets)
fig = ancestor(ax_ref, 'figure');
if isempty(fig)
    return
end
drawnow;

ref_lbl = ax_ref.YLabel;
ref_lbl.Units = 'normalized';
ref_pos = ref_lbl.Position;
ax_ref.Units = 'pixels';
ref_ax_pos_pix = ax_ref.Position;
target_x_pix = ref_ax_pos_pix(1) + ref_pos(1) * ref_ax_pos_pix(3);

if iscell(ax_targets)
    ax_list = ax_targets;
elseif isscalar(ax_targets)
    ax_list = {ax_targets};
else
    ax_list = num2cell(ax_targets);
end
for k = 1:numel(ax_list)
    ax_t = ax_list{k};
    if isequal(ax_t, ax_ref)
        continue
    end
    lbl = ax_t.YLabel;
    lbl.Units = 'normalized';
    pos = lbl.Position;
    ax_t.Units = 'pixels';
    ax_t_pos_pix = ax_t.Position;
    x_norm_new = (target_x_pix - ax_t_pos_pix(1)) / max(ax_t_pos_pix(3), eps);
    lbl.Position = [x_norm_new, pos(2), pos(3)];
end
end

function init_cbpt_report_file(report_path, meta)
if isempty(report_path)
    return
end
if isfile(report_path)
    delete(report_path);
end
lines = {};
lines{end+1} = '=== AOC Split ERS/ERD CBPT Report ===';
lines{end+1} = sprintf('Script: %s', meta.script);
lines{end+1} = sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
lines{end+1} = sprintf('Task: %s | Split: %s', meta.task_tag, meta.split_label);
lines{end+1} = sprintf('ERSD variable: %s', meta.ersd_var);
lines{end+1} = meta.split_info;
lines{end+1} = sprintf('Split groups: %s n=%d; %s n=%d', meta.group_lbl_low, meta.n_low_split, meta.group_lbl_high, meta.n_high_split);
lines{end+1} = sprintf('Time-course analysis (after completeness filter): %s n=%d; %s n=%d', ...
    meta.group_lbl_low, meta.n_low_tc, meta.group_lbl_high, meta.n_high_tc);
lines{end+1} = sprintf('Outcome metric: %s', meta.metric);
lines{end+1} = 'CBPT method: FieldTrip ft_timelockstatistics (montecarlo, cluster maxsum, indepsamplesT)';
lines{end+1} = 'Defaults per analysis block: n_perm=10000, clusteralpha=0.05, two-tailed, latency=[0 2] s post-stimulus, ds_factor=10 (20 ms bins at 500 Hz)';
lines{end+1} = '';
append_lines_to_file(report_path, lines);
end

function lines = build_cbpt_report_lines(R)
lbl_lo = sanitize_label_for_fname(R.lbl_low);
lbl_hi = sanitize_label_for_fname(R.lbl_high);
nExtreme = sum(abs(R.tvals) > R.thr.tcrit & isfinite(R.tvals));
if isempty(R.clusters)
    maxClMass = 0;
    maxClExtent = 0;
else
    maxClMass = max([0, arrayfun(@(k) R.clusters(k).mass, 1:numel(R.clusters))]);
    maxClExtent = max([0, arrayfun(@(k) R.clusters(k).extent, 1:numel(R.clusters))]);
end
tail_str = tail_label_for_report(R.tail);
lines = {};
lines{end+1} = sprintf('  [%s] n_%s=%d n_%s=%d tcrit=%.2f (%s %s vs %s) |t|>tcrit at %d timepts; max cluster mass=%.1f; max cluster extent=%d', ...
    R.tag, lbl_lo, R.nR, lbl_hi, R.nA, R.thr.tcrit, tail_str, R.lbl_low, R.lbl_high, nExtreme, maxClMass, maxClExtent);
lines{end+1} = sprintf('    Method: n_perm=%d, alpha=%.3f, %s, latency=[0 2] s, ds_factor=%d (%.0f ms bins at %d Hz)', ...
    R.n_perm, R.alpha, tail_str, R.ds_factor, R.bin_ms, R.fs);
if ~isempty(R.maxMassNull)
    lines{end+1} = sprintf('    CBPT %s: nT_ds=%d (%.0f ms bins); null mass: median=%.1f, 90th=%.1f, 95th=%.1f; null extent: median=%d, 90th=%d, 95th=%d', ...
        lower(R.modality), R.nT_ds, R.bin_ms, median(R.maxMassNull), prctile(R.maxMassNull, 90), prctile(R.maxMassNull, 95), ...
        median(R.maxExtentNull), prctile(R.maxExtentNull, 90), prctile(R.maxExtentNull, 95));
else
    lines{end+1} = sprintf('    CBPT %s: nT_ds=%d (%.0f ms bins); FieldTrip ft_timelockstatistics (cluster mass)', ...
        lower(R.modality), R.nT_ds, R.bin_ms);
end
if ~isempty(R.clusters)
    bin_width_ms = R.dt_ds * 1000;
    for k = 1:numel(R.clusters)
        idx = R.clusters(k).idx;
        t_lo = R.t_plot(idx(1));
        t_hi = R.t_plot(idx(end));
        t_start = t_lo - R.dt_ds / 2;
        t_end = t_hi + R.dt_ds / 2;
        duration_ms = R.clusters(k).extent * bin_width_ms;
        status = 'n.s.';
        if R.clusters(k).p < R.alpha
            status = 'SIGNIFICANT';
        end
        lines{end+1} = sprintf('    Cluster %d: window [%.3f, %.3f] s post-stim (duration=%d ms, %d bins @ %.0f ms); mass=%.1f; p=%.4f; %s', ...
            k, t_start, t_end, round(duration_ms), R.clusters(k).extent, bin_width_ms, R.clusters(k).mass, R.clusters(k).p, status);
    end
    if isfinite(R.thr.mass)
        lines{end+1} = sprintf('    Gap to significance: mass need +%.1f (have %.1f, need %.1f); extent need +%d (have %d, need %d)', ...
            R.thr.mass - maxClMass, maxClMass, R.thr.mass, R.thr.extent - maxClExtent, maxClExtent, R.thr.extent);
    end
else
    t_fin = R.tvals(isfinite(R.tvals));
    if ~isempty(t_fin)
        [t_max, idx_max] = max(abs(R.tvals));
        lines{end+1} = sprintf('    No clusters formed (no contiguous |t|>tcrit); largest |t|=%.2f at t=%.2f s', ...
            t_max, R.t_plot(idx_max));
    else
        lines{end+1} = '    No clusters; no valid t-values';
    end
end
lines{end+1} = '';
end

function tail_str = tail_label_for_report(tail)
if strcmpi(tail, 'onetail_pos')
    tail_str = 'one-tailed positive';
elseif strcmpi(tail, 'onetail_neg')
    tail_str = 'one-tailed negative';
else
    tail_str = 'two-tailed';
end
end

function log_cbpt_report(report_path, lines)
for i = 1:numel(lines)
    fprintf('%s\n', lines{i});
end
append_lines_to_file(report_path, lines);
end

function append_lines_to_file(file_path, lines)
if isempty(file_path) || isempty(lines)
    return
end
fid = fopen(file_path, 'a');
if fid < 0
    warning('Could not open CBPT report file: %s', file_path);
    return
end
cleanup = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end

function M = as_subjects_by_time(M, n_time)
% Force nSubjects x nTime orientation for FieldTrip inputs.
if nargin < 2 || isempty(n_time)
    error('n_time is required.');
end
if isempty(M)
    M = zeros(0, n_time);
    return
end
if isvector(M)
    v = M(:)';
    if numel(v) ~= n_time
        error('Vector length (%d) does not match n_time (%d).', numel(v), n_time);
    end
    M = v;
    return
end
if size(M, 2) == n_time
    return
elseif size(M, 1) == n_time
    M = M.';
else
    error('Matrix [%d x %d] cannot be oriented to n_time=%d.', size(M, 1), size(M, 2), n_time);
end
end

function M = group_tc_by_time(tc_grp)
% Return nSubjects x nTime matrix (avoids squeeze orientation bugs).
nG = size(tc_grp, 1);
nT = size(tc_grp, 3);
G = mean(tc_grp, 2, 'omitnan');
M = reshape(G, nG, nT);
end

function [clusters, tvals, thr, maxMassNull, maxExtentNull] = ft_cluster_permutation_1d(Rall, Aall, nPerm, alpha, tail, t_plot_ds)
if nargin < 5, tail = 'twotail'; end
if nargin < 6, t_plot_ds = []; end

nT = size(Rall, 2);
Rall = as_subjects_by_time(Rall, nT);
Aall = as_subjects_by_time(Aall, nT);
nR = size(Rall, 1);
nA = size(Aall, 1);
nT = size(Rall, 2);
if size(Aall, 2) ~= nT
    error('Rall and Aall must have the same number of time points.');
end
df = nR + nA - 2;

for t = 1:nT
    r = Rall(:, t); a = Aall(:, t);
    mr = mean(r(isfinite(r)), 'omitnan');
    ma = mean(a(isfinite(a)), 'omitnan');
    if ~isfinite(mr), mr = 0; end
    if ~isfinite(ma), ma = 0; end
    r(~isfinite(r)) = mr;
    a(~isfinite(a)) = ma;
    Rall(:, t) = r;
    Aall(:, t) = a;
end

chan_label = 'metric';
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    time_vec = t_plot_ds(:)';
else
    time_vec = (0:nT-1) / 500;
end

tl1 = build_ft_timelock_1d(Rall, time_vec, chan_label);
tl2 = build_ft_timelock_1d(Aall, time_vec, chan_label);

cfg_neigh = struct();
cfg_neigh(1).label = chan_label;
cfg_neigh(1).neighblabel = {};

cfg = struct();
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = cfg_neigh;
cfg.numrandomization = nPerm;
cfg.channel = chan_label;
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    cfg.latency = [0 2];
else
    cfg.latency = 'all';
end

if strcmpi(tail, 'onetail_pos')
    cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = alpha;
elseif strcmpi(tail, 'onetail_neg')
    cfg.tail = -1; cfg.clustertail = -1; cfg.alpha = alpha;
else
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = alpha / 2;
end
cfg.design = [ones(1, nR), 2*ones(1, nA)];
cfg.ivar = 1;

stat = ft_timelockstatistics(cfg, tl1, tl2);
tvals = stat.stat(1, :);

post_idx = find(t_plot_ds >= 0 & t_plot_ds <= 2);
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT && size(tvals, 2) < nT && ~isempty(post_idx)
    tvals_full = nan(1, nT);
    n_sel = min(size(tvals, 2), numel(post_idx));
    tvals_full(post_idx(1:n_sel)) = tvals(1:n_sel);
    tvals = tvals_full;
end

clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {}, 'p_extent', {});
has_neg = isfield(stat, 'negclusters') && ~isempty(stat.negclusters);
has_pos = isfield(stat, 'posclusters') && ~isempty(stat.posclusters);
if strcmpi(tail, 'onetail_pos')
    if has_pos, clist = stat.posclusters; lmat = stat.posclusterslabelmat; else, clist = struct('prob', {}); lmat = zeros(1, nT); end
elseif strcmpi(tail, 'onetail_neg')
    if has_neg, clist = stat.negclusters; lmat = stat.negclusterslabelmat; else, clist = struct('prob', {}); lmat = zeros(1, nT); end
else
    clist = []; lmat = zeros(1, nT);
    if has_pos, clist = stat.posclusters(:); lmat = stat.posclusterslabelmat; end
    if has_neg
        lmat_neg = stat.negclusterslabelmat;
        npos = numel(clist);
        lmat(lmat_neg > 0) = npos + lmat_neg(lmat_neg > 0);
        if isempty(clist), clist = stat.negclusters(:); else, clist = [clist; stat.negclusters(:)]; end
    end
end
for k = 1:numel(clist)
    idx = find(lmat(1, :) == k);
    if isempty(idx), continue, end
    idx = idx(1):idx(end);
    clusters(end+1).idx = idx;
    clusters(end).mass = sum(abs(tvals(idx)), 'omitnan');
    clusters(end).extent = numel(idx);
    clusters(end).p = clist(k).prob;
    clusters(end).p_extent = clist(k).prob;
end
if strcmpi(tail, 'onetail_pos') || strcmpi(tail, 'onetail_neg')
    thr.tcrit = tinv(1 - alpha, df);
else
    thr.tcrit = tinv(1 - alpha/2, df);
end
thr.mass = NaN; thr.extent = NaN;
maxMassNull = []; maxExtentNull = [];
end

function tl = build_ft_timelock_1d(data_nt, time_vec, chan_label)
[nR, nT] = size(data_nt);
if numel(time_vec) ~= nT
    error('Time vector length (%d) must match data columns (%d).', numel(time_vec), nT);
end
tl = struct();
tl.label = {chan_label};
tl.time = time_vec(:)';
tl.dimord = 'rpt_chan_time';
tl.trial = zeros(nR, 1, nT);
for k = 1:nR
    tl.trial(k, 1, :) = data_nt(k, :);
end
if ndims(tl.trial) ~= 3 || size(tl.trial, 2) ~= 1
    error('FieldTrip trial array must be nR x 1 x nT.');
end
end
