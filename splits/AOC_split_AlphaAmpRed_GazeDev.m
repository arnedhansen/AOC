%% AOC Split ERS/ERD (Subject-Level) — GazeDev (Sternberg + N-back)
% Subject-level ERSD split per task, then gaze deviation time courses
% (percent baseline) with cluster-based permutation testing.
%
% Sternberg: split0 on mean ERSD_late ([1 2]s), < 0 -> ERD, >= 0 -> ERS
% N-back: median split on pooled ERSD_full, < median -> More ERD, >= median -> Less ERD
%
% Generates (per task):
% - ERSD split inclusion figure
% - TFRs per condition for both groups (baselined, 8-14 Hz)
% - Topoplots per condition for both groups (ERSD window on baselined TFR)
% - Rainclouds for ERSD and gaze deviation
% - Time-course panels for gaze deviation (percent baseline) and ERSD: collapsed and per condition

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
pathAOC = paths.features;

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
fprintf('\n=== AOC Split ERS/ERD — Gaze Deviation (Sternberg + N-back) ===\n');
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
tasks(1).tfr_fname = 'tfr_stern.mat';
tasks(1).tfr_fname_alt = '';
tasks(1).tfr_vars = {'tfr2_bl', 'tfr4_bl', 'tfr6_bl'};
tasks(1).ersd_tc_fname = 'ersd_sternberg_timecourse.mat';
tasks(1).topo_latency = [1 2];
tasks(1).gaze_fname = 'gaze_series_sternberg_trials.mat';
tasks(1).ersd_var = 'ERSD_late';
tasks(1).fig_subdir = 'SternbergGazeDev';
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
tasks(2).tfr_fname = 'tfr_nback.mat';
tasks(2).tfr_fname_alt = 'tfr_nback_long.mat';
tasks(2).tfr_vars = {'tfr1_bl', 'tfr2_bl', 'tfr3_bl'};
tasks(2).ersd_tc_fname = 'ersd_nback_timecourse.mat';
tasks(2).topo_latency = [1 2];
tasks(2).gaze_fname = 'gaze_series_nback_trials.mat';
tasks(2).ersd_var = 'ERSD_full';
tasks(2).fig_subdir = 'NbackGazeDev';
tasks(2).group_lbl_low = 'More ERD';
tasks(2).group_lbl_high = 'Less ERD';

fontSize = 40;
tfr_winsor_cfg = struct();
tfr_winsor_cfg.enable = true;
tfr_winsor_cfg.prctile = [2 98]; % subject-level clipping per TF bin

for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
split_label = tk.split_label;
cond_vals = tk.cond_vals;
cond_codes = tk.cond_codes;
cond_labels = tk.cond_labels;
fig_prefix = sprintf('AOC_splitERSERD_GazeDev_%s_%s', task_tag, split_label);
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
    fprintf('%s (< %.4f): %d\n', tk.group_lbl_low, ersd_split_threshold, numel(erd_ids));
    fprintf('%s (>= %.4f): %d\n', tk.group_lbl_high, ersd_split_threshold, numel(ers_ids));
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
% EEG TFR (baselined; tfr_*_subj maps cell entries to subject indices)
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
metrics.ERSD = nan(nSubj, 3);
metrics.Dev = nan(nSubj, 3);          % GazeDeviationFullBL [% change]

% Time courses per subject x condition (full resolution, pre-outlier exclusion)
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
dev_tc = nan(nSubj, 3, Tf);
ersd_tc = nan(nSubj, 3, Tf);

missing_ersd = {};
missing_tfr = {};
missing_gaze = {};

%% Aggregate per-subject data
fprintf('\n=== Aggregating per-subject data (%d subjects) ===\n', nSubj);
for s = 1:nSubj
    clc; fprintf('[SPLIT ERS/ERD GAZEDEV - %s] Aggregating data for Subject %d / %d \n', upper(task_tag), s, nSubj);
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
            metrics.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
        end
    end

    % EEG TFR sources (baselined; no FOOOF)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    tfr_file = fullfile(eeg_dir, tk.tfr_fname);
    if ~isfile(tfr_file) && isfield(tk, 'tfr_fname_alt') && ~isempty(tk.tfr_fname_alt)
        tfr_file = fullfile(eeg_dir, tk.tfr_fname_alt);
    end
    try
        if ~isfile(tfr_file)
            error('Missing TFR file.');
        end
        R = load(tfr_file, tk.tfr_vars{:});
        tfr_conds = {R.(tk.tfr_vars{1}), R.(tk.tfr_vars{2}), R.(tk.tfr_vars{3})};
        if ismember(sid, erd_ids)
            for c = 1:3
                tfr_red{c}{end+1} = tfr_conds{c};
            end
            tfr_red_subj(end+1) = s;
        elseif ismember(sid, ers_ids)
            for c = 1:3
                tfr_amp{c}{end+1} = tfr_conds{c};
            end
            tfr_amp_subj(end+1) = s;
        end
    catch
        missing_tfr{end+1} = sid_str;
    end

    % ERSD time course (saved during feature extraction)
    try
        E = load(fullfile(eeg_dir, tk.ersd_tc_fname), 'ersd_timecourse');
        tcStruct = E.ersd_timecourse;
        t_ersd = tcStruct.time(:)';
        tcMat = tcStruct.ersd_occ_8_14_db;
        cond_ersd = tcStruct.condition(:);
        for c = 1:3
            row_idx = find(cond_ersd == cond_vals(c), 1);
            if isempty(row_idx)
                continue
            end
            tc_row = double(tcMat(row_idx, :));
            ersd_tc(s, c, :) = interp1(t_ersd, tc_row, t_plot, 'linear', NaN);
        end
    catch
        missing_ersd{end+1} = sid_str;
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

    % ET coverage diagnostic for the requested 0-2 s window (per subject).
    if isfield(G, 'time') && ~isempty(G.time)
        tmax_trials = nan(1, numel(G.time));
        for tr = 1:numel(G.time)
            tt_cov = G.time{tr};
            if isempty(tt_cov), continue, end
            tmax_trials(tr) = tt_cov(end);
        end
        n_reach_2s = sum(tmax_trials >= 2, 'omitnan');
        fprintf('ET coverage %s: %d/%d trials reach 2 s (tmax median=%.3f s, min=%.3f s)\n', ...
            sid_str, n_reach_2s, numel(tmax_trials), median(tmax_trials, 'omitnan'), min(tmax_trials, [], 'omitnan'));
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
                    tt = linspace(-0.5, 2, numel(x));
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
is_red = ismember(uIDs, erd_ids);
is_amp = ismember(uIDs, ers_ids);

eeg_tc = [];
%% Determine occipital channels (from first available baselined TFR)
channels = {};
for c = 1:3
    if ~isempty(tfr_red{c})
        channels = occ_channels_from_labels(tfr_red{c}{1}.label);
        break
    elseif ~isempty(tfr_amp{c})
        channels = occ_channels_from_labels(tfr_amp{c}{1}.label);
        break
    end
end
if isempty(channels)
    error('No baselined TFR data available for channel definition.');
end

%% TFRs per condition (both groups, 3x2)
fprintf('\n=== Plotting TFRs ===\n');
color_map_tfr = customcolormap_preset('red-white-blue');
plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, ...
    color_map_tfr, fig_dir_task, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg, tk.group_lbl_low, tk.group_lbl_high);

%% TFRs collapsed over conditions (both groups)
fprintf('\n=== Plotting collapsed TFRs (across conditions) ===\n');
plot_group_tfrs_collapsed(tfr_red, tfr_amp, channels, headmodel, color_map_tfr, ...
    fig_dir_task, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg, tk.group_lbl_low, tk.group_lbl_high);

%% Topoplots per condition (both groups; ERSD window on baselined TFR)
fprintf('\n=== Plotting topoplots ===\n');
plot_group_topos_ersd(tfr_red, tfr_amp, channels, headmodel, cond_labels, tk.topo_latency, fig_dir_task, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);

%% Topoplots collapsed over conditions (both groups)
fprintf('\n=== Plotting collapsed topoplots (across conditions) ===\n');
plot_group_topos_ersd_collapsed(tfr_red, tfr_amp, channels, headmodel, tk.topo_latency, fig_dir_task, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);

%% Rainclouds
fprintf('\n=== Plotting rainclouds ===\n');
plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir_task, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);

%% EEG and Gaze time courses with effect-size
close all
Tf = size(dev_tc, 3);
fprintf('\n=== Preparing ERSD time courses ===\n');
eeg_tc = ersd_tc;

% Gaze time-course representation for analysis/plotting.
fprintf('\n=== Preparing gaze time courses ===\n');
t_vec = linspace(-0.5, 2, Tf);
bl_idx = (t_vec >= -1.5) & (t_vec <= -0.5);

dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
% Percent baseline: (value/baseline - 1) * 100.
gaze_tc = (dev_tc ./ dev_bl_3d - 1) * 100;
gaze_tc(~isfinite(gaze_tc)) = NaN;
gaze_ylabel = sprintf('Gaze Deviation\nChange [%%]');
gaze_tag_base = 'GazeDev_pct';
gaze_csv_varname = 'GazeDev_pct_0_2s';
fprintf('Using baselined gaze measures (percent change) for task %s.\n', task_tag);

% Plot time courses (always save both: gaze-only and combined EEG+gaze)
close all
fontSizeTC = 40;
rng(123)
ds_factor = 10; % downsampling to ds_factor*2ms windows
tc_viz_smooth_sec = 0.05; % slight display-only smoothing for time-course plots

% Exclude subjects with incomplete gaze time courses in analysis window.
% Criterion: for each condition, signal must be finite at t=2 s and have
% >=99.5% finite samples in [0, 2] s.
tc_window_idx = (t_vec >= 0) & (t_vec <= 2);
tc_complete_min_frac = 0.995;
tc_finite_frac = squeeze(mean(isfinite(dev_tc(:, :, tc_window_idx)), 3));
tc_has_endpoint = squeeze(isfinite(dev_tc(:, :, end)));
tc_complete_by_cond = (tc_finite_frac >= tc_complete_min_frac) & tc_has_endpoint;
tc_complete_subj = all(tc_complete_by_cond, 2);
tc_excluded_subj = ~tc_complete_subj;
gaze_tc(tc_excluded_subj, :, :) = NaN;
eeg_tc(tc_excluded_subj, :, :) = NaN;
is_red_tc = is_red & tc_complete_subj;
is_amp_tc = is_amp & tc_complete_subj;
fprintf('Time-course completeness filter: finite frac >= %.3f in [0,2]s + finite endpoint at 2s\n', tc_complete_min_frac);
fprintf('Excluded incomplete gaze time-course subjects: %d\n', sum(tc_excluded_subj & (is_red | is_amp)));
if any(tc_excluded_subj & (is_red | is_amp))
    excl_ids = uIDs(tc_excluded_subj & (is_red | is_amp));
    fprintf('Excluded IDs: %s\n', sprintf('%d ', excl_ids));
end

cbpt_report_file = fullfile(stats_dir, sprintf('AOC_splitERSERD_GazeDev_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_split_AlphaAmpRed_GazeDev.m', ...
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
    'metric', 'Gaze deviation change [% baseline]'));

plot_timecourse_with_effect_CBPT(gaze_tc, is_red_tc, is_amp_tc, colors, ...
    gaze_ylabel, sprintf('%s_%s_%s', task_tag, split_label, gaze_tag_base), fig_dir_root, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc, false, 'ERSD [dB]', tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
plot_timecourse_individuals(gaze_tc, is_red_tc, is_amp_tc, colors, ...
    gaze_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_%s_%s_individuals_collapsed', task_tag, split_label, gaze_tag_base), ...
    fig_dir_task, fig_pos, fontSizeTC, fs, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high);

% Additional condition-wise time course outputs
for c = 1:numel(cond_vals)
    tc_cond = gaze_tc(:, c, :);
    eeg_tc_cond = eeg_tc(:, c, :);
    save_tag_cond = sprintf('%s_%s_%s_%s', task_tag, split_label, gaze_tag_base, sanitize_label_for_fname(cond_labels{c}));
    plot_timecourse_with_effect_CBPT(tc_cond, is_red_tc, is_amp_tc, colors, ...
        gaze_ylabel, save_tag_cond, fig_dir_task, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc_cond, false, 'ERSD [dB]', tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
    plot_timecourse_individuals(tc_cond, is_red_tc, is_amp_tc, colors, ...
        gaze_ylabel, cond_labels{c}, ...
        sprintf('%s_individuals', save_tag_cond), ...
        fig_dir_task, fig_pos, fontSizeTC, fs, tc_viz_smooth_sec, tk.group_lbl_low, tk.group_lbl_high);
end

%% Sanity checks
fprintf('\n=== Data Diagnostics [%s] ===\n', task_tag);
fprintf('Missing ERSD timecourse (%s): %d\n', tk.ersd_tc_fname, numel(unique(missing_ersd)));
fprintf('Missing TFR files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Main figures saved to: %s\n', fig_dir_root);
fprintf('Supplementary figures saved to: %s\n', fig_dir_task);

%% Export CSV for Python statistics script
fprintf('\n=== Exporting CSV for Python stats ===\n');
t_win_lo = 0;
t_win_hi = 2;
task_idx_gaze = (t_vec >= t_win_lo) & (t_vec <= t_win_hi);
dev_summary_by_load = squeeze(mean(gaze_tc(:, :, task_idx_gaze), 3, 'omitnan'));

n_rows = nSubj * numel(cond_vals);
ID_col = nan(n_rows, 1);
LoadValue_col = nan(n_rows, 1);
LoadLabel_col = strings(n_rows, 1);
Group_col = strings(n_rows, 1);
Included_col = false(n_rows, 1);
ERSD_col = nan(n_rows, 1);
GazeDev_col = nan(n_rows, 1);
GazeSummary_col = nan(n_rows, 1);

r = 0;
for s = 1:nSubj
    for c = 1:numel(cond_vals)
        r = r + 1;
        ID_col(r) = uIDs(s);
        LoadValue_col(r) = cond_vals(c);
        LoadLabel_col(r) = string(cond_labels{c});
        ERSD_col(r) = metrics.ERSD(s, c);
        GazeDev_col(r) = metrics.Dev(s, c);
        GazeSummary_col(r) = dev_summary_by_load(s, c);

        if is_red(s)
            Group_col(r) = string(tk.group_lbl_low);
            Included_col(r) = is_red_tc(s);
        elseif is_amp(s)
            Group_col(r) = string(tk.group_lbl_high);
            Included_col(r) = is_amp_tc(s);
        else
            Group_col(r) = "Excluded";
            Included_col(r) = false;
        end
    end
end

stats_tbl = table( ...
    ID_col, LoadValue_col, LoadLabel_col, Group_col, Included_col, ...
    ERSD_col, GazeDev_col, GazeSummary_col, ...
    'VariableNames', { ...
    'ID', 'LoadValue', 'LoadLabel', 'Group', 'Included', ...
    tk.ersd_var, 'GazeDeviationFullBL', gaze_csv_varname});

csv_out = fullfile(stats_dir, sprintf('AOC_splitERSERD_GazeDev_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV saved to: %s\n', csv_out);
fprintf('CBPT report saved to: %s\n', cbpt_report_file);

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
ttls = {'ERD', 'ERS'};
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
        set(eb(c).edge(1), 'Color', 'none');
        set(eb(c).edge(2), 'Color', 'none');
    end
    xlim([5 30]);
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ylim(ylim_shared);
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    title(ttl, 'FontSize', fsz + 4);
    legend([eb.mainLine], cond_labels, 'FontSize', fsz - 2, 'Location', 'best', 'Box', 'off');
    set(gca, 'FontSize', fsz);
    box off
end
pause(0.05); drawnow;
saveas(gcf, out_file);
close(gcf);
end

function plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg, group_lbl_low, group_lbl_high)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping TFR plot (incomplete condition data).');
    return
end

if nargin < 11 || isempty(winsor_cfg)
    winsor_cfg = struct('enable', false, 'prctile', [2 98]);
end
if nargin < 13 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 14 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end

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
time_idx = ga_red{1}.time >= -0.5 & ga_red{1}.time <= 2;
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
cfg.xlim = [-.5 2];
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
    title(ax, sprintf('%s - %s', group_lbl_low, cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');

    % Column 2: high group (rows 1–3)
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
    title(ax, sprintf('%s - %s', group_lbl_high, cond_labels{c}), 'FontSize', fsz_tfr, 'Interpreter', 'none');
end
pause(0.05); drawnow;
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

function plot_group_topos_ersd(tfr_red, tfr_amp, channels, headmodel, cond_labels, topo_latency, fig_dir, fig_prefix, fig_pos, fsz, group_lbl_low, group_lbl_high)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping topoplots (missing baselined TFR data).');
    return
end
if nargin < 11 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 12 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
cfg_ga = [];
cfg_ga.keepindividual = 'yes';
for c = 1:3
    ga_red{c} = ft_freqgrandaverage(cfg_ga, tfr_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage(cfg_ga, tfr_amp{c}{:});
    cfg_sel = [];
    cfg_sel.frequency = [8 14];
    cfg_sel.avgoverfreq = 'yes';
    ga_red{c} = ft_selectdata(cfg_sel, ga_red{c});
    ga_amp{c} = ft_selectdata(cfg_sel, ga_amp{c});
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
cfg.xlim = topo_latency;
cfg.colormap = rdbu_cmap(64);
cfg.zlim = 'maxabs';

figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    ax = subplot(2, 3, c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_red{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('%s - %s', group_lbl_low, cond_labels{c}), 'Interpreter', 'none');

    ax = subplot(2, 3, 3 + c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_amp{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('%s - %s', group_lbl_high, cond_labels{c}), 'Interpreter', 'none');
end
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_all.png', fig_prefix)));
close(gcf);
end

function plot_group_topos_ersd_collapsed(tfr_red, tfr_amp, channels, headmodel, topo_latency, fig_dir, fig_prefix, fig_pos, fsz, group_lbl_low, group_lbl_high)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping collapsed topoplots (missing baselined TFR data).');
    return
end
if nargin < 10 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 11 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end

tfr_red_all = [tfr_red{1}, tfr_red{2}, tfr_red{3}];
tfr_amp_all = [tfr_amp{1}, tfr_amp{2}, tfr_amp{3}];
if isempty(tfr_red_all) || isempty(tfr_amp_all)
    warning('Skipping collapsed topoplots (no pooled TFR data).');
    return
end

cfg_ga = [];
cfg_ga.keepindividual = 'yes';
ga_red = ft_freqgrandaverage(cfg_ga, tfr_red_all{:});
ga_amp = ft_freqgrandaverage(cfg_ga, tfr_amp_all{:});
cfg_sel = [];
cfg_sel.frequency = [8 14];
cfg_sel.avgoverfreq = 'yes';
ga_red = ft_selectdata(cfg_sel, ga_red);
ga_amp = ft_selectdata(cfg_sel, ga_amp);

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
cfg.xlim = topo_latency;
cfg.colormap = rdbu_cmap(64);
cfg.zlim = 'maxabs';

figure('Position', fig_pos, 'Color', 'w');
ax = subplot(1, 2, 1);
cfg.figure = ax;
ft_topoplotER(cfg, ga_red);
colorbar(ax);
set(ax, 'FontSize', fsz);
title(ax, sprintf('%s - collapsed over conditions', group_lbl_low), 'Interpreter', 'none');

ax = subplot(1, 2, 2);
cfg.figure = ax;
ft_topoplotER(cfg, ga_amp);
colorbar(ax);
set(ax, 'FontSize', fsz);
title(ax, sprintf('%s - collapsed over conditions', group_lbl_high), 'Interpreter', 'none');

pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_collapsedConditions.png', fig_prefix)));
close(gcf);
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
    title(ax, sprintf('ERD - %s', cond_labels{c}), 'Interpreter', 'none');

    % Row 2: ERS (cols 1–3)
    ax = subplot(2, 3, 3 + c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_amp{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('ERS - %s', cond_labels{c}), 'Interpreter', 'none');
end
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_all.png', fig_prefix)));
close(gcf);
end

function plot_group_tfrs_collapsed(tfr_red, tfr_amp, channels, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg, group_lbl_low, group_lbl_high)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping collapsed TFR plot (incomplete condition data).');
    return
end

if nargin < 9 || isempty(winsor_cfg)
    winsor_cfg = struct('enable', false, 'prctile', [2 98]);
end
if nargin < 11 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 12 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end

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
time_idx = ga_red.time >= -0.5 & ga_red.time <= 2;
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
cfg.xlim = [-0.5 2];
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
rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.FontSize = fsz_tfr - 4;
cbar.Label.String = 'Power [dB]';
cbar.Label.FontSize = fsz_tfr;
title(ax, group_lbl_low, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');

ax = subplot(2, 1, 2);
cfg.figure = ax;
ft_singleplotTFR(cfg, ga_amp);
colormap(ax, color_map);
set(ax, 'CLim', clim_abs);
cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr);
ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr);
rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.FontSize = fsz_tfr - 4;
cbar.Label.String = 'Power [dB]';
cbar.Label.FontSize = fsz_tfr;
title(ax, group_lbl_high, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');

pause(0.05); drawnow;
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
title(ax, 'ERD - collapsed over conditions', 'Interpreter', 'none');

ax = subplot(1, 2, 2);
cfg.figure = ax;
ft_topoplotER(cfg, ga_amp);
colorbar(ax);
set(ax, 'FontSize', fsz);
title(ax, 'ERS - collapsed over conditions', 'Interpreter', 'none');

pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_prefix, fig_pos, fsz, group_lbl_low, group_lbl_high)
if nargin < 10 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 11 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
metric_defs = { ...
    'ERSD', 'ERSD', true; ...
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
    title(group_lbl_low, 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    box off

    nexttile; hold on
    for c = 1:3
        draw_one_cloud(X(is_amp, c), c, colors(c,:), 0.3, 96, 0.45);
    end
    yline(0, '--', 'Color', [0.4 0.4 0.4]);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title(group_lbl_high, 'Interpreter', 'none');
    ylabel(ylab, 'Interpreter', 'none');
    box off

    sgtitle(varname, 'FontSize', fsz+2, 'Interpreter', 'none');
    pause(0.05); drawnow;
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

function plot_timecourse_load_compare(tc_a, tc_b, col_a, col_b, label_a, label_b, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, smooth_sec)
if nargin < 14 || isempty(smooth_sec)
    smooth_sec = 0.05;
end
if nargin < 13 || isempty(ds_factor)
    ds_factor = 10;
end
tc_a(~isfinite(tc_a)) = NaN;
tc_b(~isfinite(tc_b)) = NaN;
nT = size(tc_a, 2);
dt = 1 / fs;
t_plot = linspace(-0.5 + dt, 2, nT);
win_sm = max(1, round(smooth_sec * fs));
if win_sm > 1
    tc_a = movmean(tc_a, win_sm, 2, 'omitnan');
    tc_b = movmean(tc_b, win_sm, 2, 'omitnan');
end
mA = mean(tc_a, 1, 'omitnan');
mB = mean(tc_b, 1, 'omitnan');
nA_fin = sum(isfinite(tc_a), 1);
nB_fin = sum(isfinite(tc_b), 1);
sA = std(tc_a, 0, 1, 'omitnan') ./ max(sqrt(nA_fin), 1);
sB = std(tc_b, 0, 1, 'omitnan') ./ max(sqrt(nB_fin), 1);
sA(~isfinite(sA)) = NaN;
sB(~isfinite(sB)) = NaN;

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile([2 1]); hold on
e1 = shadedErrorBar(t_plot, mA, sA, 'lineProps', {'-'}, 'transparent', true);
e2 = shadedErrorBar(t_plot, mB, sB, 'lineProps', {'-'}, 'transparent', true);
set(e1.mainLine, 'Color', col_a, 'LineWidth', 2.5);
set(e2.mainLine, 'Color', col_b, 'LineWidth', 2.5);
set(e1.patch, 'FaceColor', col_a, 'FaceAlpha', 0.20);
set(e2.patch, 'FaceColor', col_b, 'FaceAlpha', 0.20);
set(e1.edge(1), 'Color', 'none');
set(e1.edge(2), 'Color', 'none');
set(e2.edge(1), 'Color', 'none');
set(e2.edge(2), 'Color', 'none');
xline(0, '--k');
xlim([-0.5 2]);
ylabel(ylab);
xlabel('Time [s]');
title(sprintf('%s vs %s (pooled over subjects)', label_a, label_b), 'Interpreter', 'none');
legend([e1.mainLine, e2.mainLine], ...
    {sprintf('%s (n=%d)', label_a, size(tc_a, 1)), sprintf('%s (n=%d)', label_b, size(tc_b, 1))}, ...
    'Location', 'northeast', 'FontSize', fsz-2, 'Box', 'off');
set(gca, 'FontSize', fsz-4);
box off

nexttile; hold on
n_perm = 10000;
min_per_group = 3;
d = nan(1, nT);
for t = 1:nT
    x = tc_a(:, t); y = tc_b(:, t);
    x = x(isfinite(x)); y = y(isfinite(y));
    if numel(x) < min_per_group || numel(y) < min_per_group, continue, end
    sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / max(numel(x)+numel(y)-2, 1));
    d(t) = (mean(y) - mean(x)) / max(sp, eps);
end
tc_a_ds = tc_a(:, 1:ds_factor:end);
tc_b_ds = tc_b(:, 1:ds_factor:end);
t_plot_ds = t_plot(1:ds_factor:end);
dt_ds = ds_factor * dt;
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d(tc_a_ds, tc_b_ds, n_perm, 0.05, 'twotail', t_plot_ds);
sig_cluster = false(1, numel(t_plot_ds));
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        sig_cluster(clusters(k).idx) = true;
    end
end
sig_uncorr = (abs(tvals_cl) > thr.tcrit) & isfinite(tvals_cl);
sig = sig_cluster;
if ~any(sig) && any(sig_uncorr)
    sig = sig_uncorr;
    fprintf('  [%s] (cluster n.s.; shading uncorrected |t|>tcrit for 3-back vs 1-back)\n', save_tag);
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
box off
set(gca, 'FontSize', fsz-4);
if ~any(sig_cluster) && any(sig_uncorr)
    title({'WARNING: No significant clusters; shading shows uncorrected |t| > t_{crit} (3-back vs 1-back)'}, ...
        'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'Interpreter', 'tex');
end
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_GazeDev_timecourse_%s.png', save_tag)));
close(gcf);
end

function plot_timecourse_individuals(tc, is_red, is_amp, colors, ylab, title_tag, save_tag, fig_dir, fig_pos, fsz, fs, smooth_sec, group_lbl_low, group_lbl_high)
if nargin < 13 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 14 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 2, nT);
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
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_GazeDev_timecourse_%s.png', save_tag)));
close(gcf);
end

% Cluster-based permutation test
function plot_timecourse_with_effect_CBPT(tc, is_red, is_amp, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, eeg_tc, addEEG_TC, eeg_ylab, smooth_sec, group_lbl_low, group_lbl_high, cbpt_report_file)
if nargin < 11 || isempty(ds_factor), ds_factor = 50; end
if nargin < 12, eeg_tc = []; end
if nargin < 13, addEEG_TC = ~isempty(eeg_tc); end
if nargin < 14, eeg_ylab = 'ERSD [dB]'; end
if nargin < 15 || isempty(smooth_sec), smooth_sec = 0.05; end
if nargin < 16 || isempty(group_lbl_low), group_lbl_low = 'ERD'; end
if nargin < 17 || isempty(group_lbl_high), group_lbl_high = 'ERS'; end
if nargin < 18, cbpt_report_file = ''; end
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 2, nT);
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
% Legend with colored patch boxes (clearer than thin lines)
leg_p1 = patch(NaN, NaN, colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(1,:), 'LineWidth', 1.5);
leg_p2 = patch(NaN, NaN, colors(3,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(3,:), 'LineWidth', 1.5);
legend([leg_p1 leg_p2], {[' ' group_lbl_low], [' ' group_lbl_high]}, 'Location', 'best', 'FontSize', fsz*0.75, 'Box', 'off');
if contains(save_tag, 'GazeDev_pct')
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
dt_ds = ds_factor * dt;
alpha_cbpt = 0.05;
tail_cbpt = 'twotail';
[clusters, tvals_cl, thr, maxMassNull, maxExtentNull] = ft_cluster_permutation_1d(Rall_ds, Aall_ds, n_perm, alpha_cbpt, tail_cbpt, t_plot_ds);
report_cfg = struct( ...
    'tag', save_tag, ...
    'modality', 'gaze', ...
    'nR', nR, ...
    'nA', nA, ...
    'lbl_low', group_lbl_low, ...
    'lbl_high', group_lbl_high, ...
    'n_perm', n_perm, ...
    'alpha', alpha_cbpt, ...
    'tail', tail_cbpt, ...
    'nT_ds', nT_ds, ...
    'bin_ms', ds_factor * 1000 / fs, ...
    'fs', fs, ...
    'ds_factor', ds_factor, ...
    'clusters', clusters, ...
    'tvals', tvals_cl, ...
    'thr', thr, ...
    't_plot', t_plot_ds, ...
    'dt_ds', dt_ds, ...
    'maxMassNull', maxMassNull, ...
    'maxExtentNull', maxExtentNull);
log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(report_cfg));
sig_cluster = false(1, nT_ds);
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        sig_cluster(clusters(k).idx) = true;
    end
end
sig_uncorr = (abs(tvals_cl) > thr.tcrit) & isfinite(tvals_cl);
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

% ERSD time course panel in combined figure.
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
    set(ebR.mainLine, 'Color', colors(1,:), 'LineWidth', 2.5);
    set(ebA.mainLine, 'Color', colors(3,:), 'LineWidth', 2.5);
    set(ebR.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.20);
    set(ebA.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.20);
    set(ebR.edge(1), 'Color', 'none');
    set(ebR.edge(2), 'Color', 'none');
    set(ebA.edge(1), 'Color', 'none');
    set(ebA.edge(2), 'Color', 'none');
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
    xlim([-0.5 2]);
    box off
    set(gca, 'FontSize', fsz-4);
    leg_p1 = patch(NaN, NaN, colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(1,:), 'LineWidth', 1.5);
    leg_p2 = patch(NaN, NaN, colors(3,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(3,:), 'LineWidth', 1.5);

    % EEG effect-size strip (same logic as gaze panel above)
    ax_d_eeg = nexttile; hold on
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
    [clusters_eeg, tvals_cl_eeg, thr_eeg] = ft_cluster_permutation_1d(EegR_ds, EegA_ds, n_perm, alpha_cbpt, tail_cbpt, t_plot_ds);
    report_cfg_eeg = struct( ...
        'tag', sprintf('%s_EEG', save_tag), ...
        'modality', 'EEG', ...
        'nR', size(EegR_ds, 1), ...
        'nA', size(EegA_ds, 1), ...
        'lbl_low', group_lbl_low, ...
        'lbl_high', group_lbl_high, ...
        'n_perm', n_perm, ...
        'alpha', alpha_cbpt, ...
        'tail', tail_cbpt, ...
        'nT_ds', nT_ds, ...
        'bin_ms', ds_factor * 1000 / fs, ...
        'fs', fs, ...
        'ds_factor', ds_factor, ...
        'clusters', clusters_eeg, ...
        'tvals', tvals_cl_eeg, ...
        'thr', thr_eeg, ...
        't_plot', t_plot_ds, ...
        'dt_ds', dt_ds, ...
        'maxMassNull', [], ...
        'maxExtentNull', []);
    log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(report_cfg_eeg));
    sig_cluster_eeg = false(1, nT_ds);
    for k = 1:numel(clusters_eeg)
        if clusters_eeg(k).p < 0.05
            sig_cluster_eeg(clusters_eeg(k).idx) = true;
        end
    end
    sig_uncorr_eeg = (abs(tvals_cl_eeg) > thr_eeg.tcrit) & isfinite(tvals_cl_eeg);
    sig_eeg = sig_cluster_eeg;
    if ~any(sig_eeg) && any(sig_uncorr_eeg)
        sig_eeg = sig_uncorr_eeg;
        log_cbpt_report(cbpt_report_file, {sprintf('  [%s EEG] (cluster n.s.; shading uncorrected |t|>tcrit)', save_tag)});
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
    xlim([-0.5 2]);
    box off
    set(gca, 'FontSize', fsz-4);
    if ~any(sig_cluster_eeg) && any(sig_uncorr_eeg)
        title({'WARNING: No significant clusters; shading shows uncorrected |t| > t_{crit}'}, ...
            'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'Interpreter', 'tex');
    end
end

drawnow;
align_stacked_tc_panels(ax_d, ax_gaze);
if show_eeg
    align_stacked_tc_panels(ax_d_eeg, ax_eeg);
end
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_GazeDev_timecourse_%s_CBPT.png', save_tag)));
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
