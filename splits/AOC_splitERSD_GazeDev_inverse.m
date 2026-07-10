%% AOC Split ERS/ERD (Trial-Level, Inverse) — GazeDev (Sternberg + N-back)
% Within-subject median split on gaze deviation; outcome = occipital ERSD.
% Run AOC_splitBehavior_inverse_Prep.m first.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir_root = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/splits/SplitERSERD/inverse';
stats_dir = paths.splits_stats;
if ~isfolder(fig_dir_root), mkdir(fig_dir_root); end
if ~isfolder(stats_dir), mkdir(stats_dir); end

split_file = fullfile(stats_dir, 'AOC_splitGazeDev_inverse_trialSplits.mat');
if ~isfile(split_file)
    error('Missing inverse split file. Run AOC_splitBehavior_inverse_Prep.m first:\n  %s', split_file);
end
Ssplit = load(split_file, 'split_gaze_dev_inverse');
split_data = Ssplit.split_gaze_dev_inverse;

fprintf('\n=== AOC Split ERS/ERD Inverse — Gaze Deviation (trial-level) ===\n');
fprintf('Split file: %s\n', split_file);

fig_pos = [0 0 1512 982];
fontSize = 40;
tfr_winsor_cfg = struct('enable', true, 'prctile', [2 98]);

tasks(1).tag = 'sternberg';
tasks(1).split_label = 'splitMedian_trial';
tasks(1).topo_latency = [1 2];
tasks(1).group_lbl_low = 'Low Gaze Dev';
tasks(1).group_lbl_high = 'High Gaze Dev';
tasks(1).split_var = 'GazeDeviationFullBL';
tasks(1).ersd_var = 'ERSD_late';

tasks(2).tag = 'nback';
tasks(2).split_label = 'splitMedian_trial';
tasks(2).topo_latency = [0 2];
tasks(2).group_lbl_low = 'Low Gaze Dev';
tasks(2).group_lbl_high = 'High Gaze Dev';
tasks(2).split_var = 'GazeDeviationFullBL';
tasks(2).ersd_var = 'ERSD_full';

for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
split_label = tk.split_label;
fig_prefix = sprintf('AOC_splitERSD_%s_GazeDev_inverse', task_tag);

fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));
if ~isfield(split_data, task_tag)
    warning('Skipping %s: missing split field.', task_tag);
    continue
end
task_split = split_data.(task_tag);
subj_splits = task_split.subjects;
nSubj = numel(subj_splits);
if nSubj < 2
    warning('Skipping %s: fewer than 2 subjects.', task_tag);
    continue
end

uIDs = [subj_splits.ID]';
split_info_str = sprintf('Within-subject median on %s (conditions pooled)', tk.split_var);
fprintf('%s\n', split_info_str);

split_mean_low = nan(nSubj, 1);
split_mean_high = nan(nSubj, 1);
split_thr = nan(nSubj, 1);
ersd_mean_low = nan(nSubj, 1);
ersd_mean_high = nan(nSubj, 1);
for s = 1:nSubj
    split_mean_low(s) = mean(subj_splits(s).SplitMetric(subj_splits(s).idx_low), 'omitnan');
    split_mean_high(s) = mean(subj_splits(s).SplitMetric(subj_splits(s).idx_high), 'omitnan');
    split_thr(s) = subj_splits(s).median_thr;
    ersd_mean_low(s) = mean(subj_splits(s).ERSD(subj_splits(s).idx_low), 'omitnan');
    ersd_mean_high(s) = mean(subj_splits(s).ERSD(subj_splits(s).idx_high), 'omitnan');
end

%% Inclusion (split variable)
figure('Position', fig_pos, 'Color', 'w'); hold on
x = (1:nSubj)';
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
plot(x, split_thr, 'k-', 'LineWidth', 1.2);
h_low = scatter(x, split_mean_low, 70, colors(1, :), 'filled', 'MarkerFaceAlpha', 0.85);
h_high = scatter(x, split_mean_high, 70, colors(3, :), 'filled', 'MarkerFaceAlpha', 0.85);
xlabel('Participant (index)');
ylabel(sprintf('Mean trial %s [%% BL]', tk.split_var), 'Interpreter', 'none');
title(sprintf('Trial-level median split [%s]', task_tag), 'Interpreter', 'none');
legend([h_low, h_high], {tk.group_lbl_low, tk.group_lbl_high}, 'Location', 'best', 'FontSize', fontSize-4, 'Box', 'off');
set(gca, 'FontSize', fontSize); box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir_root, sprintf('%s_inclusion.png', fig_prefix)));
close(gcf);

%% Collect TFRs and ERSD time courses
tfr_low_all = {};
tfr_high_all = {};
metrics_split = nan(nSubj, 2);
metrics_ERSD = nan(nSubj, 2);
ersd_cfg = init_ersd_tc_cfg(subj_splits(1).time);
t_vec = ersd_cfg.t_vec;
Tf = ersd_cfg.n_samp;
idx_viable = (t_vec >= 0) & (t_vec <= 2);
ersd_tc = nan(nSubj, 2, Tf);

fprintf('\n=== Aggregating ERSD within gaze-split groups ===\n');
for s = 1:nSubj
    clc; fprintf('[SPLIT INVERSE GAZEDEV - %s] Subject %d / %d\n', upper(task_tag), s, nSubj);
    metrics_split(s, 1) = split_mean_low(s);
    metrics_split(s, 2) = split_mean_high(s);
    metrics_ERSD(s, 1) = ersd_mean_low(s);
    metrics_ERSD(s, 2) = ersd_mean_high(s);

    if isfield(subj_splits(s), 'tfr_low') && ~isempty(subj_splits(s).tfr_low)
        tfr_low_all{end+1} = subj_splits(s).tfr_low; %#ok<AGROW>
        tfr_high_all{end+1} = subj_splits(s).tfr_high; %#ok<AGROW>
    end

  t_subj = subj_splits(s).time(:)';
  tcL = subj_splits(s).ersd_tc_low(:)';
  tcH = subj_splits(s).ersd_tc_high(:)';
  if numel(t_subj) == Tf
      ersd_tc(s, 1, :) = tcL;
      ersd_tc(s, 2, :) = tcH;
  else
      ersd_tc(s, 1, :) = interp1(t_subj, tcL, t_vec, 'linear', NaN);
      ersd_tc(s, 2, :) = interp1(t_subj, tcH, t_vec, 'linear', NaN);
  end
end

channels = {};
if ~isempty(tfr_low_all)
    channels = occ_channels_from_labels(tfr_low_all{1}.label);
elseif ~isempty(tfr_high_all)
    channels = occ_channels_from_labels(tfr_high_all{1}.label);
end

color_map_tfr = customcolormap_preset('red-white-blue');
if ~isempty(tfr_low_all) && ~isempty(tfr_high_all)
    fprintf('\n=== Plotting collapsed TFRs ===\n');
    plot_group_tfrs_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, color_map_tfr, ...
        fig_dir_root, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg, tk.group_lbl_low, tk.group_lbl_high);
    fprintf('\n=== Plotting collapsed topoplots ===\n');
    plot_group_topos_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, tk.topo_latency, ...
        fig_dir_root, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);
end

%% Rainclouds
fprintf('\n=== Plotting rainclouds ===\n');
plot_paired_raincloud(metrics_split(:,1), metrics_split(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'Gaze deviation [% BL]', fig_dir_root, sprintf('%s_raincloud_gazedev.png', fig_prefix), fig_pos, fontSize);
plot_paired_raincloud(metrics_ERSD(:,1), metrics_ERSD(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'ERSD [dB]', fig_dir_root, sprintf('%s_raincloud_ersd.png', fig_prefix), fig_pos, fontSize);

%% ERSD time courses
close all
fontSizeTC = 40;
rng(123)
tc_window_idx = idx_viable;
tc_complete_min_frac = 0.80;
keep_tc = true(nSubj, 1);
for g = 1:2
    Xg = reshape(ersd_tc(:, g, :), nSubj, Tf);
    [Xg, keep_g] = preprocess_ersd_subject_tc(Xg, idx_viable, ersd_cfg);
    ersd_tc(:, g, :) = reshape(Xg, [nSubj, 1, Tf]);
    frac = mean(isfinite(Xg(:, tc_window_idx)), 2);
    keep_tc = keep_tc & keep_g & (frac >= tc_complete_min_frac);
end
ersd_tc(~keep_tc, :, :) = NaN;
n_tc = sum(keep_tc);
fprintf('Included subjects for ERSD TC: %d / %d\n', n_tc, nSubj);

low_group_timecourses = reshape(ersd_tc(keep_tc, 1, :), n_tc, Tf);
high_group_timecourses = reshape(ersd_tc(keep_tc, 2, :), n_tc, Tf);

ersd_ylabel = sprintf('ERSD\nPower Change [dB]');
cbpt_report_file = fullfile(stats_dir, sprintf('AOC_splitERSD_GazeDev_inverse_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_splitERSD_GazeDev_inverse.m', ...
    'task_tag', task_tag, ...
    'split_label', split_label, ...
    'split_info', split_info_str, ...
    'split_var', tk.split_var, ...
    'ersd_var', tk.ersd_var, ...
    'group_lbl_low', tk.group_lbl_low, ...
    'group_lbl_high', tk.group_lbl_high, ...
    'n_low_split', nSubj, ...
    'n_high_split', nSubj, ...
    'n_low_tc', n_tc, ...
    'n_high_tc', n_tc, ...
    'metric', 'Occipital alpha ERSD [dB] (paired within subject)'));

tc_base = sprintf('AOC_splitERSD_%s_GazeDev_inverse_timecourse', task_tag);
report_tag = sprintf('%s_%s_GazeDev_inverse', task_tag, split_label);
plot_paired_timecourse_CBPT(low_group_timecourses, high_group_timecourses, colors, ersd_ylabel, ...
    tc_base, report_tag, ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
plot_difference_timecourse_CBPT(low_group_timecourses, high_group_timecourses, ...
    sprintf('%s_HighMinusLow', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, 'ERSD [dB]');
plot_paired_timecourse_individuals(low_group_timecourses, high_group_timecourses, colors, ersd_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_individuals', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high);

%% CSV
fprintf('\n=== Exporting CSV ===\n');
t_win = (t_vec >= 0) & (t_vec <= 2);
ersd_sum_low = mean(low_group_timecourses(:, t_win), 2, 'omitnan');
ersd_sum_high = mean(high_group_timecourses(:, t_win), 2, 'omitnan');
ids_tc = uIDs(keep_tc);
ID_col = []; Group_col = strings(0,1); Included_col = [];
Split_col = []; ERSD_col = []; ERSDSummary_col = [];
for s = 1:nSubj
    for g = 1:2
        ID_col(end+1,1) = uIDs(s);
        if g == 1
            Group_col(end+1,1) = string(tk.group_lbl_low);
            Split_col(end+1,1) = metrics_split(s,1);
            ERSD_col(end+1,1) = metrics_ERSD(s,1);
        else
            Group_col(end+1,1) = string(tk.group_lbl_high);
            Split_col(end+1,1) = metrics_split(s,2);
            ERSD_col(end+1,1) = metrics_ERSD(s,2);
        end
        incl = keep_tc(s);
        Included_col(end+1,1) = incl;
        if incl
            ix = find(ids_tc == uIDs(s), 1);
            if g == 1, ERSDSummary_col(end+1,1) = ersd_sum_low(ix);
            else, ERSDSummary_col(end+1,1) = ersd_sum_high(ix);
            end
        else
            ERSDSummary_col(end+1,1) = NaN;
        end
    end
end
stats_tbl = table(ID_col, Group_col, Included_col, Split_col, ERSD_col, ERSDSummary_col, ...
    'VariableNames', {'ID','Group','Included', tk.split_var, tk.ersd_var, 'ERSD_tc_0_2s'});
csv_out = fullfile(stats_dir, sprintf('AOC_splitERSD_GazeDev_inverse_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV: %s\n', csv_out);
fprintf('CBPT report: %s\n', cbpt_report_file);
end

function ersd_cfg = init_ersd_tc_cfg(t_vec)
ersd_cfg.t_vec = t_vec(:)';
ersd_cfg.n_samp = numel(ersd_cfg.t_vec);
dt = median(diff(ersd_cfg.t_vec), 'omitnan');
if ~isfinite(dt) || dt <= 0, dt = 0.05; end
ersd_cfg.dt = dt;
ersd_cfg.outlier_k_iqr = 1.5;
ersd_cfg.max_interp_gap_sec = 0.15;
ersd_cfg.min_subject_coverage = 0.60;
ersd_cfg.smooth_sec = 0.10;
ersd_cfg.win_sm = max(1, round(ersd_cfg.smooth_sec / dt));
end

function [X, keep_subj] = preprocess_ersd_subject_tc(X, idx_viable, ersd_cfg)
X(~isfinite(X)) = NaN;
med_metric = median(X(:, idx_viable), 2, 'omitnan');
[X, keep_subj] = exclude_outlier_trajectories(X, med_metric, ersd_cfg.outlier_k_iqr);
max_interp_gap_smp = max(1, round(ersd_cfg.max_interp_gap_sec / ersd_cfg.dt));
for s = 1:size(X, 1)
    X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
end
subj_cov = mean(isfinite(X(:, idx_viable)), 2);
keep_subj = keep_subj & (subj_cov >= ersd_cfg.min_subject_coverage);
X(~keep_subj, :) = NaN;
if ersd_cfg.win_sm > 1
    X = movmean(X, ersd_cfg.win_sm, 2, 'omitnan');
end
end

function [X_keep, keep_subj] = exclude_outlier_trajectories(X, metric, k_iqr)
med_m = median(metric, 'omitnan');
iqr_m = iqr(metric);
if ~isfinite(iqr_m) || iqr_m == 0
    low_m = -inf; high_m = inf;
else
    low_m = med_m - k_iqr * iqr_m;
    high_m = med_m + k_iqr * iqr_m;
end
keep_subj = isfinite(metric) & (metric >= low_m) & (metric <= high_m);
X_keep = X;
X_keep(~keep_subj, :) = NaN;
end

function x = fill_short_nan_gaps(x, max_gap_smp)
if isempty(x), return, end
valid = isfinite(x);
if all(~valid) || all(valid), return, end
n = numel(x); i = 1;
while i <= n
    if valid(i), i = i + 1; continue, end
    j = i;
    while j <= n && ~valid(j), j = j + 1; end
    gap_start = i; gap_end = j - 1; gap_len = gap_end - gap_start + 1;
    left = gap_start - 1; right = gap_end + 1;
    if left >= 1 && right <= n && valid(left) && valid(right) && gap_len <= max_gap_smp
        x(gap_start:gap_end) = interp1([left right], [x(left) x(right)], gap_start:gap_end);
    end
    i = j;
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
if isempty(ch), ch = labels; end
end

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm), return, end
P = freq_in.powspctrm;
if size(P, 1) < 3, return, end
P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function A = mean_over_channels_tfr(T, ch_idx)
P = T.powspctrm;
if isempty(P)
    A = nan(numel(T.freq), numel(T.time)); return
end
Psize = size(P);
chan_dim = find(Psize == numel(T.label), 1, 'first');
if isempty(chan_dim), chan_dim = 2; end
switch ndims(P)
    case 4
        if chan_dim == 2
            A = squeeze(mean(P(:, ch_idx, :, :), 2, 'omitnan'));
        else
            A = squeeze(mean(P(ch_idx, :, :, :), 1, 'omitnan'));
        end
        if ndims(A) == 3, A = squeeze(mean(A, 1, 'omitnan')); end
    case 3
        if chan_dim == 1
            A = squeeze(mean(P(ch_idx, :, :), 1, 'omitnan'));
        else
            A = squeeze(mean(P(:, ch_idx, :), 2, 'omitnan'));
            if size(A, 1) == numel(T.time), A = A'; end
        end
    otherwise
        error('Unsupported TFR dimensionality.');
end
end

function plot_group_tfrs_collapsed_paired(tfr_low, tfr_high, channels, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg, lblLow, lblHigh)
cfg = []; cfg.keepindividual = 'yes';
ga_low = ft_freqgrandaverage(cfg, tfr_low{:});
ga_high = ft_freqgrandaverage(cfg, tfr_high{:});
if winsor_cfg.enable
    ga_low = winsorize_freq_subjects(ga_low, winsor_cfg.prctile);
    ga_high = winsorize_freq_subjects(ga_high, winsor_cfg.prctile);
end
[~, ch_idx] = ismember(channels, ga_low.label);
ch_idx = ch_idx(ch_idx > 0);
freq_idx = ga_low.freq >= 5 & ga_low.freq <= 30;
time_idx = ga_low.time >= -0.5 & ga_low.time <= 2;
Ared = mean_over_channels_tfr(ga_low, ch_idx);
Aamp = mean_over_channels_tfr(ga_high, ch_idx);
mx = max([max(abs(Ared(freq_idx, time_idx)), [], 'all'), max(abs(Aamp(freq_idx, time_idx)), [], 'all')]);
if ~isfinite(mx) || mx <= 0, mx = 0.1; end
clim_abs = [-0.9*mx, 0.9*mx];
cfg = [];
cfg.channel = channels; cfg.colorbar = 'no'; cfg.zlim = 'maxabs';
cfg.xlim = [-0.5 2]; cfg.ylim = [5 30]; cfg.layout = headmodel.layANThead;
fsz_tfr = round(fsz * 0.8);
figure('Position', [0 0 1512*0.666 982], 'Color', 'w');
ax = subplot(2,1,1); cfg.figure = ax; ft_singleplotTFR(cfg, ga_low);
colormap(ax, color_map); set(ax, 'CLim', clim_abs); cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr); ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr); rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.Label.String = 'Power [dB]'; title(ax, lblLow, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');
ax = subplot(2,1,2); cfg.figure = ax; ft_singleplotTFR(cfg, ga_high);
colormap(ax, color_map); set(ax, 'CLim', clim_abs); cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr); ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr); rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.Label.String = 'Power [dB]'; title(ax, lblHigh, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_tfr_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function plot_group_topos_collapsed_paired(tfr_low, tfr_high, channels, headmodel, topo_latency, fig_dir, fig_prefix, fig_pos, fsz, lblLow, lblHigh)
cfg_ga = []; cfg_ga.keepindividual = 'yes';
ga_low = ft_freqgrandaverage(cfg_ga, tfr_low{:});
ga_high = ft_freqgrandaverage(cfg_ga, tfr_high{:});
cfg_sel = []; cfg_sel.frequency = [8 14]; cfg_sel.avgoverfreq = 'yes';
ga_low = ft_selectdata(cfg_sel, ga_low);
ga_high = ft_selectdata(cfg_sel, ga_high);
cfg = [];
cfg.layout = headmodel.layANThead;
allch = cfg.layout.label;
cfg.channel = allch(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.marker = 'off'; cfg.highlight = 'on'; cfg.highlightchannel = channels;
cfg.highlightsymbol = '.'; cfg.highlightsize = 10; cfg.figure = 'gcf';
cfg.gridscale = 300; cfg.comment = 'no'; cfg.xlim = topo_latency;
cfg.colormap = rdbu_cmap(64); cfg.zlim = 'maxabs';
figure('Position', fig_pos, 'Color', 'w');
ax = subplot(1,2,1); cfg.figure = ax; ft_topoplotER(cfg, ga_low); colorbar(ax);
set(ax, 'FontSize', fsz); title(ax, sprintf('%s - collapsed', lblLow), 'Interpreter', 'none');
ax = subplot(1,2,2); cfg.figure = ax; ft_topoplotER(cfg, ga_high); colorbar(ax);
set(ax, 'FontSize', fsz); title(ax, sprintf('%s - collapsed', lblHigh), 'Interpreter', 'none');
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function cmap = rdbu_cmap(n)
rdbu_11 = [33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; ...
    253 219 199; 244 165 130; 214 96 77; 178 24 43] / 255;
x = linspace(0, 1, size(rdbu_11, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, rdbu_11, xi, 'linear');
end

function plot_paired_raincloud(yLow, yHigh, colors, lblLow, lblHigh, ylab, fig_dir, out_name, fig_pos, fsz)
all_vals = [yLow(:); yHigh(:)]; all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals), ymax = 1; else, ymax = max(abs(all_vals))*1.15; if ymax<=0, ymax=1; end, end
figure('Position', fig_pos, 'Color', 'w'); hold on
draw_one_cloud(yLow, 1, colors(1,:), 0.35, 96, 0.45);
draw_one_cloud(yHigh, 2, colors(3,:), 0.35, 96, 0.45);
for i = 1:numel(yLow)
    if isfinite(yLow(i)) && isfinite(yHigh(i))
        plot([1 2], [yLow(i) yHigh(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim([-ymax ymax]);
set(gca, 'XTick', [1 2], 'XTickLabel', {lblLow, lblHigh}, 'FontSize', fsz-2);
ylabel(ylab, 'Interpreter', 'none');
title('Within-subject trial split (conditions collapsed)', 'Interpreter', 'none');
box off; pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, out_name)); close(gcf);
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha)
y = yvals(isfinite(yvals));
if numel(y) < 3, return, end
[f, xi] = ksdensity(y, 'NumPoints', 120);
if max(f) > 0, f = f / max(f) * 0.35; else, f = zeros(size(f)); end
x_den = xpos - 0.08; x_box = xpos + 0.03;
fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25); q3 = prctile(y, 75); med = median(y);
p5 = prctile(y, 5); p95 = prctile(y, 95);
plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [x_box-box_w/2, q1, box_w, max(q3-q1, eps)], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = (box_w/2) * 2 * (rand(numel(y),1)-0.5);
scatter(x_box + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

function plot_paired_timecourse_individuals(low_group_timecourses, high_group_timecourses, colors, ylab, title_tag, out_name, fig_dir, fig_pos, fsz, t_vec, lblLow, lblHigh)
t_plot = t_vec(:);
low_curves = low_group_timecourses; high_curves = high_group_timecourses;
low_curves(~isfinite(low_curves)) = NaN; high_curves(~isfinite(high_curves)) = NaN;
colR_light = colors(1,:)*0.35+0.65; colA_light = colors(3,:)*0.35+0.65;
figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; plot(t_plot, low_curves', 'Color', colR_light, 'LineWidth', 0.8);
plot(t_plot, mean(low_curves,1,'omitnan'), 'Color', colors(1,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblLow, size(low_curves,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
nexttile; hold on; plot(t_plot, high_curves', 'Color', colA_light, 'LineWidth', 0.8);
plot(t_plot, mean(high_curves,1,'omitnan'), 'Color', colors(3,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); xlabel('Time [s]'); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblHigh, size(high_curves,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function plot_paired_timecourse_CBPT(low_group_timecourses, high_group_timecourses, colors, ylab, out_name, report_tag, fig_dir, fig_pos, fsz, t_vec, lblLow, lblHigh, cbpt_report_file)
t_plot = t_vec(:)';
num_subjects_in = size(low_group_timecourses, 1);
[low_cb, high_cb, t_cb, n_pairs, dt_cb, drop_info] = prepare_paired_cbpt_poststim(low_group_timecourses, high_group_timecourses, t_plot);
cohens_d_cb = paired_cohens_dz_curve(low_cb, high_cb);
figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([2 1]); hold on
low_group_mean = mean(low_group_timecourses, 1, 'omitnan');
high_group_mean = mean(high_group_timecourses, 1, 'omitnan');
low_group_sem = std(low_group_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(low_group_timecourses), 1)), 1);
high_group_sem = std(high_group_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(high_group_timecourses), 1)), 1);
low_group_plot = shadedErrorBar(t_plot, low_group_mean, low_group_sem, 'lineProps', {'-'}, 'transparent', true);
high_group_plot = shadedErrorBar(t_plot, high_group_mean, high_group_sem, 'lineProps', {'-'}, 'transparent', true);
set(low_group_plot.mainLine, 'Color', colors(1,:), 'LineWidth', 2.5);
set(high_group_plot.mainLine, 'Color', colors(3,:), 'LineWidth', 2.5);
set(low_group_plot.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.20);
set(high_group_plot.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.20);
set(low_group_plot.edge(1), 'Color', 'none'); set(low_group_plot.edge(2), 'Color', 'none');
set(high_group_plot.edge(1), 'Color', 'none'); set(high_group_plot.edge(2), 'Color', 'none');
xline(0, '--k'); ylabel(ylab); xlim([-0.5 2]);
ylim(ylim_from_mean_sem(low_group_mean, low_group_sem, high_group_mean, high_group_sem));
box off; set(gca, 'FontSize', fsz-4);
leg_p1 = patch(NaN, NaN, colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(1,:), 'LineWidth', 1.5);
leg_p2 = patch(NaN, NaN, colors(3,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(3,:), 'LineWidth', 1.5);
legend([leg_p1 leg_p2], {[' ' lblLow], [' ' lblHigh]}, 'Location', 'best', 'FontSize', fsz*0.75, 'Box', 'off');
nexttile; hold on
n_perm = 10000; alpha_cbpt = 0.05; tail_cbpt = 'twotail';
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d_paired(low_cb, high_cb, n_perm, alpha_cbpt, tail_cbpt, t_cb);
log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(struct( ...
    'tag', report_tag, 'modality', 'ERSD', 'nR', n_pairs, 'n_input', num_subjects_in, ...
    'lbl_low', lblLow, 'lbl_high', lblHigh, 'n_perm', n_perm, 'alpha', alpha_cbpt, 'tail', tail_cbpt, ...
    'nT_cb', numel(t_cb), 'clusters', clusters, 'tvals', tvals_cl, 'thr', thr, ...
    't_plot', t_cb, 'dt_cb', dt_cb, 'drop_info', drop_info)));
ylims = ylim_from_effect_curve(cohens_d_cb);
shade_sig_clusters(gca, clusters, t_cb, dt_cb, ylims, alpha_cbpt);
ylim(ylims);
plot(t_cb, cohens_d_cb, 'k-', 'LineWidth', 3.5);
yline(0, '--'); xline(0, '--k');
xlabel('Time [s]'); ylabel('Cohen''s d_z');
xlim([-0.5 2]); box off; set(gca, 'FontSize', fsz-4);
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function plot_difference_timecourse_CBPT(low_group_timecourses, high_group_timecourses, out_name, fig_dir, fig_pos, fsz, t_vec, diff_ylab)
t_plot = t_vec(:)';
difference_timecourses = high_group_timecourses - low_group_timecourses;
figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([2 1]); hold on
difference_mean = mean(difference_timecourses, 1, 'omitnan');
difference_sem = std(difference_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(difference_timecourses), 1)), 1);
difference_ci = 1.96 * difference_sem;
difference_plot = shadedErrorBar(t_plot, difference_mean, difference_ci, 'lineProps', {'-'}, 'transparent', true);
set(difference_plot.mainLine, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.8);
set(difference_plot.patch, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.20);
set(difference_plot.edge(1), 'Color', 'none'); set(difference_plot.edge(2), 'Color', 'none');
yline(0, '--'); xline(0, '--k'); xlim([-0.5 2]);
ylim(ylim_from_single_mean_err(difference_mean, difference_ci));
ylabel(diff_ylab); box off; set(gca, 'FontSize', fsz-4);
nexttile; hold on
n_perm = 10000; alpha_cbpt = 0.05; tail_cbpt = 'twotail';
[low_cb, high_cb, t_cb, ~, dt_cb] = prepare_paired_cbpt_poststim(low_group_timecourses, high_group_timecourses, t_plot);
cohens_d_cb = paired_cohens_dz_curve(low_cb, high_cb);
[clusters, ~, ~] = ft_cluster_permutation_1d_paired(low_cb, high_cb, n_perm, alpha_cbpt, tail_cbpt, t_cb);
ylims = ylim_from_effect_curve(cohens_d_cb);
shade_sig_clusters(gca, clusters, t_cb, dt_cb, ylims, alpha_cbpt);
ylim(ylims);
plot(t_cb, cohens_d_cb, 'k-', 'LineWidth', 3.5);
yline(0, '--'); xline(0, '--k');
xlabel('Time [s]'); ylabel('Cohen''s d_z');
xlim([-0.5 2]); box off; set(gca, 'FontSize', fsz-4);
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function [low_cb, high_cb, t_cb, n_pairs, dt_cb, drop_info] = prepare_paired_cbpt_poststim(low, high, t_vec)
min_pairs = 3;
n_subjects_input = size(low, 1);
t_vec = t_vec(:)';
post = t_vec >= 0 & t_vec <= 2;
low_post = low(:, post); high_post = high(:, post); t_post = t_vec(post);
n_bins_post_stim = numel(t_post);
complete = all(isfinite(low_post) & isfinite(high_post), 2);
n_subjects_dropped = sum(~complete);
low_cb = low_post(complete, :); high_cb = high_post(complete, :);
n_pairs = size(low_cb, 1);
ok_t = sum(isfinite(low_cb) & isfinite(high_cb), 1) >= min_pairs;
n_bins_dropped = sum(~ok_t);
low_cb = low_cb(:, ok_t); high_cb = high_cb(:, ok_t); t_cb = t_post(ok_t);
if n_pairs < min_pairs, error('Fewer than %d subjects with complete post-stimulus data for CBPT.', min_pairs); end
if isempty(t_cb), error('No valid post-stimulus time points remain for CBPT.'); end
dt_cb = mean(diff(t_cb), 'omitnan');
drop_info = struct('n_subjects_input', n_subjects_input, 'n_subjects_cbpt', n_pairs, ...
    'n_subjects_dropped', n_subjects_dropped, 'n_bins_post_stim', n_bins_post_stim, ...
    'n_bins_cbpt', numel(t_cb), 'n_bins_dropped', n_bins_dropped, 'min_pairs', min_pairs);
fprintf('CBPT exclusions: subjects %d -> %d (dropped %d); time points %d -> %d (dropped %d)\n', ...
    drop_info.n_subjects_input, drop_info.n_subjects_cbpt, drop_info.n_subjects_dropped, ...
    drop_info.n_bins_post_stim, drop_info.n_bins_cbpt, drop_info.n_bins_dropped);
end

function dz = paired_cohens_dz_curve(low, high, min_pairs)
if nargin < 3, min_pairs = 3; end
nT = size(low, 2); dz = nan(1, nT);
for t = 1:nT
    paired_ok = isfinite(low(:, t)) & isfinite(high(:, t));
    if sum(paired_ok) < min_pairs, continue, end
    diff_t = high(paired_ok, t) - low(paired_ok, t);
    sd = std(diff_t);
    if sd > 0, dz(t) = mean(diff_t) / sd; end
end
end

function ylims = ylim_from_effect_curve(d_curve)
valid_d = d_curve(isfinite(d_curve));
mx = max(abs(valid_d), [], 'omitnan');
if isempty(valid_d) || ~isfinite(mx) || mx == 0, ylims = [-0.6 0.6];
else, ylims = [-max(mx + 0.1, 0.6), max(mx + 0.1, 0.6)]; end
end

function shade_sig_clusters(ax, clusters, t_plot_cb, dt_cb, ylims, alpha_cbpt)
axes(ax); %#ok<LAXES>
sig = false(1, numel(t_plot_cb));
for k = 1:numel(clusters), if clusters(k).p < alpha_cbpt, sig(clusters(k).idx) = true; end, end
run_start = [false, diff(sig) == 1]; run_end = [diff(sig) == -1, false];
if ~isempty(sig) && sig(1), run_start(1) = true; end
if ~isempty(sig) && sig(end), run_end(end) = true; end
starts = find(run_start); ends = find(run_end);
patch_alpha = 0.4 * ~any(sig) + 0.25 * any(sig);
for k = 1:numel(starts)
    t1 = t_plot_cb(starts(k)) - dt_cb / 2; t2 = t_plot_cb(ends(k)) + dt_cb / 2;
    patch([t1 t2 t2 t1], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.5 0.5 0.5], ...
        'FaceAlpha', patch_alpha, 'EdgeColor', 'none');
end
end

function finalize_tc_tiledlayout(tl)
drawnow; tl.OuterPosition = [0.06, 0.04, 0.90, 0.92];
end

function yl = ylim_from_single_mean_err(m, e)
lo = m - e; hi = m + e;
if isempty(lo) || all(~isfinite(lo)), yl = [-1, 1]; return, end
span = max(max(hi, [], 'omitnan') - min(lo, [], 'omitnan'), eps);
pad = 0.10 * span;
yl = [min(lo, [], 'omitnan') - pad, max(hi, [], 'omitnan') + pad];
end

function yl = ylim_from_mean_sem(m1, s1, m2, s2)
lo = min([m1 - s1, m2 - s2], [], 'omitnan');
hi = max([m1 + s1, m2 + s2], [], 'omitnan');
if isempty(lo) || all(~isfinite(lo)), yl = [-1, 1]; return, end
span = max(max(hi, [], 'omitnan') - min(lo, [], 'omitnan'), eps);
pad = 0.10 * span;
yl = [min(lo, [], 'omitnan') - pad, max(hi, [], 'omitnan') + pad];
end

function init_cbpt_report_file(report_path, meta)
if isempty(report_path), return, end
if isfile(report_path), delete(report_path); end
lines = {};
lines{end+1} = '=== AOC Split ERS/ERD Inverse CBPT Report (PAIRED / trial-level) ===';
lines{end+1} = sprintf('Script: %s', meta.script);
lines{end+1} = sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
lines{end+1} = sprintf('Task: %s | Split: %s', meta.task_tag, meta.split_label);
lines{end+1} = sprintf('Split variable: %s', meta.split_var);
lines{end+1} = sprintf('ERSD variable: %s', meta.ersd_var);
lines{end+1} = meta.split_info;
lines{end+1} = sprintf('Subjects contributing both groups: n=%d (TC n=%d)', meta.n_low_split, meta.n_low_tc);
lines{end+1} = sprintf('Outcome metric: %s', meta.metric);
lines{end+1} = 'CBPT method: FieldTrip ft_timelockstatistics (montecarlo, cluster maxsum, depsamplesT)';
lines{end+1} = 'Exclusions: no NaN imputation; subjects need complete post-stim time points; time points need >=3 paired subjects';
lines{end+1} = 'Defaults: n_perm=10000, clusteralpha=0.05, two-tailed, post-stimulus [0 2] s, no temporal binning';
lines{end+1} = sprintf('TC QC: subjects included for plotting if >=80%% finite samples in [0 2] s (n=%d of %d split subjects)', meta.n_low_tc, meta.n_low_split);
lines{end+1} = '';
append_lines_to_file(report_path, lines);
end

function lines = build_cbpt_report_lines(R)
lbl_lo = sanitize_label_for_fname(R.lbl_low);
lbl_hi = sanitize_label_for_fname(R.lbl_high);
nExtreme = sum(abs(R.tvals) > R.thr.tcrit & isfinite(R.tvals));
if isempty(R.clusters), maxClMass=0; maxClExtent=0;
else
    maxClMass = max([0, arrayfun(@(k) R.clusters(k).mass, 1:numel(R.clusters))]);
    maxClExtent = max([0, arrayfun(@(k) R.clusters(k).extent, 1:numel(R.clusters))]);
end
lines = {};
if isfield(R, 'n_input'), n_line = sprintf('n=%d pairs (of %d input)', R.nR, R.n_input);
else, n_line = sprintf('n=%d pairs', R.nR); end
lines{end+1} = sprintf('  [%s] %s tcrit=%.2f |t|>tcrit at %d timepts; max mass=%.1f; max extent=%d', ...
    R.tag, n_line, R.thr.tcrit, nExtreme, maxClMass, maxClExtent);
lines{end+1} = sprintf('    Method: n_perm=%d, alpha=%.3f, two-tailed paired (%s vs %s), post-stim [0 2] s, no temporal binning', ...
    R.n_perm, R.alpha, lbl_lo, lbl_hi);
if isfield(R, 'drop_info') && ~isempty(R.drop_info)
    d = R.drop_info;
    lines{end+1} = sprintf('    Exclusions: subjects %d -> %d (dropped %d, %.1f%% kept)', ...
        d.n_subjects_input, d.n_subjects_cbpt, d.n_subjects_dropped, ...
        100 * d.n_subjects_cbpt / max(d.n_subjects_input, 1));
    lines{end+1} = sprintf('    Exclusions: time points %d -> %d (dropped %d, %.1f%% kept; min_pairs=%d)', ...
        d.n_bins_post_stim, d.n_bins_cbpt, d.n_bins_dropped, ...
        100 * d.n_bins_cbpt / max(d.n_bins_post_stim, 1), d.min_pairs);
end
if ~isempty(R.clusters)
    for k = 1:numel(R.clusters)
        idx = R.clusters(k).idx;
        t_start = R.t_plot(idx(1)) - R.dt_cb/2;
        t_end = R.t_plot(idx(end)) + R.dt_cb/2;
        status = 'n.s.'; if R.clusters(k).p < R.alpha, status = 'SIGNIFICANT'; end
        lines{end+1} = sprintf('    Cluster %d: [%.3f, %.3f] s; mass=%.1f; p=%.4f; %s', ...
            k, t_start, t_end, R.clusters(k).mass, R.clusters(k).p, status);
    end
else
    lines{end+1} = '    No clusters formed';
end
lines{end+1} = '';
end

function txt = sanitize_label_for_fname(label)
txt = lower(label);
txt = regexprep(txt, '[^a-z0-9]+', '_');
txt = regexprep(txt, '_+', '_');
txt = regexprep(txt, '^_|_$', '');
end

function log_cbpt_report(report_path, lines)
for i = 1:numel(lines), fprintf('%s\n', lines{i}); end
append_lines_to_file(report_path, lines);
end

function append_lines_to_file(file_path, lines)
if isempty(file_path) || isempty(lines), return, end
fid = fopen(file_path, 'a');
if fid < 0, warning('Could not open CBPT report: %s', file_path); return, end
cleanup = onCleanup(@() fclose(fid));
for i = 1:numel(lines), fprintf(fid, '%s\n', lines{i}); end
end

function [clusters, tvals, thr] = ft_cluster_permutation_1d_paired(low_group_timecourses, high_group_timecourses, nPerm, alpha, tail, t_plot_ds)
if nargin < 5, tail = 'twotail'; end
if nargin < 6, t_plot_ds = []; end
nS = size(low_group_timecourses, 1);
nT = size(low_group_timecourses, 2);
if any(~isfinite(low_group_timecourses(:))) || any(~isfinite(high_group_timecourses(:)))
    error('CBPT input must not contain NaN. Use prepare_paired_cbpt_poststim first.');
end
df = nS - 1;
chan_label = 'metric';
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    time_vec = t_plot_ds(:)';
else
    time_vec = (0:nT-1)' / 20;
end
tl1 = struct('label', {{chan_label}}, 'time', time_vec, 'dimord', 'rpt_chan_time');
tl2 = tl1;
tl1.trial = reshape(low_group_timecourses, [nS, 1, nT]);
tl2.trial = reshape(high_group_timecourses, [nS, 1, nT]);
cfg_neigh = struct('label', chan_label, 'neighblabel', {{}});
cfg = struct();
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = cfg_neigh;
cfg.numrandomization = nPerm;
cfg.channel = chan_label;
cfg.latency = 'all';
if strcmpi(tail, 'onetail_pos')
    cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = alpha;
elseif strcmpi(tail, 'onetail_neg')
    cfg.tail = -1; cfg.clustertail = -1; cfg.alpha = alpha;
else
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = alpha / 2;
end
cfg.design = [1:nS, 1:nS; ones(1, nS), 2*ones(1, nS)];
cfg.uvar = 1; cfg.ivar = 2;
stat = ft_timelockstatistics(cfg, tl1, tl2);
tvals = stat.stat(1, :);
clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {});
has_neg = isfield(stat, 'negclusters') && ~isempty(stat.negclusters);
has_pos = isfield(stat, 'posclusters') && ~isempty(stat.posclusters);
clist = []; lmat = zeros(1, nT);
if has_pos, clist = stat.posclusters(:); lmat = stat.posclusterslabelmat; end
if has_neg
    lmat_neg = stat.negclusterslabelmat; npos = numel(clist);
    lmat(lmat_neg > 0) = npos + lmat_neg(lmat_neg > 0);
    if isempty(clist), clist = stat.negclusters(:); else, clist = [clist; stat.negclusters(:)]; end
end
for k = 1:numel(clist)
    idx = find(lmat(1, :) == k);
    if isempty(idx), continue, end
    idx = idx(1):idx(end);
    clusters(end+1).idx = idx; %#ok<AGROW>
    clusters(end).mass = sum(abs(tvals(idx)), 'omitnan');
    clusters(end).extent = numel(idx);
    clusters(end).p = clist(k).prob;
end
if strcmpi(tail, 'onetail_pos') || strcmpi(tail, 'onetail_neg')
    thr.tcrit = tinv(1 - alpha, df);
else
    thr.tcrit = tinv(1 - alpha/2, df);
end
end
