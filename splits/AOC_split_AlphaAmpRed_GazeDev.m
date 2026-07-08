%% AOC Split ERS/ERD (Trial-Level) — GazeDev (Sternberg + N-back)
% Loads trial indices + within-group TFRs from AOC_splits_AlphaAmpRed_Prep.m
% (within-subject median split, conditions pooled, common baseline).
% Gaze deviation time courses are rebuilt from full dataET_* (includes
% [-1.5 -0.5] s baseline), matched by Trial ID; conditions collapsed.
% CBPT uses paired (depsamplesT) design.

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir_root = fullfile(paths.figures, 'splits', 'SplitERSERD');
stats_dir = paths.splits_stats;
if ~isfolder(fig_dir_root), mkdir(fig_dir_root); end
if ~isfolder(stats_dir), mkdir(stats_dir); end

split_file = fullfile(stats_dir, 'AOC_splitAlphaAmpRed_trialSplits.mat');
if ~isfile(split_file)
    error('Missing trial split file. Run AOC_splits_AlphaAmpRed_Prep.m first:\n  %s', split_file);
end
Ssplit = load(split_file, 'split_alpha_amp_red');
split_alpha_amp_red = Ssplit.split_alpha_amp_red;

fprintf('\n=== AOC Split ERS/ERD — Gaze Deviation (trial-level) ===\n');
fprintf('Split file: %s\n', split_file);

fig_pos = [0 0 1512 982];
fontSize = 40;
tfr_winsor_cfg = struct('enable', true, 'prctile', [2 98]);

tasks(1).tag = 'sternberg';
tasks(1).split_label = 'splitMedian_trial';
tasks(1).et_fname = 'dataET_sternberg.mat';
tasks(1).gaze_trials_file = 'AOC_gaze_matrix_sternberg_trials.mat';
tasks(1).gaze_trials_var = 'gaze_data_sternberg_trials';
tasks(1).topo_latency = [1 2];
tasks(1).fig_subdir = 'SternbergGazeDev';
tasks(1).group_lbl_low = 'ERD';
tasks(1).group_lbl_high = 'ERS';
tasks(1).ersd_var = 'ERSD_late';

tasks(2).tag = 'nback';
tasks(2).split_label = 'splitMedian_trial';
tasks(2).et_fname = 'dataET_nback.mat';
tasks(2).gaze_trials_file = 'AOC_gaze_matrix_nback_trials.mat';
tasks(2).gaze_trials_var = 'gaze_data_nback_trials';
tasks(2).topo_latency = [0 2];
tasks(2).fig_subdir = 'NbackGazeDev';
tasks(2).group_lbl_low = 'More ERD';
tasks(2).group_lbl_high = 'Less ERD';
tasks(2).ersd_var = 'ERSD_full';

for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
split_label = tk.split_label;
fig_prefix = sprintf('AOC_splitERSERD_GazeDev_%s_%s', task_tag, split_label);
fig_dir_task = fullfile(fig_dir_root, tk.fig_subdir);
if ~isfolder(fig_dir_task), mkdir(fig_dir_task); end

fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));
if ~isfield(split_alpha_amp_red, task_tag)
    warning('Skipping %s: missing split field.', task_tag);
    continue
end
task_split = split_alpha_amp_red.(task_tag);
subj_splits = task_split.subjects;
nSubj = numel(subj_splits);
if nSubj < 2
    warning('Skipping %s: fewer than 2 subjects.', task_tag);
    continue
end

uIDs = [subj_splits.ID]';
split_info_str = sprintf('Within-subject median on %s (conditions pooled, common baseline)', tk.ersd_var);
fprintf('%s\n', split_info_str);

ersd_mean_low = nan(nSubj, 1);
ersd_mean_high = nan(nSubj, 1);
ersd_thr = nan(nSubj, 1);
for s = 1:nSubj
    ersd_mean_low(s) = mean(subj_splits(s).ERSD(subj_splits(s).idx_low), 'omitnan');
    ersd_mean_high(s) = mean(subj_splits(s).ERSD(subj_splits(s).idx_high), 'omitnan');
    ersd_thr(s) = subj_splits(s).median_thr;
end

%% Inclusion
figure('Position', fig_pos, 'Color', 'w'); hold on
x = (1:nSubj)';
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
plot(x, ersd_thr, 'k-', 'LineWidth', 1.2);
h_low = scatter(x, ersd_mean_low, 70, colors(1, :), 'filled', 'MarkerFaceAlpha', 0.85);
h_high = scatter(x, ersd_mean_high, 70, colors(3, :), 'filled', 'MarkerFaceAlpha', 0.85);
xlabel('Participant (index)'); ylabel('Mean trial ERSD [dB]');
title(sprintf('Trial-level median split [%s]', task_tag), 'Interpreter', 'none');
legend([h_low, h_high], {tk.group_lbl_low, tk.group_lbl_high}, 'Location', 'best', 'FontSize', fontSize-4, 'Box', 'off');
set(gca, 'FontSize', fontSize); box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir_task, sprintf('%s_inclusion.png', fig_prefix)));
close(gcf);

%% Collect TFRs and ERSD / gaze metrics
tfr_low_all = {};
tfr_high_all = {};
metrics_ERSD = nan(nSubj, 2);
metrics_Dev = nan(nSubj, 2);
fs = 500;
blink_win = 50; % 100 ms at 500 Hz
min_trial_coverage = 0.8;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
idx_viable = (t_plot >= 0) & (t_plot <= 2);
dev_tc = nan(nSubj, 2, Tf);  % group: 1=low, 2=high
eeg_tc = nan(nSubj, 2, Tf);
missing_gaze = {};

gaze_mat_file = fullfile(feat_dir, tk.gaze_trials_file);
Gmat = load(gaze_mat_file, tk.gaze_trials_var);
Tg = struct2table(Gmat.(tk.gaze_trials_var));

fprintf('\n=== Aggregating GazeDev / EEG within trial-split groups ===\n');
fprintf('Gaze TC from full %s (baseline [-1.5 -0.5] s)\n', tk.et_fname);
for s = 1:nSubj
    clc; fprintf('[SPLIT ERS/ERD GAZEDEV - %s] Subject %d / %d\n', upper(task_tag), s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder), subj_folder = sid_str; end

    ids_low = subj_splits(s).trial_ids_low(:);
    ids_high = subj_splits(s).trial_ids_high(:);
    metrics_ERSD(s, 1) = ersd_mean_low(s);
    metrics_ERSD(s, 2) = ersd_mean_high(s);

    % EEG TC from prep (already occipital alpha dB)
    t_ersd = subj_splits(s).time(:);
    eeg_tc(s, 1, :) = interp1(t_ersd, subj_splits(s).ersd_tc_low(:), t_plot(:), 'linear', NaN);
    eeg_tc(s, 2, :) = interp1(t_ersd, subj_splits(s).ersd_tc_high(:), t_plot(:), 'linear', NaN);

    if isfield(subj_splits(s), 'tfr_low') && ~isempty(subj_splits(s).tfr_low)
        tfr_low_all{end+1} = subj_splits(s).tfr_low; %#ok<AGROW>
        tfr_high_all{end+1} = subj_splits(s).tfr_high; %#ok<AGROW>
    end

    % Gaze scalar
    rows = Tg(Tg.ID == sid, :);
    if ~isempty(rows) && ismember('GazeDeviationFullBL', rows.Properties.VariableNames)
        metrics_Dev(s, 1) = mean(rows.GazeDeviationFullBL(ismember(rows.Trial, ids_low)), 'omitnan');
        metrics_Dev(s, 2) = mean(rows.GazeDeviationFullBL(ismember(rows.Trial, ids_high)), 'omitnan');
    end

    % Full ET (includes [-1.5 -0.5] baseline); match trials by Trial ID
    et_file = fullfile(feat_dir, subj_folder, 'gaze', tk.et_fname);
    if ~isfile(et_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end
    [D, load_msg] = load_dataETlong_with_tmp_fallback(et_file);
    if isempty(D) || ~isfield(D, 'dataETlong') || ~isfield(D.dataETlong, 'trial')
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        fprintf('  Warning: ET load failed for %s (%s)\n', sid_str, load_msg);
        continue
    end
    et = D.dataETlong;
    trial_ids = parse_trial_ids(et.trialinfo);
    if isempty(trial_ids)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    for g = 1:2
        if g == 1, want = ids_low; else, want = ids_high; end
        tr_idx = find(ismember(trial_ids, want));
        if isempty(tr_idx), continue, end
        dev_mat = nan(numel(tr_idx), Tf);
        for k = 1:numel(tr_idx)
            tr = tr_idx(k);
            if tr > numel(et.trial) || tr > numel(et.time)
                continue
            end
            data = double(et.trial{tr});
            t = double(et.time{tr});
            if size(data, 1) < 2 || isempty(t) || numel(t) ~= size(data, 2)
                continue
            end

            % Match AOC_gaze_dev_* / gaze fex preprocessing
            valid_idx = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
            data = data(1:min(3, size(data, 1)), valid_idx);
            t = t(valid_idx);
            if isempty(t)
                continue
            end
            data(2, :) = 600 - data(2, :);
            data = remove_blinks(data, blink_win);

            idx_base = (t >= -1.5) & (t <= -0.5);
            if ~any(idx_base)
                continue
            end
            dev = sqrt((data(1, :) - 400).^2 + (data(2, :) - 300).^2);
            gd_base = mean(dev(idx_base), 'omitnan');
            if ~isfinite(gd_base) || gd_base <= 0
                continue
            end
            dev_pct = 100 * (dev - gd_base) ./ gd_base;
            tc_interp = interp1(t, dev_pct, t_plot, 'linear', NaN);
            if mean(isfinite(tc_interp(idx_viable))) < min_trial_coverage
                continue
            end
            dev_mat(k, :) = tc_interp;
        end
        if any(isfinite(dev_mat(:)))
            dev_tc(s, g, :) = mean(dev_mat, 1, 'omitnan');
        end
    end
end

%% Channels from first TFR
channels = {};
if ~isempty(tfr_low_all)
    channels = occ_channels_from_labels(tfr_low_all{1}.label);
elseif ~isempty(tfr_high_all)
    channels = occ_channels_from_labels(tfr_high_all{1}.label);
end
if isempty(channels)
    warning('No TFR for channel definition; EEG TFR/topo panels will be skipped.');
end

%% TFR / topo (collapsed; both groups)
color_map_tfr = customcolormap_preset('red-white-blue');
if ~isempty(tfr_low_all) && ~isempty(tfr_high_all)
    fprintf('\n=== Plotting collapsed TFRs ===\n');
    plot_group_tfrs_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, color_map_tfr, ...
        fig_dir_task, fig_prefix, fig_pos, fontSize, tfr_winsor_cfg, tk.group_lbl_low, tk.group_lbl_high);
    fprintf('\n=== Plotting collapsed topoplots ===\n');
    plot_group_topos_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, tk.topo_latency, ...
        fig_dir_task, fig_prefix, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);
end

%% Rainclouds
fprintf('\n=== Plotting rainclouds ===\n');
plot_paired_raincloud(metrics_ERSD(:,1), metrics_ERSD(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'ERSD [dB]', fig_dir_task, sprintf('%s_raincloud_ersd.png', fig_prefix), fig_pos, fontSize);
plot_paired_raincloud(metrics_Dev(:,1), metrics_Dev(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'Gaze deviation', fig_dir_task, sprintf('%s_raincloud_dev.png', fig_prefix), fig_pos, fontSize);

%% Time courses
close all
fontSizeTC = 40;
rng(123)
ds_factor = 10;
tc_viz_smooth_sec = 0.05;
t_vec = t_plot;
tc_window_idx = (t_vec >= 0) & (t_vec <= 2);
tc_complete_min_frac = 0.995;
keep_tc = true(nSubj, 1);
for g = 1:2
    Xg = reshape(dev_tc(:, g, :), nSubj, Tf);
    frac = mean(isfinite(Xg(:, tc_window_idx)), 2);
    has_end = isfinite(Xg(:, end));
    keep_tc = keep_tc & (frac >= tc_complete_min_frac) & has_end;
end
dev_tc(~keep_tc, :, :) = NaN;
eeg_tc(~keep_tc, :, :) = NaN;
n_tc = sum(keep_tc);
fprintf('Included subjects for TC: %d / %d\n', n_tc, nSubj);

Rall = reshape(dev_tc(keep_tc, 1, :), n_tc, Tf);
Aall = reshape(dev_tc(keep_tc, 2, :), n_tc, Tf);
EegR = reshape(eeg_tc(keep_tc, 1, :), n_tc, Tf);
EegA = reshape(eeg_tc(keep_tc, 2, :), n_tc, Tf);

gaze_ylabel = sprintf('Gaze Deviation\nChange [%%]');
cbpt_report_file = fullfile(stats_dir, sprintf('AOC_splitERSERD_GazeDev_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_split_AlphaAmpRed_GazeDev.m', ...
    'task_tag', task_tag, ...
    'split_label', split_label, ...
    'split_info', split_info_str, ...
    'ersd_var', tk.ersd_var, ...
    'group_lbl_low', tk.group_lbl_low, ...
    'group_lbl_high', tk.group_lbl_high, ...
    'n_low_split', nSubj, ...
    'n_high_split', nSubj, ...
    'n_low_tc', n_tc, ...
    'n_high_tc', n_tc, ...
    'metric', 'Gaze deviation change [% baseline] (paired within subject)'));

plot_paired_timecourse_CBPT(Rall, Aall, colors, gaze_ylabel, ...
    sprintf('%s_%s_GazeDev_pct', task_tag, split_label), ...
    fig_dir_root, fig_pos, fontSizeTC, fs, ds_factor, t_vec, tc_viz_smooth_sec, ...
    tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file, EegR, EegA, 'ERSD [dB]');
plot_paired_timecourse_individuals(Rall, Aall, colors, gaze_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_%s_GazeDev_pct_individuals_collapsed', task_tag, split_label), ...
    fig_dir_task, fig_pos, fontSizeTC, fs, t_vec, tc_viz_smooth_sec, ...
    tk.group_lbl_low, tk.group_lbl_high);

%% CSV
fprintf('\n=== Exporting CSV ===\n');
t_win = (t_vec >= 0) & (t_vec <= 2);
dev_sum_low = mean(Rall(:, t_win), 2, 'omitnan');
dev_sum_high = mean(Aall(:, t_win), 2, 'omitnan');
ids_tc = uIDs(keep_tc);
ID_col = []; Group_col = strings(0,1); Included_col = [];
ERSD_col = []; GazeDev_col = []; GazeSummary_col = [];
for s = 1:nSubj
    for g = 1:2
        ID_col(end+1,1) = uIDs(s); %#ok<AGROW>
        if g == 1
            Group_col(end+1,1) = string(tk.group_lbl_low); %#ok<AGROW>
            ERSD_col(end+1,1) = metrics_ERSD(s,1); %#ok<AGROW>
            GazeDev_col(end+1,1) = metrics_Dev(s,1); %#ok<AGROW>
        else
            Group_col(end+1,1) = string(tk.group_lbl_high); %#ok<AGROW>
            ERSD_col(end+1,1) = metrics_ERSD(s,2); %#ok<AGROW>
            GazeDev_col(end+1,1) = metrics_Dev(s,2); %#ok<AGROW>
        end
        incl = keep_tc(s);
        Included_col(end+1,1) = incl; %#ok<AGROW>
        if incl
            ix = find(ids_tc == uIDs(s), 1);
            if g == 1, GazeSummary_col(end+1,1) = dev_sum_low(ix); %#ok<AGROW>
            else, GazeSummary_col(end+1,1) = dev_sum_high(ix); %#ok<AGROW>
            end
        else
            GazeSummary_col(end+1,1) = NaN; %#ok<AGROW>
        end
    end
end
stats_tbl = table(ID_col, Group_col, Included_col, ERSD_col, GazeDev_col, GazeSummary_col, ...
    'VariableNames', {'ID','Group','Included', tk.ersd_var, 'GazeDeviationFullBL', 'GazeDev_pct_0_2s'});
csv_out = fullfile(stats_dir, sprintf('AOC_splitERSERD_GazeDev_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV: %s\n', csv_out);
fprintf('Missing gaze: %d\n', numel(unique(missing_gaze)));
end

%% ========================= Local Functions =========================
function tids = parse_trial_ids(trialinfo)
if isempty(trialinfo)
    tids = [];
    return
end
if isvector(trialinfo)
    tids = trialinfo(:);
elseif size(trialinfo, 2) >= 2
    tids = trialinfo(:, 2);
else
    tids = trialinfo(:, 1);
end
end

function [D, err_msg] = load_dataETlong_with_tmp_fallback(et_file)
D = [];
err_msg = '';
try
    D = load(et_file, 'dataETlong');
    return
catch
    err_msg = 'direct load failed';
end
tmp_file = '';
try
    tmp_file = fullfile(tempdir, sprintf('AOC_tmp_%s.mat', char(java.util.UUID.randomUUID)));
    copyfile(et_file, tmp_file, 'f');
    D = load(tmp_file, 'dataETlong');
catch
    err_msg = 'direct and tmp load failed';
    D = [];
end
if ~isempty(tmp_file) && isfile(tmp_file)
    delete(tmp_file);
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end+1} = lab; %#ok<AGROW>
    end
end
if isempty(ch), ch = labels; end
end

function subj_folder = resolve_subject_folder(subjects, sid)
subj_folder = '';
for i = 1:numel(subjects)
    sval = str2double(subjects{i});
    if isfinite(sval) && sval == sid
        subj_folder = subjects{i}; return
    end
end
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

function plot_paired_timecourse_individuals(Rall, Aall, colors, ylab, title_tag, save_tag, fig_dir, fig_pos, fsz, fs, t_vec, smooth_sec, lblLow, lblHigh)
nT = size(Rall, 2); t_plot = t_vec(:);
win_sm = max(1, round(smooth_sec * fs));
R = Rall; A = Aall; R(~isfinite(R)) = NaN; A(~isfinite(A)) = NaN;
if win_sm > 1, R = movmean(R, win_sm, 2, 'omitnan'); A = movmean(A, win_sm, 2, 'omitnan'); end
colR_light = colors(1,:)*0.35+0.65; colA_light = colors(3,:)*0.35+0.65;
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(2,1,'TileSpacing','compact');
nexttile; hold on; plot(t_plot, R', 'Color', colR_light, 'LineWidth', 0.8);
plot(t_plot, mean(R,1,'omitnan'), 'Color', colors(1,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblLow, size(R,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
nexttile; hold on; plot(t_plot, A', 'Color', colA_light, 'LineWidth', 0.8);
plot(t_plot, mean(A,1,'omitnan'), 'Color', colors(3,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); xlabel('Time [s]'); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblHigh, size(A,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_GazeDev_timecourse_%s.png', save_tag)));
close(gcf);
end

function plot_paired_timecourse_CBPT(Rall, Aall, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, t_vec, smooth_sec, lblLow, lblHigh, cbpt_report_file, EegR, EegA, eeg_ylab)
if nargin < 16, EegR = []; EegA = []; end
if nargin < 18, eeg_ylab = 'ERSD [dB]'; end
nT = size(Rall,2); t_plot = t_vec(:); dt = mean(diff(t_plot),'omitnan');
win_sm = max(1, round(smooth_sec * fs));
if win_sm > 1
    Rall = movmean(Rall, win_sm, 2, 'omitnan');
    Aall = movmean(Aall, win_sm, 2, 'omitnan');
end
nR = size(Rall,1);
show_eeg = ~isempty(EegR) && ~isempty(EegA) && size(EegR,1)==nR;
figure('Position', fig_pos, 'Color', 'w');
if show_eeg, tiledlayout(5,1,'TileSpacing','compact'); else, tiledlayout(3,1,'TileSpacing','compact'); end
ax_gaze = nexttile([2 1]); hold on
mR = mean(Rall,1,'omitnan'); mA = mean(Aall,1,'omitnan');
sR = std(Rall,0,1,'omitnan')./max(sqrt(sum(isfinite(Rall),1)),1);
sA = std(Aall,0,1,'omitnan')./max(sqrt(sum(isfinite(Aall),1)),1);
e1 = shadedErrorBar(t_plot, mR, sR, 'lineProps', {'-'}, 'transparent', true);
e2 = shadedErrorBar(t_plot, mA, sA, 'lineProps', {'-'}, 'transparent', true);
set(e1.mainLine,'Color',colors(1,:),'LineWidth',2.5); set(e2.mainLine,'Color',colors(3,:),'LineWidth',2.5);
set(e1.patch,'FaceColor',colors(1,:),'FaceAlpha',0.20); set(e2.patch,'FaceColor',colors(3,:),'FaceAlpha',0.20);
set(e1.edge(1),'Color','none'); set(e1.edge(2),'Color','none');
set(e2.edge(1),'Color','none'); set(e2.edge(2),'Color','none');
xline(0,'--k'); ylabel(ylab); xlim([-0.5 2]); box off; set(gca,'FontSize',fsz-4);
leg_p1 = patch(NaN,NaN,colors(1,:),'FaceAlpha',0.25,'EdgeColor',colors(1,:),'LineWidth',1.5);
leg_p2 = patch(NaN,NaN,colors(3,:),'FaceAlpha',0.25,'EdgeColor',colors(3,:),'LineWidth',1.5);
legend([leg_p1 leg_p2], {[' ' lblLow], [' ' lblHigh]}, 'Location','best','FontSize',fsz*0.75,'Box','off');

ax_d = nexttile; hold on
n_perm = 10000; alpha_cbpt = 0.05; tail_cbpt = 'twotail';
d = nan(1,nT);
for t = 1:nT
    x = Rall(:,t); y = Aall(:,t); ok = isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
    if numel(x)<3, continue, end
    dd = y-x; d(t) = mean(dd)/max(std(dd),eps);
end
Rall_ds = Rall(:,1:ds_factor:end); Aall_ds = Aall(:,1:ds_factor:end);
t_plot_ds = t_plot(1:ds_factor:end); dt_ds = ds_factor*dt;
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d_paired(Rall_ds, Aall_ds, n_perm, alpha_cbpt, tail_cbpt, t_plot_ds);
log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(struct( ...
    'tag', save_tag, 'modality', 'gaze', 'nR', nR, 'nA', nR, 'lbl_low', lblLow, 'lbl_high', lblHigh, ...
    'n_perm', n_perm, 'alpha', alpha_cbpt, 'tail', tail_cbpt, 'nT_ds', numel(t_plot_ds), ...
    'bin_ms', ds_factor*1000/fs, 'fs', fs, 'ds_factor', ds_factor, 'clusters', clusters, ...
    'tvals', tvals_cl, 'thr', thr, 't_plot', t_plot_ds, 'dt_ds', dt_ds)));
sig = false(1,numel(t_plot_ds));
for k = 1:numel(clusters)
    if clusters(k).p < 0.05, sig(clusters(k).idx) = true; end
end
draw_cohen_panel(t_plot, t_plot_ds, d, sig, dt_ds, fsz);
align_ylabel_to_reference(ax_d, ax_gaze);

if show_eeg
    if win_sm > 1
        EegR = movmean(EegR, win_sm, 2, 'omitnan');
        EegA = movmean(EegA, win_sm, 2, 'omitnan');
    end
    nexttile; hold on
    mR_e = mean(EegR,1,'omitnan'); mA_e = mean(EegA,1,'omitnan');
    sR_e = std(EegR,0,1,'omitnan')./max(sqrt(sum(isfinite(EegR),1)),1);
    sA_e = std(EegA,0,1,'omitnan')./max(sqrt(sum(isfinite(EegA),1)),1);
    ebR = shadedErrorBar(t_plot, mR_e, sR_e, 'lineProps', {'-'}, 'transparent', true);
    ebA = shadedErrorBar(t_plot, mA_e, sA_e, 'lineProps', {'-'}, 'transparent', true);
    set(ebR.mainLine,'Color',colors(1,:),'LineWidth',2.5); set(ebA.mainLine,'Color',colors(3,:),'LineWidth',2.5);
    set(ebR.patch,'FaceColor',colors(1,:),'FaceAlpha',0.20); set(ebA.patch,'FaceColor',colors(3,:),'FaceAlpha',0.20);
    set(ebR.edge(1),'Color','none'); set(ebR.edge(2),'Color','none');
    set(ebA.edge(1),'Color','none'); set(ebA.edge(2),'Color','none');
    yline(0,'--'); xline(0,'--k'); ylabel(eeg_ylab); xlabel('Time [s]'); xlim([-0.5 2]); box off; set(gca,'FontSize',fsz-4);
    ax_eeg = gca;
    ax_d_eeg = nexttile; hold on
    d_eeg = nan(1,nT);
    for t = 1:nT
        x = EegR(:,t); y = EegA(:,t); ok = isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
        if numel(x)<3, continue, end
        dd = y-x; d_eeg(t) = mean(dd)/max(std(dd),eps);
    end
    [clusters_e, tvals_e, thr_e] = ft_cluster_permutation_1d_paired(EegR(:,1:ds_factor:end), EegA(:,1:ds_factor:end), n_perm, alpha_cbpt, tail_cbpt, t_plot_ds);
    log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(struct( ...
        'tag', [save_tag '_EEG'], 'modality', 'EEG', 'nR', nR, 'nA', nR, 'lbl_low', lblLow, 'lbl_high', lblHigh, ...
        'n_perm', n_perm, 'alpha', alpha_cbpt, 'tail', tail_cbpt, 'nT_ds', numel(t_plot_ds), ...
        'bin_ms', ds_factor*1000/fs, 'fs', fs, 'ds_factor', ds_factor, 'clusters', clusters_e, ...
        'tvals', tvals_e, 'thr', thr_e, 't_plot', t_plot_ds, 'dt_ds', dt_ds)));
    sig_e = false(1,numel(t_plot_ds));
    for k = 1:numel(clusters_e)
        if clusters_e(k).p < 0.05, sig_e(clusters_e(k).idx) = true; end
    end
    draw_cohen_panel(t_plot, t_plot_ds, d_eeg, sig_e, dt_ds, fsz);
    ylabel('EEG Cohen''s d_z');
    align_ylabel_to_reference(ax_d_eeg, ax_eeg);
end
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_GazeDev_timecourse_%s_CBPT.png', save_tag)));
close(gcf);
end

function draw_cohen_panel(t_plot, t_plot_ds, d, sig, dt_ds, fsz)
run_start = [false, diff(sig)==1]; run_end = [diff(sig)==-1, false];
if ~isempty(sig) && sig(1), run_start(1)=true; end
if ~isempty(sig) && sig(end), run_end(end)=true; end
starts = find(run_start); ends = find(run_end);
d_fin = d(isfinite(d)); mx = max(abs(d_fin), [], 'omitnan');
if isempty(d_fin) || ~isfinite(mx) || mx==0, ylims=[-0.6 0.6];
else, ylims=[-max(mx+0.1,0.6), max(mx+0.1,0.6)]; end
patch_alpha = 0.4*(~any(sig)) + 0.25*any(sig);
for k = 1:numel(starts)
    t1 = max(0, t_plot_ds(starts(k))-dt_ds/2);
    t2 = t_plot_ds(ends(k))+dt_ds/2;
    patch([t1 t2 t2 t1], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.5 0.5 0.5], 'FaceAlpha', patch_alpha, 'EdgeColor', 'none');
end
ylim(ylims); plot(t_plot, d, 'k-', 'LineWidth', 3.5);
yline(0,'--'); xline(0,'--k'); xlabel('Time [s]'); ylabel('Cohen''s d_z');
xlim([-0.5 2]); box off; set(gca,'FontSize',fsz-4);
end

function align_ylabel_to_reference(ax_ref, ax_targets)
fig = ancestor(ax_ref, 'figure');
if isempty(fig), return, end
drawnow;
ref_lbl = ax_ref.YLabel; ref_lbl.Units = 'normalized'; ref_pos = ref_lbl.Position;
ax_ref.Units = 'pixels'; ref_ax_pos_pix = ax_ref.Position;
target_x_pix = ref_ax_pos_pix(1) + ref_pos(1)*ref_ax_pos_pix(3);
if ~iscell(ax_targets), ax_list = {ax_targets}; else, ax_list = ax_targets; end
for k = 1:numel(ax_list)
    ax_t = ax_list{k};
    if isequal(ax_t, ax_ref), continue, end
    lbl = ax_t.YLabel; lbl.Units = 'normalized'; pos = lbl.Position;
    ax_t.Units = 'pixels'; ax_t_pos_pix = ax_t.Position;
    x_norm_new = (target_x_pix - ax_t_pos_pix(1)) / max(ax_t_pos_pix(3), eps);
    lbl.Position = [x_norm_new, pos(2), pos(3)];
end
end

function init_cbpt_report_file(report_path, meta)
if isempty(report_path), return, end
if isfile(report_path), delete(report_path); end
lines = {};
lines{end+1} = '=== AOC Split ERS/ERD CBPT Report (PAIRED / trial-level) ===';
lines{end+1} = sprintf('Script: %s', meta.script);
lines{end+1} = sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
lines{end+1} = sprintf('Task: %s | Split: %s', meta.task_tag, meta.split_label);
lines{end+1} = sprintf('ERSD variable: %s', meta.ersd_var);
lines{end+1} = meta.split_info;
lines{end+1} = sprintf('Subjects contributing both groups: n=%d (TC n=%d)', meta.n_low_split, meta.n_low_tc);
lines{end+1} = sprintf('Outcome metric: %s', meta.metric);
lines{end+1} = 'CBPT method: FieldTrip ft_timelockstatistics (montecarlo, cluster maxsum, depsamplesT)';
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
lines{end+1} = sprintf('  [%s] n=%d pairs tcrit=%.2f |t|>tcrit at %d timepts; max mass=%.1f; max extent=%d', ...
    R.tag, R.nR, R.thr.tcrit, nExtreme, maxClMass, maxClExtent);
lines{end+1} = sprintf('    Method: n_perm=%d, alpha=%.3f, paired (%s vs %s), ds_factor=%d', ...
    R.n_perm, R.alpha, lbl_lo, lbl_hi, R.ds_factor);
if ~isempty(R.clusters)
    for k = 1:numel(R.clusters)
        idx = R.clusters(k).idx;
        t_start = R.t_plot(idx(1)) - R.dt_ds/2;
        t_end = R.t_plot(idx(end)) + R.dt_ds/2;
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

function [clusters, tvals, thr] = ft_cluster_permutation_1d_paired(Rall, Aall, nPerm, alpha, tail, t_plot_ds)
if nargin < 5, tail = 'twotail'; end
if nargin < 6, t_plot_ds = []; end
nS = size(Rall,1); nT = size(Rall,2);
df = nS - 1;
for t = 1:nT
    r = Rall(:,t); a = Aall(:,t); ok = isfinite(r)&isfinite(a);
    if ~any(ok), r(:)=0; a(:)=0;
    else, mr=mean(r(ok)); ma=mean(a(ok)); r(~ok)=mr; a(~ok)=ma; end
    Rall(:,t)=r; Aall(:,t)=a;
end
chan_label = 'metric';
if ~isempty(t_plot_ds) && numel(t_plot_ds)==nT, time_vec = t_plot_ds(:); else, time_vec = (0:nT-1)'/500; end
tl1 = struct('label',{{chan_label}},'time',time_vec,'dimord','rpt_chan_time');
tl2 = tl1; tl1.trial = reshape(Rall,[nS,1,nT]); tl2.trial = reshape(Aall,[nS,1,nT]);
cfg_neigh = struct('label',chan_label,'neighblabel',{{}});
cfg = struct();
cfg.method='montecarlo'; cfg.statistic='ft_statfun_depsamplesT';
cfg.correctm='cluster'; cfg.clusteralpha=alpha; cfg.clusterstatistic='maxsum';
cfg.minnbchan=0; cfg.neighbours=cfg_neigh; cfg.numrandomization=nPerm; cfg.channel=chan_label;
if ~isempty(t_plot_ds) && numel(t_plot_ds)==nT, cfg.latency=[0 2]; else, cfg.latency='all'; end
if strcmpi(tail,'onetail_pos'), cfg.tail=1; cfg.clustertail=1; cfg.alpha=alpha;
elseif strcmpi(tail,'onetail_neg'), cfg.tail=-1; cfg.clustertail=-1; cfg.alpha=alpha;
else, cfg.tail=0; cfg.clustertail=0; cfg.alpha=alpha/2; end
cfg.design = [1:nS, 1:nS; ones(1,nS), 2*ones(1,nS)];
cfg.uvar=1; cfg.ivar=2;
stat = ft_timelockstatistics(cfg, tl1, tl2);
tvals = stat.stat(1,:);
post_idx = find(t_plot_ds >= 0 & t_plot_ds <= 2);
if ~isempty(t_plot_ds) && numel(t_plot_ds)==nT && size(tvals,2)<nT && ~isempty(post_idx)
    tvals_full = nan(1,nT); n_sel=min(size(tvals,2),numel(post_idx));
    tvals_full(post_idx(1:n_sel)) = tvals(1:n_sel); tvals = tvals_full;
end
clusters = struct('idx',{},'mass',{},'extent',{},'p',{},'p_extent',{});
has_neg = isfield(stat,'negclusters') && ~isempty(stat.negclusters);
has_pos = isfield(stat,'posclusters') && ~isempty(stat.posclusters);
clist=[]; lmat=zeros(1,nT);
if has_pos, clist=stat.posclusters(:); lmat=stat.posclusterslabelmat; end
if has_neg
    lmat_neg=stat.negclusterslabelmat; npos=numel(clist);
    lmat(lmat_neg>0)=npos+lmat_neg(lmat_neg>0);
    if isempty(clist), clist=stat.negclusters(:); else, clist=[clist; stat.negclusters(:)]; end
end
for k = 1:numel(clist)
    idx = find(lmat(1,:)==k); if isempty(idx), continue, end
    idx = idx(1):idx(end);
    clusters(end+1).idx=idx; %#ok<AGROW>
    clusters(end).mass=sum(abs(tvals(idx)),'omitnan');
    clusters(end).extent=numel(idx);
    clusters(end).p=clist(k).prob; clusters(end).p_extent=clist(k).prob;
end
if strcmpi(tail,'onetail_pos') || strcmpi(tail,'onetail_neg')
    thr.tcrit = tinv(1-alpha, df);
else
    thr.tcrit = tinv(1-alpha/2, df);
end
thr.mass=NaN; thr.extent=NaN;
end
