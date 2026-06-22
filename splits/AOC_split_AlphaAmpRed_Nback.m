%% AOC Split ERS/ERD (N-back)
% Subject-level split on pooled ERSD_full (median or fixed threshold), then
% ERD vs ERS gaze deviation time courses (percent baseline) with
% cluster-based permutation testing (FieldTrip).

%% Setup
startup
[subjects, paths, colors] = setup('AOC');
pathAOC = paths.features;

addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'splits', 'SplitERSERD', 'Nback');
stats_dir = paths.splits_stats;
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
if ~isfolder(stats_dir)
    mkdir(stats_dir);
end

fprintf('\n=== AOC Split ERS/ERD (N-back) ===\n');
fprintf('Figure directory: %s\n', fig_dir);

fig_pos = [0 0 1512 982];
fontSize = 30;

% Split definition 
nback_split_mode = 'median';
nback_split_fixed_value = -0.05;
if strcmpi(nback_split_mode, 'fixed')
    split_label = sprintf('splitFixed_%s', strrep(sprintf('%.3f', nback_split_fixed_value), '.', 'p'));
elseif strcmpi(nback_split_mode, 'median')
    split_label = 'splitMedian';
else
    error('nback_split_mode must be ''median'' or ''fixed''.');
end

task_tag = 'nback';
fig_prefix = sprintf('AOC_splitERSERD_%s_%s', task_tag, split_label);

cond_vals = [1 2 3];
cond_codes = [21 22 23];
cond_labels = {'1-back', '2-back', '3-back'};
gaze_fname = 'gaze_series_nback_trials.mat';
tfr_fname_primary = 'tfr_nback.mat';
tfr_fname_long = 'tfr_nback_long.mat';
tfr_vars = {'tfr1_bl', 'tfr2_bl', 'tfr3_bl'};
ersd_tc_fname = 'ersd_nback_timecourse.mat';

%% Load merged subject table
merged_file = fullfile(feat_dir, 'AOC_merged_data_nback.mat');
if ~isfile(merged_file)
    error('Missing merged file: %s', merged_file);
end
S = load(merged_file, 'merged_data_nback');
if ~isfield(S, 'merged_data_nback')
    error('Variable merged_data_nback not found in: %s', merged_file);
end
T = struct2table(S.merged_data_nback);
ersd_var = 'ERSD_full';
if ~ismember(ersd_var, T.Properties.VariableNames)
    error('Variable %s not found in merged_data_nback.', ersd_var);
end

uIDs = unique(T.ID);
nSubj = numel(uIDs);
ersd_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    ersd_mean(i) = mean(T.(ersd_var)(mask), 'omitnan');
end

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
switch lower(nback_split_mode)
    case 'median'
        ersd_split_threshold = median(valid_ersd_mean, 'omitnan');
        erd_ids = uIDs(split_valid & (ersd_mean < ersd_split_threshold));
        ers_ids = uIDs(split_valid & (ersd_mean >= ersd_split_threshold));
        zero_ids = [];
        invalid_ids = uIDs(~split_valid);
        split_info_str = sprintf('Median split at %.4f', ersd_split_threshold);
    case 'fixed'
        ersd_split_threshold = nback_split_fixed_value;
        erd_ids = uIDs(split_valid & (ersd_mean < ersd_split_threshold));
        ers_ids = uIDs(split_valid & (ersd_mean >= ersd_split_threshold));
        zero_ids = [];
        invalid_ids = uIDs(~split_valid);
        split_info_str = sprintf('Fixed split at %.4f', ersd_split_threshold);
    otherwise
        error('Unsupported nback_split_mode: %s', nback_split_mode);
end

fprintf('\n=== Split Summary [%s | %s] (%s, pooled loads) ===\n', task_tag, nback_split_mode, ersd_var);
fprintf('Subjects total: %d\n', nSubj);
fprintf('%s\n', split_info_str);
fprintf('Stronger ERD: %d\n', numel(erd_ids));
fprintf('Weaker ERD: %d\n', numel(ers_ids));
if ~isempty(invalid_ids)
    fprintf('Excluded (invalid/pathological ERSD): %d\n', numel(invalid_ids));
end

if numel(erd_ids) < 2 || numel(ers_ids) < 2
    error('Insufficient subjects per split group (ERD=%d, ERS=%d).', ...
        numel(erd_ids), numel(ers_ids));
end

%% ERSD split inclusion figure
fprintf('\n=== Plotting ERSD split inclusion figure ===\n');
figure('Position', fig_pos, 'Color', 'w');
hold on
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
yline(ersd_split_threshold, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
x_vals = (1:nSubj)';
idx_erd = ismember(uIDs, erd_ids);
idx_ers = ismember(uIDs, ers_ids);
idx_excl = ismember(uIDs, zero_ids) | ismember(uIDs, invalid_ids);
h_excl = scatter(x_vals(idx_excl), ersd_mean(idx_excl), 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
h_erd = scatter(x_vals(idx_erd), ersd_mean(idx_erd), 80, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.8);
h_ers = scatter(x_vals(idx_ers), ersd_mean(idx_ers), 80, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
xlabel('Participant (index)');
ylabel('ERSD [dB]');
title(sprintf('ERS/ERD Split [%s | %s]', task_tag, nback_split_mode), 'Interpreter', 'none');
legend([h_excl, h_erd, h_ers], {'Excluded', 'Stronger ERD', 'Weaker ERD'}, 'Location', 'best', 'FontSize', fontSize - 2, 'Box', 'off');
set(gca, 'FontSize', fontSize);
box off
saveas(gcf, fullfile(fig_dir, sprintf('%s_inclusion.png', fig_prefix)));
close(gcf);

%% Preallocate
tfr_red = cell(1, 3);
tfr_amp = cell(1, 3);
tfr_red_subj = [];
tfr_amp_subj = [];
for c = 1:3
    tfr_red{c} = {};
    tfr_amp{c} = {};
end

metrics = struct();
metrics.ERSD = nan(nSubj, 3);
metrics.Dev = nan(nSubj, 3);

fs = 500;
t_full = -0.5:1/fs:3;
t_plot = t_full(2:end);
Tf = numel(t_plot);
dev_tc = nan(nSubj, 3, Tf);
ersd_tc = nan(nSubj, 3, Tf);

missing_ersd = {};
missing_tfr = {};
missing_gaze = {};

%% Aggregate per-subject data
fprintf('\n=== Aggregating N-back data (%d subjects) ===\n', nSubj);
for s = 1:nSubj
    clc; fprintf('[SPLIT ERS/ERD - NBACK] Subject %d / %d\n', s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end
    subj_rows = T(T.ID == sid, :);

    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics.ERSD(s, c) = mean(subj_rows.(ersd_var)(cmask), 'omitnan');
            metrics.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
        end
    end

    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');

    tfr_file = fullfile(eeg_dir, tfr_fname_primary);
    if ~isfile(tfr_file)
        tfr_file = fullfile(eeg_dir, tfr_fname_long);
    end
    if ~isfile(tfr_file)
        missing_tfr{end+1} = sid_str; 
    else
        try
            R = load(tfr_file, tfr_vars{:});
            tfr_conds = {R.(tfr_vars{1}), R.(tfr_vars{2}), R.(tfr_vars{3})};
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
    end

    try
        E = load(fullfile(eeg_dir, ersd_tc_fname), 'ersd_timecourse');
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

    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', gaze_fname);
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

is_red = ismember(uIDs, erd_ids);
is_amp = ismember(uIDs, ers_ids);

%% Occipital channels (from baselined TFR)
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

%% Gaze percent baseline and time-course CBPT
close all
Tf = size(dev_tc, 3);
t_vec = linspace(-0.5, 3, Tf);
bl_idx = (t_vec >= -1.5) & (t_vec <= -0.5);
dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
gaze_tc = (dev_tc ./ dev_bl_3d - 1) * 100;
gaze_tc(~isfinite(gaze_tc)) = NaN;
gaze_ylabel = 'Gaze Deviation [%]';

tc_window_idx = (t_vec >= 0) & (t_vec <= 3);
tc_complete_min_frac = 0.995;
tc_finite_frac = squeeze(mean(isfinite(dev_tc(:, :, tc_window_idx)), 3));
tc_has_endpoint = squeeze(isfinite(dev_tc(:, :, end)));
tc_complete_by_cond = (tc_finite_frac >= tc_complete_min_frac) & tc_has_endpoint;
tc_complete_subj = all(tc_complete_by_cond, 2);
tc_excluded_subj = ~tc_complete_subj;
gaze_tc(tc_excluded_subj, :, :) = NaN;

is_red_tc = is_red & tc_complete_subj;
is_amp_tc = is_amp & tc_complete_subj;
fprintf('Time-course completeness filter: finite frac >= %.3f in [0,3]s + finite endpoint at 3s\n', tc_complete_min_frac);
fprintf('Excluded incomplete gaze time-course subjects (in split groups): %d\n', sum(tc_excluded_subj & (is_red | is_amp)));

fprintf('\n=== Preparing ERSD time courses ===\n');
eeg_tc = ersd_tc;
eeg_tc(tc_excluded_subj, :, :) = NaN;

fontSizeTC = 25;
rng(123);
ds_factor = 50;
tc_viz_smooth_sec = 0.10;

plot_timecourse_with_effect_CBPT(gaze_tc, is_red_tc, is_amp_tc, colors, ...
    gaze_ylabel, sprintf('%s_%s_gaze_deviation_pct', task_tag, split_label), fig_dir, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc, false, 'ERSD [dB]', tc_viz_smooth_sec);

for c = 1:numel(cond_vals)
    tc_cond = gaze_tc(:, c, :);
    eeg_tc_cond = eeg_tc(:, c, :);
    save_tag_cond = sprintf('%s_%s_gaze_deviation_pct_%s', task_tag, split_label, sanitize_label_for_fname(cond_labels{c}));
    plot_timecourse_with_effect_CBPT(tc_cond, is_red_tc, is_amp_tc, colors, ...
        gaze_ylabel, save_tag_cond, fig_dir, fig_pos, fontSizeTC, fs, ds_factor, eeg_tc_cond, false, 'ERSD [dB]', tc_viz_smooth_sec);
end

%% Diagnostics
fprintf('\n=== Data Diagnostics [nback] ===\n');
fprintf('Missing ERSD timecourse (%s): %d\n', ersd_tc_fname, numel(unique(missing_ersd)));
fprintf('Missing TFR files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

%% Export CSV for Python statistics
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
            Group_col(r) = "ERD";
            Included_col(r) = is_red_tc(s);
        elseif is_amp(s)
            Group_col(r) = "ERS";
            Included_col(r) = is_amp_tc(s);
        else
            Group_col(r) = "Excluded";
            Included_col(r) = false;
        end
    end
end

if ~isvarname(ersd_var)
    error('ERSD column name is not a valid table variable name: %s', ersd_var);
end
stats_tbl = table( ...
    ID_col, LoadValue_col, LoadLabel_col, Group_col, Included_col, ...
    ERSD_col, GazeDev_col, GazeSummary_col, ...
    'VariableNames', {'ID', 'LoadValue', 'LoadLabel', 'Group', 'Included', ...
    ersd_var, 'GazeDeviationFullBL', 'GazeDev_pct_0_2s'});

csv_out = fullfile(stats_dir, sprintf('AOC_splitERSERD_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('CSV saved to: %s\n', csv_out);

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

function plot_timecourse_with_effect_CBPT(tc, is_red, is_amp, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs, ds_factor, eeg_tc, addEEG_TC, eeg_ylab, smooth_sec)
if nargin < 11 || isempty(ds_factor), ds_factor = 50; end
if nargin < 12, eeg_tc = []; end
if nargin < 13, addEEG_TC = ~isempty(eeg_tc); end
if nargin < 14, eeg_ylab = 'ERSD [dB]'; end
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
xlim([-0.5 3]);
box off
set(gca, 'FontSize', fsz-4);
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1 leg_p2], {'Stronger ERD', 'Weaker ERD'}, 'Location', 'northeast', 'FontSize', fsz-2, 'Box', 'off');
if contains(save_tag, 'gaze_deviation_pct')
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
        t_start = t_lo - dt_ds/2;
        t_end = t_hi + dt_ds/2;
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
box off
set(gca, 'FontSize', fsz-4);
if ~any(sig_cluster) && any(sig_uncorr)
    title({'WARNING: No significant clusters; shading shows uncorrected t > t_{crit}'}, ...
        'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'Interpreter', 'tex');
end

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
    xlim([-0.5 3]);
    box off
    set(gca, 'FontSize', fsz-4);
    leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
    leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
    drawnow;
    pos_g = get(ax_gaze, 'Position');
    pos_e = get(ax_eeg, 'Position');
    left_min = min(pos_g(1), pos_e(1));
    set(ax_gaze, 'Position', [left_min, pos_g(2), pos_g(3) + (pos_g(1) - left_min), pos_g(4)]);
    set(ax_eeg,  'Position', [left_min, pos_e(2), pos_e(3) + (pos_e(1) - left_min), pos_e(4)]);

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
    box off
    set(gca, 'FontSize', fsz-4);
    if ~any(sig_cluster_eeg) && any(sig_uncorr_eeg)
        title({'WARNING: No significant clusters; shading shows uncorrected t > t_{crit}'}, ...
            'Color', [0.8 0 0], 'FontSize', max(8, fsz-6), 'Interpreter', 'tex');
    end
end

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitERSERD_timecourse_%s_CBPT.png', save_tag)));
close(gcf);
end

function [clusters, tvals, thr, maxMassNull, maxExtentNull] = ft_cluster_permutation_1d(Rall, Aall, nPerm, alpha, tail, t_plot_ds)
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
tl1.trial = reshape(Rall, [nR, 1, nT]);

tl2 = struct();
tl2.label = {chan_label};
tl2.time = time_vec;
tl2.dimord = 'rpt_chan_time';
tl2.trial = reshape(Aall, [nA, 1, nT]);

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
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    cfg.latency = [0 3];
else
    cfg.latency = 'all';
end
cfg.channel = chan_label;

if strcmpi(tail, 'onetail_pos')
    cfg.tail = 1;
    cfg.clustertail = 1;
    cfg.alpha = alpha;
elseif strcmpi(tail, 'onetail_neg')
    cfg.tail = -1;
    cfg.clustertail = -1;
    cfg.alpha = alpha;
else
    cfg.tail = 0;
    cfg.clustertail = 0;
    cfg.alpha = alpha / 2;
end

cfg.design = [ones(1, nR), 2*ones(1, nA)];
cfg.ivar = 1;

stat = ft_timelockstatistics(cfg, tl1, tl2);

tvals = stat.stat(1, :);

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

clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {}, 'p_extent', {});

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
    idx = idx(1):idx(end);
    if ~isempty(post_idx)
        valid = idx <= numel(post_idx);
        idx = post_idx(idx(valid));
        if isempty(idx), continue, end
    end
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
thr.mass = NaN;
thr.extent = NaN;

maxMassNull = [];
maxExtentNull = [];
end

function eeg_tc = extract_alpha_timecourse_tfr(tfr_red, tfr_amp, tfr_red_subj, tfr_amp_subj, nSubj, channels, Tf, t_target)
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
