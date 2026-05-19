%% AOC Split Sternberg Alpha Loads (Time Course)
% Split participants by alpha-power slope across WM loads (2,4,6), then
% visualize gaze-deviation time courses for the two slope groups.
%
% Group definition:
%   Increase: slope > +threshold
%   Decrease: slope < -threshold
%   Flat:     |slope| <= threshold (excluded from inferential comparison)
%
% Output:
% - Inclusion/split figure
% - Gaze deviation time-course figure with CBPT effect-size shading
% - CSV for Python stats

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
pathAOC = paths.features;

addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'splits', 'SplitAlphaLoadsTimeCourse');
stats_dir = paths.splits_stats;
if ~isfolder(fig_dir), mkdir(fig_dir); end
if ~isfolder(stats_dir), mkdir(stats_dir); end

fprintf('\n=== AOC Split Alpha Loads Time Course ===\n');
fprintf('Figure directory: %s\n', fig_dir);

fig_pos = [0 0 1512 982];
fontSize = 20;
cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

% Match the robust thresholding strategy from AOC_split_AlphaAmpRed.m:
% exclude near-zero slopes within slope_zero_pct of the slope_ref_percentile of |slope|.
slope_zero_pct = 5;
slope_ref_percentile = 95;

%% Load merged data
merged_file = fullfile(feat_dir, 'AOC_merged_data_sternberg.mat');
if ~isfile(merged_file)
    error('Missing file: %s', merged_file);
end
S = load(merged_file, 'merged_data_sternberg');
if ~isfield(S, 'merged_data_sternberg')
    error('Variable merged_data_sternberg not found in %s', merged_file);
end
T = struct2table(S.merged_data_sternberg);

uIDs = unique(T.ID);
nSubj = numel(uIDs);

%% Compute alpha slope per subject
alpha_by_load = nan(nSubj, 3);
alpha_slope = nan(nSubj, 1);
for s = 1:nSubj
    sid = uIDs(s);
    rows = T(T.ID == sid, :);
    for c = 1:3
        cmask = rows.Condition == cond_vals(c);
        if any(cmask)
            alpha_by_load(s, c) = mean(rows.AlphaPower_FOOOF_bl_late(cmask), 'omitnan');
        end
    end
    if all(isfinite(alpha_by_load(s, :)))
        X = [ones(3, 1), cond_vals(:)];
        b = X \ alpha_by_load(s, :)';
        alpha_slope(s) = b(2);
    end
end

valid_slope = alpha_slope(isfinite(alpha_slope));
if isempty(valid_slope)
    error('No finite slope values found for split.');
end
slope_abs_ref = prctile(abs(valid_slope), slope_ref_percentile);
slope_threshold = (slope_zero_pct / 100) * slope_abs_ref;

increase_ids = uIDs(alpha_slope > slope_threshold);
decrease_ids = uIDs(alpha_slope < -slope_threshold);
flat_ids = uIDs(abs(alpha_slope) <= slope_threshold | ~isfinite(alpha_slope));

fprintf('\n=== Split Summary (alpha slope over WM load) ===\n');
fprintf('Subjects total: %d\n', nSubj);
fprintf('Zero-band setting: %.2f%% of %dth-percentile |slope| (ref=%.6f, cutoff=%.6f)\n', ...
    slope_zero_pct, slope_ref_percentile, slope_abs_ref, slope_threshold);
fprintf('Increase (> +thr): %d\n', numel(increase_ids));
fprintf('Decrease (< -thr): %d\n', numel(decrease_ids));
fprintf('Excluded flat/invalid: %d\n', numel(flat_ids));

if numel(increase_ids) < 2 || numel(decrease_ids) < 2
    error('Insufficient subjects per slope group.');
end

%% Split inclusion figure
figure('Position', fig_pos, 'Color', 'w');
hold on
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
yline(slope_threshold, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
yline(-slope_threshold, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);

x_vals = (1:nSubj)';
idx_dec = ismember(uIDs, decrease_ids);
idx_inc = ismember(uIDs, increase_ids);
idx_excl = ismember(uIDs, flat_ids);

h_excl = scatter(x_vals(idx_excl), alpha_slope(idx_excl), 80, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
h_dec = scatter(x_vals(idx_dec), alpha_slope(idx_dec), 80, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.8);
h_inc = scatter(x_vals(idx_inc), alpha_slope(idx_inc), 80, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.8);
xlabel('Subject');
ylabel('Alpha slope [a.u./load]');
title('Alpha-load slope split: Decrease (blue), Increase (red), Excluded (grey)');
legend([h_excl, h_dec, h_inc], {'Excluded', 'Decrease', 'Increase'}, 'Location', 'best', 'FontSize', fontSize - 2, 'Box', 'off');
set(gca, 'FontSize', fontSize);
box on
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlphaLoads_inclusion.png'));
close(gcf);

%% Build metrics + gaze deviation time course
metrics_raw = struct();
metrics_raw.Alpha = alpha_by_load;
metrics_raw.Dev = nan(nSubj, 3);

fs = 500;
t_full = -0.5:1/fs:3;
t_plot = t_full(2:end);
Tf = numel(t_plot);
dev_tc_raw = nan(nSubj, 3, Tf);

missing_gaze = {};
for s = 1:nSubj
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder), subj_folder = sid_str; end

    subj_rows = T(T.ID == sid, :);
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics_raw.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
        end
    end

    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', 'gaze_series_sternberg_trials.mat');
    if ~isfile(gaze_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    G = load(gaze_file);
    if ~isfield(G, 'trialinfo') || ~isfield(G, 'gaze_x') || ~isfield(G, 'gaze_y')
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    conds = parse_trialinfo_conds(G.trialinfo);
    if isempty(conds)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    for c = 1:3
        tr_idx = find(conds == cond_codes(c));
        if isempty(tr_idx), continue, end
        dev_mat = nan(numel(tr_idx), Tf);
        for k = 1:numel(tr_idx)
            tr = tr_idx(k);
            x = double(G.gaze_x{tr});
            y = double(G.gaze_y{tr});
            if isempty(x) || isempty(y) || numel(x) ~= numel(y), continue, end
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
        dev_tc_raw(s, c, :) = nanmean(dev_mat, 1);
    end
end

is_dec = ismember(uIDs, decrease_ids);
is_inc = ismember(uIDs, increase_ids);

%% Outlier exclusion (Tukey on per-condition gaze summary)
[metrics, dev_tc] = exclude_outliers_tukey_dev(metrics_raw, dev_tc_raw);

%% Baselined gaze time course (dB)
t_vec = linspace(-0.5, 3, Tf);
bl_idx = (t_vec >= -0.5) & (t_vec <= -0.25);
dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
dev_tc_dB = 10 * log10(dev_tc ./ dev_bl_3d);
dev_tc_dB(~isfinite(dev_tc_dB)) = NaN;

%% Time course figure with CBPT shading
fontSizeTC = 25;
rng(321)
plot_timecourse_with_effect_CBPT(dev_tc_dB, is_dec, is_inc, colors, ...
    'Gaze Deviation [dB]', 'gaze_deviation_dB', fig_dir, fig_pos, fontSizeTC, fs);

%% CSV export for Python stats
split_group = repmat("Excluded", nSubj, 1);
split_group(is_dec) = "Decrease";
split_group(is_inc) = "Increase";
included = is_dec | is_inc;

t_win_lo = 0.5;
t_win_hi = 1.5;
task_idx = (t_vec >= t_win_lo) & (t_vec <= t_win_hi);
dev_dB_by_load = squeeze(mean(dev_tc_dB(:, :, task_idx), 3, 'omitnan'));

rows_n = nSubj * numel(cond_vals);
ID = nan(rows_n, 1);
LoadValue = nan(rows_n, 1);
LoadLabel = strings(rows_n, 1);
Group = strings(rows_n, 1);
Included = false(rows_n, 1);
AlphaPower_FOOOF_bl_late = nan(rows_n, 1);
Dev_dB = nan(rows_n, 1);
GazeDev_dB_0p5_1p5s = nan(rows_n, 1);
AlphaSlope = nan(rows_n, 1);
SlopeThreshold = nan(rows_n, 1);
SlopeAbsRef = nan(rows_n, 1);
SlopeZeroPct = nan(rows_n, 1);
SlopeRefPercentile = nan(rows_n, 1);

r = 1;
for s = 1:nSubj
    for c = 1:numel(cond_vals)
        ID(r) = uIDs(s);
        LoadValue(r) = cond_vals(c);
        LoadLabel(r) = string(cond_labels{c});
        Group(r) = split_group(s);
        Included(r) = included(s);
        AlphaPower_FOOOF_bl_late(r) = metrics.Alpha(s, c);
        Dev_dB(r) = metrics.Dev(s, c);
        GazeDev_dB_0p5_1p5s(r) = dev_dB_by_load(s, c);
        AlphaSlope(r) = alpha_slope(s);
        SlopeThreshold(r) = slope_threshold;
        SlopeAbsRef(r) = slope_abs_ref;
        SlopeZeroPct(r) = slope_zero_pct;
        SlopeRefPercentile(r) = slope_ref_percentile;
        r = r + 1;
    end
end

split_stats = table(ID, LoadValue, LoadLabel, Group, Included, ...
    AlphaPower_FOOOF_bl_late, Dev_dB, GazeDev_dB_0p5_1p5s, AlphaSlope, SlopeThreshold, ...
    SlopeAbsRef, SlopeZeroPct, SlopeRefPercentile);
out_csv = fullfile(stats_dir, 'AOC_splitAlphaLoads_TimeCourse_sternberg_stats_input.csv');
writetable(split_stats, out_csv);
fprintf('\nSaved Python stats input CSV: %s\n', out_csv);
fprintf('Rows: %d (included rows: %d)\n', height(split_stats), sum(split_stats.Included));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

%% ========================= Local Functions =========================
function [metrics_out, dev_tc_out] = exclude_outliers_tukey_dev(metrics_in, dev_tc_in)
metrics_out = metrics_in;
dev_tc_out = dev_tc_in;
X = metrics_out.Dev;
for c = 1:3
    vals = X(:, c);
    vals_f = vals(isfinite(vals));
    if numel(vals_f) < 4, continue, end
    q1 = prctile(vals_f, 25);
    q3 = prctile(vals_f, 75);
    iqr_val = q3 - q1;
    if iqr_val <= 0, continue, end
    lo = q1 - 1.5 * iqr_val;
    hi = q3 + 1.5 * iqr_val;
    out_mask = (vals < lo) | (vals > hi);
    X(out_mask, c) = NaN;
    dev_tc_out(out_mask, c, :) = NaN;
end
metrics_out.Dev = X;
end

function conds = parse_trialinfo_conds(trialinfo)
conds = [];
if isempty(trialinfo), return, end
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

function plot_timecourse_with_effect_CBPT(tc, is_dec, is_inc, colors, ylab, save_tag, fig_dir, fig_pos, fsz, fs)
dt = 1 / fs;
nT = size(tc, 3);
t_plot = linspace(-0.5 + dt, 3, nT);
tc(~isfinite(tc)) = NaN;

Dall = squeeze(mean(tc(is_dec, :, :), 2, 'omitnan'));
Iall = squeeze(mean(tc(is_inc, :, :), 2, 'omitnan'));
win_sm = max(1, round(0.05 * fs));
if win_sm > 1
    Dall = movmean(Dall, win_sm, 2, 'omitnan');
    Iall = movmean(Iall, win_sm, 2, 'omitnan');
end

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile([2 1]); hold on
mD = mean(Dall, 1, 'omitnan');
mI = mean(Iall, 1, 'omitnan');
nD_fin = sum(isfinite(Dall), 1);
nI_fin = sum(isfinite(Iall), 1);
sD = std(Dall, 0, 1, 'omitnan') ./ max(sqrt(nD_fin), 1);
sI = std(Iall, 0, 1, 'omitnan') ./ max(sqrt(nI_fin), 1);
sD(~isfinite(sD)) = NaN;
sI(~isfinite(sI)) = NaN;
e1 = shadedErrorBar(t_plot, mD, sD, 'lineProps', {'-'}, 'transparent', true);
e2 = shadedErrorBar(t_plot, mI, sI, 'lineProps', {'-'}, 'transparent', true);
set(e1.mainLine, 'Color', colors(1,:), 'LineWidth', 3.5);
set(e2.mainLine, 'Color', colors(3,:), 'LineWidth', 3.5);
set(e1.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.25);
set(e2.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.25);
xline(0, '--k');
ylabel(ylab);
xlim([-0.5 3]);
set(gca, 'FontSize', fsz-4);
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1 leg_p2], {'Decrease', 'Increase'}, 'Location', 'northeast', 'FontSize', fsz-2, 'Box', 'off');
box on

nexttile; hold on
n_perm = 10000;
min_per_group = 3;
d = nan(1, nT);
for t = 1:nT
    x = Dall(:, t); y = Iall(:, t);
    x = x(isfinite(x)); y = y(isfinite(y));
    if numel(x) < min_per_group || numel(y) < min_per_group, continue, end
    sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / max(numel(x)+numel(y)-2, 1));
    d(t) = (mean(y) - mean(x)) / max(sp, eps);
end

ds_factor = 10;
Dall_ds = Dall(:, 1:ds_factor:end);
Iall_ds = Iall(:, 1:ds_factor:end);
t_plot_ds = t_plot(1:ds_factor:end);
dt_ds = ds_factor * dt;
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d(Dall_ds, Iall_ds, n_perm, 0.05, 'onetail_pos', t_plot_ds);

fprintf('\n[%s] CBPT summary\n', save_tag);
fprintf('  Threshold tcrit (one-tailed, alpha=0.05): %.3f\n', thr.tcrit);

sig_cluster = false(1, numel(t_plot_ds));
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        sig_cluster(clusters(k).idx) = true;
    end
end
sig_uncorr = (tvals_cl > thr.tcrit) & isfinite(tvals_cl);
sig = sig_cluster;
if ~any(sig) && any(sig_uncorr)
    sig = sig_uncorr;
end

% Print top observed clusters (significant or not), ranked by mass.
obs_clusters = get_observed_clusters(tvals_cl, thr.tcrit, t_plot_ds, dt_ds);
for oi = 1:numel(obs_clusters)
    for fi = 1:numel(clusters)
        if isequal(obs_clusters(oi).idx, clusters(fi).idx)
            obs_clusters(oi).p = clusters(fi).p;
            break
        end
    end
end
if isempty(obs_clusters)
    fprintf('  Observed clusters above tcrit: none\n');
    cand_clusters = get_candidate_clusters(tvals_cl, t_plot_ds, dt_ds, 0);
    if isempty(cand_clusters)
        fprintf('  Candidate clusters above 0: none\n');
    else
        [~, ordc] = sort([cand_clusters.mass], 'descend');
        top_k = min(5, numel(ordc));
        fprintf('  Candidate clusters above 0 (subthreshold), showing top %d by mass\n', top_k);
        for ii = 1:top_k
            c = cand_clusters(ordc(ii));
            fprintf('    #%d onset=%.3fs offset=%.3fs mass=%.3f extent=%d peak_t=%.3f peak_p_unc=%.4g\n', ...
                ii, c.onset, c.offset, c.mass, c.extent, c.peak_t, c.peak_p_unc);
        end
    end
else
    [~, ord] = sort([obs_clusters.mass], 'descend');
    top_k = min(5, numel(ord));
    fprintf('  Observed clusters above tcrit: %d (showing top %d by mass)\n', numel(obs_clusters), top_k);
    for ii = 1:top_k
        cidx = ord(ii);
        c = obs_clusters(cidx);
        p_str = 'n/a';
        if isfinite(c.p)
            p_str = sprintf('%.4f', c.p);
        end
        sig_str = '';
        if isfinite(c.p) && c.p < 0.05
            sig_str = ' [p<0.05]';
        end
        fprintf('    #%d onset=%.3fs offset=%.3fs mass=%.3f extent=%d p=%s%s\n', ...
            ii, c.onset, c.offset, c.mass, c.extent, p_str, sig_str);
    end
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
xlim([-0.5 3]);
box on
set(gca, 'FontSize', fsz-4);

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlphaLoads_timecourse_%s_CBPT.png', save_tag)));
close(gcf);
end

function [clusters, tvals, thr] = ft_cluster_permutation_1d(Dall, Iall, nPerm, alpha, tail, t_plot_ds)
if nargin < 5, tail = 'twotail'; end
if nargin < 6, t_plot_ds = []; end

nD = size(Dall, 1);
nI = size(Iall, 1);
nT = size(Dall, 2);
df = nD + nI - 2;

for t = 1:nT
    d = Dall(:, t); i = Iall(:, t);
    md = mean(d(isfinite(d)), 'omitnan');
    mi = mean(i(isfinite(i)), 'omitnan');
    if ~isfinite(md), md = 0; end
    if ~isfinite(mi), mi = 0; end
    d(~isfinite(d)) = md;
    i(~isfinite(i)) = mi;
    Dall(:, t) = d;
    Iall(:, t) = i;
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
tl1.trial = reshape(Dall, [nD, 1, nT]);

tl2 = struct();
tl2.label = {chan_label};
tl2.time = time_vec;
tl2.dimord = 'rpt_chan_time';
tl2.trial = reshape(Iall, [nI, 1, nT]);

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
cfg.latency = [0 3];

if strcmpi(tail, 'onetail_pos')
    cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = alpha;
elseif strcmpi(tail, 'onetail_neg')
    cfg.tail = -1; cfg.clustertail = -1; cfg.alpha = alpha;
else
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = alpha / 2;
end

cfg.design = [ones(1, nD), 2*ones(1, nI)];
cfg.ivar = 1;

stat = ft_timelockstatistics(cfg, tl1, tl2);
tvals = stat.stat(1, :);

post_idx = find(t_plot_ds >= 0);
if size(tvals, 2) < nT && ~isempty(post_idx)
    tvals_full = nan(1, nT);
    n_sel = min(size(tvals, 2), numel(post_idx));
    tvals_full(post_idx(1:n_sel)) = tvals(1:n_sel);
    tvals = tvals_full;
end

clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {}, 'p_extent', {});
has_pos = isfield(stat, 'posclusters') && ~isempty(stat.posclusters);
if has_pos
    clist = stat.posclusters;
    lmat = stat.posclusterslabelmat;
else
    clist = struct('prob', {});
    lmat = zeros(1, nT);
end

for k = 1:numel(clist)
    idx = find(lmat(1, :) == k);
    if isempty(idx), continue, end
    idx = idx(1):idx(end);
    valid = idx <= numel(post_idx);
    idx = post_idx(idx(valid));
    if isempty(idx), continue, end
    clusters(end+1).idx = idx; %#ok<AGROW>
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
end

function obs = get_observed_clusters(tvals, tcrit, t_plot_ds, dt_ds)
% Build contiguous positive clusters above threshold directly from observed t-values.
sig_obs = (tvals > tcrit) & isfinite(tvals);
obs = struct('idx', {}, 'mass', {}, 'extent', {}, 'onset', {}, 'offset', {}, 'p', {});
if ~any(sig_obs)
    return
end
run_start = [false, diff(sig_obs) == 1];
run_end = [diff(sig_obs) == -1, false];
if sig_obs(1), run_start(1) = true; end
if sig_obs(end), run_end(end) = true; end
starts = find(run_start);
ends = find(run_end);
for k = 1:numel(starts)
    idx = starts(k):ends(k);
    obs(k).idx = idx; %#ok<AGROW>
    obs(k).mass = sum(tvals(idx), 'omitnan');
    obs(k).extent = numel(idx);
    obs(k).onset = max(0, t_plot_ds(idx(1)) - dt_ds/2);
    obs(k).offset = t_plot_ds(idx(end)) + dt_ds/2;
    obs(k).p = NaN; % gets filled if matched to FT cluster output
end
end

function cand = get_candidate_clusters(tvals, t_plot_ds, dt_ds, min_t)
% Candidate clusters for diagnostics, independent of cluster-forming threshold.
% For one-tailed positive tests, use contiguous runs where tvals > min_t.
if nargin < 4
    min_t = 0;
end
mask = (tvals > min_t) & isfinite(tvals);
cand = struct('idx', {}, 'mass', {}, 'extent', {}, 'onset', {}, 'offset', {}, 'peak_t', {}, 'peak_p_unc', {});
if ~any(mask)
    return
end

run_start = [false, diff(mask) == 1];
run_end = [diff(mask) == -1, false];
if mask(1), run_start(1) = true; end
if mask(end), run_end(end) = true; end
starts = find(run_start);
ends = find(run_end);
df_approx = max(sum(isfinite(tvals)) - 2, 1);

for k = 1:numel(starts)
    idx = starts(k):ends(k);
    seg = tvals(idx);
    [peak_t, ~] = max(seg, [], 'omitnan');
    if isempty(peak_t) || ~isfinite(peak_t)
        peak_t = NaN;
        p_unc = NaN;
    else
        p_unc = 1 - tcdf(peak_t, df_approx); % one-tailed, uncorrected (approximate)
    end
    cand(k).idx = idx; %#ok<AGROW>
    cand(k).mass = sum(seg, 'omitnan');
    cand(k).extent = numel(idx);
    cand(k).onset = max(0, t_plot_ds(idx(1)) - dt_ds/2);
    cand(k).offset = t_plot_ds(idx(end)) + dt_ds/2;
    cand(k).peak_t = peak_t;
    cand(k).peak_p_unc = p_unc;
end
end
