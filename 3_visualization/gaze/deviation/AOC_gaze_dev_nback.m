%% AOC Gaze Deviation - N-back conditions

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');

addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'gaze', 'deviation');
if ~isfolder(fig_dir), mkdir(fig_dir); end

fontSize = 40;
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
smooth_sec = 0.05;

cond_codes = [21 22 23];
cond_labels = {'1-back', '2-back', '3-back'};
et_fname = 'dataET_nback.mat';
blink_win = 50; % 100 ms at 500 Hz
min_trial_coverage = 0;
idx_viable = (t_plot >= 0) & (t_plot <= 2);
outlier_k_iqr = 1.5;
max_interp_gap_sec = 0.20;
min_subject_coverage = 0.85;
winsor_pct = 10;

%% Aggregate subject data 
merged_file = fullfile(feat_dir, 'AOC_merged_data_nback.mat');
S = load(merged_file, 'merged_data_nback');
T = struct2table(S.merged_data_nback);
uIDs = unique(T.ID);
nSubj = numel(uIDs);
gaze_tc = nan(nSubj, 3, Tf);
gaze_tc_raw = nan(nSubj, 3, Tf);

for s = 1:nSubj
    clc; fprintf('[GAZE DEV - NBACK] Computing gaze deviation time courses for Subject %d / %d \n', s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end

    et_file = fullfile(feat_dir, subj_folder, 'gaze', et_fname);
    if ~isfile(et_file)
        continue
    end

    [D, load_msg] = load_dataETlong_with_tmp_fallback(et_file);
    if isempty(D)
        warning('Skipping subject %s: could not load %s (%s)', subj_folder, et_file, load_msg);
        continue
    end
    if ~isfield(D, 'dataETlong') || ~isfield(D.dataETlong, 'trial') || ~isfield(D.dataETlong, 'trialinfo')
        continue
    end

    dataETlong = D.dataETlong;
    conds = parse_trialinfo_conds(dataETlong.trialinfo);
    if isempty(conds)
        continue
    end

    for c = 1:3
        tr_idx = find(conds == cond_codes(c));
        if isempty(tr_idx)
            continue
        end
        tc_trials = nan(numel(tr_idx), Tf);
        tc_trials_raw = nan(numel(tr_idx), Tf);
        for k = 1:numel(tr_idx)
            tr = tr_idx(k);
            if tr > numel(dataETlong.trial) || tr > numel(dataETlong.time)
                continue
            end

            data = double(dataETlong.trial{tr});
            t = double(dataETlong.time{tr});
            if size(data, 1) < 3 || isempty(t) || numel(t) ~= size(data, 2)
                continue
            end

            % Match feature extraction preprocessing.
            valid_idx = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
            data = data(1:3, valid_idx);
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
            tc_interp_raw = interp1(t, dev, t_plot, 'linear', NaN);
            if mean(isfinite(tc_interp(idx_viable))) < min_trial_coverage
                continue
            end
            tc_trials(k, :) = tc_interp;
            tc_trials_raw(k, :) = tc_interp_raw;
        end
        gaze_tc(s, c, :) = mean(tc_trials, 1, 'omitnan');
        gaze_tc_raw(s, c, :) = mean(tc_trials_raw, 1, 'omitnan');
    end
end

%% Plot baselined gaze deviation
close all
win_sm = max(1, round(smooth_sec * fs));
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:3
    X = squeeze(gaze_tc(:, c, :));
    X(~isfinite(X)) = NaN;

    % Exclude outlier subject trajectories before group statistics.
    med_metric = median(X(:, idx_viable), 2, 'omitnan');
    [X, ~] = exclude_outlier_trajectories(X, med_metric, outlier_k_iqr);

    % Fill short blink-related gaps to stabilize each subject trajectory.
    max_interp_gap_smp = max(1, round(max_interp_gap_sec * fs));
    for s = 1:size(X, 1)
        X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
    end

    % Keep only participants with sufficient valid coverage in the viable window.
    subj_cov = mean(isfinite(X(:, idx_viable)), 2);
    keep_cov = subj_cov >= min_subject_coverage;
    X = X(keep_cov, :);

    if win_sm > 1
        X = movmean(X, win_sm, 2, 'omitnan');
    end

    [m, se] = winsorized_nanmean_se(X, winsor_pct);

    eb = shadedErrorBar(t_plot, m, se, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Gaze Deviation [%]');
xlim([-0.5 2]);
ylim([-25 125]);
set(gca, 'FontSize', fontSize - 4);
box off
legend_handles = gobjects(1, 3);
for c = 1:3
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cond_labels, 'Location', 'northeast', 'FontSize', fontSize*0.666, 'Box', 'off');
saveas(gcf, fullfile(fig_dir, 'AOC_gaze_dev_timecourse_nback.png'));

%% Plot raw gaze deviation time courses
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:3
    X = squeeze(gaze_tc_raw(:, c, :));
    X(~isfinite(X)) = NaN;

    % Exclude outlier subject trajectories before group statistics.
    med_metric = median(X(:, idx_viable), 2, 'omitnan');
    [X, ~] = exclude_outlier_trajectories(X, med_metric, outlier_k_iqr);

    max_interp_gap_smp = max(1, round(max_interp_gap_sec * fs));
    for s = 1:size(X, 1)
        X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
    end

    subj_cov = mean(isfinite(X(:, idx_viable)), 2);
    keep_cov = subj_cov >= min_subject_coverage;
    X = X(keep_cov, :);

    if win_sm > 1
        X = movmean(X, win_sm, 2, 'omitnan');
    end

    [m, se] = winsorized_nanmean_se(X, winsor_pct);

    eb = shadedErrorBar(t_plot, m, se, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end

yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Gaze Deviation [px]');
xlim([-0.5 2]);
ylim([10 40]);
set(gca, 'FontSize', fontSize - 4);
box off
legend_handles = gobjects(1, 3);
for c = 1:3
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cond_labels, 'Location', 'northeast', 'FontSize', fontSize*0.666, 'Box', 'off');
saveas(gcf, fullfile(fig_dir, 'AOC_gaze_dev_timecourse_nback_raw.png'));

%%
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

function [D, err_msg] = load_dataETlong_with_tmp_fallback(et_file)
D = [];
err_msg = '';
try
    D = load(et_file, 'dataETlong');
    return
catch ME
    err_msg = 'direct load failed';
end

tmp_file = '';
try
    tmp_file = fullfile(tempdir, sprintf('AOC_tmp_%s.mat', char(java.util.UUID.randomUUID)));
    copyfile(et_file, tmp_file, 'f');
    D = load(tmp_file, 'dataETlong');
catch ME2
    err_msg = 'direct and tmp load failed';
    D = [];
end

if ~isempty(tmp_file) && isfile(tmp_file)
    delete(tmp_file);
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
X_keep = X(keep_subj, :);
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

function [m, se] = winsorized_nanmean_se(X, pct)
nT = size(X, 2);
m = nan(1, nT);
se = nan(1, nT);
for t = 1:nT
    xt = X(:, t);
    xt = xt(isfinite(xt));
    if isempty(xt)
        continue
    end
    xt = winsorize_values(xt, pct);
    m(t) = mean(xt);
    n = numel(xt);
    if n > 1
        se(t) = std(xt, 0) / sqrt(n);
    end
end
end

function xw = winsorize_values(x, pct)
xw = x;
if numel(xw) < 5 || pct <= 0
    return
end
lo = prctile(xw, pct);
hi = prctile(xw, 100 - pct);
xw(xw < lo) = lo;
xw(xw > hi) = hi;
end
