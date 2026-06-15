%% AOC Gaze Velocity — N-Back & Sternberg (Baselined Difference Time Courses)
% Loads dataETlong per subject, computes 2D gaze speed (computeEyeVelocity),
% baselines each trial to percent change in [-1.5, -0.5] s, and plots the
% high-minus-low load difference (3-back - 1-back; WM load 6 - WM load 2).
%
% Key outputs:
%   Baselined gaze velocity difference time courses (per task + combined figure)

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC', 0);
addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'gaze', 'velocity');
if ~isfolder(fig_dir), mkdir(fig_dir); end

fontSize = 25;
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
smooth_sec = 0.05;
blink_win = 50;
idx_viable = (t_plot >= 0) & (t_plot <= 2);
outlier_k_iqr = 1.5;
max_interp_gap_sec = 0.20;
min_subject_coverage = 0.85;
winsor_pct = 10;
vel_win = 50; % movmean window for computeEyeVelocity (samples)

task_cfg(1) = struct( ...
    'name', 'nback', ...
    'merged_file', fullfile(feat_dir, 'AOC_merged_data_nback.mat'), ...
    'merged_var', 'merged_data_nback', ...
    'et_fname', 'dataET_nback.mat', ...
    'cond_codes', [21 22 23], ...
    'low_idx', 1, ...
    'high_idx', 3, ...
    'diff_label', '3-back - 1-back', ...
    'min_trial_coverage', 0);

task_cfg(2) = struct( ...
    'name', 'sternberg', ...
    'merged_file', fullfile(feat_dir, 'AOC_merged_data_sternberg.mat'), ...
    'merged_var', 'merged_data_sternberg', ...
    'et_fname', 'dataET_sternberg.mat', ...
    'cond_codes', [22 24 26], ...
    'low_idx', 1, ...
    'high_idx', 3, ...
    'diff_label', 'WM load 6 - WM load 2', ...
    'min_trial_coverage', 0.8);

win_sm = max(1, round(smooth_sec * fs));
diff_tc = struct('name', {}, 'label', {}, 'm', {}, 'se', {}, 'n', {});

for task = 1:numel(task_cfg)
    cfg = task_cfg(task);
    fprintf('\n=== Gaze velocity: %s ===\n', cfg.name);

    S = load(cfg.merged_file, cfg.merged_var);
    T = struct2table(S.(cfg.merged_var));
    uIDs = unique(T.ID);
    nSubj = numel(uIDs);
    gaze_tc = nan(nSubj, 3, Tf);

    for s = 1:nSubj
        fprintf('Computing Subject %d / %d\n', s, nSubj);
        sid = uIDs(s);
        subj_folder = resolve_subject_folder(subjects, sid);
        if isempty(subj_folder)
            subj_folder = num2str(sid);
        end

        et_file = fullfile(feat_dir, subj_folder, 'gaze', cfg.et_fname);
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
            tr_idx = find(conds == cfg.cond_codes(c));
            if isempty(tr_idx)
                continue
            end
            tc_trials = nan(numel(tr_idx), Tf);
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

                valid_idx = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
                data = data(1:3, valid_idx);
                t = t(valid_idx);
                if isempty(t)
                    continue
                end
                data(2, :) = 600 - data(2, :);
                data = remove_blinks(data, blink_win);

                trialET = struct();
                trialET.fsample = fs;
                trialET.trial = {data(1:2, :)};
                trialET.time = {t};
                [~, velTS] = computeEyeVelocity(trialET, vel_win);
                vel = velTS.trial{1}(3, :);
                t_vel = velTS.time{1};

                idx_base = (t_vel >= -1.5) & (t_vel <= -0.5);
                if ~any(idx_base)
                    continue
                end
                vel_base = mean(vel(idx_base), 'omitnan');
                if ~isfinite(vel_base) || vel_base <= 0
                    continue
                end

                vel_pct = 100 * (vel - vel_base) ./ vel_base;
                tc_interp = interp1(t_vel, vel_pct, t_plot, 'linear', NaN);
                if mean(isfinite(tc_interp(idx_viable))) < cfg.min_trial_coverage
                    continue
                end
                tc_trials(k, :) = tc_interp;
            end
            gaze_tc(s, c, :) = mean(tc_trials, 1, 'omitnan');
        end
    end

    diff_tc(task) = summarize_velocity_diff(gaze_tc, cfg, t_plot, idx_viable, ...
        win_sm, outlier_k_iqr, max_interp_gap_sec, min_subject_coverage, winsor_pct);

    figure('Position', [0 0 1512 982], 'Color', 'w');
    plot_velocity_diff_panel(t_plot, diff_tc(task), cfg.diff_label, fontSize);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_gaze_velocity_%s_diff.png', cfg.name)));
end

%% Combined figure (both tasks)
figure('Position', [0 0 1512 982], 'Color', 'w');
for task = 1:numel(diff_tc)
    subplot(2, 1, task);
    hold on
    plot_velocity_diff_panel(t_plot, diff_tc(task), diff_tc(task).label, fontSize, false);
    title(diff_tc(task).label, 'FontSize', fontSize - 4, 'FontWeight', 'normal');
end
saveas(gcf, fullfile(fig_dir, 'AOC_gaze_velocity_diff_combined.png'));

%% Local functions
function summary = summarize_velocity_diff(gaze_tc, cfg, t_plot, idx_viable, ...
    win_sm, outlier_k_iqr, max_interp_gap_sec, min_subject_coverage, winsor_pct)

X_hi = squeeze(gaze_tc(:, cfg.high_idx, :));
X_lo = squeeze(gaze_tc(:, cfg.low_idx, :));
X_diff = X_hi - X_lo;
X_diff(~isfinite(X_diff)) = NaN;

med_metric = median(X_diff(:, idx_viable), 2, 'omitnan');
[X_diff, ~] = exclude_outlier_trajectories(X_diff, med_metric, outlier_k_iqr);

max_interp_gap_smp = max(1, round(max_interp_gap_sec * 500));
for s = 1:size(X_diff, 1)
    X_diff(s, :) = fill_short_nan_gaps(X_diff(s, :), max_interp_gap_smp);
end

subj_cov = mean(isfinite(X_diff(:, idx_viable)), 2);
X_diff = X_diff(subj_cov >= min_subject_coverage, :);

if win_sm > 1
    X_diff = movmean(X_diff, win_sm, 2, 'omitnan');
end

[m, se] = winsorized_nanmean_se(X_diff, winsor_pct);
summary.name = cfg.name;
summary.label = cfg.diff_label;
summary.m = m;
summary.se = se;
summary.n = size(X_diff, 1);
fprintf('  %s: n = %d subjects in group average\n', cfg.name, summary.n);
end

function plot_velocity_diff_panel(t_plot, summary, title_str, fontSize, do_title)
if nargin < 5
    do_title = true;
end

hold on
eb = shadedErrorBar(t_plot, summary.m, summary.se, 'lineProps', {'-'}, 'transparent', true);
set(eb.mainLine, 'Color', [0.10 0.10 0.10], 'LineWidth', 2.5);
set(eb.patch, 'FaceColor', [0.10 0.10 0.10], 'FaceAlpha', 0.20);
set(eb.edge(1), 'Color', 'none');
set(eb.edge(2), 'Color', 'none');

yline(0, '--', 'Color', [0.5 0.5 0.5]);
xline(0, '--k');

m = summary.m;
se = summary.se;
yl = [min(m - se, [], 'omitnan'), max(m + se, [], 'omitnan')];
if any(~isfinite(yl))
    yl = [-5 5];
else
    pad = 0.1 * diff(yl);
    if pad == 0
        pad = 1;
    end
    yl = yl + [-pad pad];
end
ylim(yl);
hPatch = patch([1 2 2 1], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
    'EdgeColor', 'none');
uistack(hPatch, 'bottom');

xlabel('Time [s]');
ylabel('Gaze Velocity [%]');
xlim([-0.5 2]);
set(gca, 'FontSize', fontSize - 4);
box off
if do_title
    title(sprintf('%s (n = %d)', title_str, summary.n), 'FontSize', fontSize - 2, 'FontWeight', 'normal');
end
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
    err_msg = ME.message;
end

tmp_file = '';
try
    tmp_file = fullfile(tempdir, sprintf('AOC_tmp_%s.mat', char(java.util.UUID.randomUUID)));
    copyfile(et_file, tmp_file, 'f');
    D = load(tmp_file, 'dataETlong');
catch ME2
    err_msg = sprintf('direct: %s | tmp: %s', err_msg, ME2.message);
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
