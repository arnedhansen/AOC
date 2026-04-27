%% AOC Gaze Deviation - Sternberg conditions

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');

if ispc
    seb_path = 'W:\Students\Arne\toolboxes\shadedErrorBar';
else
    seb_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar';
end
addpath(seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'gaze', 'deviation');
if ~isfolder(fig_dir), mkdir(fig_dir); end
fig_dir_ctrl = fullfile(paths.figures, 'controls', 'gaze', 'deviation');
if ~isfolder(fig_dir_ctrl), mkdir(fig_dir_ctrl); end

fontSize = 25;
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
smooth_sec = 0.05;

cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};
et_fname = 'dataET_sternberg.mat';
blink_win = 50; % 100 ms at 500 Hz
min_trial_coverage = 0;
idx_viable = (t_plot >= 0) & (t_plot <= 2);

%% Aggregate subject data using feature-extraction-comparable preprocessing
merged_file = fullfile(feat_dir, 'AOC_merged_data_sternberg.mat');
S = load(merged_file, 'merged_data_sternberg');
T = struct2table(S.merged_data_sternberg);
uIDs = unique(T.ID);
nSubj = numel(uIDs);
gaze_tc = nan(nSubj, 3, Tf);

for s = 1:nSubj
    fprintf('Computing Subject %d / %d\n', s, nSubj);
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

            idx_base = (t >= -0.5) & (t <= -0.25);
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
            tc_trials(k, :) = tc_interp;
        end
        gaze_tc(s, c, :) = mean(tc_trials, 1, 'omitnan');
    end
end

%% Plot
win_sm = max(1, round(smooth_sec * fs));
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:3
    X = squeeze(gaze_tc(:, c, :));
    X(~isfinite(X)) = NaN;
    if win_sm > 1
        X = movmean(X, win_sm, 2, 'omitnan');
    end

    m = mean(X, 1, 'omitnan');
    n_fin = sum(isfinite(X), 1);
    se = std(X, 0, 1, 'omitnan') ./ max(sqrt(n_fin), 1);
    se(~isfinite(se)) = NaN;

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
ylim([-15 375]);
set(gca, 'FontSize', fontSize - 4);
box off
legend_handles = gobjects(1, 3);
for c = 1:3
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cond_labels, 'Location', 'northeast', 'FontSize', fontSize - 2, 'Box', 'off');

saveas(gcf, fullfile(fig_dir, 'AOC_gaze_dev_timecourse_sternberg_conditions.png'));

%% Debug plots: individual subject trajectories per condition
for c = 1:3
    X = squeeze(gaze_tc(:, c, :));
    X(~isfinite(X)) = NaN;
    if win_sm > 1
        X = movmean(X, win_sm, 2, 'omitnan');
    end

    figure('Position', [0 0 1512 982], 'Color', 'w');
    hold on

    % Plot each participant with low opacity (or lightened color fallback).
    for s = 1:size(X, 1)
        xs = X(s, :);
        if ~any(isfinite(xs))
            continue
        end
        h = plot(t_plot, xs, '-', 'LineWidth', 0.8, 'Color', colors(c, :));
        try
            h.Color(4) = 0.18;
        catch
            h.Color = 0.70 * colors(c, :) + 0.30 * [1 1 1];
        end
    end

    m = mean(X, 1, 'omitnan');
    plot(t_plot, m, '-', 'LineWidth', 3.0, 'Color', colors(c, :));

    yline(0, '--');
    xline(0, '--k');
    xlabel('Time [s]');
    ylabel('Gaze Deviation [%]');
    xlim([-0.5 2]);
    %ylim([-10 375]);
    title(sprintf('Sternberg %s: individual trajectories + mean', cond_labels{c}));
    set(gca, 'FontSize', fontSize - 4);
    box off

    fname = sprintf('AOC_gaze_dev_sternberg_individual_%s.png', regexprep(lower(cond_labels{c}), '\s+', '_'));
    saveas(gcf, fullfile(fig_dir_ctrl, fname));
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
