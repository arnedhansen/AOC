%% AOC Split Alpha Amp/Red (N-back load13 only)
% Dedicated N-back script (1-back vs 3-back, pooled within subject).
% Includes gaze deviation and EEG alpha time courses in one figure.

%% Setup
startup
[subjects, paths, colors] = setup('AOC'); %#ok<ASGLU>
pathAOC = paths.features;

addpath(paths.seb_path);

feat_dir = paths.features;
fig_dir = fullfile(paths.figures, 'splits', 'SplitAlphaAmpRed', 'Nback');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
fprintf('\n=== AOC Split Alpha Amp/Red (N-back load13 only) ===\n');
fprintf('Figure directory: %s\n', fig_dir);

% Canonical figure size.
fig_pos = [0 0 1512 982];

% Configuration: load13 only.
cond_vals = [1 2 3];
cond_codes = [21 22 23];
cond_labels = {'1-back', '2-back', '3-back'}; %#ok<NASGU>
gaze_fname = 'gaze_series_nback_trials.mat';
tfr_fname_primary = 'tfr_nback.mat';
tfr_fname_long = 'tfr_nback_long.mat';
tfr_vars = {'tfr1_fooof_bl', 'tfr2_fooof_bl', 'tfr3_fooof_bl'};
split_label = 'splitLoad1vs3';

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
uIDs = unique(T.ID);
nSubj = numel(uIDs);

%% Preallocate time courses
fs = 500;
t_full = -0.5:1/fs:3;
t_plot = t_full(2:end);
Tf = numel(t_plot);

dev_tc = nan(nSubj, 3, Tf);
eeg_tc = nan(nSubj, 3, Tf);
channels = {};

missing_tfr_long = {};
missing_gaze = {};

%% Aggregate subject data
fprintf('\n=== Aggregating N-back subject data (%d subjects) ===\n', nSubj);
for s = 1:nSubj
    clc; fprintf('[SPLIT ALPHA AMP - NBACK] Aggregating N-back data for Subject %d / %d \n', s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end

    %% EEG (long TFR required for 3 s coverage)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    tfr_file = fullfile(eeg_dir, tfr_fname_primary);
    if ~isfile(tfr_file)
        % Optional fallback for setups that store an explicit long file.
        tfr_file = fullfile(eeg_dir, tfr_fname_long);
    end
    if ~isfile(tfr_file)
        missing_tfr_long{end+1} = sid_str; %#ok<AGROW>
    else
        try
            R = load(tfr_file, tfr_vars{:});
            tfr_conds = {R.(tfr_vars{1}), R.(tfr_vars{2}), R.(tfr_vars{3})};
            if isempty(channels)
                channels = occ_channels_from_labels(tfr_conds{1}.label);
            end
            for c = 1:3
                eeg_tc(s, c, :) = extract_alpha_timecourse_one_tfr(tfr_conds{c}, channels, t_plot);
            end
        catch
            missing_tfr_long{end+1} = sid_str; %#ok<AGROW>
        end
    end

    %% Gaze deviation
    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', gaze_fname);
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
            dev_mat(k, :) = interp1(tt, dev, t_plot, 'linear', NaN);
        end
        dev_tc(s, c, :) = mean(dev_mat, 1, 'omitnan');
    end
end

if isempty(channels)
    error('No N-back EEG data could be loaded from %s (or fallback %s).', tfr_fname_primary, tfr_fname_long);
end

%% Gaze representation: percent baseline
t_vec = linspace(-0.5, 3, Tf);
bl_idx = (t_vec >= -0.5) & (t_vec <= -0.25);
dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
gaze_tc = (dev_tc ./ dev_bl_3d - 1) * 100;
gaze_tc(~isfinite(gaze_tc)) = NaN;

%% Subject completeness filter (0..3 s)
tc_window_idx = (t_vec >= 0) & (t_vec <= 3);
tc_complete_min_frac = 0.995;
tc_finite_frac = squeeze(mean(isfinite(dev_tc(:, :, tc_window_idx)), 3));
tc_has_endpoint = squeeze(isfinite(dev_tc(:, :, end)));
tc_complete_by_cond = (tc_finite_frac >= tc_complete_min_frac) & tc_has_endpoint;
tc_complete_subj = all(tc_complete_by_cond, 2);

gaze_tc(~tc_complete_subj, :, :) = NaN;
eeg_tc(~tc_complete_subj, :, :) = NaN;

fprintf('Time-course completeness filter: finite frac >= %.3f in [0,3]s + finite endpoint at 3s\n', tc_complete_min_frac);
fprintf('Included subjects after gaze completeness filter: %d / %d\n', sum(tc_complete_subj), nSubj);

%% Load13 pooled comparison (1-back vs 3-back) with EEG trace
cond_1_idx = find(cond_vals == 1, 1, 'first');
cond_3_idx = find(cond_vals == 3, 1, 'first');
if isempty(cond_1_idx) || isempty(cond_3_idx)
    error('cond_vals must contain 1 and 3.');
end

tc_1back = squeeze(gaze_tc(:, cond_1_idx, :));
tc_3back = squeeze(gaze_tc(:, cond_3_idx, :));
eeg_1back = squeeze(eeg_tc(:, cond_1_idx, :));
eeg_3back = squeeze(eeg_tc(:, cond_3_idx, :));

tc_1back = tc_1back(tc_complete_subj, :);
tc_3back = tc_3back(tc_complete_subj, :);
eeg_1back = eeg_1back(tc_complete_subj, :);
eeg_3back = eeg_3back(tc_complete_subj, :);

save_tag = sprintf('nback_%s_gaze_eeg_1back_vs_3back', split_label);
plot_timecourse_load_compare_with_eeg(tc_1back, tc_3back, eeg_1back, eeg_3back, ...
    colors(1, :), colors(3, :), '1-back', '3-back', 'Gaze Deviation [%]', ...
    'Alpha Power [dB]', save_tag, fig_dir, fig_pos, 25, fs, 0.05);

%% Diagnostics
fprintf('\n=== Data Diagnostics [nback load13] ===\n');
fprintf('Missing EEG TFR (%s; fallback %s): %d\n', tfr_fname_primary, tfr_fname_long, numel(unique(missing_tfr_long)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

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
        ch{end+1} = lab; %#ok<AGROW>
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

function alpha_interp = extract_alpha_timecourse_one_tfr(T, channels, t_target)
alpha_interp = nan(1, numel(t_target));
if ~isfield(T, 'powspctrm') || ~isfield(T, 'label') || ~isfield(T, 'freq') || ~isfield(T, 'time')
    return
end
ch_idx = find(ismember(T.label, channels));
freq_idx = (T.freq >= 8) & (T.freq <= 14);
if isempty(ch_idx) || ~any(freq_idx)
    return
end
alpha_raw = squeeze(mean(mean(T.powspctrm(ch_idx, freq_idx, :), 1, 'omitnan'), 2, 'omitnan'));
if numel(alpha_raw) ~= numel(T.time)
    alpha_raw = alpha_raw(:);
end
alpha_interp = interp1(T.time(:), alpha_raw(:), t_target(:), 'linear', NaN)';
end

function plot_timecourse_load_compare_with_eeg(tc_a, tc_b, eeg_a, eeg_b, col_a, col_b, label_a, label_b, ylab_gaze, ylab_eeg, save_tag, fig_dir, fig_pos, fsz, fs, smooth_sec)
if nargin < 16 || isempty(smooth_sec)
    smooth_sec = 0.05;
end
tc_a(~isfinite(tc_a)) = NaN;
tc_b(~isfinite(tc_b)) = NaN;
eeg_a(~isfinite(eeg_a)) = NaN;
eeg_b(~isfinite(eeg_b)) = NaN;

nT = size(tc_a, 2);
t_plot = linspace(-0.5 + 1/fs, 3, nT);
win_sm = max(1, round(smooth_sec * fs));
if win_sm > 1
    tc_a = movmean(tc_a, win_sm, 2, 'omitnan');
    tc_b = movmean(tc_b, win_sm, 2, 'omitnan');
    eeg_a = movmean(eeg_a, win_sm, 2, 'omitnan');
    eeg_b = movmean(eeg_b, win_sm, 2, 'omitnan');
end

[mG_a, sG_a] = mean_sem(tc_a);
[mG_b, sG_b] = mean_sem(tc_b);
[mE_a, sE_a] = mean_sem(eeg_a);
[mE_b, sE_b] = mean_sem(eeg_b);

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile; hold on
e1 = shadedErrorBar(t_plot, mG_a, sG_a, 'lineProps', {'-'}, 'transparent', true);
e2 = shadedErrorBar(t_plot, mG_b, sG_b, 'lineProps', {'-'}, 'transparent', true);
set(e1.mainLine, 'Color', col_a, 'LineWidth', 2.5);
set(e2.mainLine, 'Color', col_b, 'LineWidth', 2.5);
set(e1.patch, 'FaceColor', col_a, 'FaceAlpha', 0.20);
set(e2.patch, 'FaceColor', col_b, 'FaceAlpha', 0.20);
set(e1.edge(1), 'Color', 'none'); set(e1.edge(2), 'Color', 'none');
set(e2.edge(1), 'Color', 'none'); set(e2.edge(2), 'Color', 'none');
xline(0, '--k');
xlim([-0.5 3]);
ylabel(ylab_gaze);
title(sprintf('%s vs %s (pooled over subjects)', label_a, label_b), 'Interpreter', 'none');
legend([e1.mainLine, e2.mainLine], ...
    {sprintf('%s gaze (n=%d)', label_a, size(tc_a, 1)), sprintf('%s gaze (n=%d)', label_b, size(tc_b, 1))}, ...
    'Location', 'northeast', 'FontSize', fsz-5, 'Box', 'off');
set(gca, 'FontSize', fsz-4);
box on

nexttile; hold on
e3 = shadedErrorBar(t_plot, mE_a, sE_a, 'lineProps', {'-'}, 'transparent', true);
e4 = shadedErrorBar(t_plot, mE_b, sE_b, 'lineProps', {'-'}, 'transparent', true);
set(e3.mainLine, 'Color', col_a, 'LineWidth', 2.5);
set(e4.mainLine, 'Color', col_b, 'LineWidth', 2.5);
set(e3.patch, 'FaceColor', col_a, 'FaceAlpha', 0.20);
set(e4.patch, 'FaceColor', col_b, 'FaceAlpha', 0.20);
set(e3.edge(1), 'Color', 'none'); set(e3.edge(2), 'Color', 'none');
set(e4.edge(1), 'Color', 'none'); set(e4.edge(2), 'Color', 'none');
xline(0, '--k');
xlim([-0.5 3]);
xlabel('Time [s]');
ylabel(ylab_eeg);
legend([e3.mainLine, e4.mainLine], ...
    {sprintf('%s EEG', label_a), sprintf('%s EEG', label_b)}, ...
    'Location', 'northeast', 'FontSize', fsz-5, 'Box', 'off');
set(gca, 'FontSize', fsz-4);
box on

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlphaAmpRed_timecourse_%s.png', save_tag)));
close(gcf);
end

function [m, sem] = mean_sem(X)
m = mean(X, 1, 'omitnan');
n_fin = sum(isfinite(X), 1);
sem = std(X, 0, 1, 'omitnan') ./ max(sqrt(n_fin), 1);
sem(~isfinite(sem)) = NaN;
end
