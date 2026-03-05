%% AOC Split Sternberg Alpha (Subject-Level)
% Subject-level split (fixed across conditions) using merged_data_sternberg:
%   mean AlphaPower_FOOOF_bl across WM2/4/6 (baselined, FOOOFed alpha)
%   < 0  -> reduction group
%   > 0  -> amplification group
%   == 0 -> excluded
%
% Uses baselined+FOOOFed alpha (power spectra, TFR, topoplots) and baselined gaze metrics (*FullBL).
%
% Outliers excluded via Tukey 1.5*IQR (per metric per condition) before visualization/analyses.
%
% Generates:
% - Power spectra (3 conditions) for both groups
% - TFRs per condition for both groups
% - Topoplots per condition for both groups + group-difference topoplots
% - Rainclouds for alpha + gaze metrics
% - Within-subject trajectories over conditions
% - Time-course panels (SPL, velocity, gaze deviation) with effect-size strips
% - Correlation panels within each group

%% Setup
startup
[subjects, pathAOC, colors, headmodel] = setup('AOC');

if ispc
    base_data = 'W:\Students\Arne\AOC';
    seb_path = 'W:\Students\Arne\toolboxes\shadedErrorBar';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
    seb_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/shadedErrorBar';
end
addpath(seb_path);

feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'tests', 'splitAlpha');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

% Keep canonical figure size requested.
fig_pos = [0 0 1512 982];

% Conditions and labels
cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

fontSize = 20;
rng(42);
alpha_zero_pct = 5; % Exclude near-zero alpha within this % of max |alpha|

%% Load subject-level merged data and define alpha split
merged_file = fullfile(feat_dir, 'merged_data_sternberg.mat');
if ~isfile(merged_file)
    error('Missing file: %s', merged_file);
end
S = load(merged_file, 'merged_data_sternberg');
if ~isfield(S, 'merged_data_sternberg')
    error('Variable merged_data_sternberg not found in %s', merged_file);
end
T = struct2table(S.merged_data_sternberg);

% Compute subject-level split value from full-window FOOOF-baselined alpha.
uIDs = unique(T.ID);
nSubj = numel(uIDs);
alpha_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    alpha_mean(i) = mean(T.AlphaPower_FOOOF_bl(mask), 'omitnan');
end

% Percentage-based exclusion band around zero.
valid_alpha_mean = alpha_mean(isfinite(alpha_mean));
if isempty(valid_alpha_mean)
    error('No finite subject-level alpha values found for split.');
end
alpha_abs_ref = max(abs(valid_alpha_mean));
alpha_zero_margin = (alpha_zero_pct / 100) * alpha_abs_ref;

reduction_ids = uIDs(alpha_mean < -alpha_zero_margin);
amplification_ids = uIDs(alpha_mean > alpha_zero_margin);
zero_ids = uIDs(abs(alpha_mean) <= alpha_zero_margin);

fprintf('\n=== Split Summary (AlphaPower_FOOOF_bl, full window) ===\n');
fprintf('Subjects total: %d\n', nSubj);
fprintf('Zero-band setting: %.2f%% of max |alpha| (max |alpha|=%.4f, cutoff=%.4f)\n', ...
    alpha_zero_pct, alpha_abs_ref, alpha_zero_margin);
fprintf('Reduction (< -%.4f): %d\n', alpha_zero_margin, numel(reduction_ids));
fprintf('Amplification (> %.4f): %d\n', alpha_zero_margin, numel(amplification_ids));
fprintf('Excluded (|alpha| <= %.4f): %d\n', alpha_zero_margin, numel(zero_ids));

if numel(reduction_ids) < 2 || numel(amplification_ids) < 2
    error('Insufficient subjects per split group.');
end

%% Preallocate containers
% EEG power spectra/topography
pow_red = cell(1, 3);
pow_amp = cell(1, 3);
for c = 1:3
    pow_red{c} = {};
    pow_amp{c} = {};
end

% EEG TFR
tfr_red = cell(1, 3);
tfr_amp = cell(1, 3);
for c = 1:3
    tfr_red{c} = {};
    tfr_amp{c} = {};
end

% Gaze summary metrics from merged subject table
metrics = struct();
metrics.Alpha = nan(nSubj, 3);
metrics.SPL = nan(nSubj, 3);           % ScanPathLengthFullBL [% change]
metrics.Dev = nan(nSubj, 3);          % GazeDeviationFullBL [% change]
metrics.BCEA = nan(nSubj, 3);         % BCEAFullBL [% change]
metrics.MS = nan(nSubj, 3);           % MSRateFullBL [% change]
metrics.Vel = nan(nSubj, 3);          % velocity [% change], computed below

% Time courses per subject x condition (full resolution)
fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
spl_tc = nan(nSubj, 3, Tf);
dev_tc = nan(nSubj, 3, Tf);
vel_tc = nan(nSubj, 3, Tf);

missing_eeg = {};
missing_tfr = {};
missing_gaze = {};

%% Aggregate per-subject data
for s = 1:nSubj
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end
    subj_rows = T(T.ID == sid, :);

    % Subject-level condition metrics from merged_data_sternberg
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics.Alpha(s, c) = mean(subj_rows.AlphaPower_FOOOF_bl(cmask), 'omitnan');
            metrics.SPL(s, c) = mean(subj_rows.ScanPathLengthFullBL(cmask), 'omitnan');
            metrics.Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
            metrics.BCEA(s, c) = mean(subj_rows.BCEAFullBL(cmask), 'omitnan');
            metrics.MS(s, c) = mean(subj_rows.MSRateFullBL(cmask), 'omitnan');
        end
    end

    % EEG power spectra and topography sources (baselined, FOOOFed)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    try
        P = load(fullfile(eeg_dir, 'power_stern_fooof.mat'), 'pow2_fooof_bl', 'pow4_fooof_bl', 'pow6_fooof_bl');
        pow_conds = {P.pow2_fooof_bl, P.pow4_fooof_bl, P.pow6_fooof_bl};
        for c = 1:3
            if ismember(sid, reduction_ids)
                pow_red{c}{end+1} = pow_conds{c}; %#ok<AGROW>
            elseif ismember(sid, amplification_ids)
                pow_amp{c}{end+1} = pow_conds{c}; %#ok<AGROW>
            end
        end
    catch
        missing_eeg{end+1} = sid_str; %#ok<AGROW>
    end

    % EEG TFR sources
    try
        % Use FOOOFed + baselined TFR to avoid raw-power scale skew.
        R = load(fullfile(eeg_dir, 'tfr_stern.mat'), 'tfr2_fooof_bl', 'tfr4_fooof_bl', 'tfr6_fooof_bl');
        tfr_conds = {R.tfr2_fooof_bl, R.tfr4_fooof_bl, R.tfr6_fooof_bl};
        for c = 1:3
            if ismember(sid, reduction_ids)
                tfr_red{c}{end+1} = tfr_conds{c}; %#ok<AGROW>
            elseif ismember(sid, amplification_ids)
                tfr_amp{c}{end+1} = tfr_conds{c}; %#ok<AGROW>
            end
        end
    catch
        missing_tfr{end+1} = sid_str; %#ok<AGROW>
    end

    % Gaze series for time courses and velocity summary
    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', 'gaze_series_sternberg_trials.mat');
    if ~isfile(gaze_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    G = load(gaze_file);
    if ~isfield(G, 'trialinfo')
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    conds = parse_trialinfo_conds(G.trialinfo);
    if isempty(conds)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    % SPL time course from ScanPathSeries
    if isfield(G, 'ScanPathSeries') && isfield(G, 'ScanPathSeriesT')
        for c = 1:3
            tr_mask = conds == cond_codes(c);
            tr_idx = find(tr_mask);
            if isempty(tr_idx)
                continue
            end
            mat = nan(numel(tr_idx), Tf);
            for k = 1:numel(tr_idx)
                tr = tr_idx(k);
                srl = G.ScanPathSeries{tr};
                tt = G.ScanPathSeriesT{tr};
                if isempty(srl) || isempty(tt)
                    continue
                end
                try
                    mat(k, :) = interp1(tt, srl, t_plot, 'linear', NaN);
                catch
                end
            end
            spl_tc(s, c, :) = nanmean(mat, 1);
        end
    end

    % Deviation and velocity time course from gaze x/y if available
    has_xy = isfield(G, 'gaze_x') && isfield(G, 'gaze_y');
    if has_xy
        for c = 1:3
            tr_mask = conds == cond_codes(c);
            tr_idx = find(tr_mask);
            if isempty(tr_idx)
                continue
            end

            dev_mat = nan(numel(tr_idx), Tf);
            vel_mat = nan(numel(tr_idx), Tf);
            vel_full_trials = nan(numel(tr_idx), 1);
            vel_bl_trials = nan(numel(tr_idx), 1);

            for k = 1:numel(tr_idx)
                tr = tr_idx(k);
                x = double(G.gaze_x{tr});
                y = double(G.gaze_y{tr});
                if isempty(x) || isempty(y) || numel(x) ~= numel(y)
                    continue
                end

                tt = linspace(-0.5, 2, numel(x));
                dev = sqrt((x - 400).^2 + (y - 300).^2);
                [vx, vy] = compute_velocity_sg(x, y, fs, 3);
                [vx, vy] = clean_velocity_components(vx, vy);
                vel = hypot(vx, vy);

                try
                    dev_mat(k, :) = interp1(tt, dev, t_plot, 'linear', NaN);
                    vel_mat(k, :) = interp1(tt, vel, t_plot, 'linear', NaN);
                catch
                end

                idx_full = tt >= 0 & tt <= 2;
                idx_bl = tt >= -0.5 & tt <= -0.25;
                v_full = mean(vel(idx_full), 'omitnan');
                v_bl = mean(vel(idx_bl), 'omitnan');
                vel_full_trials(k) = v_full;
                if isfinite(v_full) && isfinite(v_bl) && v_bl > 0
                    vel_bl_trials(k) = 100 * (v_full - v_bl) / v_bl;
                end
            end

            dev_tc(s, c, :) = nanmean(dev_mat, 1);
            vel_tc(s, c, :) = nanmean(vel_mat, 1);
            metrics.Vel(s, c) = mean(vel_bl_trials, 'omitnan');
        end
    end
end

%% Define groups in row index space
is_red = ismember(uIDs, reduction_ids);
is_amp = ismember(uIDs, amplification_ids);

%% Exclude outliers (Tukey 1.5*IQR) before visualization and analyses
fprintf('\n=== Outlier exclusion (Tukey 1.5*IQR, per metric per condition) ===\n');
[metrics, spl_tc, dev_tc, vel_tc] = exclude_outliers_tukey(metrics, spl_tc, dev_tc, vel_tc);

%% Determine occipital channels (from first available power file)
channels = {};
for c = 1:3
    if ~isempty(pow_red{c})
        p0 = pow_red{c}{1};
        channels = occ_channels_from_labels(p0.label);
        break
    elseif ~isempty(pow_amp{c})
        p0 = pow_amp{c}{1};
        channels = occ_channels_from_labels(p0.label);
        break
    end
end
if isempty(channels)
    error('No power data available for channel definition.');
end

%% -------- Power spectra (both groups) --------
plot_group_power_spectrum(pow_red, channels, colors, cond_labels, ...
    'Reduction Group Power Spectrum', ...
    fullfile(fig_dir, 'AOC_splitAlpha_powspctrm_reduction.png'), fig_pos, fontSize);

plot_group_power_spectrum(pow_amp, channels, colors, cond_labels, ...
    'Amplification Group Power Spectrum', ...
    fullfile(fig_dir, 'AOC_splitAlpha_powspctrm_amplification.png'), fig_pos, fontSize);

%% -------- TFRs per condition (both groups + diff, 3x3) --------
color_map_tfr = interp1(linspace(0,1,5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0,1,64));

plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, ...
    color_map_tfr, fig_dir, fig_pos, fontSize);

%% -------- Topoplots per condition (both groups + differences) --------
plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir, fig_pos, fontSize);

%% -------- Rainclouds --------
plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% -------- Within-subject trajectories --------
plot_metric_trajectories(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% -------- Time courses with effect-size strips --------
plot_timecourse_with_effect(spl_tc, is_red, is_amp, cond_labels, colors, ...
    'Scan Path Length [px]', 'spl', fig_dir, fig_pos, fontSize);
plot_timecourse_with_effect(vel_tc, is_red, is_amp, cond_labels, colors, ...
    'Eye Velocity [px/s]', 'velocity', fig_dir, fig_pos, fontSize);
plot_timecourse_with_effect(dev_tc, is_red, is_amp, cond_labels, colors, ...
    'Gaze Deviation [px]', 'gaze_deviation', fig_dir, fig_pos, fontSize);

%% -------- Correlation panels --------
plot_correlation_panels(metrics, is_red, is_amp, fig_dir, fig_pos, fontSize);

%% -------- Sanity checks --------
fprintf('\n=== Data Diagnostics ===\n');
fprintf('Missing EEG power (power_stern_fooof.mat): %d\n', numel(unique(missing_eeg)));
fprintf('Missing TFR files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

%% ========================= Local Functions =========================
function [metrics_out, spl_tc_out, dev_tc_out, vel_tc_out] = exclude_outliers_tukey(metrics_in, spl_tc_in, dev_tc_in, vel_tc_in)
% Apply Tukey 1.5*IQR rule per metric per condition. Set outliers to NaN.
% Also excludes corresponding time-course data for SPL, Dev, Vel.
metrics_out = metrics_in;
spl_tc_out = spl_tc_in;
dev_tc_out = dev_tc_in;
vel_tc_out = vel_tc_in;

metric_fields = {'Alpha', 'SPL', 'Vel', 'Dev', 'BCEA', 'MS'};
tc_map = containers.Map({'SPL', 'Dev', 'Vel'}, {1, 2, 3});  % index for spl/dev/vel
tc_cells = {spl_tc_out, dev_tc_out, vel_tc_out};

for m = 1:numel(metric_fields)
    key = metric_fields{m};
    X = metrics_out.(key);
    nOut = 0;
    for c = 1:3
        vals = X(:, c);
        vals_f = vals(isfinite(vals));
        if numel(vals_f) < 4
            continue
        end
        q1 = prctile(vals_f, 25);
        q3 = prctile(vals_f, 75);
        iqr_val = q3 - q1;
        if iqr_val <= 0
            continue
        end
        lo = q1 - 1.5 * iqr_val;
        hi = q3 + 1.5 * iqr_val;
        out_mask = (vals < lo) | (vals > hi);
        nOut = nOut + sum(out_mask);
        X(out_mask, c) = NaN;
        % Exclude corresponding time course for gaze metrics
        if isKey(tc_map, key)
            tc_idx = tc_map(key);
            tc = tc_cells{tc_idx};
            tc(out_mask, c, :) = NaN;
            tc_cells{tc_idx} = tc;
        end
    end
    metrics_out.(key) = X;
    if nOut > 0
        fprintf('  %s: %d outlier(s) excluded\n', key, nOut);
    end
end
spl_tc_out = tc_cells{1};
dev_tc_out = tc_cells{2};
vel_tc_out = tc_cells{3};
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

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'}) || contains(lab, {'PO'})
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

function [vx, vy] = clean_velocity_components(vx, vy)
halfwin = 11; % SG edge region for framelen=21
if numel(vx) > 2*halfwin
    vx(1:halfwin) = NaN;
    vx(end-halfwin+1:end) = NaN;
    vy(1:halfwin) = NaN;
    vy(end-halfwin+1:end) = NaN;
end
zvx = (vx - nanmean(vx)) ./ (nanstd(vx) + eps);
zvy = (vy - nanmean(vy)) ./ (nanstd(vy) + eps);
bad = abs(zvx) > 4 | abs(zvy) > 4;
vx(bad) = NaN;
vy(bad) = NaN;
end

function plot_group_power_spectrum(pow_cells, channels, colors, cond_labels, ttl, out_file, fig_pos, fsz)
if any(cellfun(@isempty, pow_cells))
    warning('Skipping %s (incomplete power condition data).', ttl);
    return
end
ga = cell(1, 3);
for c = 1:3
    ga{c} = ft_freqgrandaverage([], pow_cells{c}{:});
end

figure('Position', fig_pos, 'Color', 'w');
hold on

elecs = ismember(ga{1}.label, channels);
freqs = ga{1}.freq;

for c = 1:3
    subj_n = numel(pow_cells{c});
    subj_spec = nan(subj_n, numel(freqs));
    for s = 1:subj_n
        p = pow_cells{c}{s};
        subj_spec(s, :) = mean(p.powspctrm(elecs, :), 1, 'omitnan');
    end
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(subj_spec), 1));
    eb(c) = shadedErrorBar(freqs, m, se, 'lineProps', {'-'}, 'transparent', true); %#ok<AGROW>
    set(eb(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
    set(eb(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
    set(eb(c).edge(1), 'Color', colors(c, :));
    set(eb(c).edge(2), 'Color', colors(c, :));
end

xlim([5 20]);
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title(ttl, 'FontSize', fsz + 4);
legend([eb.mainLine], cond_labels, 'FontSize', fsz - 2, 'Location', 'best');
set(gca, 'FontSize', fsz);
box on
saveas(gcf, out_file);
close(gcf);
end

function plot_group_tfrs_all(tfr_red, tfr_amp, channels, cond_labels, cond_vals, headmodel, color_map, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, tfr_red)) || any(cellfun(@isempty, tfr_amp))
    warning('Skipping TFR plot (incomplete condition data).');
    return
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], tfr_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], tfr_amp{c}{:});
end

% Shared clim: absolute for Reduction/Amplification, separate for Diff
[~, ch_idx] = ismember(channels, ga_red{1}.label);
freq_idx = ga_red{1}.freq >= 5 & ga_red{1}.freq <= 30;
time_idx = ga_red{1}.time >= -0.5 & ga_red{1}.time <= 2;
mx = 0;
all_diff_vals = [];
for c = 1:3
    Ared = squeeze(mean(ga_red{c}.powspctrm(ch_idx, :, :), 1, 'omitnan'));
    Aamp = squeeze(mean(ga_amp{c}.powspctrm(ch_idx, :, :), 1, 'omitnan'));
    mx = max(mx, max(abs(Ared(freq_idx, time_idx)), [], 'all'));
    mx = max(mx, max(abs(Aamp(freq_idx, time_idx)), [], 'all'));
    Adiff = Aamp - Ared;
    tmp = Adiff(freq_idx, time_idx);
    all_diff_vals = [all_diff_vals; tmp(:)]; %#ok<AGROW>
end
clim_abs = [-mx mx];
all_diff_vals = all_diff_vals(isfinite(all_diff_vals));
if isempty(all_diff_vals)
    clim_diff = [-mx mx];
else
    mxdiff = prctile(abs(all_diff_vals), 99);
    clim_diff = [-max(mxdiff, eps) max(mxdiff, eps)];
end

cfg = [];
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;

% Single 3x3 figure: row 1 = Reduction, row 2 = Amplification, row 3 = Diff
figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    % Row 1: Reduction
    ax = subplot(3, 3, c);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_red{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    colorbar(ax);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'Frequency [Hz]');
    rectangle(ax, 'Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Reduction - %s', cond_labels{c}), 'FontSize', fsz + 2);

    % Row 2: Amplification
    ax = subplot(3, 3, 3 + c);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, ga_amp{c});
    colormap(ax, color_map);
    set(ax, 'CLim', clim_abs);
    colorbar(ax);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'Frequency [Hz]');
    rectangle(ax, 'Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Amplification - %s', cond_labels{c}), 'FontSize', fsz + 2);

    % Row 3: Difference (Amplification - Reduction)
    gd = ga_amp{c};
    gd.powspctrm = ga_amp{c}.powspctrm - ga_red{c}.powspctrm;
    ax = subplot(3, 3, 6 + c);
    cfg.figure = ax;
    ft_singleplotTFR(cfg, gd);
    colormap(ax, rdbu_cmap(64));
    set(ax, 'CLim', clim_diff);
    colorbar(ax);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'Frequency [Hz]');
    rectangle(ax, 'Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Amp - Red - %s', cond_labels{c}), 'FontSize', fsz + 2);
end
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlpha_tfr_all.png'));
close(gcf);
end

function plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, cond_vals, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping topoplots (missing power data).');
    return
end

ga_red = cell(1, 3);
ga_amp = cell(1, 3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], pow_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], pow_amp{c}{:});
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

cmap_abs = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
cfg.colormap = cmap_abs;
freq_idx = ga_red{1}.freq >= 8 & ga_red{1}.freq <= 14;
all_alpha = [];
all_diff = [];
for c = 1:3
    Ared = mean(ga_red{c}.powspctrm(:, freq_idx), 2, 'omitnan');
    Aamp = mean(ga_amp{c}.powspctrm(:, freq_idx), 2, 'omitnan');
    all_alpha = [all_alpha; Ared(:); Aamp(:)]; %#ok<AGROW>
    all_diff = [all_diff; (Aamp(:) - Ared(:))]; %#ok<AGROW>
end
all_alpha = all_alpha(isfinite(all_alpha));
if isempty(all_alpha)
    cfg.zlim = 'maxabs';
else
    zlo = prctile(all_alpha, 1);
    zhi = prctile(all_alpha, 99);
    if zlo == zhi
        zlo = min(all_alpha);
        zhi = max(all_alpha);
    end
    cfg.zlim = [zlo zhi];
end

% Difference map cfg (prepare before loop)
cfgd = cfg;
cmap_diff = rdbu_cmap(100);  % RdBu diverging, no cbrewer dependency
cmap_diff = flipud(max(min(cmap_diff, 1), 0));
cfgd.colormap = cmap_diff;
all_diff = all_diff(isfinite(all_diff));
if isempty(all_diff)
    cfgd.zlim = 'maxabs';
else
    mx = prctile(abs(all_diff), 99);
    cfgd.zlim = [-mx mx];
end

% Single 3x3 figure: row 1 = Reduction, row 2 = Amplification, row 3 = Diff
figure('Position', fig_pos, 'Color', 'w');
for c = 1:3
    % Row 1: Reduction
    ax = subplot(3, 3, c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_red{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Reduction - %s', cond_labels{c}));

    % Row 2: Amplification
    ax = subplot(3, 3, 3 + c);
    cfg.figure = ax;
    ft_topoplotER(cfg, ga_amp{c});
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Amplification - %s', cond_labels{c}));

    % Row 3: Difference
    gd = ga_amp{c};
    gd.powspctrm = ga_amp{c}.powspctrm - ga_red{c}.powspctrm;
    ax = subplot(3, 3, 6 + c);
    cfgd.figure = ax;
    ft_topoplotER(cfgd, gd);
    colorbar(ax);
    set(ax, 'FontSize', fsz);
    title(ax, sprintf('Amp - Red - %s', cond_labels{c}));
end
saveas(gcf, fullfile(fig_dir, 'AOC_splitAlpha_topo_all.png'));
close(gcf);
end

function plot_metric_rainclouds(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fsz)
metric_defs = { ...
    'Alpha', 'AlphaPower_FOOOF_bl', false; ...
    'SPL', 'Scan path length', true; ...
    'Vel', 'Velocity', true; ...
    'Dev', 'Gaze deviation', true; ...
    'BCEA', 'BCEA', true; ...
    'MS', 'Microsaccade rate', true; ...
    };

for m = 1:size(metric_defs, 1)
    key = metric_defs{m, 1};
    varname = metric_defs{m, 2};
    is_pct = metric_defs{m, 3};
    X = metrics.(key);

    % Shared ymax for both groups, centered around 0
    all_vals = X(isfinite(X));
    if isempty(all_vals)
        ymax = 1;
    else
        ymax = max(abs(all_vals));
        if ymax <= 0
            ymax = 1;
        end
    end
    ylim_shared = [-ymax ymax];

    if is_pct
        ylab = [varname ' [%]'];
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
    title('Reduction');
    ylabel(ylab);
    box off

    nexttile; hold on
    for c = 1:3
        draw_one_cloud(X(is_amp, c), c, colors(c,:), 0.3, 96, 0.45);
    end
    yline(0, '--', 'Color', [0.4 0.4 0.4]);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title('Amplification');
    ylabel(ylab);
    box off

    sgtitle(varname, 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_raincloud_%s.png', lower(key))));
    close(gcf);
end
end

function plot_metric_trajectories(metrics, is_red, is_amp, cond_labels, colors, fig_dir, fig_pos, fsz)
metric_defs = {'Alpha','SPL','Vel','Dev','BCEA','MS'};
ylab_defs = {'AlphaPower_FOOOF_bl', 'Scan path length [% change]', 'Velocity [% change]', ...
    'Gaze deviation [% change]', 'BCEA [% change]', 'Microsaccade rate [% change]'};
for m = 1:numel(metric_defs)
    key = metric_defs{m};
    ylab = ylab_defs{m};
    X = metrics.(key);
    % Shared ymax for both groups, centered around 0
    all_vals = X(isfinite(X));
    if isempty(all_vals)
        ymax = 1;
    else
        ymax = max(abs(all_vals));
        if ymax <= 0
            ymax = 1;
        end
    end
    ylim_shared = [-ymax ymax];

    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    nexttile; hold on
    XR = X(is_red, :);
    for i = 1:size(XR,1)
        plot(1:3, XR(i,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
    hSub = plot(1:3, XR(1,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hMean = plot(1:3, mean(XR,1,'omitnan'), '-o', 'Color', colors(1,:), 'LineWidth', 3, 'MarkerFaceColor', colors(1,:));
    yline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title('Reduction');
    ylabel(ylab);
    box on
    legend([hSub hMean], {'Subjects', 'Group mean'}, 'Location', 'best', 'FontSize', fsz-6);

    nexttile; hold on
    XA = X(is_amp, :);
    for i = 1:size(XA,1)
        plot(1:3, XA(i,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
    hSub = plot(1:3, XA(1,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hMean = plot(1:3, mean(XA,1,'omitnan'), '-o', 'Color', colors(3,:), 'LineWidth', 3, 'MarkerFaceColor', colors(3,:));
    yline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
    ylim(ylim_shared);
    set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
    title('Amplification');
    ylabel(ylab);
    box on
    legend([hSub hMean], {'Subjects', 'Group mean'}, 'Location', 'best', 'FontSize', fsz-6);

    sgtitle(sprintf('Within-Subject Trajectories: %s', key), 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trajectory_%s.png', lower(key))));
    close(gcf);
end
end

function plot_timecourse_with_effect(tc, is_red, is_amp, cond_labels, colors, ylab, save_tag, fig_dir, fig_pos, fsz)
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(4, 1, 'TileSpacing', 'compact');
t = linspace(-0.5, 2, size(tc, 3));

for c = 1:3
    nexttile(c); hold on
    R = squeeze(tc(is_red, c, :));
    A = squeeze(tc(is_amp, c, :));
    mR = mean(R, 1, 'omitnan');
    mA = mean(A, 1, 'omitnan');
    sR = std(R, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(R), 1));
    sA = std(A, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(A), 1));

    e1 = shadedErrorBar(t, mR, sR, 'lineProps', {'-'}, 'transparent', true);
    e2 = shadedErrorBar(t, mA, sA, 'lineProps', {'-'}, 'transparent', true);
    set(e1.mainLine, 'Color', colors(1,:), 'LineWidth', 2.5);
    set(e2.mainLine, 'Color', colors(3,:), 'LineWidth', 2.5);
    set(e1.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.25);
    set(e2.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.25);
    xline(0, '--k');
    ylabel(ylab);
    title(cond_labels{c});
    xlim([-0.5 2]);
    box on
    set(gca, 'FontSize', fsz-4);
    legend([e1.mainLine e2.mainLine], {'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fsz-7);
end

nexttile(4); hold on
% Effect-size strip averaged across conditions.
Rall = squeeze(mean(tc(is_red, :, :), 2, 'omitnan'));
Aall = squeeze(mean(tc(is_amp, :, :), 2, 'omitnan'));
d = nan(1, size(tc, 3));
pval = nan(1, size(tc, 3));
for i = 1:size(tc, 3)
    x = Rall(:, i);
    y = Aall(:, i);
    x = x(isfinite(x));
    y = y(isfinite(y));
    if numel(x) >= 3 && numel(y) >= 3
        sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / (numel(x)+numel(y)-2));
        d(i) = (mean(y) - mean(x)) / max(sp, eps);
        [~, pval(i)] = ttest2(x, y);
    end
end
% Grey transparent boxes for significant timepoints (contiguous runs)
dt = t(2) - t(1);
sig = pval < 0.05 & isfinite(d);
% Find contiguous runs of significant samples
run_start = [false, diff(sig) == 1];
run_end = [diff(sig) == -1, false];
if sig(1), run_start(1) = true; end
if sig(end), run_end(end) = true; end
starts = find(run_start);
ends = find(run_end);
ylims = [min(d(isfinite(d)), [], 'omitnan') - 0.1, max(d(isfinite(d)), [], 'omitnan') + 0.1];
if any(~isfinite(ylims)) || diff(ylims) < 0.1
    ylims = [-1 1];
end
for k = 1:numel(starts)
    t1 = t(starts(k)) - dt/2;
    t2 = t(ends(k)) + dt/2;
    patch([t1 t2 t2 t1], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.5 0.5 0.5], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none');
end
ylim(ylims);
plot(t, d, 'k-', 'LineWidth', 2.5);
yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Cohen''s d');
xlim([-0.5 2]);
box on
set(gca, 'FontSize', fsz-4);

saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_timecourse_%s.png', save_tag)));
close(gcf);
end

function plot_correlation_panels(metrics, is_red, is_amp, fig_dir, fig_pos, fsz)
% Use subject-mean across conditions.
A = mean(metrics.Alpha, 2, 'omitnan');
SPL = mean(metrics.SPL, 2, 'omitnan');
VEL = mean(metrics.Vel, 2, 'omitnan');
DEV = mean(metrics.Dev, 2, 'omitnan');
BCEA = mean(metrics.BCEA, 2, 'omitnan');
MS = mean(metrics.MS, 2, 'omitnan');

Xset = {SPL, VEL, DEV, BCEA, MS};
names = {'SPL [% change]', 'Velocity [% change]', 'Gaze deviation [% change]', 'BCEA [% change]', 'Microsaccades [% change]'};

groups = {'reduction', 'amplification'};
masks = {is_red, is_amp};

for g = 1:2
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(2, 3, 'TileSpacing', 'compact');
    for i = 1:numel(Xset)
        nexttile; hold on
        m = masks{g};
        x = Xset{i}(m);
        y = A(m);
        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);
        scatter(x, y, 45, 'k', 'filled', 'MarkerFaceAlpha', 0.55);
        if numel(x) >= 5
            [xfit, yfit, ylo, yhi] = fit_line_ci(x, y);
            plot(xfit, yfit, 'r-', 'LineWidth', 2.5);
            fill([xfit; flipud(xfit)], [ylo; flipud(yhi)], [1 0.6 0.6], ...
                'FaceAlpha', 0.25, 'EdgeColor', 'none');
            [r, p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
            title(sprintf('%s vs Alpha (r=%.2f, p=%.3f)', names{i}, r, p), 'FontSize', fsz-6);
        else
            title(sprintf('%s vs Alpha (n too low)', names{i}), 'FontSize', fsz-6);
        end
        xlabel(names{i});
        ylabel('AlphaPower_FOOOF_bl');
        box on
        set(gca, 'FontSize', fsz-6);
    end
    sgtitle(sprintf('Correlation Panels - %s', groups{g}), 'FontSize', fsz);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_correlations_%s.png', groups{g})));
    close(gcf);
end
end

function [xfit, yfit, ylo, yhi] = fit_line_ci(x, y)
xfit = linspace(min(x), max(x), 120)';
try
    mdl = fitlm(x, y);
    [yfit, ci] = predict(mdl, xfit);
    ylo = ci(:, 1);
    yhi = ci(:, 2);
catch
    p = polyfit(x, y, 1);
    yfit = polyval(p, xfit);
    resid = y - polyval(p, x);
    se = std(resid, 'omitnan');
    ylo = yfit - 1.96*se;
    yhi = yfit + 1.96*se;
end
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha)
y = yvals(isfinite(yvals));
if numel(y) < 3
    return
end
[f, xi] = ksdensity(y, 'NumPoints', 120);
f = f / max(f) * 0.35;
fill([xpos - f, fliplr(repmat(xpos, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25);
q3 = prctile(y, 75);
med = median(y);
p5 = prctile(y, 5);
p95 = prctile(y, 95);
plot([xpos xpos], [p5 q1], '-k', 'LineWidth', 1.2);
plot([xpos xpos], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [xpos-box_w/2, q1, box_w, q3-q1], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(xpos + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = 0.10*(rand(numel(y),1)-0.5);
scatter(xpos + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
Ts = 1 / fs;
L = numel(X);
framelen = min(21, L);
if mod(framelen, 2) == 0
    framelen = framelen - 1;
end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0
    minLegal = minLegal + 1;
end
if framelen < minLegal
    framelen = minLegal;
end
if framelen > L
    framelen = L - mod(L, 2) + 1;
end
useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / Ts) * G(:, 2)';
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
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
