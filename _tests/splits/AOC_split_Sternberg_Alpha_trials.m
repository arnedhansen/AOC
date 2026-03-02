%% AOC Split Sternberg Alpha (Trial-Level)
% Trial-level split using merged_data_sternberg_trials:
%   AlphaPowerFullBL < 0  -> reduction trials
%   AlphaPowerFullBL > 0  -> amplification trials
%   AlphaPowerFullBL == 0 -> excluded
%
% Output figures are saved to:
%   .../figures/tests/splitAlpha/trials/

%% Setup
clear
clc
close all

startup
[subjects, pathAOC, colors, headmodel] = setup('AOC');

if ispc
    base_data = 'W:\Students\Arne\AOC';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
end

feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'tests', 'splitAlpha', 'trials');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

fig_pos = [0 0 1512 982];
fontSize = 20;
rng(42);
alpha_zero_margin = 0.05; % Exclude near-zero trial alpha from amp/red split

cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

%% Load trial-level merged table and define trial groups
merged_file = fullfile(feat_dir, 'merged_data_sternberg_trials.mat');
if ~isfile(merged_file)
    error('Missing file: %s', merged_file);
end
S = load(merged_file, 'merged_data_sternberg_trials');
if ~isfield(S, 'merged_data_sternberg_trials')
    error('Variable merged_data_sternberg_trials not found in %s', merged_file);
end
T = S.merged_data_sternberg_trials;

is_red_trial = isfinite(T.AlphaPowerFullBL) & T.AlphaPowerFullBL < -alpha_zero_margin;
is_amp_trial = isfinite(T.AlphaPowerFullBL) & T.AlphaPowerFullBL > alpha_zero_margin;
is_zero_trial = isfinite(T.AlphaPowerFullBL) & abs(T.AlphaPowerFullBL) <= alpha_zero_margin;

fprintf('\n=== Trial Split Summary (AlphaPowerFullBL) ===\n');
fprintf('Total rows: %d\n', height(T));
fprintf('Reduction trials (< -%.3f): %d\n', alpha_zero_margin, nnz(is_red_trial));
fprintf('Amplification trials (> %.3f): %d\n', alpha_zero_margin, nnz(is_amp_trial));
fprintf('Excluded trials (|alpha| <= %.3f): %d\n', alpha_zero_margin, nnz(is_zero_trial));

if nnz(is_red_trial) < 20 || nnz(is_amp_trial) < 20
    error('Insufficient trial counts in one or both trial groups.');
end

uIDs = unique(T.ID);
nSubj = numel(uIDs);

%% Containers
pow_red = cell(1, 3);
pow_amp = cell(1, 3);
tfr_red = cell(1, 3);
tfr_amp = cell(1, 3);
for c = 1:3
    pow_red{c} = {};
    pow_amp{c} = {};
    tfr_red{c} = {};
    tfr_amp{c} = {};
end

metrics_red = struct('Alpha', nan(nSubj,3), 'SPL', nan(nSubj,3), 'Dev', nan(nSubj,3), ...
    'BCEA', nan(nSubj,3), 'MS', nan(nSubj,3), 'Vel', nan(nSubj,3));
metrics_amp = struct('Alpha', nan(nSubj,3), 'SPL', nan(nSubj,3), 'Dev', nan(nSubj,3), ...
    'BCEA', nan(nSubj,3), 'MS', nan(nSubj,3), 'Vel', nan(nSubj,3));

fs = 500;
t_full = -0.5:1/fs:2;
t_plot = t_full(2:end);
Tf = numel(t_plot);
spl_tc_red = nan(nSubj, 3, Tf);
spl_tc_amp = nan(nSubj, 3, Tf);
dev_tc_red = nan(nSubj, 3, Tf);
dev_tc_amp = nan(nSubj, 3, Tf);
vel_tc_red = nan(nSubj, 3, Tf);
vel_tc_amp = nan(nSubj, 3, Tf);

missing_eeg = {};
missing_tfr = {};
missing_gaze = {};

%% Aggregate per subject and condition
for s = 1:nSubj
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end

    rows_s = T(T.ID == sid, :);

    for c = 1:3
        rows_sc = rows_s(rows_s.Condition == cond_vals(c), :);
        if isempty(rows_sc)
            continue
        end
        red_mask = isfinite(rows_sc.AlphaPowerFullBL) & rows_sc.AlphaPowerFullBL < -alpha_zero_margin;
        amp_mask = isfinite(rows_sc.AlphaPowerFullBL) & rows_sc.AlphaPowerFullBL > alpha_zero_margin;

        metrics_red.Alpha(s,c) = mean(rows_sc.AlphaPowerFullBL(red_mask), 'omitnan');
        metrics_red.SPL(s,c) = mean(rows_sc.ScanPathLengthFullBL(red_mask), 'omitnan');
        metrics_red.Dev(s,c) = mean(rows_sc.GazeDeviationFullBL(red_mask), 'omitnan');
        metrics_red.BCEA(s,c) = mean(rows_sc.BCEAFullBL(red_mask), 'omitnan');
        metrics_red.MS(s,c) = mean(rows_sc.MSRateFullBL(red_mask), 'omitnan');

        metrics_amp.Alpha(s,c) = mean(rows_sc.AlphaPowerFullBL(amp_mask), 'omitnan');
        metrics_amp.SPL(s,c) = mean(rows_sc.ScanPathLengthFullBL(amp_mask), 'omitnan');
        metrics_amp.Dev(s,c) = mean(rows_sc.GazeDeviationFullBL(amp_mask), 'omitnan');
        metrics_amp.BCEA(s,c) = mean(rows_sc.BCEAFullBL(amp_mask), 'omitnan');
        metrics_amp.MS(s,c) = mean(rows_sc.MSRateFullBL(amp_mask), 'omitnan');
    end

    % Trial numbers by split and condition from merged table.
    red_trials_by_cond = cell(1,3);
    amp_trials_by_cond = cell(1,3);
    for c = 1:3
        rows_sc = rows_s(rows_s.Condition == cond_vals(c), :);
        red_trials_by_cond{c} = rows_sc.Trial(isfinite(rows_sc.AlphaPowerFullBL) & rows_sc.AlphaPowerFullBL < -alpha_zero_margin);
        amp_trials_by_cond{c} = rows_sc.Trial(isfinite(rows_sc.AlphaPowerFullBL) & rows_sc.AlphaPowerFullBL > alpha_zero_margin);
    end

    % Power spectra (trial-level source)
    eeg_dir = fullfile(pathAOC, subj_folder, 'eeg');
    try
        P = load(fullfile(eeg_dir, 'power_stern_full_trials.mat'), ...
            'powload2_full', 'powload4_full', 'powload6_full');
        pconds = {P.powload2_full, P.powload4_full, P.powload6_full};
        for c = 1:3
            [red_sel, amp_sel] = split_freq_trials_struct(pconds{c}, red_trials_by_cond{c}, amp_trials_by_cond{c});
            if ~isempty(red_sel), pow_red{c}{end+1} = red_sel; end %#ok<AGROW>
            if ~isempty(amp_sel), pow_amp{c}{end+1} = amp_sel; end %#ok<AGROW>
        end
    catch
        missing_eeg{end+1} = sid_str; %#ok<AGROW>
    end

    % TFR (FOOOFed + baselined trial-level if available; fallback = baseline-corrected raw trials)
    try
        tfr_file_fooof_trials = fullfile(eeg_dir, 'tfr_stern_fooof_trials.mat');
        has_trial_fooof = false;
        if isfile(tfr_file_fooof_trials)
            R = load(tfr_file_fooof_trials);
            if isfield(R, 'tfr_all_fooof_bl')
                tfr_src = R.tfr_all_fooof_bl;
                has_trial_fooof = true;
            elseif isfield(R, 'tfr_all_fooof')
                cfgb = [];
                cfgb.baseline = [-0.5 -0.25];
                cfgb.baselinetype = 'absolute';
                tfr_src = ft_freqbaseline(cfgb, R.tfr_all_fooof);
                has_trial_fooof = true;
            end
        end
        if ~has_trial_fooof
            R = load(fullfile(eeg_dir, 'tfr_stern_trials.mat'), 'tfr_all');
            tfr_src = R.tfr_all;
            warning('No trial-wise FOOOF TFR for subject %s. Using baseline-corrected raw trial TFR.', sid_str);
            cfgb = [];
            cfgb.baseline = [-0.5 -0.25];
            cfgb.baselinetype = 'absolute';
            tfr_src = ft_freqbaseline(cfgb, tfr_src);
        end
        if ~isfield(tfr_src, 'trialinfo')
            error('Missing trialinfo in TFR source');
        end
        trinfo = tfr_src.trialinfo;
        conds_tfr = parse_trialinfo_conds(trinfo);
        trials_tfr = parse_trialinfo_trials(trinfo);
        for c = 1:3
            cond_mask = conds_tfr == cond_codes(c);
            idx_red = find(cond_mask & ismember(trials_tfr, red_trials_by_cond{c}));
            idx_amp = find(cond_mask & ismember(trials_tfr, amp_trials_by_cond{c}));
            if ~isempty(idx_red)
                cfg = [];
                cfg.trials = idx_red;
                cfg.channel = occ_channels_from_labels(tfr_src.label);
                cfg.avgoverrpt = 'yes';
                tfr_red{c}{end+1} = ft_selectdata(cfg, tfr_src); %#ok<AGROW>
            end
            if ~isempty(idx_amp)
                cfg = [];
                cfg.trials = idx_amp;
                cfg.channel = occ_channels_from_labels(tfr_src.label);
                cfg.avgoverrpt = 'yes';
                tfr_amp{c}{end+1} = ft_selectdata(cfg, tfr_src); %#ok<AGROW>
            end
        end
    catch
        missing_tfr{end+1} = sid_str; %#ok<AGROW>
    end

    % Gaze time courses and velocity summaries
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
    conds_gz = parse_trialinfo_conds(G.trialinfo);
    trnums_gz = parse_trialinfo_trials(G.trialinfo);
    if isempty(conds_gz) || isempty(trnums_gz)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    if isfield(G, 'ScanPathSeries') && isfield(G, 'ScanPathSeriesT')
        for c = 1:3
            tr_idx_red = find(conds_gz == cond_codes(c) & ismember(trnums_gz, red_trials_by_cond{c}));
            tr_idx_amp = find(conds_gz == cond_codes(c) & ismember(trnums_gz, amp_trials_by_cond{c}));
            spl_tc_red(s,c,:) = make_scanpath_mean(G.ScanPathSeries, G.ScanPathSeriesT, tr_idx_red, t_plot);
            spl_tc_amp(s,c,:) = make_scanpath_mean(G.ScanPathSeries, G.ScanPathSeriesT, tr_idx_amp, t_plot);
        end
    end

    if isfield(G, 'gaze_x') && isfield(G, 'gaze_y')
        for c = 1:3
            tr_idx_red = find(conds_gz == cond_codes(c) & ismember(trnums_gz, red_trials_by_cond{c}));
            tr_idx_amp = find(conds_gz == cond_codes(c) & ismember(trnums_gz, amp_trials_by_cond{c}));
            [dev_tc_red(s,c,:), vel_tc_red(s,c,:), metrics_red.Vel(s,c)] = ...
                make_dev_vel_mean(G.gaze_x, G.gaze_y, tr_idx_red, t_plot, fs);
            [dev_tc_amp(s,c,:), vel_tc_amp(s,c,:), metrics_amp.Vel(s,c)] = ...
                make_dev_vel_mean(G.gaze_x, G.gaze_y, tr_idx_amp, t_plot, fs);
        end
    end
end

%% Channel set
channels = {};
for c = 1:3
    if ~isempty(pow_red{c})
        channels = occ_channels_from_labels(pow_red{c}{1}.label);
        break
    elseif ~isempty(pow_amp{c})
        channels = occ_channels_from_labels(pow_amp{c}{1}.label);
        break
    end
end
if isempty(channels)
    error('No power trial data found for plotting.');
end

%% EEG visualizations
plot_group_power_spectrum(pow_red, channels, colors, cond_labels, ...
    'Reduction Trials Power Spectrum', fullfile(fig_dir, 'AOC_splitAlpha_trials_powspctrm_reduction.png'), fig_pos, fontSize);
plot_group_power_spectrum(pow_amp, channels, colors, cond_labels, ...
    'Amplification Trials Power Spectrum', fullfile(fig_dir, 'AOC_splitAlpha_trials_powspctrm_amplification.png'), fig_pos, fontSize);

color_map_tfr = interp1(linspace(0,1,5), ...
    [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], ...
    linspace(0,1,64));
plot_group_tfrs(tfr_red, channels, cond_labels, headmodel, color_map_tfr, 'ReductionTrials', fig_dir, fig_pos, fontSize);
plot_group_tfrs(tfr_amp, channels, cond_labels, headmodel, color_map_tfr, 'AmplificationTrials', fig_dir, fig_pos, fontSize);
plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, fig_dir, fig_pos, fontSize);

%% Boxplots
plot_metric_boxplots_trial(metrics_red, metrics_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% Rainclouds
plot_metric_rainclouds_trial(metrics_red, metrics_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% Within-subject trajectories
plot_metric_trajectories_trial(metrics_red, metrics_amp, cond_labels, colors, fig_dir, fig_pos, fontSize);

%% Time courses and effect-size panels
plot_timecourse_with_effect_trial(spl_tc_red, spl_tc_amp, cond_labels, colors, ...
    'Scan Path Length [px]', 'spl', fig_dir, fig_pos, fontSize);
plot_timecourse_with_effect_trial(vel_tc_red, vel_tc_amp, cond_labels, colors, ...
    'Eye Velocity [px/s]', 'velocity', fig_dir, fig_pos, fontSize);
plot_timecourse_with_effect_trial(dev_tc_red, dev_tc_amp, cond_labels, colors, ...
    'Gaze Deviation [px]', 'gaze_deviation', fig_dir, fig_pos, fontSize);

%% Correlation panels
plot_correlation_panels_trial(metrics_red, metrics_amp, fig_dir, fig_pos, fontSize);

%% Diagnostics
fprintf('\n=== Trial-Level Diagnostics ===\n');
fprintf('Subjects in merged trials table: %d\n', nSubj);
fprintf('Missing EEG trial files: %d\n', numel(unique(missing_eeg)));
fprintf('Missing TFR trial files: %d\n', numel(unique(missing_tfr)));
fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));
fprintf('Figures saved to: %s\n', fig_dir);

%% ========================= Local Functions =========================
function subj_folder = resolve_subject_folder(subjects, sid)
subj_folder = '';
for i = 1:numel(subjects)
    s = str2double(subjects{i});
    if isfinite(s) && s == sid
        subj_folder = subjects{i};
        return
    end
end
end

function conds = parse_trialinfo_conds(trialinfo)
conds = [];
if isempty(trialinfo), return; end
if isvector(trialinfo)
    conds = trialinfo(:);
elseif size(trialinfo,2) >= 1
    if size(trialinfo,2) == 2
        conds = trialinfo(:,1);
    elseif size(trialinfo,1) == 2
        conds = trialinfo(1,:)';
    else
        conds = trialinfo(:,1);
    end
end
end

function trials = parse_trialinfo_trials(trialinfo)
trials = [];
if isempty(trialinfo), return; end
if isvector(trialinfo)
    return
elseif size(trialinfo,2) >= 2
    trials = trialinfo(:,2);
elseif size(trialinfo,1) == 2
    trials = trialinfo(2,:)';
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

function [red_sel, amp_sel] = split_freq_trials_struct(F, red_trials, amp_trials)
red_sel = [];
amp_sel = [];
if ~isfield(F, 'trialinfo') || isempty(F.trialinfo)
    return
end
trials = parse_trialinfo_trials(F.trialinfo);
if isempty(trials)
    return
end
idx_red = find(ismember(trials, red_trials));
idx_amp = find(ismember(trials, amp_trials));
if ~isempty(idx_red)
    cfg = [];
    cfg.trials = idx_red;
    red_sel = ft_selectdata(cfg, F);
end
if ~isempty(idx_amp)
    cfg = [];
    cfg.trials = idx_amp;
    amp_sel = ft_selectdata(cfg, F);
end
end

function out = make_scanpath_mean(seriesCell, timeCell, tr_idx, t_plot)
out = nan(1, numel(t_plot));
if isempty(tr_idx), return; end
mat = nan(numel(tr_idx), numel(t_plot));
for k = 1:numel(tr_idx)
    tr = tr_idx(k);
    s = seriesCell{tr};
    t = timeCell{tr};
    if isempty(s) || isempty(t), continue; end
    try
        mat(k,:) = interp1(t, s, t_plot, 'linear', NaN);
    catch
    end
end
out = nanmean(mat, 1);
end

function [dev_mean, vel_mean, vel_bl] = make_dev_vel_mean(gaze_x, gaze_y, tr_idx, t_plot, fs)
Tf = numel(t_plot);
dev_mean = nan(1, Tf);
vel_mean = nan(1, Tf);
vel_bl = NaN;
if isempty(tr_idx), return; end
dev_mat = nan(numel(tr_idx), Tf);
vel_mat = nan(numel(tr_idx), Tf);
vel_bl_trials = nan(numel(tr_idx), 1);
for k = 1:numel(tr_idx)
    tr = tr_idx(k);
    x = double(gaze_x{tr});
    y = double(gaze_y{tr});
    if isempty(x) || isempty(y) || numel(x) ~= numel(y), continue; end
    tt = linspace(-0.5, 2, numel(x));
    dev = sqrt((x - 400).^2 + (y - 300).^2);
    [vx, vy] = compute_velocity_sg(x, y, fs, 3);
    [vx, vy] = clean_velocity_components(vx, vy);
    vel = hypot(vx, vy);
    try
        dev_mat(k,:) = interp1(tt, dev, t_plot, 'linear', NaN);
        vel_mat(k,:) = interp1(tt, vel, t_plot, 'linear', NaN);
    catch
    end
    idx_full = tt >= 0 & tt <= 2;
    idx_bl = tt >= -0.5 & tt <= -0.25;
    v_full = mean(vel(idx_full), 'omitnan');
    v_bl = mean(vel(idx_bl), 'omitnan');
    if isfinite(v_full) && isfinite(v_bl) && v_bl > 0
        vel_bl_trials(k) = 10 * log10(v_full / v_bl);
    end
end

function [vx, vy] = clean_velocity_components(vx, vy)
halfwin = 11;
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
dev_mean = nanmean(dev_mat, 1);
vel_mean = nanmean(vel_mat, 1);
vel_bl = mean(vel_bl_trials, 'omitnan');
end

function plot_group_power_spectrum(pow_cells, channels, colors, cond_labels, ttl, out_file, fig_pos, fsz)
if any(cellfun(@isempty, pow_cells))
    warning('Skipping %s (incomplete power condition data).', ttl);
    return
end
figure('Position', fig_pos, 'Color', 'w');
hold on
for c = 1:3
    subj_n = numel(pow_cells{c});
    if subj_n < 2, continue; end
    p0 = pow_cells{c}{1};
    freqs = p0.freq;
    subj_spec = nan(subj_n, numel(freqs));
    for s = 1:subj_n
        P = pow_cells{c}{s};
        if contains(P.dimord, 'rpt')
            x = squeeze(mean(P.powspctrm, 1, 'omitnan')); % chan x freq
        else
            x = P.powspctrm;
        end
        elecs = ismember(P.label, channels);
        subj_spec(s,:) = mean(x(elecs,:), 1, 'omitnan');
    end
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(subj_spec), 1));
    eb(c) = shadedErrorBar(freqs, m, se, 'lineProps', {'-'}, 'transparent', true); %#ok<AGROW>
    set(eb(c).mainLine, 'Color', colors(c,:), 'LineWidth', 3);
    set(eb(c).patch, 'FaceColor', colors(c,:), 'FaceAlpha', 0.25);
end
xlim([5 20]);
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title(ttl, 'FontSize', fsz+2);
legend([eb.mainLine], cond_labels, 'Location', 'best', 'FontSize', fsz-2);
set(gca, 'FontSize', fsz);
box on
saveas(gcf, out_file);
close(gcf);
end

function plot_group_tfrs(tfr_cells, channels, cond_labels, headmodel, color_map, grp_label, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, tfr_cells))
    warning('Skipping TFR plot for %s (incomplete condition data).', grp_label);
    return
end
ga = cell(1,3);
for c = 1:3
    ga{c} = ft_freqgrandaverage([], tfr_cells{c}{:});
end
[~, ch_idx] = ismember(channels, ga{1}.label);
freq_idx = ga{1}.freq >= 5 & ga{1}.freq <= 30;
time_idx = ga{1}.time >= -0.5 & ga{1}.time <= 2;
mx = 0;
for c = 1:3
    A = squeeze(mean(ga{c}.powspctrm(ch_idx,:,:), 1, 'omitnan'));
    mx = max(mx, max(abs(A(freq_idx,time_idx)), [], 'all'));
end
clim = [-mx mx];
cfg = [];
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;
for c = 1:3
    figure('Position', fig_pos, 'Color', 'w');
    ft_singleplotTFR(cfg, ga{c});
    colormap(color_map);
    set(gca, 'CLim', clim);
    cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fsz);
    xlabel('Time [s]'); ylabel('Frequency [Hz]');
    rectangle('Position', [0,8,2,6], 'EdgeColor', 'k', 'LineWidth', 4);
    set(gca, 'FontSize', fsz);
    title(sprintf('%s - %s', grp_label, cond_labels{c}), 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_tfr_%s_cond%d.png', lower(grp_label), c)));
    close(gcf);
end
end

function plot_group_topos(pow_red, pow_amp, channels, headmodel, cond_labels, fig_dir, fig_pos, fsz)
if any(cellfun(@isempty, pow_red)) || any(cellfun(@isempty, pow_amp))
    warning('Skipping topoplots (missing power trial data).');
    return
end
ga_red = cell(1,3); ga_amp = cell(1,3);
for c = 1:3
    ga_red{c} = ft_freqgrandaverage([], pow_red{c}{:});
    ga_amp{c} = ft_freqgrandaverage([], pow_amp{c}{:});
end
cfg = [];
cfg.layout = headmodel.layANThead;
allch = cfg.layout.label;
cfg.channel = allch(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel,'M1') & ~strcmp(cfg.channel,'M2'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
cfg.colormap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
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

for c = 1:3
    figure('Position', fig_pos, 'Color', 'w');
    ft_topoplotER(cfg, ga_red{c});
    cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]');
    set(gca, 'FontSize', fsz); title(sprintf('Reduction Trials - %s', cond_labels{c}));
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_topo_reduction_cond%d.png', c)));
    close(gcf);

    figure('Position', fig_pos, 'Color', 'w');
    ft_topoplotER(cfg, ga_amp{c});
    cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]');
    set(gca, 'FontSize', fsz); title(sprintf('Amplification Trials - %s', cond_labels{c}));
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_topo_amplification_cond%d.png', c)));
    close(gcf);
end

cfgd = cfg;
cmap_diff = cbrewer('div', 'RdBu', 100);
cfgd.colormap = flipud(max(min(cmap_diff,1),0));
all_diff = all_diff(isfinite(all_diff));
if isempty(all_diff)
    cfgd.zlim = 'maxabs';
else
    mx = prctile(abs(all_diff), 99);
    cfgd.zlim = [-mx mx];
end
for c = 1:3
    gd = ga_amp{c};
    gd.powspctrm = ga_amp{c}.powspctrm - ga_red{c}.powspctrm;
    figure('Position', fig_pos, 'Color', 'w');
    ft_topoplotER(cfgd, gd);
    cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]');
    set(gca, 'FontSize', fsz); title(sprintf('Trials Topo Diff (Amp-Red) - %s', cond_labels{c}));
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_topo_diff_cond%d.png', c)));
    close(gcf);
end
end

function plot_metric_boxplots_trial(R, A, cond_labels, colors, fig_dir, fig_pos, fsz)
defs = {'SPL','ScanPathLengthFullBL [dB]'; 'Vel','VelocityFullBL [dB]'; ...
    'Dev','GazeDeviationFullBL [dB]'; 'BCEA','BCEAFullBL [dB]'; 'MS','MSRateFullBL [dB]'};
for m = 1:size(defs,1)
    key = defs{m,1}; yl = defs{m,2};
    figure('Position', fig_pos, 'Color', 'w'); tiledlayout(1,2,'TileSpacing','compact');
    nexttile; hold on
    X = R.(key); boxplot(X, 'Colors','k','Symbol','','Widths',0.5);
    for c = 1:3, scatter(c+0.08*randn(size(X,1),1), X(:,c), 35, colors(c,:), 'filled', 'MarkerFaceAlpha',0.45); end
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Reduction trials'); ylabel(yl); box on
    nexttile; hold on
    X = A.(key); boxplot(X, 'Colors','k','Symbol','','Widths',0.5);
    for c = 1:3, scatter(c+0.08*randn(size(X,1),1), X(:,c), 35, colors(c,:), 'filled', 'MarkerFaceAlpha',0.45); end
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Amplification trials'); ylabel(yl); box on
    sgtitle(sprintf('Trial-level %s by condition', key), 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_boxplot_%s.png', lower(key))));
    close(gcf);
end
end

function plot_metric_rainclouds_trial(R, A, cond_labels, colors, fig_dir, fig_pos, fsz)
defs = {'Alpha','AlphaPowerFullBL'; 'SPL','ScanPathLengthFullBL'; 'Vel','VelocityFullBL'; ...
    'Dev','GazeDeviationFullBL'; 'BCEA','BCEAFullBL'; 'MS','MSRateFullBL'};
for m = 1:size(defs,1)
    key = defs{m,1}; ylab = defs{m,2};
    figure('Position', fig_pos, 'Color', 'w'); tiledlayout(1,2,'TileSpacing','compact');
    nexttile; hold on
    X = R.(key);
    for c = 1:3, draw_one_cloud(X(:,c), c, colors(c,:), 0.3, 24, 0.45); end
    h1 = plot(nan, nan, '-', 'Color', colors(1,:), 'LineWidth', 2);
    h2 = plot(nan, nan, '-', 'Color', colors(2,:), 'LineWidth', 2);
    h3 = plot(nan, nan, '-', 'Color', colors(3,:), 'LineWidth', 2);
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Reduction trials'); ylabel(ylab); box off
    legend([h1 h2 h3], cond_labels, 'Location', 'best', 'FontSize', fsz-5);
    nexttile; hold on
    X = A.(key);
    for c = 1:3, draw_one_cloud(X(:,c), c, colors(c,:), 0.3, 24, 0.45); end
    h1 = plot(nan, nan, '-', 'Color', colors(1,:), 'LineWidth', 2);
    h2 = plot(nan, nan, '-', 'Color', colors(2,:), 'LineWidth', 2);
    h3 = plot(nan, nan, '-', 'Color', colors(3,:), 'LineWidth', 2);
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Amplification trials'); ylabel(ylab); box off
    legend([h1 h2 h3], cond_labels, 'Location', 'best', 'FontSize', fsz-5);
    sgtitle(sprintf('Trial-level raincloud: %s', key), 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_raincloud_%s.png', lower(key))));
    close(gcf);
end
end

function plot_metric_trajectories_trial(R, A, cond_labels, colors, fig_dir, fig_pos, fsz)
defs = {'Alpha','SPL','Vel','Dev','BCEA','MS'};
for m = 1:numel(defs)
    key = defs{m};
    figure('Position', fig_pos, 'Color', 'w'); tiledlayout(1,2,'TileSpacing','compact');
    nexttile; hold on
    X = R.(key);
    for i = 1:size(X,1), plot(1:3, X(i,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1); end
    hSub = plot(1:3, X(1,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hMean = plot(1:3, mean(X,1,'omitnan'), '-o', 'Color', colors(1,:), 'LineWidth', 3, 'MarkerFaceColor', colors(1,:));
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Reduction trials'); ylabel(key); box on
    legend([hSub hMean], {'Subjects', 'Group mean'}, 'Location', 'best', 'FontSize', fsz-6);
    nexttile; hold on
    X = A.(key);
    for i = 1:size(X,1), plot(1:3, X(i,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1); end
    hSub = plot(1:3, X(1,:), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hMean = plot(1:3, mean(X,1,'omitnan'), '-o', 'Color', colors(3,:), 'LineWidth', 3, 'MarkerFaceColor', colors(3,:));
    set(gca,'XTick',1:3,'XTickLabel',cond_labels,'FontSize',fsz-2); title('Amplification trials'); ylabel(key); box on
    legend([hSub hMean], {'Subjects', 'Group mean'}, 'Location', 'best', 'FontSize', fsz-6);
    sgtitle(sprintf('Trial-level trajectories: %s', key), 'FontSize', fsz+2);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_trajectory_%s.png', lower(key))));
    close(gcf);
end
end

function plot_timecourse_with_effect_trial(red_tc, amp_tc, cond_labels, colors, ylab, tag, fig_dir, fig_pos, fsz)
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(4,1,'TileSpacing','compact');
t = linspace(-0.5, 2, size(red_tc,3));
for c = 1:3
    nexttile(c); hold on
    R = squeeze(red_tc(:,c,:)); A = squeeze(amp_tc(:,c,:));
    mR = mean(R,1,'omitnan'); mA = mean(A,1,'omitnan');
    sR = std(R,0,1,'omitnan') ./ sqrt(sum(isfinite(R),1));
    sA = std(A,0,1,'omitnan') ./ sqrt(sum(isfinite(A),1));
    e1 = shadedErrorBar(t,mR,sR,'lineProps',{'-'},'transparent',true);
    e2 = shadedErrorBar(t,mA,sA,'lineProps',{'-'},'transparent',true);
    set(e1.mainLine,'Color',colors(1,:),'LineWidth',2.5); set(e2.mainLine,'Color',colors(3,:),'LineWidth',2.5);
    set(e1.patch,'FaceColor',colors(1,:),'FaceAlpha',0.25); set(e2.patch,'FaceColor',colors(3,:),'FaceAlpha',0.25);
    xline(0,'--k'); xlim([-0.5 2]); ylabel(ylab); title(cond_labels{c}); box on; set(gca,'FontSize',fsz-4);
    legend([e1.mainLine e2.mainLine], {'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fsz-7);
end
nexttile(4); hold on
Rall = squeeze(mean(red_tc,2,'omitnan')); Aall = squeeze(mean(amp_tc,2,'omitnan'));
d = nan(1,size(red_tc,3));
for i = 1:numel(d)
    x = Rall(:,i); y = Aall(:,i); x = x(isfinite(x)); y = y(isfinite(y));
    if numel(x)>=3 && numel(y)>=3
        sp = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y)) / (numel(x)+numel(y)-2));
        d(i) = (mean(y)-mean(x))/max(sp,eps);
    end
end
plot(t,d,'k-','LineWidth',2.5); yline(0,'--'); xline(0,'--k'); xlim([-0.5 2]);
xlabel('Time [s]'); ylabel('Cohen''s d'); box on; set(gca,'FontSize',fsz-4);
saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_timecourse_%s.png', tag)));
close(gcf);
end

function plot_correlation_panels_trial(R, A, fig_dir, fig_pos, fsz)
for g = 1:2
    if g == 1, M = R; gl = 'reduction_trials'; else, M = A; gl = 'amplification_trials'; end
    alpha = mean(M.Alpha,2,'omitnan');
    Xset = {mean(M.SPL,2,'omitnan'), mean(M.Vel,2,'omitnan'), mean(M.Dev,2,'omitnan'), ...
        mean(M.BCEA,2,'omitnan'), mean(M.MS,2,'omitnan')};
    names = {'SPL','Velocity','GazeDeviation','BCEA','Microsaccades'};
    figure('Position', fig_pos, 'Color', 'w'); tiledlayout(2,3,'TileSpacing','compact');
    for i = 1:numel(Xset)
        nexttile; hold on
        x = Xset{i}; y = alpha; valid = isfinite(x) & isfinite(y); x = x(valid); y = y(valid);
        scatter(x,y,45,'k','filled','MarkerFaceAlpha',0.55);
        if numel(x) >= 5
            [xfit, yfit, ylo, yhi] = fit_line_ci(x,y);
            plot(xfit,yfit,'r-','LineWidth',2.5);
            fill([xfit; flipud(xfit)], [ylo; flipud(yhi)], [1 0.6 0.6], 'FaceAlpha',0.25,'EdgeColor','none');
            [r,p] = corr(x,y,'Type','Spearman','Rows','complete');
            title(sprintf('%s vs Alpha (r=%.2f, p=%.3f)', names{i}, r, p), 'FontSize', fsz-6);
        else
            title(sprintf('%s vs Alpha (n too low)', names{i}), 'FontSize', fsz-6);
        end
        xlabel(names{i}); ylabel('AlphaPowerFullBL'); box on; set(gca,'FontSize',fsz-6);
    end
    sgtitle(sprintf('Trial-level correlations - %s', strrep(gl,'_',' ')), 'FontSize', fsz);
    saveas(gcf, fullfile(fig_dir, sprintf('AOC_splitAlpha_trials_correlations_%s.png', gl)));
    close(gcf);
end
end

function [xfit, yfit, ylo, yhi] = fit_line_ci(x, y)
xfit = linspace(min(x), max(x), 120)';
try
    mdl = fitlm(x, y);
    [yfit, ci] = predict(mdl, xfit);
    ylo = ci(:,1); yhi = ci(:,2);
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
if numel(y) < 3, return; end
[f, xi] = ksdensity(y, 'NumPoints', 120);
f = f / max(f) * 0.35;
fill([xpos - f, fliplr(repmat(xpos,1,numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y,25); q3 = prctile(y,75); med = median(y); p5 = prctile(y,5); p95 = prctile(y,95);
plot([xpos xpos], [p5 q1], '-k', 'LineWidth', 1.2);
plot([xpos xpos], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [xpos-box_w/2, q1, box_w, q3-q1], 'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(xpos + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = 0.10*(rand(numel(y),1)-0.5);
scatter(xpos + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', 'none');
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
Ts = 1 / fs;
L = numel(X);
framelen = min(21, L);
if mod(framelen,2)==0, framelen = framelen - 1; end
minLegal = polyOrd + 3;
if mod(minLegal,2)==0, minLegal = minLegal + 1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L, framelen = L - mod(L,2) + 1; end
useFallback = framelen < 5;
if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / Ts) * G(:,2)';
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end
