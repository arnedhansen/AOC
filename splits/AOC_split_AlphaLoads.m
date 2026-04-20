%% AOC Split Sternberg Alpha Loads
% Stratifies participants by alpha power slope across WM load (2,4,6).
% Alpha increase vs decrease subgroups; TFRs, power spectra, behavioral (RT, ACC),
% optional gaze density. Uses 1-3 s retention window for alpha extraction.
%%
% Groups:
%       inc: positive alpha–load slope (increasing alpha across WM loads)
%       dec: negative alpha–load slope (decreasing alpha across WM loads)
%       flat: near-zero slope

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
path = paths.features;

if ispc
    base_data = 'W:\Students\Arne\AOC';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
end
feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'splits', 'SplitAlphaLoads');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
cond_vals = [2 4 6];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

% Load merged subject-level table for gaze deviation summaries.
merged_file = fullfile(feat_dir, 'AOC_merged_data_sternberg.mat');
if ~isfile(merged_file)
    error('Missing file: %s', merged_file);
end
S_merged = load(merged_file, 'merged_data_sternberg');
if ~isfield(S_merged, 'merged_data_sternberg')
    error('Variable merged_data_sternberg not found in %s', merged_file);
end
T_merged = struct2table(S_merged.merged_data_sternberg);

log_dir = fullfile(base_data, 'data', 'controls', 'logs');
if ~isfolder(log_dir)
    mkdir(log_dir);
end
cmdlog_file = fullfile(log_dir, sprintf('AOC_split_AlphaLoads_commandwindow_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
diary('off');
diary(cmdlog_file);
cleanup_diary = onCleanup(@() diary('off'));
fprintf('Command window log file: %s\n', cmdlog_file);

%% Figure setup
fig_pos = [0 0 1512 982];
fontSize = 15;
color_map = customcolormap_preset('red-white-blue');
winsor_cfg = struct();
winsor_cfg.enable = true;
winsor_cfg.prctile = [2 98]; % subject-level clipping per TF bin
split_cfg = struct();
split_cfg.use_threshold_override = true; % set true to override symmetric tail split
split_cfg.threshold = 0.015 % default absolute slope threshold for override

% Gaze heatmap color limits (diverging maps): robust scale for sparse, heavy-tailed fields
gaze_zlim_cfg = struct();
gaze_zlim_cfg.min_abs = 1e-18;   % values with |x| <= this treated as empty-bin zeros for scaling
gaze_zlim_cfg.prctile = 90;      % percentile of |x| (above min_abs)
gaze_zlim_cfg.k_mad = 4;         % symmetric MAD fence: median(|x|)+k*MAD(|x|)
gaze_zlim_cfg.fallback_abs = 3;  % if scale degenerates

%% Loop over subjects - load EEG TFR (specParam, baselined)
clc
cfg_bl = [];
cfg_bl.baseline     = [-.5 -.25];
cfg_bl.baselinetype = 'absolute';
for subj = 1:length(subjects)
    clc
    fprintf('LOADING Subject %d / %d', subj, length(subjects))
    eeg_dir = fullfile(path, subjects{subj}, 'eeg');
    f = fullfile(eeg_dir, 'tfr_stern.mat');
    if ~isfile(f)
        error('Missing: %s', f);
    end
    datTFR = load(f);
    if isfield(datTFR, 'tfr2_fooof_bl')
        load2{subj} = datTFR.tfr2_fooof_bl;
        load4{subj} = datTFR.tfr4_fooof_bl;
        load6{subj} = datTFR.tfr6_fooof_bl;
    else
        load2{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr2_fooof);
        load4{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr4_fooof);
        load6{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr6_fooof);
    end
end
clc; disp('LOADING FINISHED')

%% Determine occipital channels
% Occipital channels
occ_channels = {};
labels = load2{1, 1}.label;
for i = 1:length(labels)
    label = labels{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Compute grand average
cfg = [];
cfg.keepindividual = 'yes';
ga2 = ft_freqgrandaverage(cfg, load2{:});
ga4 = ft_freqgrandaverage(cfg, load4{:});
ga6 = ft_freqgrandaverage(cfg, load6{:});
if winsor_cfg.enable
    ga2 = winsorize_freq_subjects(ga2, winsor_cfg.prctile);
    ga4 = winsorize_freq_subjects(ga4, winsor_cfg.prctile);
    ga6 = winsorize_freq_subjects(ga6, winsor_cfg.prctile);
end

%% Select alpha power (retention 1-3 s)
cfg = [];
cfg.frequency = [8 14];
cfg.latency = [1 3];
cfg.channel = channels;
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
val2 = ft_selectdata(cfg,ga2);
val4 = ft_selectdata(cfg,ga4);
val6 = ft_selectdata(cfg,ga6);
alpha2 = val2.powspctrm;
alpha4 = val4.powspctrm;
alpha6 = val6.powspctrm;

%% Use subject-level regression to identify slopes
nSubj = length(subjects);

for subj = 1:nSubj
    y = [alpha2(subj), alpha4(subj), alpha6(subj)];
    Xmat = [ones(3,1), [2;4;6]];
    b = Xmat \ y';
    slope(subj,1) = b(2);
end

%% Split and Plot slope distribution (inclusion)
% Default: zero-centered symmetric tails.
% Optional override: absolute slope threshold (|slope| > threshold).
nSubj = numel(slope);
if split_cfg.use_threshold_override
    thr = split_cfg.threshold;
    if ~isscalar(thr) || ~isfinite(thr) || thr < 0
        error('split_cfg.threshold must be a finite, non-negative scalar.');
    end
    idx_inc = slope > thr;   % amplification
    idx_dec = slope < -thr;  % reduction
    idx_flat = ~(idx_inc | idx_dec); % intermediate
else
    % Zero-centered symmetric tails:
    % - amplification: largest positive slopes
    % - reduction: most negative slopes
    % - intermediate: remaining participants
    k = floor(nSubj/3); % requested equal tail size
    if k < 1
        error('Not enough subjects for a zero-centered tail split.');
    end

    pos_idx = find(slope > 0);
    [~, ord_pos] = sort(slope(pos_idx), 'descend');
    pos_sel = pos_idx(ord_pos(1:min(k, numel(pos_idx))));

    neg_idx = find(slope < 0);
    [~, ord_neg] = sort(slope(neg_idx), 'ascend');
    neg_sel = neg_idx(ord_neg(1:min(k, numel(neg_idx))));

    % Enforce equal N in both tails; handle sign-imbalanced cohorts gracefully.
    k_eff = min(numel(pos_sel), numel(neg_sel));
    if k_eff < k
        warning('Only %d participants per tail possible (requested %d).', k_eff, k);
    end
    pos_sel = pos_sel(1:k_eff);
    neg_sel = neg_sel(1:k_eff);

    idx_inc = false(nSubj,1); % amplification
    idx_dec  = false(nSubj,1); % reduction
    idx_flat   = true(nSubj,1);  % intermediate
    idx_inc(pos_sel) = true;
    idx_dec(neg_sel)  = true;
    idx_flat(idx_inc | idx_dec) = false;
end

% counts
n_inc = sum(idx_inc);
n_dec = sum(idx_dec);
n_f = sum(idx_flat);

% Group-defining cutoffs (for visualization only)
if any(idx_dec)
    t1 = max(slope(idx_dec));
else
    t1 = NaN;
end
if any(idx_inc)
    t2 = min(slope(idx_inc));
else
    t2 = NaN;
end

% Plot
close all
figure('Position', fig_pos, 'Color', 'w');
hold on
% Use one shared, constant bin width and align edges with group cutoffs.
n_bins_middle = 8;
bin_w = (t2 - t1) / n_bins_middle;
if ~isfinite(bin_w) || bin_w <= 0
    all_s = slope(isfinite(slope));
    if isempty(all_s)
        bin_w = 1;
    else
        bin_w = (max(all_s) - min(all_s)) / 12;
        if ~isfinite(bin_w) || bin_w <= 0
            bin_w = 1e-3;
        end
    end
end
min_s = min(slope);
max_s = max(slope);
k_left = ceil((t1 - min_s) / bin_w);
k_right = ceil((max_s - t1) / bin_w);
bin_edges = t1 + (-k_left:k_right) * bin_w;
histogram(slope(idx_inc), 'BinEdges', bin_edges, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6);
histogram(slope(idx_dec),  'BinEdges', bin_edges, 'FaceColor', [0 0 0.8], 'FaceAlpha', 0.6);
histogram(slope(idx_flat),   'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);
xline(t1, 'k--', 'LineWidth', 2);
xline(t2, 'k--', 'LineWidth', 2);
xlabel('Alpha power slope')
ylabel('Participants')
title('Linear Slope of Alpha Power across WM Load (2, 4, 6 items)', 'FontSize', 20)
legend({sprintf('Increase (N=%d)', n_inc), ...
    sprintf('Decrease (N=%d)', n_dec), ...
    sprintf('intermediate (N=%d)', n_f), ...
    'symmetric-tail cutoffs'}, 'Box', 'off')
box on
set(gca, 'FontSize', 15)
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_inclusion.png'));

fprintf('\nZero-centered symmetric-tail classification summary:\n');
fprintf('Increase (positive tail):     %d\n', n_inc);
fprintf('Decrease (negative tail):     %d\n', n_dec);
fprintf('Intermediate:     %d\n', n_f);

%% Grand averages per subgroup
cfg = [];
cfg.keepindividual = 'yes';
ga2_inc = ft_freqgrandaverage(cfg,load2{idx_inc});
ga4_inc = ft_freqgrandaverage(cfg,load4{idx_inc});
ga6_inc = ft_freqgrandaverage(cfg,load6{idx_inc});
ga2_dec = ft_freqgrandaverage(cfg,load2{idx_dec});
ga4_dec = ft_freqgrandaverage(cfg,load4{idx_dec});
ga6_dec = ft_freqgrandaverage(cfg,load6{idx_dec});
if winsor_cfg.enable
    ga2_inc = winsorize_freq_subjects(ga2_inc, winsor_cfg.prctile);
    ga4_inc = winsorize_freq_subjects(ga4_inc, winsor_cfg.prctile);
    ga6_inc = winsorize_freq_subjects(ga6_inc, winsor_cfg.prctile);
    ga2_dec = winsorize_freq_subjects(ga2_dec, winsor_cfg.prctile);
    ga4_dec = winsorize_freq_subjects(ga4_dec, winsor_cfg.prctile);
    ga6_dec = winsorize_freq_subjects(ga6_dec, winsor_cfg.prctile);
end

%% TFR EEG
close all
disp(upper('Plotting EEG TFRs...'))
mx_tfr = 0.25;
target_ticks = 11;
tick_step = 0.05 * ceil(((2 * mx_tfr) / max(target_ticks - 1, 1)) / 0.05);
tick_step = max(tick_step, 0.05);
tick_max = tick_step * ceil(mx_tfr / tick_step);
cb_ticks = -tick_max:tick_step:tick_max;

cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [5 30];
cfg.xlim = [-.5 3];
cfg.zlim = [-tick_max tick_max];
cfg.channel = channels;
cfg.layout = headmodel.layANThead;
figure('Position', fig_pos, 'Color', 'w');
subplot(3,2,1); ft_singleplotTFR(cfg, ga2_inc); title('Increase: WM load 2', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,3); ft_singleplotTFR(cfg, ga4_inc); title('Increase: WM load 4', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,5); ft_singleplotTFR(cfg, ga6_inc); title('Increase: WM load 6', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;

subplot(3,2,2); ft_singleplotTFR(cfg, ga2_dec); title('Decrease: WM load 2', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,4); ft_singleplotTFR(cfg, ga4_dec); title('Decrease: WM load 4', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,6); ft_singleplotTFR(cfg, ga6_dec); title('Decrease: WM load 6', 'FontSize', fontSize);
hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 4; c.Ticks = cb_ticks; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
colormap(gcf, color_map);
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_TFR_OVERVIEW.png'));

% Save each EEG TFR panel as an individual figure
tfr_data = {ga2_inc, ga4_inc, ga6_inc, ga2_dec, ga4_dec, ga6_dec};
tfr_titles = {'Increase: WM load 2', 'Increase: WM load 4', 'Increase: WM load 6', ...
    'Decrease: WM load 2', 'Decrease: WM load 4', 'Decrease: WM load 6'};
tfr_filenames = {'AOC_split_AlphaLoads_TFR_Increase2.png', ...
    'AOC_split_AlphaLoads_TFR_Increase4.png', ...
    'AOC_split_AlphaLoads_TFR_Increase6.png', ...
    'AOC_split_AlphaLoads_TFR_Decrease2.png', ...
    'AOC_split_AlphaLoads_TFR_Decrease4.png', ...
    'AOC_split_AlphaLoads_TFR_Decrease6.png'};

for ii = 1:numel(tfr_data)
    figure('Position', fig_pos, 'Color', 'w');
    ft_singleplotTFR(cfg, tfr_data{ii});
    hold on; rectangle('Position', [1 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
    title(tfr_titles{ii}, 'FontSize', fontSize);
    ax = gca;
    c = colorbar;
    xlabel(ax, 'Time [s]', 'FontSize', fontSize);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
    set(ax, 'FontSize', fontSize);
    c.FontSize = fontSize - 4;
    c.Ticks = cb_ticks;
    c.Label.String = 'Power [dB]';
    c.Label.FontSize = fontSize;
    colormap(gcf, color_map);
    drawnow;
    saveas(gcf, fullfile(fig_dir, tfr_filenames{ii}));
    close(gcf);
end

%% Select powspctrm (retention 1-3 s)
cfg = [];
cfg.latency = [1 3];
cfg.avgovertime = 'yes';
ga2_inc_powspctrm = ft_selectdata(cfg,ga2_inc);
ga2_inc_powspctrm = rmfield(ga2_inc_powspctrm,'time');

ga4_inc_powspctrm = ft_selectdata(cfg,ga4_inc);
ga4_inc_powspctrm = rmfield(ga4_inc_powspctrm,'time');

ga6_inc_powspctrm = ft_selectdata(cfg,ga6_inc);
ga6_inc_powspctrm = rmfield(ga6_inc_powspctrm,'time');

ga2_dec_powspctrm = ft_selectdata(cfg,ga2_dec);
ga2_dec_powspctrm = rmfield(ga2_dec_powspctrm,'time');

ga4_dec_powspctrm = ft_selectdata(cfg,ga4_dec);
ga4_dec_powspctrm = rmfield(ga4_dec_powspctrm,'time');

ga6_dec_powspctrm = ft_selectdata(cfg,ga6_dec);
ga6_dec_powspctrm = rmfield(ga6_dec_powspctrm,'time');

%% Plot Powerspectra
close all
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

[~, ch_idx_inc] = ismember(channels, ga2_inc_powspctrm.label);
ch_idx_inc = ch_idx_inc(ch_idx_inc > 0);
freqs_inc = ga2_inc_powspctrm.freq;

addpath('W:\Students\Arne\toolboxes\shadedErrorBar\')
nexttile; hold on
ga_inc = {ga2_inc_powspctrm, ga4_inc_powspctrm, ga6_inc_powspctrm};
for c = 1:3
    subj_spec = squeeze(mean(ga_inc{c}.powspctrm(:, ch_idx_inc, :), 2, 'omitnan'));
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
    eb_inc(c) = shadedErrorBar(freqs_inc, m, se, 'lineProps', {'-'}, 'transparent', true); 
    set(eb_inc(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
    set(eb_inc(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
    set(eb_inc(c).edge(1), 'Color', colors(c, :));
    set(eb_inc(c).edge(2), 'Color', colors(c, :));
end
xlim([2 30]);
xlabel('Frequency [Hz]', 'FontSize', fontSize-4);
yline(0, '--')
ylim([-0.2 0.2])
ylabel('Power [dB]', 'FontSize', fontSize);
title('Increasing slope', 'FontSize', fontSize);
% Legend with colored patch boxes
leg_p_inc = gobjects(1, 3);
for c = 1:3
    leg_p_inc(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_inc, {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', fontSize - 2);
set(gca, 'FontSize', fontSize);
box on

[~, ch_idx_dec] = ismember(channels, ga2_dec_powspctrm.label);
ch_idx_dec = ch_idx_dec(ch_idx_dec > 0);
freqs_dec = ga2_dec_powspctrm.freq;

nexttile; hold on
ga_dec = {ga2_dec_powspctrm, ga4_dec_powspctrm, ga6_dec_powspctrm};
for c = 1:3
    subj_spec = squeeze(mean(ga_dec{c}.powspctrm(:, ch_idx_dec, :), 2, 'omitnan'));
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
    eb_dec(c) = shadedErrorBar(freqs_dec, m, se, 'lineProps', {'-'}, 'transparent', true); 
    set(eb_dec(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
    set(eb_dec(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
    set(eb_dec(c).edge(1), 'Color', colors(c, :));
    set(eb_dec(c).edge(2), 'Color', colors(c, :));
end
xlim([2 30]);
xlabel('Frequency [Hz]', 'FontSize', fontSize-4);
yline(0, '--')
ylim([-0.25 0.25])
ylabel('Power [dB]', 'FontSize', fontSize);
title('Decreasing slope', 'FontSize', fontSize);
% Legend with colored patch boxes
leg_p_dec = gobjects(1, 3);
for c = 1:3
    leg_p_dec(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_dec, {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', fontSize - 2);
set(gca, 'FontSize', fontSize);
box on

drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_powspctrm.png'));

%% Gaze density
disp(upper('Building gaze heatmaps from raw ET (no cache)...'))
gaze_dwell_time = build_sternberg_split_alphaLoads_gaze_from_raw_et(subjects, path);

% Unpack raw
allgazebase2 = gaze_dwell_time.allgazebase2;
allgazebase4 = gaze_dwell_time.allgazebase4;
allgazebase6 = gaze_dwell_time.allgazebase6;
allgazetasklate2 = gaze_dwell_time.allgazetasklate2;
allgazetasklate4 = gaze_dwell_time.allgazetasklate4;
allgazetasklate6 = gaze_dwell_time.allgazetasklate6;

% CBPT inputs spatially downsample (memory issues on server)
% gaze_cbpt_bins: scalar N → N×N; or [width height] screen px (e.g. [800 600]) → freq×time = [height width].
gaze_cbpt_bins = [200, 150];
allgazebase2_cbpt = downsample_gaze_cells_powspctrm(allgazebase2, gaze_cbpt_bins);
allgazebase4_cbpt = downsample_gaze_cells_powspctrm(allgazebase4, gaze_cbpt_bins);
allgazebase6_cbpt = downsample_gaze_cells_powspctrm(allgazebase6, gaze_cbpt_bins);
allgazetasklate2_cbpt = downsample_gaze_cells_powspctrm(allgazetasklate2, gaze_cbpt_bins);
allgazetasklate4_cbpt = downsample_gaze_cells_powspctrm(allgazetasklate4, gaze_cbpt_bins);
allgazetasklate6_cbpt = downsample_gaze_cells_powspctrm(allgazetasklate6, gaze_cbpt_bins);

%% Baseline (subtraction: tasklate - baseline)
disp('BASELINING GAZE DATA (SUBTRACTION)')
load2_gaze = cell(size(allgazebase2));
load4_gaze = cell(size(allgazebase4));
load6_gaze = cell(size(allgazebase6));
for subj = 1:length(allgazebase6)
    load2_gaze{subj} = allgazebase2{subj};
    load4_gaze{subj} = allgazebase4{subj};
    load6_gaze{subj} = allgazebase6{subj};
    load2_gaze{subj}.powspctrm = allgazetasklate2{subj}.powspctrm - allgazebase2{subj}.powspctrm;
    load4_gaze{subj}.powspctrm = allgazetasklate4{subj}.powspctrm - allgazebase4{subj}.powspctrm;
    load6_gaze{subj}.powspctrm = allgazetasklate6{subj}.powspctrm - allgazebase6{subj}.powspctrm;
end
load2_gaze_cbpt = downsample_gaze_cells_powspctrm(load2_gaze, gaze_cbpt_bins);
load4_gaze_cbpt = downsample_gaze_cells_powspctrm(load4_gaze, gaze_cbpt_bins);
load6_gaze_cbpt = downsample_gaze_cells_powspctrm(load6_gaze, gaze_cbpt_bins);

%% Compute gaze grand average
disp(upper('Computing gaze grand average...'))
cfg = [];
cfg.keepindividual = 'yes';
ga2_inc_gaze = ft_freqgrandaverage(cfg, load2_gaze{idx_inc});
ga4_inc_gaze = ft_freqgrandaverage(cfg, load4_gaze{idx_inc});
ga6_inc_gaze = ft_freqgrandaverage(cfg, load6_gaze{idx_inc});
ga2_dec_gaze = ft_freqgrandaverage(cfg, load2_gaze{idx_dec});
ga4_dec_gaze = ft_freqgrandaverage(cfg, load4_gaze{idx_dec});
ga6_dec_gaze = ft_freqgrandaverage(cfg, load6_gaze{idx_dec});

%% Plot Gaze TFRs
close all
clc
disp(upper('Plotting Gaze TFRs...'))
cfg = [];
cfg.figure = 'gcf';
% Strict symmetric max-abs z-limits for subtraction baseline maps
allGazeDiff = [ ...
    ga2_inc_gaze.powspctrm(:); ...
    ga4_inc_gaze.powspctrm(:); ...
    ga6_inc_gaze.powspctrm(:); ...
    ga2_dec_gaze.powspctrm(:); ...
    ga4_dec_gaze.powspctrm(:); ...
    ga6_dec_gaze.powspctrm(:) ];
allGazeDiff = allGazeDiff(isfinite(allGazeDiff));
zlimAbs = robust_diverging_zlim_abs(allGazeDiff, gaze_zlim_cfg);
cfg.zlim = [-zlimAbs zlimAbs];
figure('Position', fig_pos, 'Color', 'w');

% Gaze density relative to baseline (subtraction; tasklate - baseline)
subplot(3,2,1);
ft_singleplotTFR(cfg, ga2_inc_gaze);
title('Increasing slope: WM load 2', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,3);
ft_singleplotTFR(cfg, ga4_inc_gaze);
title('Increasing slope: WM load 4', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,5);
ft_singleplotTFR(cfg, ga6_inc_gaze);
title('Increasing slope: WM load 6', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,2);
ft_singleplotTFR(cfg, ga2_dec_gaze);
title('Decreasing slope: WM load 2', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,4);
ft_singleplotTFR(cfg, ga4_dec_gaze);
title('Decreasing slope: WM load 4', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,6);
ft_singleplotTFR(cfg, ga6_dec_gaze);
title('Decreasing slope: WM load 6', 'FontSize', fontSize);
ax = gca;
xlim(ax, [0 800]);
ylim(ax, [0 600]);
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
hold(ax, 'on');
plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

colormap(gcf, color_map);
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR.png'));

%% Compute CBPT stats: paired baseline-vs-task tests separately for each load (2/4/6) and each subgroup (increase vs decrease)
clc
% Six CBPT maps (3 loads x 2 groups); alpha = 0.05 per map (same as occ_split_stern_up_down_with_load.m).
alpha_cbpt_taskVsBase = 0.05;

cfg                  = [];
cfg.spmversion       = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_diff';  % paired mean difference, task vs baseline (same as supervisor script)
cfg.clusterthreshold = 'nonparametric_common';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = alpha_cbpt_taskVsBase;
cfg.correcttail      = 'alpha';
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
subj = numel(ga2_inc_gaze.powspctrm(:,1,1,1));
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
cfg_cbpt_gaze_taskVsBase_inc = cfg;
design_cbpt_gaze_taskVsBase_inc = design;
clc; disp(upper('[stat_inc_2 INC]')); [stat_inc_2] = ft_freqstatistics(cfg, allgazetasklate2_cbpt{idx_inc},allgazebase2_cbpt{idx_inc});
clc; disp(upper('[stat_inc_4 INC]')); [stat_inc_4] = ft_freqstatistics(cfg, allgazetasklate4_cbpt{idx_inc},allgazebase4_cbpt{idx_inc});
clc; disp(upper('[stat_inc_6 INC]')); [stat_inc_6] = ft_freqstatistics(cfg, allgazetasklate6_cbpt{idx_inc},allgazebase6_cbpt{idx_inc});
% Hide non-significant pixels in visualization
stat_inc_2.stat(stat_inc_2.mask==0) = NaN;
stat_inc_4.stat(stat_inc_4.mask==0) = NaN;
stat_inc_6.stat(stat_inc_6.mask==0) = NaN;

subj = numel(ga2_dec_gaze.powspctrm(:,1,1,1));
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
cfg_cbpt_gaze_taskVsBase_dec = cfg;
design_cbpt_gaze_taskVsBase_dec = design;
clc; disp(upper('[stat_dec_2 DEC]')); [stat_dec_2] = ft_freqstatistics(cfg, allgazetasklate2_cbpt{idx_dec},allgazebase2_cbpt{idx_dec});
clc; disp(upper('[stat_dec_4 DEC]')); [stat_dec_4] = ft_freqstatistics(cfg, allgazetasklate4_cbpt{idx_dec},allgazebase4_cbpt{idx_dec});
clc; disp(upper('[stat_dec_6 DEC]')); [stat_dec_6] = ft_freqstatistics(cfg, allgazetasklate6_cbpt{idx_dec},allgazebase6_cbpt{idx_dec});
% Hide non-significant pixels in visualization
stat_dec_2.stat(stat_dec_2.mask==0) = NaN;
stat_dec_4.stat(stat_dec_4.mask==0) = NaN;
stat_dec_6.stat(stat_dec_6.mask==0) = NaN;

%% Plot CBPT Raw Stats
stat_dec_2.cfg=[];
stat_dec_4.cfg=[];
stat_dec_6.cfg=[];
close all
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
allStatVals = [stat_inc_2.stat(:); stat_inc_4.stat(:); stat_inc_6.stat(:); ...
    stat_dec_2.stat(:); stat_dec_4.stat(:); stat_dec_6.stat(:)];
allStatVals = allStatVals(isfinite(allStatVals));
zlimAbs = robust_diverging_zlim_abs(allStatVals, gaze_zlim_cfg);
cfg.zlim = [-zlimAbs zlimAbs];
cfg.figure = 'gcf';

figure('Position', fig_pos, 'Color', 'w');
subplot(3,2,1);
ft_singleplotTFR(cfg, stat_inc_2);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Increasing slope: WM load 2', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,3);
ft_singleplotTFR(cfg,stat_inc_4);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Increasing slope: WM load 4', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,5);
ft_singleplotTFR(cfg,stat_inc_6);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Increasing slope: WM load 6', 'FontSize', fontSize, 'Interpreter', 'none')

% plot decreasing-slope group
subplot(3,2,2);
ft_singleplotTFR(cfg,stat_dec_2);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Decreasing slope: WM load 2', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,4);
ft_singleplotTFR(cfg,stat_dec_4);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Decreasing slope: WM load 4', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,6);
ft_singleplotTFR(cfg,stat_dec_6);
img = findobj(gca, 'Type', 'image');
if ~isempty(img) && isprop(img(1), 'CData')
    set(img(1), 'AlphaData', ~isnan(img(1).CData));
end
set(gca, 'Color', [1 1 1]);
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
hold on
plot(400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = 18;   % optional
title('Decreasing slope: WM load 6', 'FontSize', fontSize, 'Interpreter', 'none')
colormap(customcolormap_preset('red-white-blue'));

drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_cbpt_taskVsBase_incDec.png'));

%% CBPT Stats: repeated-measures (dependent-samples) F test across the three loads
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesFunivariate';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 1;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

subj = numel(ga2_inc_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3];
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_omnibus_inc = cfg;
[statF_gaze_inc] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_inc}, load4_gaze_cbpt{idx_inc}, load6_gaze_cbpt{idx_inc});
statF_gaze_inc.stat(statF_gaze_inc.mask==0)=0;% set everything not relevant to zero

cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesFunivariate';

cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;

cfg.clusterstatistic = 'maxsum';

cfg.tail             = 1;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

subj = numel(ga2_dec_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3];
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_omnibus_dec = cfg;
[statF_gaze_dec] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_dec}, load4_gaze_cbpt{idx_dec}, load6_gaze_cbpt{idx_dec});
statF_gaze_dec.stat(statF_gaze_dec.mask==0)=0;% set everything not relevant to zero

%% Plot F-stats
close all
cfg = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.zlim = [0 14];
cfg.figure = 'gcf';
cfg.colormap = 'YlOrRd';

figure('Position', [0 0 1512*0.5 982], 'Color', 'w');

% plot
ax1 = subplot(2,1,1);
ft_singleplotTFR(cfg, statF_gaze_inc);
title('Increasing slope: Omnibus load effect', 'FontSize', fontSize, 'Interpreter', 'none');
xlim(ax1, [0 800]);
ylim(ax1, [0 600]);
hold(ax1, 'on');
plot(ax1, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');

ax2 = subplot(2,1,2);
ft_singleplotTFR(cfg, statF_gaze_dec);
title('Decreasing slope: Omnibus load effect', 'FontSize', fontSize, 'Interpreter', 'none');
xlim(ax2, [0 800]);
ylim(ax2, [0 600]);
hold(ax2, 'on');
plot(ax2, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');

% fix colormaps
cmap = colormap(ax1);
cmap(1,:) = [1 1 1];
colormap(ax1, cmap);
colormap(ax2, cmap);

% format
set(ax1, 'FontSize', fontSize);
xlabel(ax1, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax1, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax1);
c.LineWidth = 1; c.FontSize = fontSize-2;
c.Ticks = [0 10]; c.Label.String = 'F-value'; c.Label.FontSize = fontSize-2;

set(ax2, 'FontSize', fontSize);
xlabel(ax2, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax2, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax2);
c.LineWidth = 1; c.FontSize = fontSize-2;
c.Ticks = [0 10]; c.Label.String = 'F-value'; c.Label.FontSize = fontSize-2;

drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_statF_omnibus.png'));

%% Test linear (and quadratic) trends across WM load
addpath(fileparts(mfilename('fullpath')));

% ---------------- Linear load trend ----------------
disp(upper('Computing linear load trend...'))
cfg                  = [];
cfg.method           = 'montecarlo';
% Contrast weights for ft_statfun_loadtrend are [-1 0 1] for loads [2 4 6].
% Positive t means increase from low to high load; negative t means decrease.
cfg.statistic        = 'ft_statfun_loadtrend';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

subj = numel(ga2_inc_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadtrend_inc = cfg;
[statT_gaze_inc] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_inc}, load4_gaze_cbpt{idx_inc}, load6_gaze_cbpt{idx_inc});
statT_gaze_inc.stat(statT_gaze_inc.mask==0)=0;

cfg                  = [];
cfg.method           = 'montecarlo';
% Contrast weights for ft_statfun_loadtrend are [-1 0 1] for loads [2 4 6].
% Positive t means increase from low to high load; negative t means decrease.
cfg.statistic        = 'ft_statfun_loadtrend';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
subj = numel(ga2_dec_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadtrend_dec = cfg;
[statT_gaze_dec] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_dec}, load4_gaze_cbpt{idx_dec}, load6_gaze_cbpt{idx_dec});
statT_gaze_dec.stat(statT_gaze_dec.mask==0)=0;
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.zlim = [-4 4];
cfg.xlim = [0 800];
cfg.ylim = [0 600];
cfg.figure = 'gcf';
figure('Position', [0 0 1512/2 982], 'Color', 'w');
subplot(2,1,1); ft_singleplotTFR(cfg,statT_gaze_inc); title('Increasing slope: Linear load trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;

subplot(2,1,2); ft_singleplotTFR(cfg,statT_gaze_dec); title('Decreasing slope: Linear load trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;
colormap(gcf, color_map); drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadTrend.png'));

%% ---------------- Quadratic load trend ----------------
disp(upper('Computing quadratic load trend...'))
cfg                  = [];
cfg.method           = 'montecarlo';
% Contrast weights for ft_statfun_loadquadratic are [1 -2 1] for loads [2 4 6].
% Positive t means U-shape; negative t means inverted U-shape.
cfg.statistic        = 'ft_statfun_loadquadratic';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

subj = numel(ga2_inc_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadquadratic_inc = cfg;
[statT_gaze_quad_inc] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_inc}, load4_gaze_cbpt{idx_inc}, load6_gaze_cbpt{idx_inc});
statT_gaze_quad_inc.stat(statT_gaze_quad_inc.mask==0)=0;

cfg                  = [];
cfg.method           = 'montecarlo';
% Contrast weights for ft_statfun_loadquadratic are [1 -2 1] for loads [2 4 6].
% Positive t means U-shape; negative t means inverted U-shape.
cfg.statistic        = 'ft_statfun_loadquadratic';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

subj = numel(ga2_dec_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadquadratic_dec = cfg;
[statT_gaze_quad_dec] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_dec}, load4_gaze_cbpt{idx_dec}, load6_gaze_cbpt{idx_dec});
statT_gaze_quad_dec.stat(statT_gaze_quad_dec.mask==0)=0;
%% plot
cfg = [];
cfg.zlim = [-0.5 0.5];
cfg.parameter = 'stat'; cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
cfg.xlim = [0 800];
cfg.ylim = [0 600];
cfg.figure = 'gcf';
figure('Position', [0 0 1512/2 982], 'Color', 'w');
subplot(2,1,1); ft_singleplotTFR(cfg,statT_gaze_quad_inc); title('Increasing slope: Quadratic load trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;

subplot(2,1,2); ft_singleplotTFR(cfg,statT_gaze_quad_dec); title('Decreasing slope: Quadratic load trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;
colormap(gcf, color_map); drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadQuadratic.png'));

%% Behavioral data
behav_file = fullfile(feat_dir, 'AOC_behavioral_matrix_sternberg.mat');
if ~isfile(behav_file)
    error('Missing: %s', behav_file);
end
load(behav_file);

nSubj = length(subjects);

% Preallocate
RT2 = nan(nSubj,1); RT4 = nan(nSubj,1); RT6 = nan(nSubj,1);
ACC2 = nan(nSubj,1); ACC4 = nan(nSubj,1); ACC6 = nan(nSubj,1);

for i = 1:nSubj

    subj_id = str2double(subjects{i}); % assumes subjects like '301'

    idx = [behav_data_sternberg.ID] == subj_id;

    cond = [behav_data_sternberg(idx).Condition];
    RT   = [behav_data_sternberg(idx).ReactionTime];
    ACC  = [behav_data_sternberg(idx).Accuracy];

    RT2(i) = RT(cond==2);
    RT4(i) = RT(cond==4);
    RT6(i) = RT(cond==6);

    ACC2(i) = ACC(cond==2);
    ACC4(i) = ACC(cond==4);
    ACC6(i) = ACC(cond==6);
end
%%
% Reaction Time example
rt_inc_2 = RT2(idx_inc);
rt_inc_4 = RT4(idx_inc);
rt_inc_6 = RT6(idx_inc);

rt_dec_2 = RT2(idx_dec);
rt_dec_4 = RT4(idx_dec);
rt_dec_6 = RT6(idx_dec);

%%
figure('Position', [0 0 1512 982], 'Color', 'w');
clf;
cond_cols = colors(1:3, :);

% LEFT: alpha increase group
subplot(2,2,1);
plot_split_raincloud_triplet({rt_inc_2, rt_inc_4, rt_inc_6}, cond_labels, cond_cols, [0 1.5], ...
    'Reaction Time [s]', sprintf('Increase (N=%d)', sum(idx_inc)), fontSize);

% RIGHT: alpha decrease group
subplot(2,2,2);
plot_split_raincloud_triplet({rt_dec_2, rt_dec_4, rt_dec_6}, cond_labels, cond_cols, [0 1.5], ...
    'Reaction Time [s]', sprintf('Decrease (N=%d)', sum(idx_dec)), fontSize);
set(gcf,'color','w');

%% Accuracy example
acc_inc_2 = ACC2(idx_inc);
acc_inc_4 = ACC4(idx_inc);
acc_inc_6 = ACC6(idx_inc);

acc_dec_2 = ACC2(idx_dec);
acc_dec_4 = ACC4(idx_dec);
acc_dec_6 = ACC6(idx_dec);

% remove NaNs (important!)
acc_inc_2 = acc_inc_2(~isnan(acc_inc_2)); acc_inc_4 = acc_inc_4(~isnan(acc_inc_4)); acc_inc_6 = acc_inc_6(~isnan(acc_inc_6));
acc_dec_2 = acc_dec_2(~isnan(acc_dec_2)); acc_dec_4 = acc_dec_4(~isnan(acc_dec_4)); acc_dec_6 = acc_dec_6(~isnan(acc_dec_6));

%%
% ================= LEFT: alpha increase =================
subplot(2,2,3);
plot_split_raincloud_triplet({acc_inc_2, acc_inc_4, acc_inc_6}, cond_labels, cond_cols, [50 110], ...
    'Accuracy [%]', sprintf('Increase (N=%d)', sum(idx_inc)), fontSize);

% ================= RIGHT: alpha decrease =================
subplot(2,2,4);
plot_split_raincloud_triplet({acc_dec_2, acc_dec_4, acc_dec_6}, cond_labels, cond_cols, [50 110], ...
    'Accuracy [%]', sprintf('Decrease (N=%d)', sum(idx_dec)), fontSize);

set(gcf,'color','w');
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_RT_ACC.png'));

%% Gaze deviation rainclouds (Increase vs Decrease; matched AlphaAmpRed aesthetics)
DEV2 = nan(nSubj,1); DEV4 = nan(nSubj,1); DEV6 = nan(nSubj,1);
for i = 1:nSubj
    subj_id = str2double(subjects{i});
    subj_rows = T_merged(T_merged.ID == subj_id, :);
    if isempty(subj_rows)
        continue
    end
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            dev_val = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
            if c == 1
                DEV2(i) = dev_val;
            elseif c == 2
                DEV4(i) = dev_val;
            else
                DEV6(i) = dev_val;
            end
        end
    end
end

dev_inc = [DEV2(idx_inc), DEV4(idx_inc), DEV6(idx_inc)];
dev_dec = [DEV2(idx_dec), DEV4(idx_dec), DEV6(idx_dec)];
plot_gaze_deviation_raincloud_split(dev_inc, dev_dec, cond_labels, colors, ...
    sum(idx_inc), sum(idx_dec), fig_dir, fig_pos, fontSize);

%% Subject-level alpha–gaze coupling over WM load (scalar + LME; no spatial CBPT)
% Tests whether alpha power vs load (slope) covaries inversely with gaze deviation vs load
% on the same [2,4,6] design. Complements voxel-wise regression CBPT above.
disp(upper('Alpha–gaze coupling over WM load (subject-level models)...'))

gaze_dev_slope = nan(nSubj, 1);
for ii = 1:nSubj
    ydev = [DEV2(ii), DEV4(ii), DEV6(ii)];
    if all(isfinite(ydev))
        Xmat = [ones(3, 1), [2; 4; 6]];
        b = Xmat \ ydev(:);
        gaze_dev_slope(ii) = b(2);
    end
end

mask_couple = isfinite(slope(:)) & isfinite(gaze_dev_slope(:));
% Figure and coupling stats here: alpha increase/decrease groups only (exclude intermediate).
mask_couple_fig = mask_couple & (idx_inc | idx_dec);
n_couple = sum(mask_couple);
n_couple_fig = sum(mask_couple_fig);
mI = mask_couple_fig & idx_inc;
mD = mask_couple_fig & idx_dec;
n_inc_fig = sum(mI);
n_dec_fig = sum(mD);
r_p_inc = NaN; p_p_inc = NaN; r_s_inc = NaN; p_s_inc = NaN; p_pearson_inv_inc = NaN;
r_p_dec = NaN; p_p_dec = NaN; r_s_dec = NaN; p_s_dec = NaN; p_pearson_inv_dec = NaN;

fprintf('\n========== SUBJECT-LEVEL ALPHA–GAZE COUPLING (WM LOAD 2–6) ==========\n');
fprintf('Gaze deviation slope: OLS of GazeDeviationFullBL on load [2,4,6] (slope = b_load).\n');
fprintf('Alpha slope: same regression on occipital alpha power (same loads).\n');

mask_gaze_fin = isfinite(gaze_dev_slope);
n_gaze_fin = sum(mask_gaze_fin);
fprintf('\n--- Descriptive statistics ---\n');
fprintf('N with finite gaze deviation slope: %d\n', n_gaze_fin);
if n_gaze_fin > 0
    gs_all = gaze_dev_slope(mask_gaze_fin);
    fprintf('  Gaze dev slope: min=%.5g, max=%.5g, mean=%.5g, SD=%.5g, median=%.5g\n', ...
        min(gs_all), max(gs_all), mean(gs_all), std(gs_all), median(gs_all));
end
fprintf('N with finite alpha and gaze slopes (any): %d; increase/decrease groups (figure / coupling): %d\n', ...
    n_couple, n_couple_fig);
if n_couple_fig > 0
    sa_d = slope(mask_couple_fig);
    sg_d = gaze_dev_slope(mask_couple_fig);
    fprintf('  Alpha slope (coupling N): min=%.5g, max=%.5g, mean=%.5g, SD=%.5g\n', ...
        min(sa_d), max(sa_d), mean(sa_d), std(sa_d));
    fprintf('  Gaze dev slope (coupling N): min=%.5g, max=%.5g, mean=%.5g, SD=%.5g\n', ...
        min(sg_d), max(sg_d), mean(sg_d), std(sg_d));
    n_inc_c = sum(mask_couple_fig & idx_inc);
    n_dec_c = sum(mask_couple_fig & idx_dec);
    fprintf('  Counts in this sample: alpha increase=%d, alpha decrease=%d\n', n_inc_c, n_dec_c);
end

fprintf('\n--- Within-group correlation (alpha slope vs gaze dev slope), per group ---\n');
fprintf('(Matches separate regression lines on the scatter; need N >= 3 per group for r and p.)\n');
if n_inc_fig >= 3
    sa_i = slope(mI); sg_i = gaze_dev_slope(mI);
    [r_p_inc, p_p_inc] = corr(sa_i, sg_i, 'Type', 'Pearson', 'rows', 'complete');
    [r_s_inc, p_s_inc] = corr(sa_i, sg_i, 'Type', 'Spearman', 'rows', 'complete');
    df_i = n_inc_fig - 2;
    if df_i > 0 && abs(r_p_inc) < 1 - eps
        t_ri = r_p_inc * sqrt(df_i / max(eps, 1 - r_p_inc^2));
        p_pearson_inv_inc = tcdf(t_ri, df_i);
    end
    fprintf('Alpha increase (N=%d): Pearson r=%.5f, p=%.6g; H1 inv. coupling (rho<0): p=%.6g; Spearman rho=%.5f, p=%.6g\n', ...
        n_inc_fig, r_p_inc, p_p_inc, p_pearson_inv_inc, r_s_inc, p_s_inc);
else
    fprintf('Alpha increase: N=%d (skipped correlation; need N>=3)\n', n_inc_fig);
end
if n_dec_fig >= 3
    sa_dv = slope(mD); sg_dv = gaze_dev_slope(mD);
    [r_p_dec, p_p_dec] = corr(sa_dv, sg_dv, 'Type', 'Pearson', 'rows', 'complete');
    [r_s_dec, p_s_dec] = corr(sa_dv, sg_dv, 'Type', 'Spearman', 'rows', 'complete');
    df_d = n_dec_fig - 2;
    if df_d > 0 && abs(r_p_dec) < 1 - eps
        t_rd = r_p_dec * sqrt(df_d / max(eps, 1 - r_p_dec^2));
        p_pearson_inv_dec = tcdf(t_rd, df_d);
    end
    fprintf('Alpha decrease (N=%d): Pearson r=%.5f, p=%.6g; H1 inv. coupling (rho<0): p=%.6g; Spearman rho=%.5f, p=%.6g\n', ...
        n_dec_fig, r_p_dec, p_p_dec, p_pearson_inv_dec, r_s_dec, p_s_dec);
else
    fprintf('Alpha decrease: N=%d (skipped correlation; need N>=3)\n', n_dec_fig);
end
if n_inc_fig >= 2
    p_ols_inc = polyfit(slope(mI), gaze_dev_slope(mI), 1);
    fprintf('OLS, alpha increase: gaze_dev_slope = %.5g + %.5g * alpha_slope\n', p_ols_inc(2), p_ols_inc(1));
end
if n_dec_fig >= 2
    p_ols_dec = polyfit(slope(mD), gaze_dev_slope(mD), 1);
    fprintf('OLS, alpha decrease: gaze_dev_slope = %.5g + %.5g * alpha_slope\n', p_ols_dec(2), p_ols_dec(1));
end
if n_couple_fig < 3
    fprintf('Note: total N in coupling sample = %d (very small combined sample).\n', n_couple_fig);
end

% --- Figures: top = alpha vs gaze scatter; bottom = alpha slope histogram (inclusion-style) + gaze slopes (one bar per participant) ---
% (mask_gaze_fin / n_gaze_fin defined in descriptive block above)
col_inc = [0.8 0 0];
col_dec = [0 0 0.8];
col_flat = [0.5 0.5 0.5];
mu_gs = nan(1, 2);
sem_gs = nan(1, 2);
n_gs = zeros(1, 2);
for gi = 1:2
    if gi == 1
        gmask = mask_couple_fig & idx_inc;
    else
        gmask = mask_couple_fig & idx_dec;
    end
    vals = gaze_dev_slope(gmask);
    vals = vals(isfinite(vals));
    n_gs(gi) = numel(vals);
    if n_gs(gi) >= 1
        mu_gs(gi) = mean(vals);
        if n_gs(gi) >= 2
            sem_gs(gi) = std(vals) / sqrt(n_gs(gi));
        else
            sem_gs(gi) = 0;
        end
    end
end

fprintf('\n--- Building alpha–gaze coupling figures ---\n');
figure('Position', fig_pos, 'Color', 'w');
axS = subplot(2, 2, [1 2]);
hold(axS, 'on');
h1 = scatter(axS, NaN, NaN, 40, col_inc, 'filled', 'MarkerEdgeColor', col_inc * 0.35);
h2 = scatter(axS, NaN, NaN, 40, col_dec, 'filled', 'MarkerEdgeColor', col_dec * 0.35);
idx_plot = find(mask_couple_fig);
for ii = idx_plot(:)'
    if idx_inc(ii)
        scatter(axS, slope(ii), gaze_dev_slope(ii), 36, col_inc, 'filled', 'MarkerFaceAlpha', 0.78, ...
            'MarkerEdgeColor', col_inc * 0.35, 'LineWidth', 0.45);
    elseif idx_dec(ii)
        scatter(axS, slope(ii), gaze_dev_slope(ii), 36, col_dec, 'filled', 'MarkerFaceAlpha', 0.78, ...
            'MarkerEdgeColor', col_dec * 0.35, 'LineWidth', 0.45);
    end
end
if sum(mI) >= 2
    sa_i = slope(mI);
    p_inc = polyfit(sa_i, gaze_dev_slope(mI), 1);
    xl_i = [min(sa_i), max(sa_i)];
    pad_i = 0.05 * max(abs(diff(xl_i)), eps);
    xl_i = [xl_i(1) - pad_i, xl_i(2) + pad_i];
    plot(axS, xl_i, polyval(p_inc, xl_i), '-', 'Color', col_inc * 0.55, 'LineWidth', 1.65);
end
if sum(mD) >= 2
    sa_d = slope(mD);
    p_dec = polyfit(sa_d, gaze_dev_slope(mD), 1);
    xl_d = [min(sa_d), max(sa_d)];
    pad_d = 0.05 * max(abs(diff(xl_d)), eps);
    xl_d = [xl_d(1) - pad_d, xl_d(2) + pad_d];
    plot(axS, xl_d, polyval(p_dec, xl_d), '-', 'Color', col_dec * 0.55, 'LineWidth', 1.65);
end
xline(axS, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
yline(axS, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
hold(axS, 'off');
sa_sc = slope(mask_couple_fig);
sg_sc = gaze_dev_slope(mask_couple_fig);
sx_sc = max(abs(sa_sc(isfinite(sa_sc))));
sy_sc = max(abs(sg_sc(isfinite(sg_sc))));
if ~isfinite(sx_sc) || sx_sc <= 0, sx_sc = eps; end
if ~isfinite(sy_sc) || sy_sc <= 0, sy_sc = eps; end
xlim(axS, [-1 1] * sx_sc * 1.05);
ylim(axS, [-1 1] * sy_sc * 1.05);
xlabel(axS, 'Alpha power slope (a.u. per item)', 'FontSize', fontSize - 1);
ylabel(axS, 'Gaze deviation slope (a.u. per item)', 'FontSize', fontSize - 1);
title(axS, ' ', 'FontSize', fontSize, 'Interpreter', 'tex');
set(axS, 'FontSize', fontSize - 2);
legend(axS, [h1, h2], {'Alpha increase', 'Alpha decrease'}, ...
    'Location', 'best', 'Box', 'off', 'Interpreter', 'tex', 'FontSize', fontSize - 3);
txt = {};
if n_inc_fig >= 3
    txt{end+1} = sprintf(['Alpha increase (N=%d): Pearson \\it r\\rm = %.3f, \\it p\\rm = %.4g; ', ...
        'H_1: \\rho < 0: \\it p\\rm = %.4g; Spearman \\rho = %.3f, \\it p\\rm = %.4g'], ...
        n_inc_fig, r_p_inc, p_p_inc, p_pearson_inv_inc, r_s_inc, p_s_inc);
else
    txt{end+1} = sprintf('Alpha increase: N = %d (correlation requires N \\geq 3)', n_inc_fig);
end
if n_dec_fig >= 3
    txt{end+1} = sprintf(['Alpha decrease (N=%d): Pearson \\it r\\rm = %.3f, \\it p\\rm = %.4g; ', ...
        'H_1: \\rho < 0: \\it p\\rm = %.4g; Spearman \\rho = %.3f, \\it p\\rm = %.4g'], ...
        n_dec_fig, r_p_dec, p_p_dec, p_pearson_inv_dec, r_s_dec, p_s_dec);
else
    txt{end+1} = sprintf('Alpha decrease: N = %d (correlation requires N \\geq 3)', n_dec_fig);
end
text(axS, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', fontSize - 4, 'Interpreter', 'tex');

fprintf('\n--- Gaze deviation slope by group (same sample as scatter) ---\n');
fprintf('  group                n    mean(gaze slope)       SEM\n');
tail_lbl = {'alpha increase', 'alpha decrease'};
for gi = 1:2
    if n_gs(gi) < 1
        fprintf('  %-20s   --          --              --\n', tail_lbl{gi});
    else
        fprintf('  %-20s  %3d  %18.5g  %12.5g\n', tail_lbl{gi}, n_gs(gi), mu_gs(gi), sem_gs(gi));
    end
end

% Bottom left: alpha slope histogram (same binning as inclusion figure at start of script)
axAlpha = subplot(2, 2, 3);
hold(axAlpha, 'on');
if numel(bin_edges) >= 2 && all(isfinite(bin_edges(:)))
    histogram(axAlpha, slope(idx_inc), 'BinEdges', bin_edges, 'FaceColor', col_inc, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_dec), 'BinEdges', bin_edges, 'FaceColor', col_dec, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_flat), 'BinEdges', bin_edges, 'FaceColor', col_flat, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    lim_a = max(abs([bin_edges(1), bin_edges(end)]));
else
    histogram(axAlpha, slope(idx_inc), 'FaceColor', col_inc, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_dec), 'FaceColor', col_dec, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_flat), 'FaceColor', col_flat, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    all_s = slope(isfinite(slope));
    lim_a = max(abs([min(all_s), max(all_s)]));
end
xline(axAlpha, 0, 'k-', 'LineWidth', 1);
if isfinite(t1), xline(axAlpha, t1, 'k--', 'LineWidth', 1.5); end
if isfinite(t2), xline(axAlpha, t2, 'k--', 'LineWidth', 1.5); end
hold(axAlpha, 'off');
if isfinite(lim_a) && lim_a > 0
    xlim(axAlpha, [-lim_a lim_a]);
end
xlabel(axAlpha, 'Alpha power slope', 'FontSize', fontSize - 1);
ylabel(axAlpha, 'Participants', 'FontSize', fontSize - 1);
title(axAlpha, 'Alpha power slope across WM load', 'FontSize', fontSize - 1, 'Interpreter', 'none');
set(axAlpha, 'FontSize', fontSize - 2);

% Bottom right: one bar per participant at x = gaze deviation slope (sorted); uniform bar height; color = alpha split
axGaze = subplot(2, 2, 4);
idx_g = find(mask_couple_fig);
gv_raw = gaze_dev_slope(idx_g);
[~, ord] = sort(gv_raw);
idx_ord = idx_g(ord);
gv_sorted = gv_raw(ord);
n_b = numel(idx_ord);
C_b = zeros(n_b, 3);
for k = 1:n_b
    ii = idx_ord(k);
    if idx_inc(ii)
        C_b(k, :) = col_inc;
    else
        C_b(k, :) = col_dec;
    end
end
if n_b >= 1
    h_bar = ones(n_b, 1);
    b = bar(axGaze, gv_sorted, h_bar, 'BarWidth', 0.92, 'EdgeColor', 'none');
    b.FaceColor = 'flat';
    b.CData = C_b;
    xline(axGaze, 0, 'k-', 'LineWidth', 1);
    mg = max(abs(gv_sorted(isfinite(gv_sorted))));
    if isfinite(mg) && mg > 0
        xlim(axGaze, [-mg mg] * 1.05);
    end
    ylim(axGaze, [0 1.12]);
    set(axGaze, 'YTick', []);
    hold(axGaze, 'on');
    plot(axGaze, NaN, NaN, 's', 'MarkerFaceColor', col_inc, 'MarkerEdgeColor', col_inc * 0.35, 'MarkerSize', 9);
    plot(axGaze, NaN, NaN, 's', 'MarkerFaceColor', col_dec, 'MarkerEdgeColor', col_dec * 0.35, 'MarkerSize', 9);
    legend(axGaze, {'Alpha increase', 'Alpha decrease'}, 'Location', 'northeast', ...
        'Box', 'off', 'FontSize', fontSize - 4, 'Interpreter', 'none');
    hold(axGaze, 'off');
end
xlabel(axGaze, 'Gaze deviation slope (a.u. per item)', 'FontSize', fontSize - 1);
ylabel(axGaze, '', 'FontSize', fontSize - 1);
title(axGaze, 'Gaze deviation slopes', 'FontSize', fontSize - 1, 'Interpreter', 'none');
set(axGaze, 'FontSize', fontSize - 2);

sgtitle('Subject-level Alpha–gaze coupling over WM load', 'FontSize', fontSize + 2, ...
    'FontWeight', 'bold', 'Interpreter', 'tex');
drawnow;
fig_coupling_path = fullfile(fig_dir, 'AOC_split_AlphaLoads_alpha_gaze_slope_coupling.png');
saveas(gcf, fig_coupling_path);
fprintf('Saved figure: %s\n', fig_coupling_path);

% Mixed model: does the within-subject gaze~load slope depend on alpha slope?
Subject_ag = [];
Load_ag = [];
Dev_ag = [];
Alpha_ag = [];
idx_ag = find(all(isfinite([DEV2, DEV4, DEV6]), 2));
mu_alpha = mean(slope(idx_ag));
for k = 1:numel(idx_ag)
    ii = idx_ag(k);
    dev_vals = [DEV2(ii), DEV4(ii), DEV6(ii)];
    for c = 1:3
        Subject_ag = [Subject_ag; ii];
        Load_ag = [Load_ag; cond_vals(c)];
        Dev_ag = [Dev_ag; dev_vals(c)];
        Alpha_ag = [Alpha_ag; slope(ii) - mu_alpha];
    end
end
tbl_ag = table(Subject_ag, Load_ag, Dev_ag, Alpha_ag, ...
    'VariableNames', {'Subject', 'Load', 'Dev', 'AlphaSlope_c'});
tbl_ag.Subject = categorical(tbl_ag.Subject);
tbl_ag.Load = (tbl_ag.Load - 4) / 2;

disp('--- LME: Dev ~ Load * AlphaSlope (between-subject) + (Load|Subject) ---')
fprintf('Load coded as (Load-4)/2 (0 = load 4; +1 = load 6).\n');
fprintf('AlphaSlope_c is mean-centered across subjects with complete gaze (3 loads / subject).\n');
try
    n_subj_ag = numel(categories(categorical(tbl_ag.Subject)));
    n_obs_ag = height(tbl_ag);
    fprintf('LME data: %d subjects with complete gaze deviation at loads 2/4/6, %d rows total.\n', ...
        n_subj_ag, n_obs_ag);
    fprintf('Grand mean alpha slope (subtracted as AlphaSlope_c): %.6g\n', mu_alpha);

    [lme_ag, lme_ag_formula] = fitlme_with_random_slope_fallback(tbl_ag, 'Dev ~ Load * AlphaSlope_c');
    fprintf('Fitted formula: %s\n', lme_ag_formula);
    fprintf('\nANOVA (marginal tests):\n');
    disp(anova(lme_ag))

    fprintf('Fixed-effects coefficients:\n');
    disp(lme_ag.Coefficients)

    coef_ag = lme_ag.Coefficients;
    idx_int = contains(coef_ag.Name, ':') & (contains(coef_ag.Name, 'Load') & contains(coef_ag.Name, 'AlphaSlope'));
    ix_int = find(idx_int, 1, 'first');
    if ~isempty(ix_int)
        row = coef_ag(ix_int, :);
        fprintf('Primary test (Load:AlphaSlope_c): inverse coupling if estimate < 0.\n');
        fprintf('  Estimate=%.6g, SE=%.6g, t=%.4g, DF=%.4g, p=%.6g\n', ...
            row.Estimate(1), row.SE(1), row.tStat(1), row.DF(1), row.pValue(1));
    else
        fprintf('Note: Load:AlphaSlope interaction row not found in coefficient table (check names).\n');
    end
    try
        fprintf('\n95%% CI for fixed effects (coefCI):\n');
        disp(coefCI(lme_ag));
    catch %#ok<CTCH>
    end
    fprintf('\n========== END SUBJECT-LEVEL ALPHA–GAZE COUPLING OUTPUT ==========\n\n');
end

%% LME (RT / ACC / gaze deviation; categorical load)
nSubj = length(subjects);


RT2 = nan(nSubj,1); RT4 = nan(nSubj,1); RT6 = nan(nSubj,1);
ACC2 = nan(nSubj,1); ACC4 = nan(nSubj,1); ACC6 = nan(nSubj,1);

for i = 1:nSubj

    subj_id = str2double(subjects{i});

    idx = [behav_data_sternberg.ID] == subj_id;

    if sum(idx)==0
        continue
    end

    cond = [behav_data_sternberg(idx).Condition];
    RT   = [behav_data_sternberg(idx).ReactionTime];
    ACC  = [behav_data_sternberg(idx).Accuracy];

    if sum(cond==2)==0 || sum(cond==4)==0 || sum(cond==6)==0
        continue
    end

    RT2(i) = RT(cond==2);
    RT4(i) = RT(cond==4);
    RT6(i) = RT(cond==6);

    ACC2(i) = ACC(cond==2);
    ACC4(i) = ACC(cond==4);
    ACC6(i) = ACC(cond==6);
end

%% ================= BUILD TABLE =================
Subject = [];
Load = [];
Group = [];
RT = [];
ACC = [];

for i = 1:nSubj

    if idx_inc(i)
        g = 1;
    elseif idx_dec(i)
        g = -1;
    else
        continue
    end

    if any(isnan([RT2(i), RT4(i), RT6(i), ACC2(i), ACC4(i), ACC6(i)]))
        continue
    end

    Subject = [Subject; i; i; i];
    Load = [Load; 2; 4; 6];
    Group = [Group; g; g; g];

    RT = [RT; RT2(i); RT4(i); RT6(i)];
    ACC = [ACC; ACC2(i); ACC4(i); ACC6(i)];
end

tbl = table(Subject, Load, Group, RT, ACC);

tbl.Subject = categorical(tbl.Subject);
tbl.Group   = categorical(tbl.Group, [1 -1], {'inc' 'dec'});

%% ================= LME ANALYSIS =================
disp('--- RT model ---')
% Treat Load as categorical for a full 3-level repeated-measures effect.
tbl.Load = categorical(tbl.Load);
lme_RT = fitlme(tbl, 'RT ~ Load * Group + (1|Subject)');
disp(anova(lme_RT))

disp('--- ACC model ---')
lme_ACC = fitlme(tbl, 'ACC ~ Load * Group + (1|Subject)');
disp(anova(lme_ACC))

disp('--- Gaze deviation model ---')
Subject_dev = [];
Load_dev = [];
Group_dev = [];
Dev = [];
for i = 1:nSubj
    if idx_inc(i)
        g = 1;
    elseif idx_dec(i)
        g = -1;
    else
        continue
    end
    dev_vals = [DEV2(i), DEV4(i), DEV6(i)];
    for c = 1:3
        if ~isfinite(dev_vals(c))
            continue
        end
        Subject_dev = [Subject_dev; i];
        Load_dev = [Load_dev; cond_vals(c)];
        Group_dev = [Group_dev; g];
        Dev = [Dev; dev_vals(c)];
    end
end
tbl_dev = table(Subject_dev, Load_dev, Group_dev, Dev, ...
    'VariableNames', {'Subject', 'Load', 'Group', 'Dev'});
tbl_dev.Subject = categorical(tbl_dev.Subject);
tbl_dev.Load = categorical(tbl_dev.Load);
tbl_dev.Group = categorical(tbl_dev.Group, [1 -1], {'inc' 'dec'});
[lme_DEV, dev_formula] = fitlme_with_random_slope_fallback(tbl_dev, 'Dev ~ Load * Group');
fprintf('Gaze deviation model formula: %s\n', dev_formula);
disp(anova(lme_DEV))

%{
tbl_dev_num = tbl_dev; % keep numeric load coding for linear-trend tests

disp('--- Gaze linear-trend hypothesis model ---')
tbl_dev_trend = tbl_dev_num;
tbl_dev_trend.Subject = categorical(tbl_dev_trend.Subject);
tbl_dev_trend.Group = categorical(tbl_dev_trend.Group, [1 -1], {'inc' 'dec'});
% Center and scale load so one unit corresponds to one WM load step (2 items).
tbl_dev_trend.Load = (tbl_dev_trend.Load - 4) / 2;
[lme_DEV_trend, dev_trend_formula] = fitlme_with_random_slope_fallback(tbl_dev_trend, 'Dev ~ Load * Group');
fprintf('Gaze linear-trend model formula: %s\n', dev_trend_formula);
disp(anova(lme_DEV_trend))

% Directional follow-up tests for the preregistered slope-direction hypothesis:
% inc group slope < 0; dec group slope > 0.
tbl_dev_inc = tbl_dev_trend(tbl_dev_trend.Group == 'inc', :);
tbl_dev_dec = tbl_dev_trend(tbl_dev_trend.Group == 'dec', :);
[lme_DEV_inc, dev_inc_formula] = fitlme_with_random_slope_fallback(tbl_dev_inc, 'Dev ~ Load');
[lme_DEV_dec, dev_dec_formula] = fitlme_with_random_slope_fallback(tbl_dev_dec, 'Dev ~ Load');

coef_inc = lme_DEV_inc.Coefficients;
idx_load_inc = strcmp(coef_inc.Name, 'Load');
if ~any(idx_load_inc)
    error('Load coefficient not found in increase-group gaze model.');
end
est_inc = coef_inc.Estimate(idx_load_inc);
t_inc = coef_inc.tStat(idx_load_inc);
df_inc = lme_DEV_inc.DFE;
p_inc_two = coef_inc.pValue(idx_load_inc);
p_inc_one = tcdf(t_inc, df_inc); % H1: slope < 0

coef_dec = lme_DEV_dec.Coefficients;
idx_load_dec = strcmp(coef_dec.Name, 'Load');
if ~any(idx_load_dec)
    error('Load coefficient not found in decrease-group gaze model.');
end
est_dec = coef_dec.Estimate(idx_load_dec);
t_dec = coef_dec.tStat(idx_load_dec);
df_dec = lme_DEV_dec.DFE;
p_dec_two = coef_dec.pValue(idx_load_dec);
p_dec_one = 1 - tcdf(t_dec, df_dec); % H1: slope > 0

fprintf('Increase-group trend model: %s\n', dev_inc_formula);
fprintf('Increase group load slope (per +2 items): %.4g, t(%g)=%.4g, p(two-sided)=%.4g, p(one-sided, <0)=%.4g\n', ...
    est_inc, df_inc, t_inc, p_inc_two, p_inc_one);
fprintf('Decrease-group trend model: %s\n', dev_dec_formula);
fprintf('Decrease group load slope (per +2 items): %.4g, t(%g)=%.4g, p(two-sided)=%.4g, p(one-sided, >0)=%.4g\n', ...
    est_dec, df_dec, t_dec, p_dec_two, p_dec_one);


%% ================= TOST ANALYSIS =================
delta_RT = 0.05;   % 50 ms
delta_ACC = 5;     % 5%

% Use one independent value per subject (mean across loads) for between-group TOST.
RT_inc = mean([RT2(idx_inc), RT4(idx_inc), RT6(idx_inc)], 2, 'omitnan');
RT_dec  = mean([RT2(idx_dec),  RT4(idx_dec),  RT6(idx_dec)], 2, 'omitnan');

ACC_inc = mean([ACC2(idx_inc), ACC4(idx_inc), ACC6(idx_inc)], 2, 'omitnan');
ACC_dec  = mean([ACC2(idx_dec),  ACC4(idx_dec),  ACC6(idx_dec)], 2, 'omitnan');

fprintf('\n--- TOST (strict) ---\n')

[p1, p2, eq_RT] = tost_welch(RT_inc, RT_dec, delta_RT);
fprintf('RT equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_RT);

[p1, p2, eq_ACC] = tost_welch(ACC_inc, ACC_dec, delta_ACC);
fprintf('ACC equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_ACC);

%% ================= RT SENSITIVITY =================
fprintf('\n--- RT equivalence sensitivity ---\n')

[p1, p2, eq] = tost_welch(RT_inc, RT_dec, 0.05);
fprintf('RT delta=0.05: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(RT_inc, RT_dec, 0.1);
fprintf('RT delta=0.10: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(RT_inc, RT_dec, 0.15);
fprintf('RT delta=0.15: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);

%% ================= ACC ROBUSTNESS =================
fprintf('\n--- ACC robustness ---\n')

[p1, p2, eq] = tost_welch(ACC_inc, ACC_dec, 3);
fprintf('ACC delta=3: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(ACC_inc, ACC_dec, 5);
fprintf('ACC delta=5: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(ACC_inc, ACC_dec, 10);
fprintf('ACC delta=10: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
%}

%% Helper functions
function zlimAbs = robust_diverging_zlim_abs(v, cfg)
% Symmetric z-limit for diverging maps: min( prctile(|x|,p), median(|x|)+k*MAD(|x|) )
% on values above cfg.min_abs (sparse gaze grids: most bins ~0). Reduces inflation
% from rare extreme pixels vs a single high global percentile.
a = abs(v(isfinite(v)));
a_pos = a(a > cfg.min_abs);
if ~isempty(a_pos)
    a = a_pos;
end
if isempty(a)
    zlimAbs = cfg.fallback_abs;
    return
end
p_lim = prctile(a, cfg.prctile);
med_a = median(a);
% Default mad(): median abs deviation from median, scaled (~sigma for Gaussian)
mad_lim = med_a + cfg.k_mad * mad(a);
zlimAbs = min(p_lim, mad_lim);
if ~isfinite(zlimAbs) || zlimAbs <= 0
    zlimAbs = cfg.fallback_abs;
end
end

function plot_tfr_matrix_panel(subplot_idx, ga_data, cfg_sel, clim, fsz)
freq = ft_selectdata(cfg_sel, ga_data);
meanpow = squeeze(mean(freq.powspctrm, 1));
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(1, 40, 500);
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, ...
    tim_grid_interp, freq_grid_interp, 'spline');
subplot(3, 2, subplot_idx);
ft_plot_matrix(flip(pow_interp));
ax = gca; hold(ax, 'on');
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);
xticks(round(interp1(tim_interp, 1:numel(tim_interp), [-0.5 0 1 2], 'linear', 'extrap')));
xticklabels({'-0.5', '0', '1', '2'});
yticks([1 125 250 375]);
yticklabels({'40', '30', '20', '10'});
set(gca, 'FontSize', fsz);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
caxis(clim);
c = colorbar;
c.LineWidth = 1;
c.FontSize = fsz - 2;
c.Ticks = clim;
ylabel(c, 'dB');
end

function sig_label = getSigLabel(p)
if p < 0.001
    sig_label = '***';
elseif p < 0.01
    sig_label = '**';
elseif p < 0.05
    sig_label = '*';
else
    sig_label = '';
end
end

function plot_gaze_deviation_raincloud_split(dev_inc, dev_dec, cond_labels, colors, n_inc, n_dec, fig_dir, fig_pos, fsz)
all_vals = [dev_inc(:); dev_dec(:)];
all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals)
    ymax = 1;
else
    ymax = max(abs(all_vals));
    if ymax <= 0
        ymax = 1;
    end
    ymax = ymax * 1.15;
end
ylim_shared = [-ymax ymax];

figure('Position', fig_pos, 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile;
hold on
for c = 1:3
    draw_one_cloud(dev_inc(:, c), c, colors(c, :), 0.3, 96, 0.45);
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim(ylim_shared);
set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
title(sprintf('Increase (N=%d)', n_inc), 'Interpreter', 'none');
ylabel('Gaze deviation [%]', 'Interpreter', 'none');
box off

nexttile;
hold on
for c = 1:3
    draw_one_cloud(dev_dec(:, c), c, colors(c, :), 0.3, 96, 0.45);
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim(ylim_shared);
set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels, 'FontSize', fsz-2);
title(sprintf('Decrease (N=%d)', n_dec), 'Interpreter', 'none');
ylabel('Gaze deviation [%]', 'Interpreter', 'none');
box off

sgtitle('Gaze deviation', 'FontSize', fsz+2, 'Interpreter', 'none');
drawnow;
saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_raincloud_gazedeviation.png'));
close(gcf);
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha, density_offset, box_offset, dot_jitter_halfwidth)
y = yvals(isfinite(yvals));
if numel(y) < 3
    return
end
[f, xi] = ksdensity(y, 'NumPoints', 120);
if nargin < 7 || isempty(density_offset), density_offset = 0.08; end
if nargin < 8 || isempty(box_offset), box_offset = 0.03; end
if nargin < 9 || isempty(dot_jitter_halfwidth), dot_jitter_halfwidth = box_w / 2; end
if max(f) > 0
    f = f / max(f) * 0.35;
else
    f = zeros(size(f));
end
x_den = xpos - density_offset;
x_box = xpos + box_offset;
fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25);
q3 = prctile(y, 75);
med = median(y);
p5 = prctile(y, 5);
p95 = prctile(y, 95);
plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [x_box-box_w/2, q1, box_w, q3-q1], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = dot_jitter_halfwidth * 2 * (rand(numel(y),1)-0.5);
scatter(x_box + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

function plot_split_raincloud_triplet(data_cells, cond_labels, cond_cols, y_limits, y_lab, ttl, fsz)
hold on
for cc = 1:3
    draw_one_cloud(data_cells{cc}, cc, cond_cols(cc, :), 0.3, 96, 0.45);
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylabel(y_lab);
title(ttl);
set(gca, 'FontSize', fsz-2);
set(gca, 'XTick', 1:3, 'XTickLabel', cond_labels);
ylim(y_limits);
xlim([0.7 3.3]);
box off
end
%{
function [p1, p2, equivalent] = tost(x1, x2, delta, alpha)

if nargin < 4
    alpha = 0.05;
end

% remove NaNs
x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

n1 = length(x1);
n2 = length(x2);

m1 = mean(x1);
m2 = mean(x2);

s1 = var(x1);
s2 = var(x2);

% pooled SE
SE = sqrt(s1/n1 + s2/n2);

df = n1 + n2 - 2;

% TOST tests
t1 = (m1 - m2 + delta) / SE; % lower bound
t2 = (m1 - m2 - delta) / SE; % upper bound

p1 = 1 - tcdf(t1, df);
p2 = tcdf(t2, df);

equivalent = (p1 < alpha) && (p2 < alpha);
end
%}

% ----- Helper functions for gaze dwell heatmaps -----
function conds = parse_trialinfo_conds_for_sternberg_heatmap(trialinfo)
% Returns per-trial condition codes (22/24/26) from the trialinfo array.
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

function pos_cells = extractPosCellsForGazeWindow( ...
    gaze_x, gaze_y, trial_idx, t_series, t_win, min_valid_samples)
% Extracts per-trial [x;y] matrices for a given time window.
% NaN positions are allowed and ignored by `histcounts2`.
pos_cells = {};
if isempty(trial_idx)
    return
end

for ii = 1:numel(trial_idx)
    trl = trial_idx(ii);
    x = double(gaze_x{trl});
    y = double(gaze_y{trl});
    if isempty(x) || isempty(y) || numel(x) ~= numel(y)
        continue
    end
    x = x(:)'; y = y(:)';

    tt = linspace(t_series(1), t_series(2), numel(x));
    idx_win = tt >= t_win(1) & tt <= t_win(2);
    if ~any(idx_win)
        continue
    end

    xw = x(idx_win);
    yw = y(idx_win);

    valid = isfinite(xw) & isfinite(yw);
    if sum(valid) < min_valid_samples
        continue
    end

    xw(~valid) = NaN;
    yw(~valid) = NaN;

    pos_cells{end+1} = [xw; yw]; 
end
end

function freq_raw = computeGazeHeatmapFromPosCells( ...
    pos_cells, x_grid, y_grid, fs, smoothing)
% Implements the supervisor's `computeGazeHeatmap()` logic (raw + normalized dwell maps).
nTime = numel(x_grid) - 1;
nFreq = numel(y_grid) - 1;

% Empty fallback: return all-zero maps with correct FieldTrip structure.
if isempty(pos_cells)
    freq_raw = struct();
    freq_raw.powspctrm = zeros(1, nFreq, nTime);
    freq_raw.time = x_grid(2:end);
    freq_raw.freq = y_grid(2:end);
    freq_raw.label = {'et'};
    freq_raw.dimord = 'chan_freq_time';
    return
end

pos = horzcat(pos_cells{:});
binned = histcounts2(pos(1, :), pos(2, :), x_grid, y_grid);
dwell_time = binned / fs;
smoothed = imgaussfilt(dwell_time, smoothing);

freq_raw = struct();
freq_raw.powspctrm(1, :, :) = flipud(smoothed');
freq_raw.time = x_grid(2:end);
freq_raw.freq = y_grid(2:end);
freq_raw.label = {'et'};
freq_raw.dimord = 'chan_freq_time';

denom = sum(dwell_time(:));
if denom > 0
    norm_time = dwell_time / denom;
    norm_smooth = imgaussfilt(norm_time, smoothing);
else
end
end

function gaze_dwell_time = build_sternberg_split_alphaLoads_gaze_from_raw_et(subjects, base_path)

n_subj = length(subjects);
num_bins = 1000;
smoothing = 5;
fsample = 500;
x_grid = linspace(0, 800, num_bins + 1);
y_grid = linspace(0, 600, num_bins + 1);

allgazebase2 = cell(1, n_subj); allgazetasklate2 = cell(1, n_subj);
allgazebase4 = cell(1, n_subj); allgazetasklate4 = cell(1, n_subj);
allgazebase6 = cell(1, n_subj); allgazetasklate6 = cell(1, n_subj);

for subj = 1:n_subj
    clc
    fprintf('Building dwell time for Subject %d/%d (%s)\n', subj, n_subj, subjects{subj});
    et_file = fullfile(base_path, subjects{subj}, 'gaze', 'dataET_sternberg.mat');
    if ~isfile(et_file)
        et_file = fullfile(base_path, subjects{subj}, 'gaze', 'dataET_sternberg');
    end
    if ~isfile(et_file)
        error('Missing ET file for subject %s: %s', subjects{subj}, et_file);
    end

    tmp = load(et_file);
    et = select_et_struct(tmp, subjects{subj});

    n_trials = numel(et.trial);
    trialinfo_vec = align_trialinfo_to_trials(et.trialinfo, n_trials, subjects{subj});
    if isempty(trialinfo_vec)
        error('No trials available after trialinfo alignment for subject %s.', subjects{subj});
    end
    idx2 = find(trialinfo_vec == 22);
    idx4 = find(trialinfo_vec == 24);
    idx6 = find(trialinfo_vec == 26);

    allgazebase2{subj} = extract_gaze_window(et, idx2, [-0.5 -0.25], x_grid, y_grid, fsample, smoothing);
    allgazebase4{subj} = extract_gaze_window(et, idx4, [-0.5 -0.25], x_grid, y_grid, fsample, smoothing);
    allgazebase6{subj} = extract_gaze_window(et, idx6, [-0.5 -0.25], x_grid, y_grid, fsample, smoothing);

    trial_tmin = nan(1, n_trials);
    trial_tmax = nan(1, n_trials);
    for tr = 1:n_trials
        tt = et.time{tr};
        if isempty(tt)
            continue
        end
        trial_tmin(tr) = tt(1);
        trial_tmax(tr) = tt(end);
    end
    allgazetasklate2{subj} = extract_gaze_window(et, idx2, [1 3], x_grid, y_grid, fsample, smoothing);
    allgazetasklate4{subj} = extract_gaze_window(et, idx4, [1 3], x_grid, y_grid, fsample, smoothing);
    allgazetasklate6{subj} = extract_gaze_window(et, idx6, [1 3], x_grid, y_grid, fsample, smoothing);
end

gaze_dwell_time = struct();
gaze_dwell_time.allgazebase2 = allgazebase2;
gaze_dwell_time.allgazetasklate2 = allgazetasklate2;
gaze_dwell_time.allgazebase4 = allgazebase4;
gaze_dwell_time.allgazetasklate4 = allgazetasklate4;
gaze_dwell_time.allgazebase6 = allgazebase6;
gaze_dwell_time.allgazetasklate6 = allgazetasklate6;
end

function et = select_et_struct(tmp, subj_label)
if isfield(tmp, 'dataETlong')
    et = tmp.dataETlong;
elseif isfield(tmp, 'dataet')
    et = tmp.dataet;
elseif isfield(tmp, 'dataET')
    et = tmp.dataET;
else
    error('No ET data structure found in dataET_sternberg for subject %s.', subj_label);
end

if ~isfield(et, 'trial') || ~isfield(et, 'trialinfo') || ~isfield(et, 'time')
    error('ET structure for subject %s is missing trial/trialinfo/time fields.', subj_label);
end
end

function freq_raw = extract_gaze_window(et, trial_idx, latency, x_grid, y_grid, fsample, smoothing)
cfg = [];
cfg.channel = {'L-GAZE-X', 'L-GAZE-Y'};
cfg.latency = latency;
cfg.trials = trial_idx;
sel = ft_selectdata(cfg, et);
sel = preprocess_et_trials_remove_blinks(sel);
freq_raw = compute_gaze_heatmap(sel, x_grid, y_grid, fsample, smoothing);
end

function sel = preprocess_et_trials_remove_blinks(sel)
% Same blink pipeline as AOC_gaze_fex_sternberg.m: in-screen mask, Y inversion (600 - y),
% then remove_blinks (100 ms window = 50 samples at 500 Hz). Column count is preserved.
if isempty(sel.trial)
    return
end
if exist('remove_blinks', 'file') ~= 2
    error('remove_blinks.m not found on MATLAB path (add /functions).');
end
win_size = 50;
for i = 1:numel(sel.trial)
    data = sel.trial{i};
    if isempty(data)
        continue
    end
    data = data(1:2, :);
    nS = size(data, 2);
    if isfield(sel, 'time') && numel(sel.time) >= i && ~isempty(sel.time{i})
        t = sel.time{i}(:)';
    else
        t = 1:nS;
    end
    if numel(t) ~= nS
        L = min(nS, numel(t));
        data = data(:, 1:L);
        t = t(1:L);
        nS = L;
    end
    valid = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
    data = data(:, valid);
    t = t(valid);
    if isempty(data)
        sel.trial{i} = zeros(2, 0);
        sel.time{i} = zeros(1, 0);
        continue
    end
    data(2, :) = 600 - data(2, :);
    data3 = [data; zeros(1, size(data, 2))];
    data3 = remove_blinks(data3, win_size);
    sel.trial{i} = data3(1:2, :);
    sel.time{i} = t;
end
end

function freq_raw = compute_gaze_heatmap(data, x_grid, y_grid, fsample, smoothing)
% Always use canonical computeGazeHeatmap from /functions.
if exist('computeGazeHeatmap', 'file') ~= 2
    error('computeGazeHeatmap.m not found on MATLAB path.');
end

n_bins = numel(x_grid) - 1;

if isempty(data.trial)
    freq_raw = struct();
    freq_raw.powspctrm = zeros(1, n_bins, n_bins);
    freq_raw.time = x_grid(2:end);
    freq_raw.freq = y_grid(2:end);
    freq_raw.label = {'et'};
    freq_raw.dimord = 'chan_freq_time';
    return
end

pos = horzcat(data.trial{:});
data3 = [pos; zeros(1, size(pos, 2))]; % canonical function expects >= 3 rows
freq_raw = computeGazeHeatmap(data3, n_bins, smoothing);

% Match the desired orientation for downstream plotting.
freq_raw.powspctrm(1, :, :) = flipud(squeeze(freq_raw.powspctrm(1, :, :)));
end

function trialinfo_vec = align_trialinfo_to_trials(trialinfo, n_trials, subj_label)
% FieldTrip-style trialinfo is commonly either:
% - nTrials x 1 vector (one value per trial)
% - nTrials x 2 matrix (two per-trial values/columns)
% - (2*nTrials) x 1 vector (duplicated rows)
%
% IMPORTANT: trialinfo(:) linearizes column-wise; this is *not* equivalent to
% "duplicated rows" when trialinfo is nTrials x 2.

if isempty(trialinfo)
    trialinfo_vec = [];
    return
end

ti = trialinfo;
ti_sz = size(ti);
is_2d = ndims(ti) == 2;

% Case A: nTrials x 2 (or more) matrix -> keep first column (condition),
% silently ignore the remaining columns (e.g., trial index).
if is_2d && ti_sz(1) == n_trials && ti_sz(2) >= 2
    trialinfo_vec = ti(:, 1);
    return
end

trialinfo_vec = ti(:);
n_ti = numel(trialinfo_vec);

if n_trials == 0 || n_ti == 0
    trialinfo_vec = [];
    return
end

if n_ti == n_trials
    return
end

% Common in these files: every trialinfo entry appears twice.
if n_ti == 2 * n_trials
    ti2 = reshape(trialinfo_vec, 2, n_trials)';
    if all(ti2(:, 1) == ti2(:, 2))
        trialinfo_vec = ti2(:, 1);
        fprintf('Warning: subject %s had duplicated trialinfo rows (2x). Collapsing to %d trials.\n', subj_label, n_trials);
    else
        trialinfo_vec = ti2(:, 1);
        fprintf('Warning: subject %s trialinfo is 2x trial count but pair values differ. Using first entry per pair.\n', subj_label);
    end
    return
end

n_common = min(n_trials, n_ti);
trialinfo_vec = trialinfo_vec(1:n_common);
fprintf('Warning: trial/trialinfo length mismatch for subject %s (trial=%d, trialinfo=%d). Truncating to %d.\n', ...
    subj_label, n_trials, n_ti, n_common);
end

function cells_out = downsample_gaze_cells_powspctrm(cells_in, target_bins)
% target_bins: scalar N → output N×N; or [screen_width, screen_height] in pixels.
% Gaze powspctrm is 1×freq×time with freq = vertical (y), time = horizontal (x); see compute_gaze_heatmap.
cells_out = cells_in;
if isscalar(target_bins)
    sz_ft = [target_bins, target_bins];
elseif numel(target_bins) == 2
    sz_ft = [target_bins(2), target_bins(1)]; % [nFreq, nTime] = [height, width]
else
    error('target_bins must be a scalar or a two-element [width height] vector.');
end
for ii = 1:numel(cells_out)
    if isempty(cells_out{ii}) || ~isfield(cells_out{ii}, 'powspctrm')
        continue
    end

    P = cells_out{ii}.powspctrm;
    if isempty(P)
        continue
    end

    % Expected shape for gaze maps in this script: 1 x freq x time.
    if ndims(P) ~= 3 || size(P, 1) ~= 1
        continue
    end

    P2 = squeeze(P(1, :, :));
    if ~isequal(size(P2), sz_ft)
        P2 = imresize(P2, sz_ft, 'bilinear');
    end
    cells_out{ii}.powspctrm = reshape(P2, [1, sz_ft(1), sz_ft(2)]);

    if isfield(cells_out{ii}, 'freq') && ~isempty(cells_out{ii}.freq)
        cells_out{ii}.freq = linspace(min(cells_out{ii}.freq), max(cells_out{ii}.freq), sz_ft(1));
    end
    if isfield(cells_out{ii}, 'time') && ~isempty(cells_out{ii}.time)
        cells_out{ii}.time = linspace(min(cells_out{ii}.time), max(cells_out{ii}.time), sz_ft(2));
    end
end
end

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm)
    return
end

P = freq_in.powspctrm;
if size(P, 1) < 3
    % Too few subjects for stable percentile clipping.
    return
end

P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

%{
%% ================= TOST FUNCTION =================
function [p1, p2, equivalent] = tost_welch(x1, x2, delta, alpha)

if nargin < 4
    alpha = 0.05;
end

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

n1 = length(x1);
n2 = length(x2);

m1 = mean(x1);
m2 = mean(x2);

v1 = var(x1);
v2 = var(x2);

SE = sqrt(v1/n1 + v2/n2);

df = (v1/n1 + v2/n2)^2 / ...
    ((v1^2)/(n1^2*(n1-1)) + (v2^2)/(n2^2*(n2-1)));

t1 = (m1 - m2 + delta) / SE;
t2 = (m1 - m2 - delta) / SE;

p1 = 1 - tcdf(t1, df);
p2 = tcdf(t2, df);

equivalent = (p1 < alpha) && (p2 < alpha);
end
%}

function [lme, used_formula] = fitlme_with_random_slope_fallback(tbl_in, fixed_formula)
formula_rs = sprintf('%s + (Load|Subject)', fixed_formula);
formula_ri = sprintf('%s + (1|Subject)', fixed_formula);
try
    lme = fitlme(tbl_in, formula_rs);
    used_formula = formula_rs;
catch
    lme = fitlme(tbl_in, formula_ri);
    used_formula = formula_ri;
end
end

function [est, ci, coeff_name] = get_group_effect_ci(lme, alpha_ci)
if nargin < 2
    alpha_ci = 0.10;
end
coef_names = lme.CoefficientNames;
% Main effect only: interaction terms are named Group_*:Load_* or Load_*:Group_*;
% fixedEffects names may reorder Load:Group, so the first "Group_" match can be an
% interaction and fail strcmp — or pick the wrong term for equivalence.
idx = startsWith(coef_names, 'Group_') & cellfun(@(s) isempty(strfind(s, ':')), coef_names);
if ~any(idx)
    error('No Group main-effect coefficient found in model.');
end
coeff_name = coef_names{find(idx, 1, 'first')};

[fe, fe_betanames] = fixedEffects(lme);
% R2013b+ returns betanames as a table with variable Name; older/docs cell paths differ.
if istable(fe_betanames) && ismember('Name', fe_betanames.Properties.VariableNames)
    fe_name_list = fe_betanames.Name;
elseif iscell(fe_betanames)
    fe_name_list = fe_betanames;
else
    fe_name_list = fe_betanames;
end
k_fe = find(ismember(fe_name_list(:), coeff_name), 1);
if isempty(k_fe)
    error('Group coefficient not found in fixed effects output.');
end
est = fe(k_fe);
ci_tbl = coefCI(lme, 'Alpha', alpha_ci);
ci = ci_tbl(strcmp(lme.CoefficientNames, coeff_name), :);
end

function tf = ci_within_bounds(ci, bounds)
tf = ci(1) >= bounds(1) && ci(2) <= bounds(2);
end