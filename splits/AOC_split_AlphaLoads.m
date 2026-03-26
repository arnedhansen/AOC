%% AOC Split Sternberg Alpha Loads
% Stratifies participants by alpha power slope across WM load (2,4,6).
% Alpha increase vs decrease subgroups; TFRs, power spectra, behavioral (RT, ACC),
% optional gaze density. Uses 1-2 s retention window for alpha extraction.
%%
% Groups:
%       Jensen: amplification of alpha over WM loads
%       N-back: reduction of alpha over WM loads
%       Flat: flat slope

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
log_dir = fullfile(base_data, 'results', 'logs');
if ~isfolder(log_dir)
    mkdir(log_dir);
end
cmdlog_file = fullfile(log_dir, sprintf('AOC_split_AlphaLoads_commandwindow_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
diary('off');
diary(cmdlog_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Command window log file: %s\n', cmdlog_file);

%% Figure setup
fig_pos = [0 0 1512 982];
fontSize = 15;
color_map = customcolormap_preset('red-white-blue');

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
disp('LOADING FINISHED')

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

%% Select alpha power (retention 1-2 s)
cfg = [];
cfg.frequency = [8 14];
cfg.latency = [1 2];
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

%% Use bootstrap regression to identify slope uncertainty per subject
nSubj = length(subjects);
n_boot = 2000;
rng(1, 'twister'); % reproducible bootstrap labels

slope = nan(nSubj, 1);
slope_ci_low = nan(nSubj, 1);
slope_ci_high = nan(nSubj, 1);
ci_valid = false(nSubj, 1);

for subj = 1:nSubj
    a2_trials = extract_alpha_trials(load2{subj}, channels, [8 14], [1 2]);
    a4_trials = extract_alpha_trials(load4{subj}, channels, [8 14], [1 2]);
    a6_trials = extract_alpha_trials(load6{subj}, channels, [8 14], [1 2]);

    if isempty(a2_trials) || isempty(a4_trials) || isempty(a6_trials)
        continue
    end

    y_obs = [mean(a2_trials, 'omitnan'), mean(a4_trials, 'omitnan'), mean(a6_trials, 'omitnan')];
    b_obs = [ones(3,1), [2;4;6]] \ y_obs';
    slope(subj) = b_obs(2);

    boot_slopes = bootstrap_alpha_slope(a2_trials, a4_trials, a6_trials, n_boot);
    boot_slopes = boot_slopes(isfinite(boot_slopes));
    if isempty(boot_slopes)
        continue
    end
    ci = prctile(boot_slopes, [2.5 97.5]);
    slope_ci_low(subj) = ci(1);
    slope_ci_high(subj) = ci(2);
    ci_valid(subj) = true;
end

%% Split and Plot slope distribution (inclusion)
idx_jensen = ci_valid & (slope_ci_low > 0);
idx_nback  = ci_valid & (slope_ci_high < 0);
idx_flat   = ~(idx_jensen | idx_nback);

% counts
n_j = sum(idx_jensen);
n_n = sum(idx_nback);
n_f = sum(idx_flat);

% Plot
close all
figure('Position', fig_pos, 'Color', 'w');
hold on
histogram(slope(idx_jensen), 20, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6);
histogram(slope(idx_nback),  20, 'FaceColor', [0 0 0.8], 'FaceAlpha', 0.6);
histogram(slope(idx_flat),   11, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);
xline(0, 'k--', 'LineWidth', 2);
xlabel('Alpha power slope')
ylabel('Participants')
title('Linear Slope of Alpha Power across WM Load (2, 4, 6 items)', 'FontSize', 20)
legend({sprintf('Alpha increase with load (N=%d)', n_j), ...
    sprintf('Alpha decrease with load (N=%d)', n_n), ...
    sprintf('intermediate (N=%d)', n_f), ...
    'zero-slope boundary'}, 'Box', 'off')
box on
set(gca, 'FontSize', 15)
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_histogram_inclusion.png'));

fprintf('\nBootstrap slope classification summary (95%% CI rule):\n');
fprintf('Increase (CI > 0):      %d\n', n_j);
fprintf('Decrease (CI < 0):      %d\n', n_n);
fprintf('Indeterminate (CI ~ 0): %d\n', n_f);
if n_j == 0 || n_n == 0
    error(['CI-based classification produced an empty increase/decrease subgroup. ' ...
        'With the current data representation, downstream subgroup comparisons are not defined.']);
end

%% Grand averages per subgroup
cfg = [];
cfg.keepindividual = 'yes';
ga2jensen = ft_freqgrandaverage(cfg,load2{idx_jensen});
ga4jensen = ft_freqgrandaverage(cfg,load4{idx_jensen});
ga6jensen = ft_freqgrandaverage(cfg,load6{idx_jensen});
ga2nback = ft_freqgrandaverage(cfg,load2{idx_nback});
ga4nback = ft_freqgrandaverage(cfg,load4{idx_nback});
ga6nback = ft_freqgrandaverage(cfg,load6{idx_nback});

%% TFR
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [5 30];
cfg.xlim = [-.5 2];
cfg.zlim = [-.25 .25];
cfg.channel = channels;
cfg.layout = headmodel.layANThead;
figure('Position', fig_pos, 'Color', 'w');
subplot(3,2,1); ft_singleplotTFR(cfg, ga2jensen); title('Increase: WM load 2', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,3); ft_singleplotTFR(cfg, ga4jensen); title('Increase: WM load 4', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,5); ft_singleplotTFR(cfg, ga6jensen); title('Increase: WM load 6', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;

subplot(3,2,2); ft_singleplotTFR(cfg, ga2nback); title('Decrease: WM load 2', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,4); ft_singleplotTFR(cfg, ga4nback); title('Decrease: WM load 4', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
subplot(3,2,6); ft_singleplotTFR(cfg, ga6nback); title('Decrease: WM load 6', 'FontSize', fontSize);
ax = gca; c = colorbar; xlabel(ax, 'Time [s]', 'FontSize', fontSize); ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
set(ax, 'FontSize', fontSize); c.FontSize = fontSize - 2; c.Label.String = 'Power [dB]'; c.Label.FontSize = fontSize;
colormap(gcf, color_map);
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_TFR_OVERVIEW.png'));

% Save each EEG TFR panel as an individual figure
tfr_data = {ga2jensen, ga4jensen, ga6jensen, ga2nback, ga4nback, ga6nback};
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
    title(tfr_titles{ii}, 'FontSize', fontSize);
    ax = gca;
    c = colorbar;
    xlabel(ax, 'Time [s]', 'FontSize', fontSize);
    ylabel(ax, 'Frequency [Hz]', 'FontSize', fontSize-4);
    set(ax, 'FontSize', fontSize);
    c.FontSize = fontSize - 2;
    c.Label.String = 'Power [dB]';
    c.Label.FontSize = fontSize;
    colormap(gcf, color_map);
    drawnow;
    saveas(gcf, fullfile(fig_dir, tfr_filenames{ii}));
    close(gcf);
end

%% Select powspctrm (retention 1-2 s)
cfg = [];
cfg.latency = [1 2];
cfg.avgovertime = 'yes';
ga2jensen_powspctrm = ft_selectdata(cfg,ga2jensen);
ga2jensen_powspctrm = rmfield(ga2jensen_powspctrm,'time');

ga4jensen_powspctrm = ft_selectdata(cfg,ga4jensen);
ga4jensen_powspctrm = rmfield(ga4jensen_powspctrm,'time');

ga6jensen_powspctrm = ft_selectdata(cfg,ga6jensen);
ga6jensen_powspctrm = rmfield(ga6jensen_powspctrm,'time');

ga2nback_powspctrm = ft_selectdata(cfg,ga2nback);
ga2nback_powspctrm = rmfield(ga2nback_powspctrm,'time');

ga4nback_powspctrm = ft_selectdata(cfg,ga4nback);
ga4nback_powspctrm = rmfield(ga4nback_powspctrm,'time');

ga6nback_powspctrm = ft_selectdata(cfg,ga6nback);
ga6nback_powspctrm = rmfield(ga6nback_powspctrm,'time');

%% Plot Powerspectra
close all
figure('Position', fig_pos, 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

[~, ch_idx_j] = ismember(channels, ga2jensen_powspctrm.label);
ch_idx_j = ch_idx_j(ch_idx_j > 0);
freqs_j = ga2jensen_powspctrm.freq;

addpath('W:\Students\Arne\toolboxes\shadedErrorBar\')
nexttile; hold on
ga_j = {ga2jensen_powspctrm, ga4jensen_powspctrm, ga6jensen_powspctrm};
for c = 1:3
    subj_spec = squeeze(mean(ga_j{c}.powspctrm(:, ch_idx_j, :), 2, 'omitnan'));
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
    eb_j(c) = shadedErrorBar(freqs_j, m, se, 'lineProps', {'-'}, 'transparent', true); %#ok<AGROW>
    set(eb_j(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
    set(eb_j(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
    set(eb_j(c).edge(1), 'Color', colors(c, :));
    set(eb_j(c).edge(2), 'Color', colors(c, :));
end
xlim([2 30]);
xlabel('Frequency [Hz]', 'FontSize', fontSize-4);
yline(0, '--')
ylim([-0.075 0.225])
ylabel('Power [dB]', 'FontSize', fontSize);
title('Increase group', 'FontSize', fontSize);
% Legend with colored patch boxes
leg_p_j = gobjects(1, 3);
for c = 1:3
    leg_p_j(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_j, {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', fontSize - 2);
set(gca, 'FontSize', fontSize);
box on

[~, ch_idx_n] = ismember(channels, ga2nback_powspctrm.label);
ch_idx_n = ch_idx_n(ch_idx_n > 0);
freqs_n = ga2nback_powspctrm.freq;

nexttile; hold on
ga_n = {ga2nback_powspctrm, ga4nback_powspctrm, ga6nback_powspctrm};
for c = 1:3
    subj_spec = squeeze(mean(ga_n{c}.powspctrm(:, ch_idx_n, :), 2, 'omitnan'));
    m = mean(subj_spec, 1, 'omitnan');
    se = std(subj_spec, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(subj_spec), 1)), 1);
    eb_n(c) = shadedErrorBar(freqs_n, m, se, 'lineProps', {'-'}, 'transparent', true); %#ok<AGROW>
    set(eb_n(c).mainLine, 'Color', colors(c, :), 'LineWidth', 3);
    set(eb_n(c).patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
    set(eb_n(c).edge(1), 'Color', colors(c, :));
    set(eb_n(c).edge(2), 'Color', colors(c, :));
end
xlim([2 30]);
xlabel('Frequency [Hz]', 'FontSize', fontSize-4);
yline(0, '--')
ylim([-0.075 0.225])
ylabel('Power [dB]', 'FontSize', fontSize);
title('Decrease group', 'FontSize', fontSize);
% Legend with colored patch boxes
leg_p_n = gobjects(1, 3);
for c = 1:3
    leg_p_n(c) = patch(NaN, NaN, colors(c, :), 'EdgeColor', 'none');
end
legend(leg_p_n, {'WM load 2', 'WM load 4', 'WM load 6'}, ...
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

% Keep full-resolution maps for visualization.
% Use downsampled copies only for permutation stats to avoid OOM.
gaze_cbpt_bins = 500;
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
cfg = [];
cfg.keepindividual = 'yes';
ga2jensen_gaze = ft_freqgrandaverage(cfg, load2_gaze{idx_jensen});
ga4jensen_gaze = ft_freqgrandaverage(cfg, load4_gaze{idx_jensen});
ga6jensen_gaze = ft_freqgrandaverage(cfg, load6_gaze{idx_jensen});
ga2nback_gaze = ft_freqgrandaverage(cfg, load2_gaze{idx_nback});
ga4nback_gaze = ft_freqgrandaverage(cfg, load4_gaze{idx_nback});
ga6nback_gaze = ft_freqgrandaverage(cfg, load6_gaze{idx_nback});

%% Plot Gaze TFRs
close all
cfg = [];
cfg.figure = 'gcf';
% Strict symmetric max-abs z-limits for subtraction baseline maps
allGazeDiff = [ ...
    ga2jensen_gaze.powspctrm(:); ...
    ga4jensen_gaze.powspctrm(:); ...
    ga6jensen_gaze.powspctrm(:); ...
    ga2nback_gaze.powspctrm(:); ...
    ga4nback_gaze.powspctrm(:); ...
    ga6nback_gaze.powspctrm(:) ];
allGazeDiff = allGazeDiff(isfinite(allGazeDiff));
zlimAbs = max(abs(allGazeDiff));
if ~isfinite(zlimAbs) || zlimAbs == 0
    zlimAbs = 3;
end
cfg.zlim = [-zlimAbs zlimAbs];
figure('Position', fig_pos, 'Color', 'w');

% Gaze density relative to baseline (subtraction; tasklate - baseline)
subplot(3,2,1);
ft_singleplotTFR(cfg, ga2jensen_gaze);
title('Amplification: WM Load 2', 'FontSize', fontSize);
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
ft_singleplotTFR(cfg, ga4jensen_gaze);
title('Amplification: WM Load 4', 'FontSize', fontSize);
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
ft_singleplotTFR(cfg, ga6jensen_gaze);
title('Amplification: WM Load 6', 'FontSize', fontSize);
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
ft_singleplotTFR(cfg, ga2nback_gaze);
title('Reduction: WM Load 2', 'FontSize', fontSize);
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
ft_singleplotTFR(cfg, ga4nback_gaze);
title('Reduction: WM Load 4', 'FontSize', fontSize);
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
ft_singleplotTFR(cfg, ga6nback_gaze);
title('Reduction: WM Load 6', 'FontSize', fontSize);
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
% Six CBPT maps are tested here (3 loads x 2 groups).
% Use a family-wise corrected alpha across this test family.
family_alpha = 0.05;
n_cbpt_tests_taskVsBase = 6;
alpha_cbpt_taskVsBase = family_alpha / n_cbpt_tests_taskVsBase;

    cfg                  = [];
    cfg.spmversion       = 'spm12';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';  % paired: baseline vs task
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
    subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
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
    cfg_cbpt_gaze_taskVsBase_jensen = cfg;
    design_cbpt_gaze_taskVsBase_jensen = design;
    clc; disp(upper('[stat_inc_2 JENSEN]')); [stat_inc_2] = ft_freqstatistics(cfg, allgazetasklate2_cbpt{idx_jensen},allgazebase2_cbpt{idx_jensen});
    clc; disp(upper('[stat_inc_4 JENSEN]')); [stat_inc_4] = ft_freqstatistics(cfg, allgazetasklate4_cbpt{idx_jensen},allgazebase4_cbpt{idx_jensen});
    clc; disp(upper('[stat_inc_6 JENSEN]')); [stat_inc_6] = ft_freqstatistics(cfg, allgazetasklate6_cbpt{idx_jensen},allgazebase6_cbpt{idx_jensen});
    % Hide non-significant pixels in visualization
    stat_inc_2.stat(stat_inc_2.mask==0) = NaN;
    stat_inc_4.stat(stat_inc_4.mask==0) = NaN;
    stat_inc_6.stat(stat_inc_6.mask==0) = NaN;

    subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
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
    cfg_cbpt_gaze_taskVsBase_nback = cfg;
    design_cbpt_gaze_taskVsBase_nback = design;
    clc; disp(upper('[stat_inc_n_2 N-back]')); [stat_inc_n_2] = ft_freqstatistics(cfg, allgazetasklate2_cbpt{idx_nback},allgazebase2_cbpt{idx_nback});
    clc; disp(upper('[stat_inc_n_4 N-back]')); [stat_inc_n_4] = ft_freqstatistics(cfg, allgazetasklate4_cbpt{idx_nback},allgazebase4_cbpt{idx_nback});
    clc; disp(upper('[stat_inc_n_6 N-back]')); [stat_inc_n_6] = ft_freqstatistics(cfg, allgazetasklate6_cbpt{idx_nback},allgazebase6_cbpt{idx_nback});
    % Hide non-significant pixels in visualization
    stat_inc_n_2.stat(stat_inc_n_2.mask==0) = NaN;
    stat_inc_n_4.stat(stat_inc_n_4.mask==0) = NaN;
    stat_inc_n_6.stat(stat_inc_n_6.mask==0) = NaN;


% cohensd=((stat_inc_2.stat)./sqrt(subj));
% stat_inc_2.stat=cohensd;
% stat_inc_2.stat(stat_inc_2.mask==0)=0;% set everything not relevant to zero
%
% cohensd=((stat_inc_4.stat)./sqrt(subj));
% stat_inc_4.stat=cohensd;
% stat_inc_4.stat(stat_inc_4.mask==0)=0;% set everything not relevant to zero
%
% cohensd=((stat_inc_6.stat)./sqrt(subj));
% stat_inc_6.stat=cohensd;
% stat_inc_6.stat(stat_inc_6.mask==0)=0;% set everything not relevant to zero

%% Plot CBPT Raw Stats
stat_inc_n_2.cfg=[];
stat_inc_n_4.cfg=[];
stat_inc_n_6.cfg=[];
close all
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'maxabs';
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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Amplification: WM Load 2', 'FontSize', fontSize, 'Interpreter', 'none')

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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Amplification: WM Load 4', 'FontSize', fontSize, 'Interpreter', 'none')

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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Amplification: WM Load 6', 'FontSize', fontSize, 'Interpreter', 'none')

% plot decreaseing
subplot(3,2,2);
ft_singleplotTFR(cfg,stat_inc_n_2);
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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Reduction: WM Load 2', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,4);
ft_singleplotTFR(cfg,stat_inc_n_4);
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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Reduction: WM Load 4', 'FontSize', fontSize, 'Interpreter', 'none')

subplot(3,2,6);
ft_singleplotTFR(cfg,stat_inc_n_6);
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
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Reduction: WM Load 6', 'FontSize', fontSize, 'Interpreter', 'none')
colormap(customcolormap_preset('red-white-blue'));

drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_statInc.png'));

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

subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3];
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_omnibus_jensen = cfg;
[statF_gaze] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_jensen}, load4_gaze_cbpt{idx_jensen}, load6_gaze_cbpt{idx_jensen});
statF_gaze.stat(statF_gaze.mask==0)=0;% set everything not relevant to zero

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

    subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
    n_2  = subj;
    n_4  = subj;
    n_6 =  subj;

    cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3];
    cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6];
    cfg.ivar                  = 1;
    cfg.uvar                  = 2;
cfg_cbpt_gaze_omnibus_nback = cfg;
[statF_gaze_n] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_nback}, load4_gaze_cbpt{idx_nback}, load6_gaze_cbpt{idx_nback});
statF_gaze_n.stat(statF_gaze_n.mask==0)=0;% set everything not relevant to zero

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
ft_singleplotTFR(cfg, statF_gaze);
title('Amplification: Omnibus Load Effect', 'FontSize', fontSize, 'Interpreter', 'none');
xlim(ax1, [0 800]);
ylim(ax1, [0 600]);
hold(ax1, 'on');
plot(ax1, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');

ax2 = subplot(2,1,2);
ft_singleplotTFR(cfg, statF_gaze_n);
title('Reduction: Omnibus Load Effect', 'FontSize', fontSize, 'Interpreter', 'none');
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

subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadtrend_jensen = cfg;
[statT_gaze] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_jensen}, load4_gaze_cbpt{idx_jensen}, load6_gaze_cbpt{idx_jensen});
statT_gaze.stat(statT_gaze.mask==0)=0;

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
subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadtrend_nback = cfg;
[statT_gaze_n] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_nback}, load4_gaze_cbpt{idx_nback}, load6_gaze_cbpt{idx_nback});
statT_gaze_n.stat(statT_gaze_n.mask==0)=0;

cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.zlim = [-5 5];
cfg.xlim = [0 800];
cfg.ylim = [0 600];
cfg.figure = 'gcf';
figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1,2,1); ft_singleplotTFR(cfg,statT_gaze); title('Amplification: Linear Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;

subplot(1,2,2); ft_singleplotTFR(cfg,statT_gaze_n); title('Reduction: Linear Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadTrend.png'));
close(gcf);

% ---------------- Quadratic load trend ----------------
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

subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadquadratic_jensen = cfg;
[statT_gaze_quad] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_jensen}, load4_gaze_cbpt{idx_jensen}, load6_gaze_cbpt{idx_jensen});
statT_gaze_quad.stat(statT_gaze_quad.mask==0)=0;

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

subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
cfg.design(1,:)           = [ones(1,subj), ones(1,subj)*2, ones(1,subj)*3];
cfg.design(2,:)           = [1:subj, 1:subj, 1:subj];
cfg.ivar                  = 1;
cfg.uvar                  = 2;
cfg_cbpt_gaze_loadquadratic_nback = cfg;
[statT_gaze_quad_n] = ft_freqstatistics(cfg, load2_gaze_cbpt{idx_nback}, load4_gaze_cbpt{idx_nback}, load6_gaze_cbpt{idx_nback});
statT_gaze_quad_n.stat(statT_gaze_quad_n.mask==0)=0;

cfg = []; cfg.parameter = 'stat'; cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline'; cfg.zlim = [-5 5]; cfg.xlim = [0 800]; cfg.ylim = [0 600]; cfg.figure = 'gcf';
figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1,2,1); ft_singleplotTFR(cfg,statT_gaze_quad); title('Amplification: Quadratic Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;

subplot(1,2,2); ft_singleplotTFR(cfg,statT_gaze_quad_n); title('Reduction: Quadratic Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');
ax = gca; set(ax, 'FontSize', fontSize); xlim(ax, [0 800]); ylim(ax, [0 600]); hold(ax, 'on'); plot(ax, 400, 300, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k'); xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize); ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
c = colorbar(ax); c.LineWidth = 1; c.FontSize = fontSize - 2; c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)]; c.Label.String = 't-value'; c.Label.FontSize = fontSize - 2;
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadQuadratic.png'));

clear statT_gaze_quad statT_gaze_quad_n

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
sb2 = RT2(idx_jensen);
sb4 = RT4(idx_jensen);
sb6 = RT6(idx_jensen);

nb1 = RT2(idx_nback);
nb2 = RT4(idx_nback);
nb3 = RT6(idx_nback);

%%
figure('Position', [0 0 1512 982], 'Color', 'w');
clf;

positions = [.3, .6, .9];
box_w = 0.05;
density_scale = 0.045;
density_offset = 0.06;
box_offset = 0.03;
cond_cols = [0 0 0; 0 0 1; 1 0 0];

% LEFT: alpha increase group
subplot(2,2,1);
plot_split_raincloud_triplet({sb2, sb4, sb6}, positions, cond_cols, [0 1.5], ...
    'Reaction Time [s]', sprintf('Alpha increase with load (N=%d)', sum(idx_jensen)), ...
    18, density_scale, density_offset, box_offset, box_w);

% RIGHT: alpha decrease group
subplot(2,2,2);
plot_split_raincloud_triplet({nb1, nb2, nb3}, positions, cond_cols, [0 1.5], ...
    'Reaction Time [s]', sprintf('Alpha decrease with load (N=%d)', sum(idx_nback)), ...
    18, density_scale, density_offset, box_offset, box_w);
set(gcf,'color','w');

%% Accuracy example
sb2 = ACC2(idx_jensen);
sb4 = ACC4(idx_jensen);
sb6 = ACC6(idx_jensen);

nb1 = ACC2(idx_nback);
nb2 = ACC4(idx_nback);
nb3 = ACC6(idx_nback);

% remove NaNs (important!)
sb2 = sb2(~isnan(sb2)); sb4 = sb4(~isnan(sb4)); sb6 = sb6(~isnan(sb6));
nb1 = nb1(~isnan(nb1)); nb2 = nb2(~isnan(nb2)); nb3 = nb3(~isnan(nb3));

%%
% ================= LEFT: alpha increase =================
subplot(2,2,3);
plot_split_raincloud_triplet({sb2, sb4, sb6}, positions, cond_cols, [50 110], ...
    'Accuracy [%]', sprintf('Alpha increase with load (N=%d)', sum(idx_jensen)), ...
    18, density_scale, density_offset, box_offset, box_w);

% ================= RIGHT: alpha decrease =================
subplot(2,2,4);
plot_split_raincloud_triplet({nb1, nb2, nb3}, positions, cond_cols, [50 110], ...
    'Accuracy [%]', sprintf('Alpha decrease with load (N=%d)', sum(idx_nback)), ...
    18, density_scale, density_offset, box_offset, box_w);

set(gcf,'color','w');
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_RT_ACC.png'));

%% LME / TOST (test effects)
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

%% ================= GROUPING =================
% Group labels are already defined above using bootstrap slope CIs.

%% ================= BUILD TABLE =================
Subject = [];
Load = [];
Group = [];
RT = [];
ACC = [];

for i = 1:nSubj

    if idx_jensen(i)
        g = 1;
    elseif idx_nback(i)
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
tbl.Group   = categorical(tbl.Group);

%% ================= LME ANALYSIS =================
disp('--- RT model ---')
% Treat Load as categorical for a full 3-level repeated-measures effect.
tbl.Load = categorical(tbl.Load);
lme_RT = fitlme(tbl, 'RT ~ Load * Group + (1|Subject)');
disp(anova(lme_RT))

disp('--- ACC model ---')
lme_ACC = fitlme(tbl, 'ACC ~ Load * Group + (1|Subject)');
disp(anova(lme_ACC))


%% ================= TOST ANALYSIS =================
delta_RT = 0.05;   % 50 ms
delta_ACC = 5;     % 5%

% Use one independent value per subject (mean across loads) for between-group TOST.
RT_jensen = mean([RT2(idx_jensen), RT4(idx_jensen), RT6(idx_jensen)], 2, 'omitnan');
RT_nback  = mean([RT2(idx_nback),  RT4(idx_nback),  RT6(idx_nback)], 2, 'omitnan');

ACC_jensen = mean([ACC2(idx_jensen), ACC4(idx_jensen), ACC6(idx_jensen)], 2, 'omitnan');
ACC_nback  = mean([ACC2(idx_nback),  ACC4(idx_nback),  ACC6(idx_nback)], 2, 'omitnan');

fprintf('\n--- TOST (strict) ---\n')

[p1, p2, eq_RT] = tost_welch(RT_jensen, RT_nback, delta_RT);
fprintf('RT equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_RT);

[p1, p2, eq_ACC] = tost_welch(ACC_jensen, ACC_nback, delta_ACC);
fprintf('ACC equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_ACC);

%% ================= RT SENSITIVITY =================
fprintf('\n--- RT equivalence sensitivity ---\n')

[p1, p2, eq] = tost_welch(RT_jensen, RT_nback, 0.05);
fprintf('RT delta=0.05: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(RT_jensen, RT_nback, 0.1);
fprintf('RT delta=0.10: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(RT_jensen, RT_nback, 0.15);
fprintf('RT delta=0.15: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);

%% ================= ACC ROBUSTNESS =================
fprintf('\n--- ACC robustness ---\n')

[p1, p2, eq] = tost_welch(ACC_jensen, ACC_nback, 3);
fprintf('ACC delta=3: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(ACC_jensen, ACC_nback, 5);
fprintf('ACC delta=5: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);
[p1, p2, eq] = tost_welch(ACC_jensen, ACC_nback, 10);
fprintf('ACC delta=10: p1=%.4g, p2=%.4g, equivalent=%d\n', p1, p2, eq);

%% ================= MIXED-EFFECTS EQUIVALENCE =================
fprintf('\n--- Mixed-effects equivalence (CI-in-bounds) ---\n')

tbl_me = tbl;
idx_finite = isfinite(tbl_me.RT) & isfinite(tbl_me.ACC);
tbl_me = tbl_me(idx_finite, :);

% RT model is fit on log scale to improve residual behavior and interpretability.
idx_rt_valid = tbl_me.RT > 0;
tbl_rt = tbl_me(idx_rt_valid, :);
tbl_rt.logRT = log(tbl_rt.RT);

[lme_rt_eq, rt_formula] = fitlme_with_random_slope_fallback(tbl_rt, 'logRT ~ Load * Group');
[rt_est, rt_ci90, rt_group_term] = get_group_effect_ci(lme_rt_eq, 0.10);

rt_bounds_ratio = [0.95, 1.05];
rt_bounds_log = log(rt_bounds_ratio);
rt_eq = ci_within_bounds(rt_ci90, rt_bounds_log);

fprintf('RT model formula: %s\n', rt_formula);
fprintf('RT group effect term: %s\n', rt_group_term);
fprintf('RT log-scale estimate=%.4g, 90%% CI=[%.4g, %.4g], bounds=[%.4g, %.4g], equivalent=%d\n', ...
    rt_est, rt_ci90(1), rt_ci90(2), rt_bounds_log(1), rt_bounds_log(2), rt_eq);
fprintf('RT ratio estimate=%.4g, 90%% CI=[%.4g, %.4g], ratio bounds=[%.2f, %.2f]\n', ...
    exp(rt_est), exp(rt_ci90(1)), exp(rt_ci90(2)), rt_bounds_ratio(1), rt_bounds_ratio(2));

[lme_acc_eq, acc_formula] = fitlme_with_random_slope_fallback(tbl_me, 'ACC ~ Load * Group');
[acc_est, acc_ci90, acc_group_term] = get_group_effect_ci(lme_acc_eq, 0.10);

acc_bounds = [-5, 5];
acc_eq = ci_within_bounds(acc_ci90, acc_bounds);

fprintf('ACC model formula: %s\n', acc_formula);
fprintf('ACC group effect term: %s\n', acc_group_term);
fprintf('ACC estimate=%.4g, 90%% CI=[%.4g, %.4g], bounds=[%.4g, %.4g], equivalent=%d\n', ...
    acc_est, acc_ci90(1), acc_ci90(2), acc_bounds(1), acc_bounds(2), acc_eq);

fprintf('\n--- Mixed-effects RT sensitivity (ratio bounds) ---\n')
rt_sens_ratio = [0.95 1.05; 0.90 1.10; 0.85 1.15];
for ii = 1:size(rt_sens_ratio, 1)
    b = log(rt_sens_ratio(ii, :));
    eq = ci_within_bounds(rt_ci90, b);
    fprintf('RT ratio bounds=[%.2f, %.2f]: equivalent=%d\n', rt_sens_ratio(ii,1), rt_sens_ratio(ii,2), eq);
end

fprintf('\n--- Mixed-effects ACC sensitivity (point bounds) ---\n')
acc_sens = [3 5 10];
for ii = 1:numel(acc_sens)
    b = [-acc_sens(ii), acc_sens(ii)];
    eq = ci_within_bounds(acc_ci90, b);
    fprintf('ACC bounds=[%d, %d]: equivalent=%d\n', b(1), b(2), eq);
end
%% Helper functions
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

function alpha_trials = extract_alpha_trials(tfr_in, channels, freq_win, time_win)
% Returns resampling units for alpha power within selected channel/frequency/time.
% Priority:
%   1) trial-level units when `rpt` dimension exists
%   2) otherwise flattened channel/frequency/time samples as fallback units
cfg = [];
cfg.channel = channels;
cfg.frequency = freq_win;
cfg.latency = time_win;
sel = ft_selectdata(cfg, tfr_in);
vals = sel.powspctrm;
if isempty(vals)
    alpha_trials = [];
    return
end

dimord = '';
if isfield(sel, 'dimord') && ischar(sel.dimord)
    dimord = sel.dimord;
end

if contains(dimord, 'rpt')
    n_rpt = size(vals, 1);
    vals2d = reshape(vals, n_rpt, []);
    alpha_trials = mean(vals2d, 2, 'omitnan');
else
    alpha_trials = vals(:);
end
alpha_trials = alpha_trials(isfinite(alpha_trials));
end

function boot_slopes = bootstrap_alpha_slope(a2_trials, a4_trials, a6_trials, n_boot)
% Bootstraps slope estimates by resampling trials per load with replacement.
boot_slopes = nan(n_boot, 1);
Xmat = [ones(3,1), [2;4;6]];
n2 = numel(a2_trials);
n4 = numel(a4_trials);
n6 = numel(a6_trials);
if min([n2 n4 n6]) < 2
    return
end
for bb = 1:n_boot
    m2 = mean(a2_trials(randi(n2, n2, 1)), 'omitnan');
    m4 = mean(a4_trials(randi(n4, n4, 1)), 'omitnan');
    m6 = mean(a6_trials(randi(n6, n6, 1)), 'omitnan');
    if any(~isfinite([m2 m4 m6]))
        continue
    end
    b = Xmat \ [m2; m4; m6];
    boot_slopes(bb) = b(2);
end
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

function plot_split_raincloud_triplet(data_cells, positions, cond_cols, y_limits, y_lab, ttl, fsz, density_scale, density_offset, box_offset, box_w)
hold on
for cc = 1:3
    y = data_cells{cc};
    y = y(isfinite(y));
    if numel(y) < 3
        continue
    end

    [f, xi] = ksdensity(y, 'NumPoints', 120);
    if max(f) > 0
        f = (f ./ max(f)) * density_scale;
    else
        f = zeros(size(f));
    end

    x_den = positions(cc) - density_offset;
    x_box = positions(cc) + box_offset;

    fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], cond_cols(cc, :), ...
        'FaceAlpha', 0.30, 'EdgeColor', cond_cols(cc, :), 'LineWidth', 1);

    q1 = prctile(y, 25);
    q3 = prctile(y, 75);
    med = median(y);
    p5 = prctile(y, 5);
    p95 = prctile(y, 95);

    plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
    plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
    rectangle('Position', [x_box-box_w/2, q1, box_w, q3-q1], ...
        'FaceColor', [cond_cols(cc, :) 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);

    jit = box_w * (rand(numel(y),1)-0.5);
    scatter(x_box + jit, y, 24, cond_cols(cc,:), 'filled', 'MarkerFaceAlpha', 0.50, ...
        'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

ylabel(y_lab);
title(ttl);
set(gca, 'FontSize', fsz);
set(gca, 'XTick', positions, 'XTickLabel', {'load 2','load 4','load 6'});
ylim(y_limits);
xlim([positions(1)-0.12 positions(end)+0.12]);
box on
end
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

    pos_cells{end+1} = [xw; yw]; %#ok<AGROW>
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

    allgazetasklate2{subj} = extract_gaze_window(et, idx2, [1 2], x_grid, y_grid, fsample, smoothing);
    allgazetasklate4{subj} = extract_gaze_window(et, idx4, [1 2], x_grid, y_grid, fsample, smoothing);
    allgazetasklate6{subj} = extract_gaze_window(et, idx6, [1 2], x_grid, y_grid, fsample, smoothing);
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
freq_raw = compute_gaze_heatmap(sel, x_grid, y_grid, fsample, smoothing);
end

function freq_raw = compute_gaze_heatmap(data, x_grid, y_grid, fsample, smoothing)
% RAW: use the canonical `computeGazeHeatmap` from the functions repo.
if isempty(data.trial)
    pos = zeros(2, 0);
else
    pos = horzcat(data.trial{:});
end
data3 = [pos; zeros(1, size(pos, 2))]; % computeGazeHeatmap expects at least 3 rows
freq_raw = computeGazeHeatmap(data3, numel(x_grid) - 1, smoothing);
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
cells_out = cells_in;
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
    if ~isequal(size(P2), [target_bins, target_bins])
        P2 = imresize(P2, [target_bins, target_bins], 'bilinear');
    end
    cells_out{ii}.powspctrm = reshape(P2, [1, target_bins, target_bins]);

    if isfield(cells_out{ii}, 'freq') && ~isempty(cells_out{ii}.freq)
        cells_out{ii}.freq = linspace(min(cells_out{ii}.freq), max(cells_out{ii}.freq), target_bins);
    end
    if isfield(cells_out{ii}, 'time') && ~isempty(cells_out{ii}.time)
        cells_out{ii}.time = linspace(min(cells_out{ii}.time), max(cells_out{ii}.time), target_bins);
    end
end
end

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
idx = startsWith(coef_names, 'Group_');
if ~any(idx)
    error('No Group main-effect coefficient found in model.');
end
coeff_name = coef_names{find(idx, 1, 'first')};

[fe, fe_names] = fixedEffects(lme);
fe_idx = strcmp(fe_names, coeff_name);
if ~any(fe_idx)
    error('Group coefficient not found in fixed effects output.');
end
est = fe(fe_idx);
ci_tbl = coefCI(lme, alpha_ci);
ci = ci_tbl(strcmp(lme.CoefficientNames, coeff_name), :);
end

function tf = ci_within_bounds(ci, bounds)
tf = ci(1) >= bounds(1) && ci(2) <= bounds(2);
end