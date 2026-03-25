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

%% Figure setup
fig_pos = [0 0 1512 982];
fontSize = 15;
color_map = interp1(linspace(0,1,5), [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], linspace(0,1,64));

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

%% Use regression to identify slopes
for subj = 1:length(alpha2)
    y = [alpha2(subj), alpha4(subj), alpha6(subj)];
    Xmat = [ones(3,1), [2;4;6]];
    b = Xmat \ y';
    slope(subj) = b(2);
end

%% Split and Plot slope distribution (inclusion)
thr = 0.01;
idx_jensen = slope > thr;
idx_nback  = slope < -thr;
idx_flat   = abs(slope) <= thr;

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
xline(thr,  'k--', 'LineWidth', 2);
xline(-thr, 'k--', 'LineWidth', 2);
xlabel('Alpha power slope')
ylabel('Participants')
title('Linear Slope of Alpha Power across WM Load (2, 4, 6 items)', 'FontSize', 20)
legend({sprintf('Alpha increase with load (N=%d)', n_j), ...
    sprintf('Alpha decrease with load (N=%d)', n_n), ...
    sprintf('intermediate (N=%d)', n_f), ...
    'threshold'})
box on
set(gca, 'FontSize', 15)
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_histogram_inclusion.png'));

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
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_TFR.png'));

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
gaze_ready = false;
gaze_cache_file = fullfile(feat_dir, 'sternberg_split_alphaLoads_gaze.mat');

if ~isfile(gaze_cache_file)
    fprintf('Gaze cache missing. Building gaze heatmaps from dataET_sternberg and saving to %s.\n', gaze_cache_file);
    gaze_dwell_time = build_sternberg_split_alphaLoads_gaze_from_raw_et(subjects, path);
    save(gaze_cache_file, 'gaze_dwell_time', 'subjects', '-v7.3');
end

if isfile(gaze_cache_file)
    tic
    disp(upper('Loading gaze cache...'))
    S = load(gaze_cache_file, 'gaze_dwell_time');
    toc
    if ~isfield(S, 'gaze_dwell_time')
        error('Gaze cache found at %s but does not contain `gaze_dwell_time`.', gaze_cache_file);
    end
    gaze_dwell_time = S.gaze_dwell_time;

    % Unpack raw
    allgazebase2 = gaze_dwell_time.allgazebase2;
    allgazebase4 = gaze_dwell_time.allgazebase4;
    allgazebase6 = gaze_dwell_time.allgazebase6;
    allgazetasklate2 = gaze_dwell_time.allgazetasklate2;
    allgazetasklate4 = gaze_dwell_time.allgazetasklate4;
    allgazetasklate6 = gaze_dwell_time.allgazetasklate6;

    % Convert from dwell proportion (computeGazeHeatmap output) to rough dwell time [s/bin]
    % NOTE: computeGazeHeatmap normalizes by window duration; multiplying by the nominal
    % window length yields an approximate dwell-time map (still smoothed).
    t_base_len = 0.25; % [-0.5 -0.25]
    t_late_len = 1.0;  % [1 2]

    allgazebase2 = scale_gaze_cells_powspctrm(allgazebase2, t_base_len);
    allgazebase4 = scale_gaze_cells_powspctrm(allgazebase4, t_base_len);
    allgazebase6 = scale_gaze_cells_powspctrm(allgazebase6, t_base_len);

    allgazetasklate2 = scale_gaze_cells_powspctrm(allgazetasklate2, t_late_len);
    allgazetasklate4 = scale_gaze_cells_powspctrm(allgazetasklate4, t_late_len);
    allgazetasklate6 = scale_gaze_cells_powspctrm(allgazetasklate6, t_late_len);

    gaze_ready = true;
else
    fprintf('Skipping gaze density: %s not found.\n', gaze_cache_file);
end

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
% Robust symmetric z-limits for subtraction baseline maps
allGazeDiff = [ ...
    ga2jensen_gaze.powspctrm(:); ...
    ga4jensen_gaze.powspctrm(:); ...
    ga6jensen_gaze.powspctrm(:); ...
    ga2nback_gaze.powspctrm(:); ...
    ga4nback_gaze.powspctrm(:); ...
    ga6nback_gaze.powspctrm(:) ];
zlimAbs = max(abs(prctile(allGazeDiff, [1 99])));
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
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,3);
ft_singleplotTFR(cfg, ga4jensen_gaze);
title('Amplification: WM Load 4', 'FontSize', fontSize);
ax = gca;
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,5);
ft_singleplotTFR(cfg, ga6jensen_gaze);
title('Amplification: WM Load 6', 'FontSize', fontSize);
ax = gca;
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,2);
ft_singleplotTFR(cfg, ga2nback_gaze);
title('Reduction: WM Load 2', 'FontSize', fontSize);
ax = gca;
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,4);
ft_singleplotTFR(cfg, ga4nback_gaze);
title('Reduction: WM Load 4', 'FontSize', fontSize);
ax = gca;
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

subplot(3,2,6);
ft_singleplotTFR(cfg, ga6nback_gaze);
title('Reduction: WM Load 6', 'FontSize', fontSize);
ax = gca;
xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
set(ax, 'FontSize', fontSize);
c = colorbar(ax);
c.Label.String = 'Gaze Density [a.u.]';
c.Label.FontSize = fontSize - 2;
c.FontSize = fontSize - 2;

colormap(gcf, color_map);
drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_subtractionBaseline.png'));

%% Compute CBPT stats: paired baseline-vs-task tests separately for each load (2/4/6) and each subgroup (increase vs decrease)
clc
cbpt_file_gaze_taskVsBase = fullfile(feat_dir, 'AOC_split_AlphaLoads_CBPT_gaze_tasklate_minus_base.mat');
% if isfile(cbpt_file_gaze_taskVsBase)
%     disp('Loading cached CBPT: gaze tasklate vs baseline...')
%     load(cbpt_file_gaze_taskVsBase);
% else
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
    cfg.alpha            = 0.05;
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
    clc; disp(upper('[stat_inc_2 JENSEN]')); [stat_inc_2] = ft_freqstatistics(cfg, allgazetasklate2{idx_jensen},allgazebase2{idx_jensen});
    clc; disp(upper('[stat_inc_4 JENSEN]')); [stat_inc_4] = ft_freqstatistics(cfg, allgazetasklate4{idx_jensen},allgazebase4{idx_jensen});
    clc; disp(upper('[stat_inc_6 JENSEN]')); [stat_inc_6] = ft_freqstatistics(cfg, allgazetasklate6{idx_jensen},allgazebase6{idx_jensen});
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
    clc; disp(upper('[stat_inc_n_2 N-back]')); [stat_inc_n_2] = ft_freqstatistics(cfg, allgazetasklate2{idx_nback},allgazebase2{idx_nback});
    clc; disp(upper('[stat_inc_n_4 N-back]')); [stat_inc_n_4] = ft_freqstatistics(cfg, allgazetasklate4{idx_nback},allgazebase4{idx_nback});
    clc; disp(upper('[stat_inc_n_6 N-back]')); [stat_inc_n_6] = ft_freqstatistics(cfg, allgazetasklate6{idx_nback},allgazebase6{idx_nback});
    % Hide non-significant pixels in visualization
    stat_inc_n_2.stat(stat_inc_n_2.mask==0) = NaN;
    stat_inc_n_4.stat(stat_inc_n_4.mask==0) = NaN;
    stat_inc_n_6.stat(stat_inc_n_6.mask==0) = NaN;

    save(cbpt_file_gaze_taskVsBase, ...
        'stat_inc_2', 'stat_inc_4', 'stat_inc_6', ...
        'stat_inc_n_2', 'stat_inc_n_4', 'stat_inc_n_6', ...
        'cfg_cbpt_gaze_taskVsBase_jensen', 'design_cbpt_gaze_taskVsBase_jensen', ...
        'cfg_cbpt_gaze_taskVsBase_nback', 'design_cbpt_gaze_taskVsBase_nback', '-v7.3');
% end

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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
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
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Label.String = 't-value';
c.Label.FontSize = 18;   % optional
title('Reduction: WM Load 6', 'FontSize', fontSize, 'Interpreter', 'none')
colormap(color_map);

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
cbpt_file_gaze_omnibus = fullfile(feat_dir, 'AOC_split_AlphaLoads_CBPT_gaze_omnibus_loadEffect.mat');

% if isfile(cbpt_file_gaze_omnibus)
%     disp('Loading cached CBPT: gaze omnibus load effect...')
%     load(cbpt_file_gaze_omnibus);
% else
    cfg_cbpt_gaze_omnibus_jensen = cfg;
    [statF_gaze] = ft_freqstatistics(cfg, load2_gaze{idx_jensen}, load4_gaze{idx_jensen}, load6_gaze{idx_jensen});
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
    [statF_gaze_n] = ft_freqstatistics(cfg, load2_gaze{idx_nback}, load4_gaze{idx_nback}, load6_gaze{idx_nback});
    statF_gaze_n.stat(statF_gaze_n.mask==0)=0;% set everything not relevant to zero

    save(cbpt_file_gaze_omnibus, ...
        'statF_gaze', 'statF_gaze_n', ...
        'cfg_cbpt_gaze_omnibus_jensen', 'cfg_cbpt_gaze_omnibus_nback', '-v7.3');
%end

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

ax2 = subplot(2,1,2);
ft_singleplotTFR(cfg, statF_gaze_n);
title('Reduction: Omnibus Load Effect', 'FontSize', fontSize, 'Interpreter', 'none');

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
funcs_paths = {fullfile(fileparts(mfilename('fullpath')), '..', 'funcs'), '/Volumes/Homestore/OCC/arne/funcs'};
has_loadtrend = false;
has_loadquadratic = false;
for pp = 1:numel(funcs_paths)
    if isfolder(funcs_paths{pp})
        addpath(funcs_paths{pp});
        if exist('ft_statfun_loadtrend', 'file')
            has_loadtrend = true;
        end
        if exist('ft_statfun_loadquadratic', 'file')
            has_loadquadratic = true;
        end
        if has_loadtrend && has_loadquadratic
            break;
        end
    end
end

% ---------------- Linear load trend ----------------
cbpt_file_gaze_loadtrend = fullfile(feat_dir, 'AOC_split_AlphaLoads_CBPT_gaze_linearLoadTrend.mat');
do_plot_loadtrend = false;
if isfile(cbpt_file_gaze_loadtrend)
    disp('Loading cached CBPT: gaze linear load trend...')
    load(cbpt_file_gaze_loadtrend);
    do_plot_loadtrend = true;
elseif has_loadtrend
    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_loadtrend'; % positive gaze increase with load at that location, negative gaze decrease with load at that location
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;

    cfg.clusterstatistic = 'maxsum';

    cfg.tail             = 0;
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
    cfg_cbpt_gaze_loadtrend_jensen = cfg;
    [statF_gaze] = ft_freqstatistics(cfg, load2_gaze{idx_jensen}, load4_gaze{idx_jensen}, load6_gaze{idx_jensen});
    statF_gaze.stat(statF_gaze.mask==0)=0;% set everything not relevant to zero


    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_loadtrend'; % positive gaze increase with load at that location, negative gaze decrease with load at that location
    % cfg.statistic        = 'ft_statfun_loadquadratic'; % positive U shape lowest at the middle load, negative inverted U shape peaks at middle load

    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;

    cfg.clusterstatistic = 'maxsum';

    cfg.tail             = 0;
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
    cfg_cbpt_gaze_loadtrend_nback = cfg;
    [statF_gaze_n] = ft_freqstatistics(cfg, load2_gaze{idx_nback}, load4_gaze{idx_nback}, load6_gaze{idx_nback});
    statF_gaze_n.stat(statF_gaze_n.mask==0)=0;% set everything not relevant to zero
    save(cbpt_file_gaze_loadtrend, ...
        'statF_gaze', 'statF_gaze_n', ...
        'cfg_cbpt_gaze_loadtrend_jensen', 'cfg_cbpt_gaze_loadtrend_nback', '-v7.3');
    do_plot_loadtrend = true;
end

if do_plot_loadtrend
    %%
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle        = 'outline';
    cfg.zlim = 'maxabs';
    % cfg.colormap = 'YlOrRd';
    cfg.zlim = [-5 5];
    % cfg.xlim =[300 500];
    % cfg.ylim =[200 400];
    cfg.figure = 'gcf';
    figure('Position', [0 0 1512 982], 'Color', 'w');
    subplot(3,2,1);
    ft_singleplotTFR(cfg,statF_gaze);
    title('Amplification: Linear Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    subplot(3,2,2);
    ft_singleplotTFR(cfg,statF_gaze_n);
    title('Reduction: Linear Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;
    %% zoom
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle        = 'outline';
    cfg.zlim = 'maxabs';
    % cfg.colormap = 'YlOrRd';
    cfg.zlim = [-5 5];
    cfg.xlim =[300 500];
    cfg.ylim =[200 400];
    cfg.figure = 'gcf';

    subplot(3,2,3);
    ft_singleplotTFR(cfg,statF_gaze);
    title('Amplification: Linear Load Trend (Zoom)', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    subplot(3,2,4);
    ft_singleplotTFR(cfg,statF_gaze_n);
    title('Reduction: Linear Load Trend (Zoom)', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;
    drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadTrend.png'));
    close(gcf);
end

if ~has_loadtrend
    fprintf('Skipping linear trend: ft_statfun_loadtrend not found.\n');
end

% ---------------- Quadratic load trend ----------------
cbpt_file_gaze_loadquadratic = fullfile(feat_dir, 'AOC_split_AlphaLoads_CBPT_gaze_quadraticLoadTrend.mat');
do_plot_loadquadratic = false;
if isfile(cbpt_file_gaze_loadquadratic)
    disp('Loading cached CBPT: gaze quadratic load trend...')
    load(cbpt_file_gaze_loadquadratic);
    do_plot_loadquadratic = true;
elseif has_loadquadratic
    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_loadquadratic'; % positive U-shape lowest at mid load, negative inverted U-shape peaks at mid load

    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;

    cfg.clusterstatistic = 'maxsum';

    cfg.tail             = 0;
    cfg.clustertail      = cfg.tail;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;

    subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
    n_2  = subj;
    n_4  = subj;
    n_6  =  subj;

    cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2, ones(1,n_6)*3];
    cfg.design(2,:)           = [1:n_2, 1:n_4, 1:n_6];
    cfg.ivar                  = 1;
    cfg.uvar                  = 2;
    cfg_cbpt_gaze_loadquadratic_jensen = cfg;
    [statF_gaze_quad] = ft_freqstatistics(cfg, load2_gaze{idx_jensen}, load4_gaze{idx_jensen}, load6_gaze{idx_jensen});
    statF_gaze_quad.stat(statF_gaze_quad.mask==0)=0; % set everything not relevant to zero

    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_loadquadratic';

    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;

    cfg.clusterstatistic = 'maxsum';

    cfg.tail             = 0;
    cfg.clustertail      = cfg.tail;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;

    subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
    n_2  = subj;
    n_4  = subj;
    n_6  = subj;

    cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2, ones(1,n_6)*3];
    cfg.design(2,:)           = [1:n_2, 1:n_4, 1:n_6];
    cfg.ivar                  = 1;
    cfg.uvar                  = 2;
    cfg_cbpt_gaze_loadquadratic_nback = cfg;
    [statF_gaze_quad_n] = ft_freqstatistics(cfg, load2_gaze{idx_nback}, load4_gaze{idx_nback}, load6_gaze{idx_nback});
    statF_gaze_quad_n.stat(statF_gaze_quad_n.mask==0)=0; % set everything not relevant to zero

    save(cbpt_file_gaze_loadquadratic, ...
        'statF_gaze_quad', 'statF_gaze_quad_n', ...
        'cfg_cbpt_gaze_loadquadratic_jensen', 'cfg_cbpt_gaze_loadquadratic_nback', '-v7.3');
    do_plot_loadquadratic = true;
else
    fprintf('Skipping quadratic trend: ft_statfun_loadquadratic not found.\n');
end

if do_plot_loadquadratic
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle        = 'outline';
    cfg.zlim = 'maxabs';
    cfg.zlim = [-5 5];
    cfg.figure = 'gcf';

    figure('Position', [0 0 1512 982], 'Color', 'w');

    subplot(3,2,1);
    ft_singleplotTFR(cfg,statF_gaze_quad);
    title('Amplification: Quadratic Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    subplot(3,2,2);
    ft_singleplotTFR(cfg,statF_gaze_quad_n);
    title('Reduction: Quadratic Load Trend', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    % ---------------- zoom ----------------
    cfg         = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle        = 'outline';
    cfg.zlim = 'maxabs';
    cfg.zlim = [-5 5];
    cfg.xlim =[300 500];
    cfg.ylim =[200 400];
    cfg.figure = 'gcf';

    subplot(3,2,3);
    ft_singleplotTFR(cfg,statF_gaze_quad);
    title('Amplification: Quadratic Load Trend (Zoom)', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    subplot(3,2,4);
    ft_singleplotTFR(cfg,statF_gaze_quad_n);
    title('Reduction: Quadratic Load Trend (Zoom)', 'FontSize', fontSize, 'Interpreter', 'none');

    ax = gca;
    set(ax, 'FontSize', fontSize);
    xlabel(ax, 'Screen Width [px]', 'FontSize', fontSize);
    ylabel(ax, 'Screen Height [px]', 'FontSize', fontSize-2);
    c = colorbar(ax);
    c.LineWidth = 1;
    c.FontSize = fontSize - 2;
    c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
    c.Label.String = 't-value';
    c.Label.FontSize = fontSize - 2;

    drawnow; saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_gaze_TFR_loadQuadratic.png'));
end

clear statF_gaze_quad statF_gaze_quad_n

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
figure('Position', fig_pos, 'Color', 'w');
clf;

positions = [.3, .6, .9];
jitter = 0.05;

% LEFT: alpha increase group
subplot(2,2,1); hold on

[f,x] = ksdensity(sb2); fill(positions(1)+f*0.05,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(sb4); fill(positions(2)+f*0.05,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(sb6); fill(positions(3)+f*0.05,x,'r','FaceAlpha',0.5);

boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2)+(rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3)+(rand(size(sb6))-0.5)*jitter, sb6, 'r.');

ylabel('Reaction Time [s]')
title(sprintf('Alpha increase with load (N=%d)', sum(idx_jensen)))
set(gca,'FontSize',18); box on;
ylim([0 1.5])

% RIGHT: alpha decrease group
subplot(2,2,2); hold on

[f,x] = ksdensity(nb1); fill(positions(1)+f*0.03,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(nb2); fill(positions(2)+f*0.03,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(nb3); fill(positions(3)+f*0.03,x,'r','FaceAlpha',0.5);

boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(nb1))-0.5)*jitter, nb1, 'k.');
scatter(positions(2)+(rand(size(nb2))-0.5)*jitter, nb2, 'b.');
scatter(positions(3)+(rand(size(nb3))-0.5)*jitter, nb3, 'r.');

ylabel('Reaction Time [s]')
title(sprintf('Alpha decrease with load (N=%d)', sum(idx_nback)))
set(gca,'FontSize',18); box on;
ylim([0 1.5])
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
% figure(31); clf;

positions = [.3, .6, .9];
jitter = 0.05;

% ================= LEFT: alpha increase =================
subplot(2,2,3); hold on

[f,x] = ksdensity(sb2); fill(positions(1)+f*3,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(sb4); fill(positions(2)+f*2,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(sb6); fill(positions(3)+f*2,x,'r','FaceAlpha',0.5);

boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2)+(rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3)+(rand(size(sb6))-0.5)*jitter, sb6, 'r.');

ylabel('Accuracy [%]')
title(sprintf('Alpha increase with load (N=%d)', sum(idx_jensen)))
set(gca,'FontSize',18); box on;

ylim([50 110])   % adjust if needed


% ================= RIGHT: alpha decrease =================
subplot(2,2,4); hold on

[f,x] = ksdensity(nb1); fill(positions(1)+f*3,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(nb2); fill(positions(2)+f*3,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(nb3); fill(positions(3)+f*3,x,'r','FaceAlpha',0.5);

boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(nb1))-0.5)*jitter, nb1, 'k.');
scatter(positions(2)+(rand(size(nb2))-0.5)*jitter, nb2, 'b.');
scatter(positions(3)+(rand(size(nb3))-0.5)*jitter, nb3, 'r.');

ylabel('Accuracy [%]')
title(sprintf('Alpha decrease with load (N=%d)', sum(idx_nback)))
set(gca,'FontSize',18); box on;

ylim([50 110])   % same scale for comparison

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
thr = 0.015;

idx_jensen = slope > thr;
idx_nback  = slope < -thr;

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
lme_RT = fitlme(tbl, 'RT ~ Load * Group + (1|Subject)');
disp(anova(lme_RT))

disp('--- ACC model ---')
lme_ACC = fitlme(tbl, 'ACC ~ Load * Group + (1|Subject)');
disp(anova(lme_ACC))


%% ================= TOST ANALYSIS =================
delta_RT = 0.05;   % 50 ms
delta_ACC = 5;     % 5%

RT_jensen = [RT2(idx_jensen); RT4(idx_jensen); RT6(idx_jensen)];
RT_nback  = [RT2(idx_nback);  RT4(idx_nback);  RT6(idx_nback)];

ACC_jensen = [ACC2(idx_jensen); ACC4(idx_jensen); ACC6(idx_jensen)];
ACC_nback  = [ACC2(idx_nback);  ACC4(idx_nback);  ACC6(idx_nback)];

fprintf('\n--- TOST (strict) ---\n')

[p1, p2, eq_RT] = tost_welch(RT_jensen, RT_nback, delta_RT);
fprintf('RT equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_RT);

[p1, p2, eq_ACC] = tost_welch(ACC_jensen, ACC_nback, delta_ACC);
fprintf('ACC equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_ACC);

%% ================= RT SENSITIVITY =================
fprintf('\n--- RT equivalence sensitivity ---\n')

tost_welch(RT_jensen, RT_nback, 0.05)
tost_welch(RT_jensen, RT_nback, 0.1)
tost_welch(RT_jensen, RT_nback, 0.15)

%% ================= ACC ROBUSTNESS (optional) =================
fprintf('\n--- ACC robustness ---\n')

tost_welch(ACC_jensen, ACC_nback, 3)
tost_welch(ACC_jensen, ACC_nback, 5)
tost_welch(ACC_jensen, ACC_nback, 10)
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
smoothing = 10;
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
freq_raw = computeGazeHeatmap(data3, numel(x_grid), smoothing);
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

function cells_out = scale_gaze_cells_powspctrm(cells_in, scale_factor)
cells_out = cells_in;
for ii = 1:numel(cells_out)
    if isempty(cells_out{ii}) || ~isfield(cells_out{ii}, 'powspctrm')
        continue
    end
    cells_out{ii}.powspctrm = cells_out{ii}.powspctrm .* scale_factor;
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