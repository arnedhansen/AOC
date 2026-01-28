%% AOC Omnibus — 2x2 Figures for N-back and Sternberg
% Creates 2x2 figures with:
%   Subplot 1: TFR visualization of F-test results for EEG data
%   Subplot 2: TFR visualization of F-test results for ET data
%   Subplot 3: Raincloud plots for EEG data condition differences (highest vs. lowest load)
%   Subplot 4: Raincloud plots for ET data condition differences (highest vs. lowest load)
%
% Key outputs:
%   AOC_omnibus_sternberg_2x2_figure.png
%   AOC_omnibus_nback_2x2_figure.png

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

% Set up save directories
figures_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/omnibus';
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end 

%% Load variables
tic
disp('Loading omnibus data...')
if ispc
    load W:\Students\Arne\AOC\data\features\omnibus_data.mat
else
    load /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data.mat
end 
toc

%% Prepare TFR data for F-tests (full time-frequency for visualization)
% Using consistent [-.5 2] second window for unbiased data-driven analysis
disp('Preparing TFR data for F-tests...')
cfg = [];
cfg.latency = [-.5 2];  % Consistent window including baseline for both tasks
for subj = 1:length(subjects)
    load2_tfr{subj} = ft_selectdata(cfg, load2{subj});
    load4_tfr{subj} = ft_selectdata(cfg, load4{subj});
    load6_tfr{subj} = ft_selectdata(cfg, load6{subj});
end

for subj = 1:length(subjects)
    load1nb_tfr{subj} = ft_selectdata(cfg, load1nb{subj});
    load2nb_tfr{subj} = ft_selectdata(cfg, load2nb{subj});
    load3nb_tfr{subj} = ft_selectdata(cfg, load3nb{subj});
end

%% Compute grand averages for TFR data (with keepindividual for F-tests)
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2tfr = ft_freqgrandaverage(cfg, load2_tfr{:});
ga_sb_4tfr = ft_freqgrandaverage(cfg, load4_tfr{:});
ga_sb_6tfr = ft_freqgrandaverage(cfg, load6_tfr{:});
ga_nb_1tfr = ft_freqgrandaverage(cfg, load1nb_tfr{:});
ga_nb_2tfr = ft_freqgrandaverage(cfg, load2nb_tfr{:});
ga_nb_3tfr = ft_freqgrandaverage(cfg, load3nb_tfr{:});

%% Prepare neighbours for cluster-based statistics
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/elec_aligned.mat');
cfg = [];
cfg.method = 'distance';
cfg.elec = elec_aligned;
cfg.layout = headmodel.layANThead;
cfg.feedback = 'yes';
cfg.neighbourdist = 40;
neighbours = ft_prepare_neighbours(cfg);

%% Compute F-tests for EEG data
disp('Computing F-tests for EEG data...')
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesFunivariate';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;
cfg.tail = 1;
cfg.clustertail = cfg.tail;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;

n_U = numel(subjects);
n_P = numel(subjects);
n_N = numel(subjects);

clear design
design = zeros(2, 3*n_U);
cfg.design(1,:) = [ones(1,n_U), ones(1,n_P)*2, ones(1,n_N)*3];
cfg.design(2,:) = [1:n_U, 1:n_P, 1:n_N];
cfg.ivar = 1;
cfg.uvar = 2;

[statFnb] = ft_freqstatistics(cfg, ga_nb_1tfr, ga_nb_2tfr, ga_nb_3tfr);
[statFsb] = ft_freqstatistics(cfg, ga_sb_2tfr, ga_sb_4tfr, ga_sb_6tfr);

%% Identify significant electrodes from F-test results
% Using alpha band (8-14 Hz) and [0 2] time window
disp('Identifying significant electrodes from F-test results...')

% For Sternberg: find electrodes with significant F-values
cfg_sig = [];
cfg_sig.parameter = 'stat';
cfg_sig.maskparameter = 'mask';
cfg_sig.frequency = [8 14];  % Alpha band
cfg_sig.latency = [0 2];  % Consistent window for both tasks
cfg_sig.avgoverfreq = 'yes';
cfg_sig.avgovertime = 'yes';
statFsb_avg = ft_selectdata(cfg_sig, statFsb);

% Find channels with significant clusters
sb_sig_channels = {};
sb_sig_mask = squeeze(statFsb_avg.mask);
sb_sig_stat = squeeze(statFsb_avg.stat);
for ch = 1:length(statFsb_avg.label)
    if sb_sig_mask(ch) > 0 && sb_sig_stat(ch) > 0
        sb_sig_channels{end+1} = statFsb_avg.label{ch};
    end
end
fprintf('Sternberg significant channels (%d): %s\n', length(sb_sig_channels), strjoin(sb_sig_channels, ', '));

% For N-back: find electrodes with significant F-values
statFnb_avg = ft_selectdata(cfg_sig, statFnb);

nb_sig_channels = {};
nb_sig_mask = squeeze(statFnb_avg.mask);
nb_sig_stat = squeeze(statFnb_avg.stat);
for ch = 1:length(statFnb_avg.label)
    if nb_sig_mask(ch) > 0 && nb_sig_stat(ch) > 0
        nb_sig_channels{end+1} = statFnb_avg.label{ch};
    end
end
fprintf('N-back significant channels (%d): %s\n', length(nb_sig_channels), strjoin(nb_sig_channels, ', '));

% Error if no significant channels found
if isempty(sb_sig_channels)
    error('No significant channels found for Sternberg. Something is severely wrong with the analysis.');
end
if isempty(nb_sig_channels)
    error('No significant channels found for N-back. Something is severely wrong with the analysis.');
end

%% Prepare gaze data for F-test analysis
disp('Preparing gaze data for F-test analysis...')

% Gaze heatmap parameters
sampling_rate = 500;
threshold = 20;
pad_ms = 150;
num_bins = 1000;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);

% Initialize gaze data structures
allgazetask_sb2 = {}; allgazetask_sb4 = {}; allgazetask_sb6 = {};
allgazebase_sb2 = {}; allgazebase_sb4 = {}; allgazebase_sb6 = {};
allgazetask_nb1 = {}; allgazetask_nb2 = {}; allgazetask_nb3 = {};
allgazebase_nb1 = {}; allgazebase_nb2 = {}; allgazebase_nb3 = {};

% Process Sternberg gaze data
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(path, subj, 'gaze');
    
    if ~exist(subj_dir, 'dir')
        warning('Subject directory %s not found. Skipping...', subj_dir);
        continue;
    end
    
    cd(subj_dir);
    
    % Load dataET_sternberg
    if ~exist('dataET_sternberg.mat', 'file')
        warning('dataET_sternberg.mat not found for subject %s. Skipping...', subj);
        continue;
    end
    load('dataET_sternberg.mat');
    
    if ~exist('dataETlong', 'var')
        warning('dataETlong variable not found for subject %s. Skipping...', subj);
        continue;
    end
    
    % Trial types
    if size(dataETlong.trialinfo, 2) == 1
        trl2 = find(dataETlong.trialinfo == 22);
        trl4 = find(dataETlong.trialinfo == 24);
        trl6 = find(dataETlong.trialinfo == 26);
    else
        trl2 = find(dataETlong.trialinfo(:,1) == 22);
        trl4 = find(dataETlong.trialinfo(:,1) == 24);
        trl6 = find(dataETlong.trialinfo(:,1) == 26);
    end
    
    % Select eye tracking channels
    if ~all(ismember({'L-GAZE-X','L-GAZE-Y'}, dataETlong.label))
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        dataETlong = ft_selectdata(cfg, dataETlong);
    end
    nTrials = numel(dataETlong.trial);
    
    % Blink correction
    dataetnan = dataETlong;
    valid_trials = true(nTrials, 1);
    
    for i = 1:nTrials
        x = dataETlong.trial{i}(1, :);
        y = dataETlong.trial{i}(2, :);
        t = dataETlong.time{i};
        
        [x_nan, y_nan, ~, ~, ~, is_valid] = removeAndInterpolateBlinks_checktrials(x, y, t, sampling_rate, threshold, pad_ms);
        
        if ~is_valid
            valid_trials(i) = false;
            continue;
        end
        
        dataetnan.trial{i}(1,:) = x_nan;
        dataetnan.trial{i}(2,:) = y_nan;
    end
    
    % Remove bad trials
    cfg = []; cfg.trials = find(valid_trials);
    dataETlong = ft_selectdata(cfg, dataETlong);
    dataetnan = ft_selectdata(cfg, dataetnan);
    
    % Recompute condition indices after cleaning
    if size(dataETlong.trialinfo, 2) == 1
        trl2 = find(dataETlong.trialinfo == 22);
        trl4 = find(dataETlong.trialinfo == 24);
        trl6 = find(dataETlong.trialinfo == 26);
    else
        trl2 = find(dataETlong.trialinfo(:,1) == 22);
        trl4 = find(dataETlong.trialinfo(:,1) == 24);
        trl6 = find(dataETlong.trialinfo(:,1) == 26);
    end
    
    % Process each condition
    for cond_idx = 1:3
        if cond_idx == 1, trl = trl2; cond = 2;
        elseif cond_idx == 2, trl = trl4; cond = 4;
        else, trl = trl6; cond = 6;
        end
        
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        cfg.latency = [0 2]; dat_task = ft_selectdata(cfg, dataetnan);  % Consistent [0 2] window
        cfg.latency = [-.5 -.25]; dat_base = ft_selectdata(cfg, dataetnan);  % Consistent baseline window
        
        % Compute heatmaps
        [freq_task, ~] = computeGazeHeatmap(dat_task, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_base, ~] = computeGazeHeatmap(dat_base, x_grid, y_grid, sampling_rate, smooth_val);
        
        eval(sprintf('allgazetask_sb%d{s} = freq_task;', cond));
        eval(sprintf('allgazebase_sb%d{s} = freq_base;', cond));
    end
end

% Process N-back gaze data
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(path, subj, 'gaze');
    
    if ~exist(subj_dir, 'dir')
        continue;
    end
    
    cd(subj_dir);
    
    if ~exist('dataET_nback.mat', 'file')
        warning('dataET_nback.mat not found for subject %s. Skipping...', subj);
        continue;
    end
    load('dataET_nback.mat');
    
    % Select eye tracking channels
    if ~all(ismember({'L-GAZE-X','L-GAZE-Y'}, dataet.label))
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        dataet = ft_selectdata(cfg, dataet);
    end
    nTrials = numel(dataet.trial);
    
    % Blink correction
    dataetnan = dataet;
    valid_trials = true(nTrials, 1);
    
    for i = 1:nTrials
        x = dataet.trial{i}(1, :);
        y = dataet.trial{i}(2, :);
        t = dataet.time{i};
        
        [x_nan, y_nan, ~, ~, ~, is_valid] = removeAndInterpolateBlinks_checktrials(x, y, t, sampling_rate, threshold, pad_ms);
        
        if ~is_valid
            valid_trials(i) = false;
            continue;
        end
        
        dataetnan.trial{i}(1,:) = x_nan;
        dataetnan.trial{i}(2,:) = y_nan;
    end
    
    cfg = []; cfg.trials = find(valid_trials);
    dataet = ft_selectdata(cfg, dataet);
    dataetnan = ft_selectdata(cfg, dataetnan);
    
    % Trial types for N-back
    if size(dataet.trialinfo, 2) == 1
        trl1 = find(dataet.trialinfo == 1);
        trl2 = find(dataet.trialinfo == 2);
        trl3 = find(dataet.trialinfo == 3);
    else
        trl1 = find(dataet.trialinfo(:,1) == 1);
        trl2 = find(dataet.trialinfo(:,1) == 2);
        trl3 = find(dataet.trialinfo(:,1) == 3);
    end
    
    % Process each condition
    for cond_idx = 1:3
        if cond_idx == 1, trl = trl1; cond = 1;
        elseif cond_idx == 2, trl = trl2; cond = 2;
        else, trl = trl3; cond = 3;
        end
        
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        cfg.latency = [0 2]; dat_task = ft_selectdata(cfg, dataetnan);
        cfg.latency = [-.5 -.25]; dat_base = ft_selectdata(cfg, dataetnan);  % Consistent baseline window
        
        % Compute heatmaps
        [freq_task, ~] = computeGazeHeatmap(dat_task, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_base, ~] = computeGazeHeatmap(dat_base, x_grid, y_grid, sampling_rate, smooth_val);
        
        eval(sprintf('allgazetask_nb%d{s} = freq_task;', cond));
        eval(sprintf('allgazebase_nb%d{s} = freq_base;', cond));
    end
end

% Remove empty cells and create baseline-corrected data
sb2_gaze = {}; sb4_gaze = {}; sb6_gaze = {};
nb1_gaze = {}; nb2_gaze = {}; nb3_gaze = {};

for subj = 1:length(subjects)
    if ~isempty(allgazebase_sb2) && length(allgazebase_sb2) >= subj && ~isempty(allgazebase_sb2{subj})
        sb2_gaze{subj} = allgazebase_sb2{subj};
        sb2_gaze{subj}.powspctrm = allgazetask_sb2{subj}.powspctrm - allgazebase_sb2{subj}.powspctrm;
        sb4_gaze{subj} = allgazebase_sb4{subj};
        sb4_gaze{subj}.powspctrm = allgazetask_sb4{subj}.powspctrm - allgazebase_sb4{subj}.powspctrm;
        sb6_gaze{subj} = allgazebase_sb6{subj};
        sb6_gaze{subj}.powspctrm = allgazetask_sb6{subj}.powspctrm - allgazebase_sb6{subj}.powspctrm;
    end
    if ~isempty(allgazebase_nb1) && length(allgazebase_nb1) >= subj && ~isempty(allgazebase_nb1{subj})
        nb1_gaze{subj} = allgazebase_nb1{subj};
        nb1_gaze{subj}.powspctrm = allgazetask_nb1{subj}.powspctrm - allgazebase_nb1{subj}.powspctrm;
        nb2_gaze{subj} = allgazebase_nb2{subj};
        nb2_gaze{subj}.powspctrm = allgazetask_nb2{subj}.powspctrm - allgazebase_nb2{subj}.powspctrm;
        nb3_gaze{subj} = allgazebase_nb3{subj};
        nb3_gaze{subj}.powspctrm = allgazetask_nb3{subj}.powspctrm - allgazebase_nb3{subj}.powspctrm;
    end
end

% Remove empty cells
sb2_gaze = sb2_gaze(~cellfun(@isempty, sb2_gaze));
sb4_gaze = sb4_gaze(~cellfun(@isempty, sb4_gaze));
sb6_gaze = sb6_gaze(~cellfun(@isempty, sb6_gaze));
nb1_gaze = nb1_gaze(~cellfun(@isempty, nb1_gaze));
nb2_gaze = nb2_gaze(~cellfun(@isempty, nb2_gaze));
nb3_gaze = nb3_gaze(~cellfun(@isempty, nb3_gaze));

% Compute F-test for gaze data (Sternberg)
if ~isempty(sb2_gaze) && ~isempty(sb4_gaze) && ~isempty(sb6_gaze)
    cfg_gaze = [];
    cfg_gaze.method = 'montecarlo';
    cfg_gaze.statistic = 'ft_statfun_depsamplesFunivariate';
    cfg_gaze.correctm = 'cluster';
    cfg_gaze.clusteralpha = 0.05;
    cfg_gaze.clusterstatistic = 'maxsum';
    cfg_gaze.neighbours = [];
    cfg_gaze.tail = 1;
    cfg_gaze.clustertail = cfg_gaze.tail;
    cfg_gaze.alpha = 0.05;
    cfg_gaze.numrandomization = 1000;
    
    n_U = length(sb2_gaze);
    n_P = length(sb4_gaze);
    n_N = length(sb6_gaze);
    
    design_gaze = zeros(2, 3*n_U);
    design_gaze(1,:) = [ones(1,n_U), ones(1,n_P)*2, ones(1,n_N)*3];
    design_gaze(2,:) = [1:n_U, 1:n_P, 1:n_N];
    cfg_gaze.design = design_gaze;
    cfg_gaze.ivar = 1;
    cfg_gaze.uvar = 2;
    
    [statFgaze_sb] = ft_freqstatistics(cfg_gaze, sb2_gaze{:}, sb4_gaze{:}, sb6_gaze{:});
    statFgaze_sb.stat(statFgaze_sb.mask==0) = 0;
else
    warning('Insufficient gaze data for Sternberg F-test');
    statFgaze_sb = [];
end

% Compute F-test for gaze data (N-back)
if ~isempty(nb1_gaze) && ~isempty(nb2_gaze) && ~isempty(nb3_gaze)
    n_U = length(nb1_gaze);
    n_P = length(nb2_gaze);
    n_N = length(nb3_gaze);
    
    design_gaze = zeros(2, 3*n_U);
    design_gaze(1,:) = [ones(1,n_U), ones(1,n_P)*2, ones(1,n_N)*3];
    design_gaze(2,:) = [1:n_U, 1:n_P, 1:n_N];
    cfg_gaze.design = design_gaze;
    
    [statFgaze_nb] = ft_freqstatistics(cfg_gaze, nb1_gaze{:}, nb2_gaze{:}, nb3_gaze{:});
    statFgaze_nb.stat(statFgaze_nb.mask==0) = 0;
else
    warning('Insufficient gaze data for N-back F-test');
    statFgaze_nb = [];
end

%% Create final 2x2 figures for both tasks
disp('Creating final 2x2 figures...');

% STERNBERG FIGURE
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');

% Subplot 1: TFR visualization of F-test results for EEG
subplot(2,2,1);
cfg = [];
cfg.channel = sb_sig_channels;
cfg.avgoverchan = 'yes';
cfg.frequency = [1 30];  % Visualization frequency range
cfg.latency = [-.5 2];  % Consistent visualization window
freq_sb = ft_selectdata(cfg, statFsb);
meanpow = squeeze(mean(freq_sb.stat, 1));  % Average over channels

tim_interp = linspace(freq_sb.time(1), freq_sb.time(end), 500);
freq_interp = linspace(1, 30, 500);  % Match visualization frequency range
[tim_grid_orig, freq_grid_orig] = meshgrid(freq_sb.time, freq_sb.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, ...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, double(squeeze(freq_sb.mask)), ...
    tim_grid_interp, freq_grid_interp, 'nearest', 0);

ft_plot_matrix(flip(pow_interp), 'highlightstyle', 'outline', 'highlight', flip(abs(round(mask_interp))));
ax = gca; hold(ax, 'on');
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);
xticks(round(interp1(tim_interp, 1:numel(tim_interp), [-.5 0 1 2])));
xticklabels({'-0.5','0','1','2'});
yticks([1 125 250 375]);
yticklabels({'30','20','10','1'});
set(gca, 'Fontsize', 18);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([0 max(pow_interp(:))]);
colormap(gca, 'YlOrRd');
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 16;
title(cb, 'F-values');
title('EEG F-test (Sternberg)');

% Subplot 2: TFR visualization of F-test results for ET
subplot(2,2,2);
if ~isempty(statFgaze_sb)
    cfg = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    cfg.colormap = 'YlOrRd';
    ft_singleplotTFR(cfg, statFgaze_sb);
    set(gca, 'Fontsize', 18);
    xlabel('x [px]');
    ylabel('y [px]');
    grid on;
    c = colorbar;
    c.LineWidth = 1;
    c.FontSize = 16;
    title(c, 'F-statistic');
    title('ET F-test (Sternberg)');
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET F-test (Sternberg)');
end

% Subplot 3: Raincloud plots for EEG (highest vs lowest load)
subplot(2,2,3);
% Extract highest vs lowest load for EEG (using significant channels)
cfg_rain = [];
cfg_rain.channel = sb_sig_channels;
cfg_rain.avgoverchan = 'yes';
cfg_rain.frequency = [8 14];  % Alpha band
cfg_rain.latency = [0 2];  % Consistent window
cfg_rain.avgoverfreq = 'yes';
cfg_rain.avgovertime = 'yes';

eeg_low_sb = [];
eeg_high_sb = [];
for subj = 1:length(load2)
    if length(load2) >= subj && length(load6) >= subj
        tmp_low = ft_selectdata(cfg_rain, load2{subj});
        tmp_high = ft_selectdata(cfg_rain, load6{subj});
        if ~isempty(tmp_low) && ~isempty(tmp_high)
            eeg_low_sb(end+1) = tmp_low.powspctrm;
            eeg_high_sb(end+1) = tmp_high.powspctrm;
        end
    end
end
% Raincloud plot for EEG
positions = [0.3, 0.7];
[f_low, xi_low] = ksdensity(eeg_low_sb);
fill(positions(1) + f_low*0.05, xi_low, [0.30, 0.75, 0.93], 'FaceAlpha', 0.5); hold on;
[f_high, xi_high] = ksdensity(eeg_high_sb);
fill(positions(2) + f_high*0.05, xi_high, [0.97, 0.26, 0.26], 'FaceAlpha', 0.5);
box_h = boxplot([eeg_low_sb(:), eeg_high_sb(:)], 'Labels', {'Low', 'High'}, 'Widths', 0.05, 'Positions', positions);
set(box_h, 'LineWidth', 2);
jitter = 0.05;
scatter(positions(1) + (rand(size(eeg_low_sb))-0.5)*jitter, eeg_low_sb, 'k.', 'SizeData', 50);
scatter(positions(2) + (rand(size(eeg_high_sb))-0.5)*jitter, eeg_high_sb, 'k.', 'SizeData', 50);
[~, p] = ttest(eeg_low_sb, eeg_high_sb);
sig_label = getSigLabel(p);
y_max = max([eeg_low_sb(:); eeg_high_sb(:)]) * 1.25;
if ~isempty(sig_label)
    line(positions, [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean(positions), y_max + 0.02, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'center');
end
title('EEG: Load 2 vs Load 6');
ylabel('Change from baseline');
xlim([0 1]);
box on;
set(gca, 'FontSize', 16);

% Subplot 4: Raincloud plots for ET (highest vs lowest load)
subplot(2,2,4);
% Extract gaze data for highest vs lowest load
et_low_sb = [];
et_high_sb = [];
for subj = 1:min(length(sb2_gaze), length(sb6_gaze))
    if ~isempty(sb2_gaze{subj}) && ~isempty(sb6_gaze{subj})
        % Average over spatial dimension
        et_low_sb(end+1) = mean(sb2_gaze{subj}.powspctrm(:), 'omitnan');
        et_high_sb(end+1) = mean(sb6_gaze{subj}.powspctrm(:), 'omitnan');
    end
end
if ~isempty(et_low_sb) && ~isempty(et_high_sb)
    % Raincloud plot for ET
    positions = [0.3, 0.7];
    [f_low, xi_low] = ksdensity(et_low_sb);
    fill(positions(1) + f_low*0.05, xi_low, [0.30, 0.75, 0.93], 'FaceAlpha', 0.5); hold on;
    [f_high, xi_high] = ksdensity(et_high_sb);
    fill(positions(2) + f_high*0.05, xi_high, [0.97, 0.26, 0.26], 'FaceAlpha', 0.5);
    box_h = boxplot([et_low_sb(:), et_high_sb(:)], 'Labels', {'Low', 'High'}, 'Widths', 0.05, 'Positions', positions);
    set(box_h, 'LineWidth', 2);
    jitter = 0.05;
    scatter(positions(1) + (rand(size(et_low_sb))-0.5)*jitter, et_low_sb, 'k.', 'SizeData', 50);
    scatter(positions(2) + (rand(size(et_high_sb))-0.5)*jitter, et_high_sb, 'k.', 'SizeData', 50);
    [~, p] = ttest(et_low_sb, et_high_sb);
    sig_label = getSigLabel(p);
    y_max = max([et_low_sb(:); et_high_sb(:)]) * 1.25;
    if ~isempty(sig_label)
        line(positions, [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
        text(mean(positions), y_max + 0.02, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
    title('ET: Load 2 vs Load 6');
    ylabel('Change from baseline');
    xlim([0 1]);
    box on;
    set(gca, 'FontSize', 16);
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET: Load 2 vs Load 6');
end

sgtitle('Sternberg: F-tests and Condition Differences', 'FontSize', 20, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_2x2_figure.png'));

% N-BACK FIGURE
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');

% Subplot 1: TFR visualization of F-test results for EEG
subplot(2,2,1);
cfg = [];
cfg.channel = nb_sig_channels;
cfg.avgoverchan = 'yes';
cfg.frequency = [1 30];  % Visualization frequency range
cfg.latency = [-.5 2];  % Consistent visualization window
freq_nb = ft_selectdata(cfg, statFnb);
meanpow = squeeze(mean(freq_nb.stat, 1));  % Average over channels

tim_interp = linspace(freq_nb.time(1), freq_nb.time(end), 500);
freq_interp = linspace(1, 30, 500);  % Match visualization frequency range
[tim_grid_orig, freq_grid_orig] = meshgrid(freq_nb.time, freq_nb.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, ...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, double(squeeze(freq_nb.mask)), ...
    tim_grid_interp, freq_grid_interp, 'nearest', 0);

ft_plot_matrix(flip(pow_interp), 'highlightstyle', 'outline', 'highlight', flip(abs(round(mask_interp))));
ax = gca; hold(ax, 'on');
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);
xticks(round(interp1(tim_interp, 1:numel(tim_interp), [-.5 0 1 2])));
xticklabels({'-0.5','0','1','2'});
yticks([1 125 250 375]);
yticklabels({'30','20','10','1'});
set(gca, 'Fontsize', 18);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([0 max(pow_interp(:))]);
colormap(gca, 'YlOrRd');
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 16;
title(cb, 'F-values');
title('EEG F-test (N-back)');

% Subplot 2: TFR visualization of F-test results for ET
subplot(2,2,2);
if ~isempty(statFgaze_nb)
    cfg = [];
    cfg.parameter = 'stat';
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    cfg.colormap = 'YlOrRd';
    ft_singleplotTFR(cfg, statFgaze_nb);
    set(gca, 'Fontsize', 18);
    xlabel('x [px]');
    ylabel('y [px]');
    grid on;
    c = colorbar;
    c.LineWidth = 1;
    c.FontSize = 16;
    title(c, 'F-statistic');
    title('ET F-test (N-back)');
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET F-test (N-back)');
end

% Subplot 3: Raincloud plots for EEG (highest vs lowest load)
subplot(2,2,3);
cfg_rain = [];
cfg_rain.channel = nb_sig_channels;
cfg_rain.avgoverchan = 'yes';
cfg_rain.frequency = [8 14];  % Alpha band
cfg_rain.latency = [0 2];  % Consistent window
cfg_rain.avgoverfreq = 'yes';
cfg_rain.avgovertime = 'yes';
eeg_low_nb = [];
eeg_high_nb = [];
for subj = 1:length(load1nb)
    if length(load1nb) >= subj && length(load3nb) >= subj
        tmp_low = ft_selectdata(cfg_rain, load1nb{subj});
        tmp_high = ft_selectdata(cfg_rain, load3nb{subj});
        if ~isempty(tmp_low) && ~isempty(tmp_high)
            eeg_low_nb(end+1) = tmp_low.powspctrm;
            eeg_high_nb(end+1) = tmp_high.powspctrm;
        end
    end
end
% Raincloud plot for EEG
positions = [0.3, 0.7];
[f_low, xi_low] = ksdensity(eeg_low_nb);
fill(positions(1) + f_low*0.05, xi_low, [0.30, 0.75, 0.93], 'FaceAlpha', 0.5); hold on;
[f_high, xi_high] = ksdensity(eeg_high_nb);
fill(positions(2) + f_high*0.05, xi_high, [0.97, 0.26, 0.26], 'FaceAlpha', 0.5);
box_h = boxplot([eeg_low_nb(:), eeg_high_nb(:)], 'Labels', {'Low', 'High'}, 'Widths', 0.05, 'Positions', positions);
set(box_h, 'LineWidth', 2);
jitter = 0.05;
scatter(positions(1) + (rand(size(eeg_low_nb))-0.5)*jitter, eeg_low_nb, 'k.', 'SizeData', 50);
scatter(positions(2) + (rand(size(eeg_high_nb))-0.5)*jitter, eeg_high_nb, 'k.', 'SizeData', 50);
[~, p] = ttest(eeg_low_nb, eeg_high_nb);
sig_label = getSigLabel(p);
y_max = max([eeg_low_nb(:); eeg_high_nb(:)]) * 1.25;
if ~isempty(sig_label)
    line(positions, [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean(positions), y_max + 0.02, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'center');
end
title('EEG: Load 1 vs Load 3');
ylabel('Change from baseline');
xlim([0 1]);
box on;
set(gca, 'FontSize', 16);

% Subplot 4: Raincloud plots for ET (highest vs lowest load)
subplot(2,2,4);
et_low_nb = [];
et_high_nb = [];
for subj = 1:min(length(nb1_gaze), length(nb3_gaze))
    if ~isempty(nb1_gaze{subj}) && ~isempty(nb3_gaze{subj})
        et_low_nb(end+1) = mean(nb1_gaze{subj}.powspctrm(:), 'omitnan');
        et_high_nb(end+1) = mean(nb3_gaze{subj}.powspctrm(:), 'omitnan');
    end
end
if ~isempty(et_low_nb) && ~isempty(et_high_nb)
    % Raincloud plot for ET
    positions = [0.3, 0.7];
    [f_low, xi_low] = ksdensity(et_low_nb);
    fill(positions(1) + f_low*0.05, xi_low, [0.30, 0.75, 0.93], 'FaceAlpha', 0.5); hold on;
    [f_high, xi_high] = ksdensity(et_high_nb);
    fill(positions(2) + f_high*0.05, xi_high, [0.97, 0.26, 0.26], 'FaceAlpha', 0.5);
    box_h = boxplot([et_low_nb(:), et_high_nb(:)], 'Labels', {'Low', 'High'}, 'Widths', 0.05, 'Positions', positions);
    set(box_h, 'LineWidth', 2);
    jitter = 0.05;
    scatter(positions(1) + (rand(size(et_low_nb))-0.5)*jitter, et_low_nb, 'k.', 'SizeData', 50);
    scatter(positions(2) + (rand(size(et_high_nb))-0.5)*jitter, et_high_nb, 'k.', 'SizeData', 50);
    [~, p] = ttest(et_low_nb, et_high_nb);
    sig_label = getSigLabel(p);
    y_max = max([et_low_nb(:); et_high_nb(:)]) * 1.25;
    if ~isempty(sig_label)
        line(positions, [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
        text(mean(positions), y_max + 0.02, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
    title('ET: Load 1 vs Load 3');
    ylabel('Change from baseline');
    xlim([0 1]);
    box on;
    set(gca, 'FontSize', 16);
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET: Load 1 vs Load 3');
end

sgtitle('N-back: F-tests and Condition Differences', 'FontSize', 20, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_2x2_figure.png'));

disp('2x2 figures created successfully!');

%% Helper functions

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

function [x_nan, y_nan, x_interp, y_interp, blink_mask, is_valid] = removeAndInterpolateBlinks_checktrials(x, y, t, fs, threshold, pad_ms)
    % Remove and interpolate blinks based on velocity threshold
    % Inputs:
    %   x, y: gaze coordinates
    %   t: time vector
    %   fs: sampling rate (Hz)
    %   threshold: velocity threshold for blink detection (pixels/sample)
    %   pad_ms: padding around blinks in milliseconds
    % Outputs:
    %   x_nan, y_nan: data with blinks set to NaN
    %   x_interp, y_interp: data with blinks interpolated
    %   blink_mask: logical mask indicating blink periods
    %   is_valid: whether trial has sufficient valid data
    
    % Initialize outputs
    x_nan = x;
    y_nan = y;
    x_interp = x;
    y_interp = y;
    
    % Detect missing/invalid data (zeros or NaN)
    invalid = ~isfinite(x) | ~isfinite(y) | (x == 0 & y == 0);
    
    % Compute velocity
    dx = diff([x(1), x]);
    dy = diff([y(1), y]);
    velocity = sqrt(dx.^2 + dy.^2);
    
    % Detect blinks based on velocity threshold
    blink_velocity = velocity > threshold;
    
    % Combine invalid data and high velocity as blink indicators
    blink_mask = invalid | blink_velocity;
    
    % Add padding around blinks
    pad_samples = round(pad_ms * fs / 1000);
    if pad_samples > 0
        blink_idx = find(blink_mask);
        for i = 1:length(blink_idx)
            start_idx = max(1, blink_idx(i) - pad_samples);
            end_idx = min(length(blink_mask), blink_idx(i) + pad_samples);
            blink_mask(start_idx:end_idx) = true;
        end
    end
    
    % Set blinks to NaN
    x_nan(blink_mask) = NaN;
    y_nan(blink_mask) = NaN;
    
    % Interpolate blinks
    valid_idx = ~blink_mask & isfinite(x) & isfinite(y);
    if sum(valid_idx) > 2
        x_interp = fillmissing(x_nan, 'linear');
        y_interp = fillmissing(y_nan, 'linear');
        
        % Edge handling: if first/last samples are NaN, use nearest valid
        if isnan(x_interp(1))
            first_valid = find(isfinite(x_interp), 1, 'first');
            if ~isempty(first_valid)
                x_interp(1:first_valid-1) = x_interp(first_valid);
                y_interp(1:first_valid-1) = y_interp(first_valid);
            end
        end
        if isnan(x_interp(end))
            last_valid = find(isfinite(x_interp), 1, 'last');
            if ~isempty(last_valid)
                x_interp(last_valid+1:end) = x_interp(last_valid);
                y_interp(last_valid+1:end) = y_interp(last_valid);
            end
        end
    else
        x_interp = x_nan;
        y_interp = y_nan;
    end
    
    % Check if trial is valid (at least 50% valid samples)
    is_valid = sum(~blink_mask) / length(blink_mask) >= 0.5;
end

function [freq_raw, freq_norm] = computeGazeHeatmap(data, x_grid, y_grid, fs, smoothing)
    pos = horzcat(data.trial{:});
    binned = histcounts2(pos(1,:), pos(2,:), x_grid, y_grid);
    dwell_time = binned / fs;
    smoothed = imgaussfilt(dwell_time, smoothing);

    freq_raw = [];
    freq_raw.powspctrm(1,:,:) = flipud(smoothed');
    freq_raw.time = x_grid(2:end);
    freq_raw.freq = y_grid(2:end);
    freq_raw.label = {'et'};
    freq_raw.dimord = 'chan_freq_time';

    norm_time = dwell_time / sum(dwell_time(:));
    norm_smooth = imgaussfilt(norm_time, smoothing);

    freq_norm = [];
    freq_norm.powspctrm(1,:,:) = flipud(norm_smooth');
    freq_norm.time = x_grid(2:end);
    freq_norm.freq = y_grid(2:end);
    freq_norm.label = {'et'};
    freq_norm.dimord = 'chan_freq_time';
end
