%% AOC Omnibus — Figures for N-back and Sternberg (plain: raw TFR, no FOOOF)
%
% WORKFLOW:
%   1. Load TFR data per condition (load2, load4, load6; load1nb, load2nb, load3nb)
%   2. Time-average [0-2s] → power spectra (chan × freq, no time dimension)
%   3. F-test on power spectra → identify significant electrodes via mask
%   4. Extract electrodes where mask==1 in alpha band [8-14 Hz]
%   5. Visualize: topoplot of F-values, TFR of grand avg at sig electrodes, rainclouds
%
% Uses omnibus_data.mat from AOC_omnibus_prep.m (raw TFR, dB baseline, 3–30 Hz).
% Stats: AOC_omnibus_statFnb_statFsb.mat

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

% Set up save directories (cross-platform)
if ispc
    figures_dir = 'W:\Students\Arne\AOC\figures\stats\omnibus';
    data_dir = 'W:\Students\Arne\AOC\data\features';
else
    figures_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/omnibus';
    data_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features';
end

%% Load variables
tic
disp('Loading omnibus data...')
if ispc
    load W:\Students\Arne\AOC\data\features\omnibus_data.mat
else
    load /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data.mat
end 
toc

%% Prepare power spectra for F-tests (time-averaged)
% F-test on time-averaged power spectra (freq × electrode), not full TFR
% This concentrates the effect and reduces the search space
disp('Preparing power spectra for F-tests (time-averaged)...')

% Sternberg: average over task window [0 2]s
cfg = [];
cfg.latency = [0 2];
cfg.frequency = [3 30];
cfg.avgovertime = 'yes';
for subj = 1:length(subjects)
    load2_pow{subj} = ft_selectdata(cfg, load2{subj});
    load4_pow{subj} = ft_selectdata(cfg, load4{subj});
    load6_pow{subj} = ft_selectdata(cfg, load6{subj});
    % Update dimord and remove time field for power spectra
    load2_pow{subj}.dimord = 'chan_freq';
    load4_pow{subj}.dimord = 'chan_freq';
    load6_pow{subj}.dimord = 'chan_freq';
    if isfield(load2_pow{subj}, 'time'), load2_pow{subj} = rmfield(load2_pow{subj}, 'time'); end
    if isfield(load4_pow{subj}, 'time'), load4_pow{subj} = rmfield(load4_pow{subj}, 'time'); end
    if isfield(load6_pow{subj}, 'time'), load6_pow{subj} = rmfield(load6_pow{subj}, 'time'); end
end

% N-back: average over task window [0 2]s
cfg.latency = [0 2];
for subj = 1:length(subjects)
    load1nb_pow{subj} = ft_selectdata(cfg, load1nb{subj});
    load2nb_pow{subj} = ft_selectdata(cfg, load2nb{subj});
    load3nb_pow{subj} = ft_selectdata(cfg, load3nb{subj});
    % Update dimord and remove time field for power spectra
    load1nb_pow{subj}.dimord = 'chan_freq';
    load2nb_pow{subj}.dimord = 'chan_freq';
    load3nb_pow{subj}.dimord = 'chan_freq';
    if isfield(load1nb_pow{subj}, 'time'), load1nb_pow{subj} = rmfield(load1nb_pow{subj}, 'time'); end
    if isfield(load2nb_pow{subj}, 'time'), load2nb_pow{subj} = rmfield(load2nb_pow{subj}, 'time'); end
    if isfield(load3nb_pow{subj}, 'time'), load3nb_pow{subj} = rmfield(load3nb_pow{subj}, 'time'); end
end

%% Compute grand averages for power spectra (with keepindividual for F-tests)
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2pow = ft_freqgrandaverage(cfg, load2_pow{:});
ga_sb_4pow = ft_freqgrandaverage(cfg, load4_pow{:});
ga_sb_6pow = ft_freqgrandaverage(cfg, load6_pow{:});
ga_nb_1pow = ft_freqgrandaverage(cfg, load1nb_pow{:});
ga_nb_2pow = ft_freqgrandaverage(cfg, load2nb_pow{:});
ga_nb_3pow = ft_freqgrandaverage(cfg, load3nb_pow{:});

%% Prepare full TFR grand averages for visualization (keep time dimension)
disp('Preparing full TFR grand averages for visualization...')
cfg = [];
cfg.latency = [-.5 3];
cfg.frequency = [3 30];
for subj = 1:length(subjects)
    load2_tfr{subj} = ft_selectdata(cfg, load2{subj});
    load4_tfr{subj} = ft_selectdata(cfg, load4{subj});
    load6_tfr{subj} = ft_selectdata(cfg, load6{subj});
end
cfg.latency = [-.5 2];
for subj = 1:length(subjects)
    load1nb_tfr{subj} = ft_selectdata(cfg, load1nb{subj});
    load2nb_tfr{subj} = ft_selectdata(cfg, load2nb{subj});
    load3nb_tfr{subj} = ft_selectdata(cfg, load3nb{subj});
end

% Grand averages for TFR visualization (without keepindividual)
cfg = [];
ga_sb_2tfr = ft_freqgrandaverage(cfg, load2_tfr{:});
ga_sb_4tfr = ft_freqgrandaverage(cfg, load4_tfr{:});
ga_sb_6tfr = ft_freqgrandaverage(cfg, load6_tfr{:});
ga_nb_1tfr = ft_freqgrandaverage(cfg, load1nb_tfr{:});
ga_nb_2tfr = ft_freqgrandaverage(cfg, load2nb_tfr{:});
ga_nb_3tfr = ft_freqgrandaverage(cfg, load3nb_tfr{:});

%% Prepare neighbours for cluster-based statistics
if ispc
    load('W:\Students\Arne\toolboxes\headmodel\elec_aligned.mat');
else
    load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/elec_aligned.mat');
end
cfg = [];
cfg.method = 'distance';
cfg.elec = elec_aligned;
cfg.layout = headmodel.layANThead;
cfg.feedback = 'yes';
cfg.neighbourdist = 0.0675; % 5 neighbours per channel
close all
neighbours = ft_prepare_neighbours(cfg);

%% Compute F-tests for EEG data (on power spectra, not full TFR)
disp('Computing F-tests for EEG power spectra...');
n_subj = numel(subjects);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesFunivariate';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
cfg.design(1,:)      = [ones(1,n_subj), ones(1,n_subj)*2, ones(1,n_subj)*3];
cfg.design(2,:)      = [1:n_subj, 1:n_subj, 1:n_subj];
cfg.ivar             = 1;
cfg.uvar             = 2;

disp('  N-back F-test...');
[statFnb] = ft_freqstatistics(cfg, ga_nb_1pow, ga_nb_2pow, ga_nb_3pow);

disp('  Sternberg F-test...');
[statFsb] = ft_freqstatistics(cfg, ga_sb_2pow, ga_sb_4pow, ga_sb_6pow);

save(fullfile(data_dir, 'AOC_omnibus_statFnb_statFsb.mat'), 'statFnb', 'statFsb');
disp('Stats saved.')

%% Identify significant electrodes from F-test results
% Extract electrodes directly from mask in alpha frequency range [8-14 Hz]

% Load
% load(fullfile(data_dir, 'AOC_omnibus_statFnb_statFsb.mat'), 'statFnb', 'statFsb');

disp(upper('Identifying significant electrodes from F-test results (alpha band 8-14 Hz)...'));

% Extract channels where mask == 1 in alpha band (any frequency bin significant)
% For power spectra, mask is chan x freq (no time dimension)
sb_sig_channels = {};
mask_sb = statFsb.mask;
freq_idx_sb = statFsb.freq >= 8 & statFsb.freq <= 14;
for ch = 1:length(statFsb.label)
    m = mask_sb(ch, freq_idx_sb);
    if any(m(:) > 0.5)
        sb_sig_channels{end+1} = statFsb.label{ch};
    end
end

nb_sig_channels = {};
mask_nb = statFnb.mask;
freq_idx_nb = statFnb.freq >= 8 & statFnb.freq <= 14;
for ch = 1:length(statFnb.label)
    m = mask_nb(ch, freq_idx_nb);
    if any(m(:) > 0.5)
        nb_sig_channels{end+1} = statFnb.label{ch};
    end
end

% Output chosen channels and counts for each task
fprintf('\nSternberg: %d significant channel(s): %s\n', length(sb_sig_channels), strjoin(sb_sig_channels, ', '));
fprintf('N-back:    %d significant channel(s): %s\n\n', length(nb_sig_channels), strjoin(nb_sig_channels, ', '));

% Topoplots (same F-test colormap as TFR: white→yellow→orange→red)
close all
key = [0 0.33 0.66 1];
key_rgb = [1 1 1; 1 1 0; 1 0.65 0; 1 0 0];
cmap_f_topo = interp1(key, key_rgb, linspace(0, 1, 64));
cfg_avg = [];
cfg_avg.parameter = 'stat';
cfg_avg.frequency = [8 14];
cfg_avg.avgoverfreq = 'yes';
statFnb_avg = ft_selectdata(cfg_avg, statFnb);
statFsb_avg = ft_selectdata(cfg_avg, statFsb);
statFnb_avg.stat = abs(statFnb_avg.stat);
statFsb_avg.stat = abs(statFsb_avg.stat);
cfg_topo = [];
cfg_topo.layout = headmodel.layANThead;
cfg_topo.parameter = 'stat';
cfg_topo.colormap = cmap_f_topo;
cfg_topo.marker = 'off';
cfg_topo.highlight = 'on';
cfg_topo.highlightsymbol = '.';
cfg_topo.highlightsize = 14;
cfg_topo.comment = 'no';
figure('Color', 'w', 'Position', [0, 0, 600, 600]);
cfg_topo.highlightchannel = nb_sig_channels;
ft_topoplotER(cfg_topo, statFnb_avg);
clim([0 max(statFnb_avg.stat(:))]);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16; title(cb, 'F-values');
text(-0.45, 0.1, nb_sig_channels)
title('N-back Significant Electrodes');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_topo_f-stats_nback.png'));
figure('Color', 'w', 'Position', [600, 0, 600, 600]);
cfg_topo.highlightchannel = sb_sig_channels;
ft_topoplotER(cfg_topo, statFsb_avg);
clim([0 max(statFnb_avg.stat(:))]);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16; title(cb, 'F-values');
text(-0.45, 0.1, sb_sig_channels)
title('Sternberg Significant Electrodes');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_topo_f-stats_stern.png'));

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
        % Select trials first, then check time range
        dataetnan_trl = ft_selectdata(cfg, dataetnan);
        % Skip if no trials for this condition (avoids ft_selectdata/ft_checkdata error on empty)
        if isempty(dataetnan_trl.trial) || numel(dataetnan_trl.trial) == 0
            eval(sprintf('allgazetask_sb%d{s} = [];', cond));
            eval(sprintf('allgazebase_sb%d{s} = [];', cond));
            continue;
        end
        % Check available time range for task window
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            all_times = cellfun(@(x) [min(x), max(x)], dataetnan_trl.time, 'UniformOutput', false);
            all_times_mat = cell2mat(all_times(:));
            min_time = min(all_times_mat(:,1));
            max_time = max(all_times_mat(:,2));
            task_time_range = [max(0, min_time), min(2, max_time)];
        else
            task_time_range = [0 2];  % Default range
        end
        cfg.latency = task_time_range;
        cfg = rmfield(cfg, 'trials');  % selection by time only; trials already subset in dataetnan_trl
        dat_task = ft_selectdata(cfg, dataetnan_trl);  % Consistent [0 2] window
        % Check available time range for baseline window
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            base_time_range = [max(-.5, min_time), min(-.25, max_time)];
            if base_time_range(1) < base_time_range(2)
                cfg.latency = base_time_range;
                dat_base = ft_selectdata(cfg, dataetnan_trl);  % Consistent baseline window
            else
                % If baseline window not available, create empty structure
                dat_base = dat_task;
                dat_base.trial = cell(1, numel(dat_task.trial));
                dat_base.time = cell(1, numel(dat_task.trial));
                for i = 1:numel(dat_task.trial)
                    dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                    dat_base.time{i} = [];
                end
            end
        else
            % If no time data, create empty baseline structure
            dat_base = dat_task;
            dat_base.trial = cell(1, numel(dat_task.trial));
            dat_base.time = cell(1, numel(dat_task.trial));
            for i = 1:numel(dat_task.trial)
                dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                dat_base.time{i} = [];
            end
        end
        
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
    
    % Trial types for N-back (trigger codes 21/22/23 = 1/2/3-back)
    if size(dataet.trialinfo, 2) == 1
        trl1 = find(dataet.trialinfo == 21);
        trl2 = find(dataet.trialinfo == 22);
        trl3 = find(dataet.trialinfo == 23);
    else
        trl1 = find(dataet.trialinfo(:,1) == 21);
        trl2 = find(dataet.trialinfo(:,1) == 22);
        trl3 = find(dataet.trialinfo(:,1) == 23);
    end
    
    % Process each condition
    for cond_idx = 1:3
        if cond_idx == 1, trl = trl1; cond = 1;
        elseif cond_idx == 2, trl = trl2; cond = 2;
        else, trl = trl3; cond = 3;
        end
        
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        % Select trials first, then check time range
        dataetnan_trl = ft_selectdata(cfg, dataetnan);
        % Skip if no trials for this condition (avoids ft_selectdata/ft_checkdata error on empty)
        if isempty(dataetnan_trl.trial) || numel(dataetnan_trl.trial) == 0
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        % Check available time range for task window
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            all_times = cellfun(@(x) [min(x), max(x)], dataetnan_trl.time, 'UniformOutput', false);
            all_times_mat = cell2mat(all_times(:));
            min_time = min(all_times_mat(:,1));
            max_time = max(all_times_mat(:,2));
            task_time_range = [max(0, min_time), min(2, max_time)];
        else
            task_time_range = [0 2];  % Default range
        end
        % Skip if no overlap between [0 2] and data (avoids 0 trials -> ft_checkdata error)
        if task_time_range(1) >= task_time_range(2)
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        cfg.latency = task_time_range;
        cfg = rmfield(cfg, 'trials');  % selection by time only; trials already subset in dataetnan_trl
        dat_task = ft_selectdata(cfg, dataetnan_trl);  % Consistent [0 2] window
        % Skip if latency selection left no trials (robustness)
        if numel(dat_task.trial) == 0
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        % Check available time range for baseline window
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            base_time_range = [max(-.5, min_time), min(-.25, max_time)];
            if base_time_range(1) < base_time_range(2)
                cfg.latency = base_time_range;
                dat_base = ft_selectdata(cfg, dataetnan_trl);  % Consistent baseline window
            else
                % If baseline window not available, create empty structure
                dat_base = dat_task;
                dat_base.trial = cell(1, numel(dat_task.trial));
                dat_base.time = cell(1, numel(dat_task.trial));
                for i = 1:numel(dat_task.trial)
                    dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                    dat_base.time{i} = [];
                end
            end
        else
            % If no time data, create empty baseline structure
            dat_base = dat_task;
            dat_base.trial = cell(1, numel(dat_task.trial));
            dat_base.time = cell(1, numel(dat_task.trial));
            for i = 1:numel(dat_task.trial)
                dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                dat_base.time{i} = [];
            end
        end
        
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
    disp(upper('Computing F-test for gaze data - Sternberg (cluster-based permutation analysis)...'));
    tic
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
    toc
    disp(' ');
else
    warning('Insufficient gaze data for Sternberg F-test');
    statFgaze_sb = [];
end

% Compute F-test for gaze data (N-back)
if ~isempty(nb1_gaze) && ~isempty(nb2_gaze) && ~isempty(nb3_gaze)
    disp(upper('Computing F-test for gaze data - N-back (cluster-based permutation analysis)...'));
    tic
    n_U = length(nb1_gaze);
    n_P = length(nb2_gaze);
    n_N = length(nb3_gaze);
    
    design_gaze = zeros(2, 3*n_U);
    design_gaze(1,:) = [ones(1,n_U), ones(1,n_P)*2, ones(1,n_N)*3];
    design_gaze(2,:) = [1:n_U, 1:n_P, 1:n_N];
    cfg_gaze.design = design_gaze;
    
    [statFgaze_nb] = ft_freqstatistics(cfg_gaze, nb1_gaze{:}, nb2_gaze{:}, nb3_gaze{:});
    statFgaze_nb.stat(statFgaze_nb.mask==0) = 0;
    toc
    disp(' ');
else
    warning('Insufficient gaze data for N-back F-test');
    statFgaze_nb = [];
end

%% Create figures
clc
close all
disp('Creating figures...');
% F-test colormap
key = [0 0.33 0.66 1];
key_rgb = [1 1 1; 1 1 0; 1 0.65 0; 1 0 0];
cmap_f = interp1(key, key_rgb, linspace(0, 1, 64));

% --- Sternberg: EEG TFR (using full TFR grand averages for visualization) ---
% Compute high-low difference from TFR grand averages
cfg_diff = [];
cfg_diff.operation = 'subtract';
cfg_diff.parameter = 'powspctrm';
ga_sb_diff = ft_math(cfg_diff, ga_sb_6tfr, ga_sb_2tfr);

figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
cfg = [];
cfg.channel = sb_sig_channels;
cfg.avgoverchan = 'yes';
cfg.frequency = [4 30];
actual_time_range = [max(-.5, ga_sb_diff.time(1)), min(2, ga_sb_diff.time(end))];
cfg.latency = actual_time_range;
freq_sb = ft_selectdata(cfg, ga_sb_diff);
meanpow = squeeze(mean(freq_sb.powspctrm, 1));
tim_interp = linspace(freq_sb.time(1), freq_sb.time(end), 500);
freq_interp = linspace(4, 30, 500);
[tim_grid_orig, freq_grid_orig] = meshgrid(freq_sb.time, freq_sb.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');
ft_plot_matrix(flip(pow_interp));
ax = gca; hold(ax, 'on');
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k--', 'LineWidth', 1);
xticks(round(interp1(tim_interp, 1:numel(tim_interp), [-.5 0 .5 1 1.5 2])));
xticklabels({'-0.5','0','0.5','1','1.5','2'});
yticks(round(interp1([5 30], [1 500], [5 10 15 20 25 30])));
yticklabels({'30','25','20','15','10','5'});  % labels only: top=30 Hz, bottom=5 Hz (data unchanged)
ylim([1 500]);
set(gca, 'Fontsize', 18);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
max_abs = max(abs(pow_interp(:)));
clim([-max_abs max_abs]);
colormap(cmap_f);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16;
title(cb, 'Power (dB)');
title('EEG TFR: Sternberg High-Low');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_EEG_TFR.png'));

% --- Sternberg: ET TFR (white→red colormap) ---
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
if ~isempty(statFgaze_sb)
    cfg = []; cfg.parameter = 'stat'; cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
    cfg.colormap = cmap_f;
    ft_singleplotTFR(cfg, statFgaze_sb);
    set(gca, 'Fontsize', 18);
    xlabel('x [px]'); ylabel('y [px]');
    c = colorbar; c.LineWidth = 1; c.FontSize = 16; title(c, 'F-statistic');
    title('ET F-test (Sternberg)');
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET F-test (Sternberg)');
end
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_ET_TFR.png'));

% --- Sternberg: EEG raincloud (all loads 2, 4, 6) ---
cfg_rain = []; cfg_rain.channel = sb_sig_channels; cfg_rain.avgoverchan = 'yes';
cfg_rain.frequency = [8 14]; cfg_rain.latency = [0 2]; cfg_rain.avgoverfreq = 'yes'; cfg_rain.avgovertime = 'yes';
eeg_sb2 = []; eeg_sb4 = []; eeg_sb6 = [];
for subj = 1:length(load2)
    if length(load2) >= subj && length(load4) >= subj && length(load6) >= subj
        t2 = ft_selectdata(cfg_rain, load2{subj}); t4 = ft_selectdata(cfg_rain, load4{subj}); t6 = ft_selectdata(cfg_rain, load6{subj});
        if ~isempty(t2) && ~isempty(t4) && ~isempty(t6)
            eeg_sb2(end+1) = t2.powspctrm; eeg_sb4(end+1) = t4.powspctrm; eeg_sb6(end+1) = t6.powspctrm;
        end
    end
end
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
positions = [0.2, 0.5, 0.8];
for k = 1:3
    tmp = {eeg_sb2, eeg_sb4, eeg_sb6}; dat = tmp{k};
    if isempty(dat), continue; end
    [f, xi] = ksdensity(dat);
    % Close polygon at position line so density is a proper half-violin
    x_fill = [positions(k), positions(k) + f*0.04, positions(k)];
    y_fill = [xi(1), xi, xi(end)];
    fill(x_fill, y_fill, colors(k,:), 'FaceAlpha', 0.5); hold on;
end
box_h = boxplot([eeg_sb2(:), eeg_sb4(:), eeg_sb6(:)], 'Labels', {'Load 2', 'Load 4', 'Load 6'}, 'Widths', 0.04, 'Positions', positions);
set(box_h, 'LineWidth', 2, 'Color', 'k');
jitter = 0.04;
scatter(positions(1) + (rand(size(eeg_sb2))-0.5)*jitter, eeg_sb2, 'k.', 'SizeData', 100);
scatter(positions(2) + (rand(size(eeg_sb4))-0.5)*jitter, eeg_sb4, 'k.', 'SizeData', 100);
scatter(positions(3) + (rand(size(eeg_sb6))-0.5)*jitter, eeg_sb6, 'k.', 'SizeData', 100);
yline(0, 'k--');
max_abs = max(abs([eeg_sb2(:); eeg_sb4(:); eeg_sb6(:)]));
ylim([-max_abs, max_abs]*1.05);
title('EEG: Sternberg (Load 2, 4, 6)');
xlabel('WM load'); ylabel('Change from baseline (dB)');
xlim([0 1]); box on; set(gca, 'FontSize', 16);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_EEG_raincloud_matlab.png'));

% --- Sternberg: ET raincloud (all loads 2, 4, 6) ---
et_sb2 = []; et_sb4 = []; et_sb6 = [];
for subj = 1:min([length(sb2_gaze), length(sb4_gaze), length(sb6_gaze)])
    if ~isempty(sb2_gaze{subj}) && ~isempty(sb4_gaze{subj}) && ~isempty(sb6_gaze{subj})
        et_sb2(end+1) = mean(sb2_gaze{subj}.powspctrm(:), 'omitnan');
        et_sb4(end+1) = mean(sb4_gaze{subj}.powspctrm(:), 'omitnan');
        et_sb6(end+1) = mean(sb6_gaze{subj}.powspctrm(:), 'omitnan');
    end
end
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
if ~isempty(et_sb2) && ~isempty(et_sb4) && ~isempty(et_sb6)
    positions = [0.2, 0.5, 0.8];
    for k = 1:3
        tmp = {et_sb2, et_sb4, et_sb6}; dat = tmp{k};
        [f, xi] = ksdensity(dat);
        x_fill = [positions(k), positions(k) + f*0.04, positions(k)];
        y_fill = [xi(1), xi, xi(end)];
        fill(x_fill, y_fill, colors(k,:), 'FaceAlpha', 0.5); hold on;
    end
    box_h = boxplot([et_sb2(:), et_sb4(:), et_sb6(:)], 'Labels', {'Load 2', 'Load 4', 'Load 6'}, 'Widths', 0.04, 'Positions', positions);
    set(box_h, 'LineWidth', 2, 'Color', 'k');
    jitter = 0.04;
    scatter(positions(1) + (rand(size(et_sb2))-0.5)*jitter, et_sb2, 'k.', 'SizeData', 100);
    scatter(positions(2) + (rand(size(et_sb4))-0.5)*jitter, et_sb4, 'k.', 'SizeData', 100);
    scatter(positions(3) + (rand(size(et_sb6))-0.5)*jitter, et_sb6, 'k.', 'SizeData', 100);
    yline(0, 'k--');
    max_abs = max(abs([et_sb2(:); et_sb4(:); et_sb6(:)]));
    ylim([-max_abs, max_abs]*1.05);
    title('ET: Sternberg (Load 2, 4, 6)');
    xlabel('WM load'); ylabel('Change from baseline (a.u.)');
    xlim([0 1]); box on; set(gca, 'FontSize', 16);
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET: Sternberg (Load 2, 4, 6)');
end
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_ET_raincloud_matlab.png'));

% --- N-back: EEG TFR (using full TFR grand averages for visualization) ---
% Compute high-low difference from TFR grand averages
cfg_diff = [];
cfg_diff.operation = 'subtract';
cfg_diff.parameter = 'powspctrm';
ga_nb_diff = ft_math(cfg_diff, ga_nb_3tfr, ga_nb_1tfr);

figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
cfg = []; cfg.channel = nb_sig_channels; cfg.avgoverchan = 'yes'; cfg.frequency = [4 30];
actual_time_range = [max(-.5, ga_nb_diff.time(1)), min(2, ga_nb_diff.time(end))]; cfg.latency = actual_time_range;
freq_nb = ft_selectdata(cfg, ga_nb_diff);
meanpow = squeeze(mean(freq_nb.powspctrm, 1));
tim_interp = linspace(freq_nb.time(1), freq_nb.time(end), 500);
freq_interp = linspace(4, 30, 500);
[tim_grid_orig, freq_grid_orig] = meshgrid(freq_nb.time, freq_nb.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');
ft_plot_matrix(flip(pow_interp));
ax = gca; hold(ax, 'on');
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k--', 'LineWidth', 1);
xticks(round(interp1(tim_interp, 1:numel(tim_interp), [-.5 0 .5 1 1.5 2]))); xticklabels({'-0.5','0','0.5','1','1.5','2'});
yticks(round(interp1([5 30], [1 500], [5 10 15 20 25 30]))); yticklabels({'30','25','20','15','10','5'});  % labels only
ylim([1 500]);
set(gca, 'Fontsize', 18); xlabel('Time [sec]'); ylabel('Frequency [Hz]');
max_abs = max(abs(pow_interp(:)));
clim([-max_abs max_abs]);
colormap(cmap_f);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16; title(cb, 'Power (dB)');
title('EEG TFR: N-back High-Low');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_EEG_TFR.png'));

% --- N-back: ET TFR (white→red colormap) ---
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
if ~isempty(statFgaze_nb)
    cfg = []; cfg.parameter = 'stat'; cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
    cfg.colormap = cmap_f;
    ft_singleplotTFR(cfg, statFgaze_nb);
    set(gca, 'Fontsize', 18); xlabel('x [px]'); ylabel('y [px]');
    c = colorbar; c.LineWidth = 1; c.FontSize = 16; title(c, 'F-statistic');
    title('ET F-test (N-back)');
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET F-test (N-back)');
end
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_ET_TFR.png'));

% --- N-back: EEG raincloud (all loads 1-, 2-, 3-back) ---
cfg_rain = []; cfg_rain.channel = nb_sig_channels; cfg_rain.avgoverchan = 'yes';
cfg_rain.frequency = [8 14]; cfg_rain.latency = [0 2]; cfg_rain.avgoverfreq = 'yes'; cfg_rain.avgovertime = 'yes';
eeg_nb1 = []; eeg_nb2 = []; eeg_nb3 = [];
for subj = 1:length(load1nb)
    if length(load1nb) >= subj && length(load2nb) >= subj && length(load3nb) >= subj
        t1 = ft_selectdata(cfg_rain, load1nb{subj}); t2 = ft_selectdata(cfg_rain, load2nb{subj}); t3 = ft_selectdata(cfg_rain, load3nb{subj});
        if ~isempty(t1) && ~isempty(t2) && ~isempty(t3)
            eeg_nb1(end+1) = t1.powspctrm; eeg_nb2(end+1) = t2.powspctrm; eeg_nb3(end+1) = t3.powspctrm;
        end
    end
end
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
positions = [0.2, 0.5, 0.8];
for k = 1:3
    tmp = {eeg_nb1, eeg_nb2, eeg_nb3}; dat = tmp{k};
    if isempty(dat), continue; end
    [f, xi] = ksdensity(dat);
    x_fill = [positions(k), positions(k) + f*0.04, positions(k)];
    y_fill = [xi(1), xi, xi(end)];
    fill(x_fill, y_fill, colors(k,:), 'FaceAlpha', 0.5); hold on;
end
box_h = boxplot([eeg_nb1(:), eeg_nb2(:), eeg_nb3(:)], 'Labels', {'1-back', '2-back', '3-back'}, 'Widths', 0.04, 'Positions', positions);
set(box_h, 'LineWidth', 2, 'Color', 'k');
jitter = 0.04;
scatter(positions(1) + (rand(size(eeg_nb1))-0.5)*jitter, eeg_nb1, 'k.', 'SizeData', 100);
scatter(positions(2) + (rand(size(eeg_nb2))-0.5)*jitter, eeg_nb2, 'k.', 'SizeData', 100);
scatter(positions(3) + (rand(size(eeg_nb3))-0.5)*jitter, eeg_nb3, 'k.', 'SizeData', 100);
yline(0, 'k--');
max_abs = max(abs([eeg_nb1(:); eeg_nb2(:); eeg_nb3(:)]));
ylim([-max_abs, max_abs]*1.05);
title('EEG: N-back (1-, 2-, 3-back)');
xlabel('N-back load'); ylabel('Change from baseline (dB)');
xlim([0 1]); box on; set(gca, 'FontSize', 16);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_EEG_raincloud_matlab.png'));

% --- N-back: ET raincloud (all loads 1-, 2-, 3-back) ---
et_nb1 = []; et_nb2 = []; et_nb3 = [];
for subj = 1:min([length(nb1_gaze), length(nb2_gaze), length(nb3_gaze)])
    if ~isempty(nb1_gaze{subj}) && ~isempty(nb2_gaze{subj}) && ~isempty(nb3_gaze{subj})
        et_nb1(end+1) = mean(nb1_gaze{subj}.powspctrm(:), 'omitnan');
        et_nb2(end+1) = mean(nb2_gaze{subj}.powspctrm(:), 'omitnan');
        et_nb3(end+1) = mean(nb3_gaze{subj}.powspctrm(:), 'omitnan');
    end
end
figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
if ~isempty(et_nb1) && ~isempty(et_nb2) && ~isempty(et_nb3)
    positions = [0.2, 0.5, 0.8];
    for k = 1:3
        tmp = {et_nb1, et_nb2, et_nb3}; dat = tmp{k};
        [f, xi] = ksdensity(dat);
        x_fill = [positions(k), positions(k) + f*0.04, positions(k)];
        y_fill = [xi(1), xi, xi(end)];
        fill(x_fill, y_fill, colors(k,:), 'FaceAlpha', 0.5); hold on;
    end
    box_h = boxplot([et_nb1(:), et_nb2(:), et_nb3(:)], 'Labels', {'1-back', '2-back', '3-back'}, 'Widths', 0.04, 'Positions', positions);
    set(box_h, 'LineWidth', 2, 'Color', 'k');
    jitter = 0.04;
    scatter(positions(1) + (rand(size(et_nb1))-0.5)*jitter, et_nb1, 'k.', 'SizeData', 100);
    scatter(positions(2) + (rand(size(et_nb2))-0.5)*jitter, et_nb2, 'k.', 'SizeData', 100);
    scatter(positions(3) + (rand(size(et_nb3))-0.5)*jitter, et_nb3, 'k.', 'SizeData', 100);
    yline(0, 'k--');
    max_abs = max(abs([et_nb1(:); et_nb2(:); et_nb3(:)]));
    ylim([-max_abs, max_abs]*1.05);
    title('ET: N-back (1-, 2-, 3-back)');
    xlabel('N-back load'); ylabel('Change from baseline (a.u.)');
    xlim([0 1]); box on; set(gca, 'FontSize', 16);
else
    text(0.5, 0.5, 'Insufficient gaze data', 'HorizontalAlignment', 'center', 'FontSize', 16);
    title('ET: N-back (1-, 2-, 3-back)');
end
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_ET_raincloud_matlab.png'));

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
