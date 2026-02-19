%% AOC Omnibus — Figures for N-back and Sternberg (FOOOF, 3D CBPT)
%
% WORKFLOW:
%   1. Load FOOOF TFR data per condition (load2, load4, load6; load1nb, load2nb, load3nb)
%   2. Compute grand averages with keepindividual (full TFR: chan x freq x time)
%   3. F-test on full TFR (3D CBPT: clusters over chan x freq x time)
%   4. Extract significant electrodes from mask in alpha band [8-14 Hz]
%   5. Visualize: topoplot of F-values, TFR at sig electrodes, rainclouds
%
% Uses omnibus_data_FOOOF.mat (FOOOF TFR, absolute baseline).
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
disp('Loading omnibus data (FOOOF)...')
if ispc
    load W:\Students\Arne\AOC\data\features\omnibus_data_FOOOF.mat
else
    load /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data_FOOOF.mat
end
toc

%% Prepare TFR data for F-tests
disp('Preparing TFR data for F-tests...')

% Select task window [0 2]s and freq [3 30] Hz
% Sternberg 
cfg = [];
cfg.latency = [0 2];
cfg.frequency = [3 30];
for subj = 1:length(subjects)
    load2_tfr{subj} = ft_selectdata(cfg, load2{subj});
    load4_tfr{subj} = ft_selectdata(cfg, load4{subj});
    load6_tfr{subj} = ft_selectdata(cfg, load6{subj});
end

% N-back
for subj = 1:length(subjects)
    load1nb_tfr{subj} = ft_selectdata(cfg, load1nb{subj});
    load2nb_tfr{subj} = ft_selectdata(cfg, load2nb{subj});
    load3nb_tfr{subj} = ft_selectdata(cfg, load3nb{subj});
end

%% Compute grand averages (with keepindividual for F-tests)
disp('Computing grand averages for F-tests...')
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2 = ft_freqgrandaverage(cfg, load2_tfr{:});
ga_sb_4 = ft_freqgrandaverage(cfg, load4_tfr{:});
ga_sb_6 = ft_freqgrandaverage(cfg, load6_tfr{:});
ga_nb_1 = ft_freqgrandaverage(cfg, load1nb_tfr{:});
ga_nb_2 = ft_freqgrandaverage(cfg, load2nb_tfr{:});
ga_nb_3 = ft_freqgrandaverage(cfg, load3nb_tfr{:});

%% Compute grand averages for visualization (without keepindividual, wider time window)
disp('Computing grand averages for visualization...')
cfg = [];
cfg.latency = [-.5 2];
cfg.frequency = [3 30];
for subj = 1:length(subjects)
    load2_vis{subj} = ft_selectdata(cfg, load2{subj});
    load4_vis{subj} = ft_selectdata(cfg, load4{subj});
    load6_vis{subj} = ft_selectdata(cfg, load6{subj});
    load1nb_vis{subj} = ft_selectdata(cfg, load1nb{subj});
    load2nb_vis{subj} = ft_selectdata(cfg, load2nb{subj});
    load3nb_vis{subj} = ft_selectdata(cfg, load3nb{subj});
end

cfg = [];
ga_sb_2vis = ft_freqgrandaverage(cfg, load2_vis{:});
ga_sb_4vis = ft_freqgrandaverage(cfg, load4_vis{:});
ga_sb_6vis = ft_freqgrandaverage(cfg, load6_vis{:});
ga_nb_1vis = ft_freqgrandaverage(cfg, load1nb_vis{:});
ga_nb_2vis = ft_freqgrandaverage(cfg, load2nb_vis{:});
ga_nb_3vis = ft_freqgrandaverage(cfg, load3nb_vis{:});

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
cfg.neighbourdist = 0.0675; % ca. 5 neighbours per channel
close all
neighbours = ft_prepare_neighbours(cfg);

%% Compute F-tests for EEG data (3D CBPT: chan x freq x time)
disp('Computing F-tests for EEG TFR (3D cluster-based permutation)...');
n_subj = numel(subjects);
all_labels = ga_nb_1.label;

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

disp('  N-back...');
[statFnb] = ft_freqstatistics(cfg, ga_nb_1, ga_nb_2, ga_nb_3);

disp('  Sternberg...');
[statFsb] = ft_freqstatistics(cfg, ga_sb_2, ga_sb_4, ga_sb_6);

save(fullfile(data_dir, 'AOC_omnibus_statFnb_statFsb_FOOOF.mat'), 'statFnb', 'statFsb');
disp('Stats saved.')

%% Identify significant electrodes from F-test results
clc
disp(upper('Identifying significant electrodes from F-test results (alpha band 8-14 Hz)...'));

% Define posterior ROI for electrode selection: channels with P, O, or I in name
posterior_idx = cellfun(@(x) contains(x, {'PO', 'O', 'I'}), statFnb.label);
posterior_chans = statFnb.label(posterior_idx);
fprintf('Posterior ROI: %d channels\n', length(posterior_chans));

% Extract posterior channels where mask == 1 in alpha band (any freq bin, any time point)
% For 3D CBPT, mask is chan x freq x time
sb_sig_channels = {};
mask_sb = statFsb.mask;
freq_idx_sb = statFsb.freq >= 8 & statFsb.freq <= 14;
for ch = 1:length(statFsb.label)
    if ~posterior_idx(ch), continue; end
    m = mask_sb(ch, freq_idx_sb, :);
    if any(m(:) > 0.5)
        sb_sig_channels{end+1} = statFsb.label{ch};
    end
end

nb_sig_channels = {};
mask_nb = statFnb.mask;
freq_idx_nb = statFnb.freq >= 8 & statFnb.freq <= 14;
for ch = 1:length(statFnb.label)
    if ~posterior_idx(ch), continue; end
    m = mask_nb(ch, freq_idx_nb, :);
    if any(m(:) > 0.5)
        nb_sig_channels{end+1} = statFnb.label{ch};
    end
end

% Fallback: if no significant channels, use top 10 posterior by max alpha F-stat
if isempty(sb_sig_channels)
    stat_sb_alpha = statFsb.stat(:, freq_idx_sb, :);
    max_f_sb = max(max(stat_sb_alpha, [], 3), [], 2);  % max over freq and time
    max_f_sb(~posterior_idx) = -Inf;  % exclude non-posterior
    [~, sort_idx] = sort(max_f_sb, 'descend');
    top10 = sort_idx(1:min(10, sum(posterior_idx)));
    sb_sig_channels = statFsb.label(top10)';
    fprintf('\nSternberg: no significant clusters — using top 10 posterior channels by F-stat:\n');
    for i = 1:length(top10)
        fprintf('  %2d. %s (F=%.3f)\n', i, statFsb.label{top10(i)}, max_f_sb(top10(i)));
    end
else
    fprintf('\nSternberg: %d significant posterior channel(s): %s\n', length(sb_sig_channels), strjoin(sb_sig_channels, ', '));
end

if isempty(nb_sig_channels)
    stat_nb_alpha = statFnb.stat(:, freq_idx_nb, :);
    max_f_nb = max(max(stat_nb_alpha, [], 3), [], 2);
    max_f_nb(~posterior_idx) = -Inf;
    [~, sort_idx] = sort(max_f_nb, 'descend');
    top10 = sort_idx(1:min(10, sum(posterior_idx)));
    nb_sig_channels = statFnb.label(top10)';
    fprintf('N-back:    no significant clusters — using top 10 posterior channels by F-stat:\n');
    for i = 1:length(top10)
        fprintf('  %2d. %s (F=%.3f)\n', i, statFnb.label{top10(i)}, max_f_nb(top10(i)));
    end
else
    fprintf('N-back:    %d significant posterior channel(s): %s\n', length(nb_sig_channels), strjoin(nb_sig_channels, ', '));
end
fprintf('\n');

%% Print cluster summary
fprintf('=== CLUSTER SUMMARY ===\n');
if isfield(statFnb, 'posclusters') && ~isempty(statFnb.posclusters)
    fprintf('N-back clusters: %d\n', length(statFnb.posclusters));
    for i = 1:min(5, length(statFnb.posclusters))
        fprintf('  Cluster %d: p=%.4f, stat=%.2f\n', i, statFnb.posclusters(i).prob, statFnb.posclusters(i).clusterstat);
    end
else
    fprintf('N-back: NO clusters found\n');
end
if isfield(statFsb, 'posclusters') && ~isempty(statFsb.posclusters)
    fprintf('Sternberg clusters: %d\n', length(statFsb.posclusters));
    for i = 1:min(5, length(statFsb.posclusters))
        fprintf('  Cluster %d: p=%.4f, stat=%.2f\n', i, statFsb.posclusters(i).prob, statFsb.posclusters(i).clusterstat);
    end
else
    fprintf('Sternberg: NO clusters found\n');
end
fprintf('\n');

%% Topoplots of F-values (averaged over alpha freq and time)
close all
key = [0 0.33 0.66 1];
key_rgb = [1 1 1; 1 1 0; 1 0.65 0; 1 0 0];
cmap_f_topo = interp1(key, key_rgb, linspace(0, 1, 64));

cfg_avg = [];
cfg_avg.parameter = 'stat';
cfg_avg.frequency = [8 14];
cfg_avg.avgoverfreq = 'yes';
cfg_avg.avgovertime = 'yes';
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
title('N-back Significant Electrodes', 'FontSize', 15);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_topo_f-stats_nback_FOOOF.png'));

figure('Color', 'w', 'Position', [600, 0, 600, 600]);
cfg_topo.highlightchannel = sb_sig_channels;
ft_topoplotER(cfg_topo, statFsb_avg);
clim([0 max(statFsb_avg.stat(:))]);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16; title(cb, 'F-values');
text(-0.45, 0.1, sb_sig_channels)
title('Sternberg Significant Electrodes', 'FontSize', 15);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_topo_f-stats_stern_FOOOF.png'));

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
        dataetnan_trl = ft_selectdata(cfg, dataetnan);
        if isempty(dataetnan_trl.trial) || numel(dataetnan_trl.trial) == 0
            eval(sprintf('allgazetask_sb%d{s} = [];', cond));
            eval(sprintf('allgazebase_sb%d{s} = [];', cond));
            continue;
        end
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            all_times = cellfun(@(x) [min(x), max(x)], dataetnan_trl.time, 'UniformOutput', false);
            all_times_mat = cell2mat(all_times(:));
            min_time = min(all_times_mat(:,1));
            max_time = max(all_times_mat(:,2));
            task_time_range = [max(0, min_time), min(2, max_time)];
        else
            task_time_range = [0 2];
        end
        cfg.latency = task_time_range;
        cfg = rmfield(cfg, 'trials');
        dat_task = ft_selectdata(cfg, dataetnan_trl);
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            base_time_range = [max(-.5, min_time), min(-.25, max_time)];
            if base_time_range(1) < base_time_range(2)
                cfg.latency = base_time_range;
                dat_base = ft_selectdata(cfg, dataetnan_trl);
            else
                dat_base = dat_task;
                dat_base.trial = cell(1, numel(dat_task.trial));
                dat_base.time = cell(1, numel(dat_task.trial));
                for i = 1:numel(dat_task.trial)
                    dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                    dat_base.time{i} = [];
                end
            end
        else
            dat_base = dat_task;
            dat_base.trial = cell(1, numel(dat_task.trial));
            dat_base.time = cell(1, numel(dat_task.trial));
            for i = 1:numel(dat_task.trial)
                dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                dat_base.time{i} = [];
            end
        end
        
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
    
    if ~all(ismember({'L-GAZE-X','L-GAZE-Y'}, dataet.label))
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        dataet = ft_selectdata(cfg, dataet);
    end
    nTrials = numel(dataet.trial);
    
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
    
    if size(dataet.trialinfo, 2) == 1
        trl1 = find(dataet.trialinfo == 21);
        trl2 = find(dataet.trialinfo == 22);
        trl3 = find(dataet.trialinfo == 23);
    else
        trl1 = find(dataet.trialinfo(:,1) == 21);
        trl2 = find(dataet.trialinfo(:,1) == 22);
        trl3 = find(dataet.trialinfo(:,1) == 23);
    end
    
    for cond_idx = 1:3
        if cond_idx == 1, trl = trl1; cond = 1;
        elseif cond_idx == 2, trl = trl2; cond = 2;
        else, trl = trl3; cond = 3;
        end
        
        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        dataetnan_trl = ft_selectdata(cfg, dataetnan);
        if isempty(dataetnan_trl.trial) || numel(dataetnan_trl.trial) == 0
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            all_times = cellfun(@(x) [min(x), max(x)], dataetnan_trl.time, 'UniformOutput', false);
            all_times_mat = cell2mat(all_times(:));
            min_time = min(all_times_mat(:,1));
            max_time = max(all_times_mat(:,2));
            task_time_range = [max(0, min_time), min(2, max_time)];
        else
            task_time_range = [0 2];
        end
        if task_time_range(1) >= task_time_range(2)
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        cfg.latency = task_time_range;
        cfg = rmfield(cfg, 'trials');
        dat_task = ft_selectdata(cfg, dataetnan_trl);
        if numel(dat_task.trial) == 0
            eval(sprintf('allgazetask_nb%d{s} = [];', cond));
            eval(sprintf('allgazebase_nb%d{s} = [];', cond));
            continue;
        end
        if ~isempty(dataetnan_trl.time) && numel(dataetnan_trl.time) > 0
            base_time_range = [max(-.5, min_time), min(-.25, max_time)];
            if base_time_range(1) < base_time_range(2)
                cfg.latency = base_time_range;
                dat_base = ft_selectdata(cfg, dataetnan_trl);
            else
                dat_base = dat_task;
                dat_base.trial = cell(1, numel(dat_task.trial));
                dat_base.time = cell(1, numel(dat_task.trial));
                for i = 1:numel(dat_task.trial)
                    dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                    dat_base.time{i} = [];
                end
            end
        else
            dat_base = dat_task;
            dat_base.trial = cell(1, numel(dat_task.trial));
            dat_base.time = cell(1, numel(dat_task.trial));
            for i = 1:numel(dat_task.trial)
                dat_base.trial{i} = zeros(size(dat_task.trial{i}, 1), 0);
                dat_base.time{i} = [];
            end
        end
        
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

sb2_gaze = sb2_gaze(~cellfun(@isempty, sb2_gaze));
sb4_gaze = sb4_gaze(~cellfun(@isempty, sb4_gaze));
sb6_gaze = sb6_gaze(~cellfun(@isempty, sb6_gaze));
nb1_gaze = nb1_gaze(~cellfun(@isempty, nb1_gaze));
nb2_gaze = nb2_gaze(~cellfun(@isempty, nb2_gaze));
nb3_gaze = nb3_gaze(~cellfun(@isempty, nb3_gaze));

% Compute F-test for gaze data (Sternberg)
if ~isempty(sb2_gaze) && ~isempty(sb4_gaze) && ~isempty(sb6_gaze)
    disp('Computing F-test for gaze data - Sternberg...');
    tic
    cfg_gaze = [];
    cfg_gaze.method = 'montecarlo';
    cfg_gaze.statistic = 'ft_statfun_depsamplesFunivariate';
    cfg_gaze.correctm = 'cluster';
    cfg_gaze.clusteralpha = 0.05;
    cfg_gaze.clusterstatistic = 'maxsum';
    cfg_gaze.neighbours = [];
    cfg_gaze.tail = 1;
    cfg_gaze.clustertail = 1;
    cfg_gaze.alpha = 0.05;
    cfg_gaze.numrandomization = 1000;
    
    n_gaze = length(sb2_gaze);
    cfg_gaze.design(1,:) = [ones(1,n_gaze), ones(1,n_gaze)*2, ones(1,n_gaze)*3];
    cfg_gaze.design(2,:) = [1:n_gaze, 1:n_gaze, 1:n_gaze];
    cfg_gaze.ivar = 1;
    cfg_gaze.uvar = 2;
    
    [statFgaze_sb] = ft_freqstatistics(cfg_gaze, sb2_gaze{:}, sb4_gaze{:}, sb6_gaze{:});
    statFgaze_sb.stat(statFgaze_sb.mask==0) = 0;
    toc
else
    warning('Insufficient gaze data for Sternberg F-test');
    statFgaze_sb = [];
end

% Compute F-test for gaze data (N-back)
if ~isempty(nb1_gaze) && ~isempty(nb2_gaze) && ~isempty(nb3_gaze)
    disp('Computing F-test for gaze data - N-back...');
    tic
    n_gaze = length(nb1_gaze);
    cfg_gaze.design(1,:) = [ones(1,n_gaze), ones(1,n_gaze)*2, ones(1,n_gaze)*3];
    cfg_gaze.design(2,:) = [1:n_gaze, 1:n_gaze, 1:n_gaze];
    
    [statFgaze_nb] = ft_freqstatistics(cfg_gaze, nb1_gaze{:}, nb2_gaze{:}, nb3_gaze{:});
    statFgaze_nb.stat(statFgaze_nb.mask==0) = 0;
    toc
else
    warning('Insufficient gaze data for N-back F-test');
    statFgaze_nb = [];
end

%% Create figures
disp('Creating figures...');
% Colormaps
key = [0 0.33 0.66 1];
key_rgb = [1 1 1; 1 1 0; 1 0.65 0; 1 0 0];
cmap_f = interp1(key, key_rgb, linspace(0, 1, 64));
cmap_bwr = flipud(cbrewer('div', 'RdBu', 64));

% --- Sternberg: EEG TFR (high-low difference at sig electrodes) ---
cfg_diff = [];
cfg_diff.operation = 'subtract';
cfg_diff.parameter = 'powspctrm';
ga_sb_diff = ft_math(cfg_diff, ga_sb_6vis, ga_sb_2vis);

figure('Color', 'w', 'Position', [0, 0, 1512, 982]);
cfg = []; cfg.channel = sb_sig_channels; cfg.avgoverchan = 'yes'; cfg.frequency = [4 30];
actual_time_range = [max(-.5, ga_sb_diff.time(1)), min(2, ga_sb_diff.time(end))]; cfg.latency = actual_time_range;
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
yticklabels({'30','25','20','15','10','5'});
ylim([1 500]);
set(gca, 'Fontsize', 18);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
max_abs = max(abs(pow_interp(:)));
clim([-max_abs max_abs]);
colormap(cmap_bwr);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16;
title(cb, 'Power');
title('EEG TFR: Sternberg High-Low');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_EEG_TFR_FOOOF.png'));

% --- Sternberg: ET TFR ---
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
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_ET_TFR_FOOOF.png'));

% --- Sternberg: EEG raincloud ---
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
xlabel('WM load'); ylabel('Change from baseline (log_{10} power)');
xlim([0 1]); box on; set(gca, 'FontSize', 16);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_EEG_raincloud_matlab_FOOOF.png'));

% --- Sternberg: ET raincloud ---
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
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternberg_ET_raincloud_matlab_FOOOF.png'));

% --- N-back: EEG TFR (high-low difference at sig electrodes) ---
cfg_diff = [];
cfg_diff.operation = 'subtract';
cfg_diff.parameter = 'powspctrm';
ga_nb_diff = ft_math(cfg_diff, ga_nb_3vis, ga_nb_1vis);

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
yticks(round(interp1([5 30], [1 500], [5 10 15 20 25 30]))); yticklabels({'30','25','20','15','10','5'});
ylim([1 500]);
set(gca, 'Fontsize', 18); xlabel('Time [sec]'); ylabel('Frequency [Hz]');
max_abs = max(abs(pow_interp(:)));
clim([-max_abs max_abs]);
colormap(cmap_bwr);
cb = colorbar; cb.LineWidth = 1; cb.FontSize = 16; title(cb, 'Power');
title('EEG TFR: N-back High-Low');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_EEG_TFR_FOOOF.png'));

% --- N-back: ET TFR ---
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
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_ET_TFR_FOOOF.png'));

% --- N-back: EEG raincloud ---
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
xlabel('N-back load'); ylabel('Change from baseline (log_{10} power)');
xlim([0 1]); box on; set(gca, 'FontSize', 16);
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_EEG_raincloud_matlab_FOOOF.png'));

% --- N-back: ET raincloud ---
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
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_ET_raincloud_matlab_FOOOF.png'));

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
    x_nan = x;
    y_nan = y;
    x_interp = x;
    y_interp = y;
    
    invalid = ~isfinite(x) | ~isfinite(y) | (x == 0 & y == 0);
    
    dx = diff([x(1), x]);
    dy = diff([y(1), y]);
    velocity = sqrt(dx.^2 + dy.^2);
    
    blink_velocity = velocity > threshold;
    blink_mask = invalid | blink_velocity;
    
    pad_samples = round(pad_ms * fs / 1000);
    if pad_samples > 0
        blink_idx = find(blink_mask);
        for i = 1:length(blink_idx)
            start_idx = max(1, blink_idx(i) - pad_samples);
            end_idx = min(length(blink_mask), blink_idx(i) + pad_samples);
            blink_mask(start_idx:end_idx) = true;
        end
    end
    
    x_nan(blink_mask) = NaN;
    y_nan(blink_mask) = NaN;
    
    valid_idx = ~blink_mask & isfinite(x) & isfinite(y);
    if sum(valid_idx) > 2
        x_interp = fillmissing(x_nan, 'linear');
        y_interp = fillmissing(y_nan, 'linear');
        
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
