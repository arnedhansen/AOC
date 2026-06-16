clear all; close all;
load('/Volumes/Homestore/OCC/arne/subjects.mat');
base_dir = '/Volumes/Homestore/OCC/arne/merged';
addpath('/Volumes/Homestore/OCC/arne/funcs');
% Gaze heatmap parameters
sampling_rate = 500;
thresh = 20;
pad_ms = 150;
num_bins = 500;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);
%% loop over subj
for s = 1:length(subjects)
    subj = subjects{s};
    cd(fullfile(base_dir, subj));
    load([subj '_Nback_all.mat']);
    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
    dataet = ft_selectdata(cfg, dataNback);
  
    nTrials = numel(dataet.trial);

    % Initialize outputs
    dataetnan = dataet;
    dataetinterp = dataet;
    dataetinterp.trialinfo_blink = false(nTrials, 1);
    dataetinterp.blink_mask = cell(nTrials, 1);
    valid_trials = true(nTrials, 1);

    for i = 1:nTrials
        x = dataet.trial{i}(1,:);
        y = dataet.trial{i}(2,:);
        t = dataet.time{i};

        [x_nan, y_nan, x_interp, y_interp, blink_mask, is_valid] = ...
            removeAndInterpolateBlinks_checktrials(x, y, t, sampling_rate, thresh, pad_ms);

        if ~is_valid
            valid_trials(i) = false;
            continue;
        end

        dataetnan.trial{i}(1,:) = x_nan;
        dataetnan.trial{i}(2,:) = y_nan;
        dataetinterp.trial{i}(1,:) = x_interp;
        dataetinterp.trial{i}(2,:) = y_interp;
        dataetinterp.blink_mask{i} = blink_mask;
        dataetinterp.trialinfo_blink(i) = any(blink_mask);
    end

    % Remove invalid trials from structures
    cfg = []; cfg.trials = find(valid_trials);
    dataetnan = ft_selectdata(cfg, dataetnan);
    dataetinterp = ft_selectdata(cfg, dataetinterp);
    dataNback = ft_selectdata(cfg, dataNback);
%%
    % Recompute trial indices after trimming
    trl1 = find(dataetnan.trialinfo == 21);
    trl2 = find(dataetnan.trialinfo == 22);
    trl3 = find(dataetnan.trialinfo == 23);

    % Loop over conditions (1-back, 2-back, 3-back)
    for cond = 1:3
        eval(sprintf('trl = trl%d;', cond));

        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        cfg.latency = [0 2]; dat_early = ft_selectdata(cfg, dataetnan);
        % Compute heatmaps (raw and normalized)
        [freq_early, freq_early_norm] = computeGazeHeatmap(dat_early, x_grid, y_grid, sampling_rate, smooth_val);
 
        eval(sprintf('allgazetaskearly%d{s} = freq_early;', cond));

        eval(sprintf('allgazetaskearly%d_norm{s} = freq_early_norm;', cond));

    end
end

%% grand average of gaze heatmaps for baseline, early interval and late interval

gagaze_e1=ft_freqgrandaverage([],allgazetaskearly1{:});
gagaze_e2=ft_freqgrandaverage([],allgazetaskearly2{:});
gagaze_e3=ft_freqgrandaverage([],allgazetaskearly3{:});
%%

% close all
figure;
cfg = [];
cfg.zlim = [-.001  .001];% .5 is 50% of the time the gaze is at that position
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,gagaze_e1);
title('Load 1');
subplot(2,2,2);
ft_singleplotTFR(cfg,gagaze_e2);
title('Load 2');
subplot(2,2,3);
ft_singleplotTFR(cfg,gagaze_e3);
title('Load 3');

%% apply stats
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
subj = numel(subjects);
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
[stat1vs3] = ft_freqstatistics(cfg, allgazetaskearly3{:},allgazetaskearly1{:});
cohensd=((stat1vs3.stat)./sqrt(numel(subjects)));
stat1vs3.stat=cohensd;
stat1vs3.stat(stat1vs3.mask==0)=0;% set everything not relevant to zero
%% plot stat
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
cfg.zlim = [-.5 .5];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,stat1vs3);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('Load 3 vs. Load 1')
%% now do the same for sternberg
% Gaze heatmap parameters
sampling_rate = 500;
threshold = 20;
pad_ms = 150;
num_bins = 500;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);

for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir);

    load([subj '_Sternberg_all.mat']);

    % Trial types
    trl2 = find(dataSternberg.trialinfo == 22);
    trl4 = find(dataSternberg.trialinfo == 24);
    trl6 = find(dataSternberg.trialinfo == 26);

    % Select eye tracking channels
    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
    dataet = ft_selectdata(cfg, dataSternberg);
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

    % Remove bad trials
    cfg = []; cfg.trials = find(valid_trials);
    dataSternberg = ft_selectdata(cfg, dataSternberg);
    dataetnan = ft_selectdata(cfg, dataetnan);

    % Recompute condition indices after cleaning
    trl2 = find(dataSternberg.trialinfo == 22);
    trl4 = find(dataSternberg.trialinfo == 24);
    trl6 = find(dataSternberg.trialinfo == 26);

    % Define conditions
    cond_list = [2, 4, 6];
    for c = 1:numel(cond_list)
        cond = cond_list(c);
        eval(sprintf('trl = trl%d;', cond));

        cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
        cfg.trials = trl;
        cfg.latency = [1 2]; dat_early = ft_selectdata(cfg, dataetnan);

        % Compute both raw and normalized heatmaps
        [freq_early, freq_early_norm] = computeGazeHeatmap(dat_early, x_grid, y_grid, sampling_rate, smooth_val);
        eval(sprintf('allgazetaskearly%d{s} = freq_early;', cond));
        eval(sprintf('allgazetaskearly%d_norm{s} = freq_early_norm;', cond));

    end
end
%% grand average of gaze heatmaps for baseline, early interval and late interval

gagaze_e2=ft_freqgrandaverage([],allgazetaskearly2{:});
gagaze_e4=ft_freqgrandaverage([],allgazetaskearly4{:});
gagaze_e6=ft_freqgrandaverage([],allgazetaskearly6{:});
%% plot
% close all
figure;
cfg = [];
% cfg.zlim = [-10 10];
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,gagaze_2);
title('Load 2');
subplot(2,2,2);
ft_singleplotTFR(cfg,gagaze_4);
title('Load 4');
subplot(2,2,3);
ft_singleplotTFR(cfg,gagaze_6);
title('Load 6');
%% do stat
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
subj = numel(subjects);
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
[stat2vs4] = ft_freqstatistics(cfg, allgazetaskearly4{:},allgazetaskearly2{:});
cohensd=((stat2vs4.stat)./sqrt(numel(subjects)));
stat2vs4.stat=cohensd;
stat2vs4.stat(stat2vs4.mask==0)=0;% set everything not relevant to zero

%% plot stat
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
cfg.zlim = [-.5 .5];
cfg.figure = 'gcf';
% figure;
subplot(3,2,2);
ft_singleplotTFR(cfg,stat2vs4);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');

grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')
title('Load 4 vs. Load 2')
%% helper
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