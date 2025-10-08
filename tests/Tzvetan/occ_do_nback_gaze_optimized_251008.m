% Updated script for computing time-based gaze heatmaps with blink correction and trial exclusion
clear all; close all;
% Subject IDs
load('/Volumes/Homestore/OCC/arne/subjects.mat');

base_dir = '/Volumes/Homestore/OCC/arne/merged';

addpath('/Volumes/Homestore/OCC/arne/funcs');
% Gaze heatmap parameters
sampling_rate = 500;
thresh = 20;
pad_ms = 150;
num_bins = 1000;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);
%%
for s = 1:length(subjects)
    subj = subjects{s};
    cd(fullfile(base_dir, subj));
    load([subj '_Nback_all.mat']);

    % Select only eye tracking channels
    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
    dataet = ft_selectdata(cfg, dataNback);
    %%
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

%         cfg.latency = [0 0.5]; dat_early = ft_selectdata(cfg, dataNback);
%         cfg.latency = [0.5 2]; dat_late  = ft_selectdata(cfg, dataNback);
%         cfg.latency = [-1 0];  dat_base  = ft_selectdata(cfg, dataNback);
        cfg.latency = [0 1]; dat_early = ft_selectdata(cfg, dataetnan);
        cfg.latency = [1 2]; dat_late  = ft_selectdata(cfg, dataetnan);
        cfg.latency = [-1.25 -.25];  dat_base  = ft_selectdata(cfg, dataetnan);

        % Compute heatmaps (raw and normalized)
        [freq_early, freq_early_norm] = computeGazeHeatmap(dat_early, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_late,  freq_late_norm]  = computeGazeHeatmap(dat_late,  x_grid, y_grid, sampling_rate, smooth_val);
        [freq_base,  freq_base_norm]  = computeGazeHeatmap(dat_base,  x_grid, y_grid, sampling_rate, smooth_val);

        eval(sprintf('allgazetaskearly%d{s} = freq_early;', cond));
        eval(sprintf('allgazetasklate%d{s}  = freq_late;', cond));
        eval(sprintf('allgazebase%d{s}      = freq_base;', cond));

        eval(sprintf('allgazetaskearly%d_norm{s} = freq_early_norm;', cond));
        eval(sprintf('allgazetasklate%d_norm{s}  = freq_late_norm;', cond));
        eval(sprintf('allgazebase%d_norm{s}      = freq_base_norm;', cond));
    end
end

%% grand average of gaze heatmaps for baseline, early interval and late interval
allgazebase2nb=allgazebase2;
allgazetasklate2nb=allgazetasklate2;
allgazetaskearly2nb=allgazetaskearly2;
allgazebase2nb_norm=allgazebase2_norm;
allgazetasklate2nb_norm=allgazetasklate2_norm;
allgazetaskearly2nb_norm=allgazetaskearly2_norm;

gagaze_b1=ft_freqgrandaverage([],allgazebase1{:});
gagaze_l1=ft_freqgrandaverage([],allgazetasklate1{:});
gagaze_e1=ft_freqgrandaverage([],allgazetaskearly1{:});

gagaze_b2=ft_freqgrandaverage([],allgazebase2nb{:});
gagaze_l2=ft_freqgrandaverage([],allgazetasklate2nb{:});
gagaze_e2=ft_freqgrandaverage([],allgazetaskearly2nb{:});

gagaze_b3=ft_freqgrandaverage([],allgazebase3{:});
gagaze_l3 =ft_freqgrandaverage([],allgazetasklate3{:});
gagaze_e3=ft_freqgrandaverage([],allgazetaskearly3{:});
%%
for subj=1:length(allgazebase3)
    load1early{subj}=allgazebase1{subj};
    load2early{subj}=allgazebase1{subj};
    load3early{subj}=allgazebase3{subj};
    
    load1late{subj}=allgazebase1{subj};
    load2late{subj}=allgazebase1{subj};
    load3late{subj}=allgazebase3{subj};
    
    load1late{subj}.powspctrm=allgazebase1{subj}.powspctrm-allgazetasklate1{subj}.powspctrm;
    load2late{subj}.powspctrm=allgazebase2{subj}.powspctrm-allgazetasklate2{subj}.powspctrm;
    load3late{subj}.powspctrm=allgazebase3{subj}.powspctrm-allgazetasklate3{subj}.powspctrm;
    
    load1early{subj}.powspctrm=allgazebase1{subj}.powspctrm-allgazetaskearly1{subj}.powspctrm;
    load2early{subj}.powspctrm=allgazebase2{subj}.powspctrm-allgazetaskearly2{subj}.powspctrm;
    load3early{subj}.powspctrm=allgazebase3{subj}.powspctrm-allgazetaskearly3{subj}.powspctrm;
    
end
%%
gagaze_1=ft_freqgrandaverage([],load1early{:});
gagaze_2 =ft_freqgrandaverage([],load2early{:});
gagaze_3=ft_freqgrandaverage([],load3early{:});
%%

% close all
figure;
cfg = [];
cfg.zlim = [-.001  .001];% .5 is 50% of the time the gaze is at that position
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,gagaze_1);
title('baseline');
subplot(2,2,2);
ft_singleplotTFR(cfg,gagaze_2);
title('early');
subplot(2,2,3);
ft_singleplotTFR(cfg,gagaze_3);
title('late');

%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.latency =[300 500];
% cfg.frequency =[200 400];
% cfg.statistic        = 'ft_statfun_diff';
% cfg.clusterthreshold ='nonparametric_common';
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
%%
[stat1early] = ft_freqstatistics(cfg, allgazetaskearly1{:},allgazebase1{:});
cohensd=((stat1early.stat)./sqrt(numel(subjects)));
stat1early.stat=cohensd;
stat1early.stat(stat1early.mask==0)=0;% set everything not relevant to zero
%%
[stat2early] = ft_freqstatistics(cfg, allgazetaskearly2{:},allgazebase2{:});
cohensd=((stat2early.stat)./sqrt(numel(subjects)));
stat2early.stat=cohensd;
stat2early.stat(stat2early.mask==0)=0;% set everything not relevant to zero


[stat3early] = ft_freqstatistics(cfg, allgazetaskearly3{:},allgazebase3{:});
cohensd=((stat3early.stat)./sqrt(numel(subjects)));
stat3early.stat=cohensd;
stat3early.stat(stat3early.mask==0)=0;% set everything not relevant to zero

[stat1late] = ft_freqstatistics(cfg, allgazetasklate1{:},allgazebase1{:});
cohensd=((stat1late.stat)./sqrt(numel(subjects)));
stat1late.stat=cohensd;
stat1late.stat(stat1late.mask==0)=0;

[stat2late] = ft_freqstatistics(cfg, allgazetasklate2{:}, allgazebase2{:});
cohensd=((stat2late.stat)./sqrt(numel(subjects)));
stat2late.stat=cohensd;
stat2late.stat(stat2late.mask==0)=0;

[stat3late] = ft_freqstatistics(cfg, allgazetasklate3{:},allgazebase3{:});
cohensd=((stat3late.stat)./sqrt(numel(subjects)));
stat3late.stat=cohensd;
stat3late.stat(stat3late.mask==0)=0;
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
cfg.zlim = [-.5 .5];
% cfg.xlim =[300 500];
% cfg.ylim =[200 400];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,stat1early);

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
title('0-1 sec')

subplot(3,2,2);
ft_singleplotTFR(cfg,stat1late);

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
title('1-2sec')

subplot(3,2,3);
ft_singleplotTFR(cfg,stat2early);

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
title('0-1sec')

subplot(3,2,4);
ft_singleplotTFR(cfg,stat2late);

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
title('1-2sec')

subplot(3,2,5);
ft_singleplotTFR(cfg,stat3early);

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
title('0-1sec')

subplot(3,2,6);
ft_singleplotTFR(cfg,stat3late);

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
title('1-2sec')
% hgsave(gcf, 'gaze_sternberg.fig', '-v7.3');   % legacy equivalent
%%
save statsnback stat1early stat2early stat3early stat1late stat2late stat3late -v7.3
save nbackgaze allgazebase1 allgazetasklate1 allgazetaskearly1 allgazebase2nb allgazetasklate2nb allgazetaskearly2nb allgazebase3 allgazetasklate3 allgazetaskearly3 -v7.3
save nbackgaze_norm allgazebase1_norm allgazetasklate1_norm allgazetaskearly1_norm allgazebase2nb_norm allgazetasklate2nb_norm allgazetaskearly2nb_norm allgazebase3_norm allgazetasklate3_norm allgazetaskearly3_norm -v7.3
%%
% ----- Helper function -----
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