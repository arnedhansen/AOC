% Updated script for computing both raw and normalized gaze heatmaps for Sternberg task
clear all; close all;
% Subject IDs
load('/Volumes/Homestore/OCC/arne/subjects.mat');

base_dir = '/Volumes/Homestore/OCC/arne/merged';

addpath('/Volumes/Homestore/OCC/arne/funcs');
% Gaze heatmap parameters
sampling_rate = 500;
threshold = 20;
pad_ms = 150;
num_bins = 1000;
smooth_val = 5;
x_grid = linspace(0, 800, num_bins);
y_grid = linspace(0, 600, num_bins);
%%
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

%         cfg.latency = [0 1.5]; dat_early = ft_selectdata(cfg, dataSternberg);
%         cfg.latency = [1.5 3]; dat_late  = ft_selectdata(cfg, dataSternberg);
%         cfg.latency = [-1.5 0]; dat_base = ft_selectdata(cfg, dataSternberg);
        cfg.latency = [0 1]; dat_early = ft_selectdata(cfg, dataetnan);
        cfg.latency = [1 3]; dat_late  = ft_selectdata(cfg, dataetnan);
        cfg.latency = [-1.25 -.25]; dat_base = ft_selectdata(cfg, dataetnan);

        % Compute both raw and normalized heatmaps
        [freq_early, freq_early_norm] = computeGazeHeatmap(dat_early, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_late, freq_late_norm] = computeGazeHeatmap(dat_late, x_grid, y_grid, sampling_rate, smooth_val);
        [freq_base, freq_base_norm] = computeGazeHeatmap(dat_base, x_grid, y_grid, sampling_rate, smooth_val);

        eval(sprintf('allgazetaskearly%d{s} = freq_early;', cond));
        eval(sprintf('allgazetasklate%d{s} = freq_late;', cond));
        eval(sprintf('allgazebase%d{s} = freq_base;', cond));

        eval(sprintf('allgazetaskearly%d_norm{s} = freq_early_norm;', cond));
        eval(sprintf('allgazetasklate%d_norm{s} = freq_late_norm;', cond));
        eval(sprintf('allgazebase%d_norm{s} = freq_base_norm;', cond));
    end
end
%% grand average of gaze heatmaps for baseline, early interval and late interval
gagaze_b2=ft_freqgrandaverage([],allgazebase2{:});
gagaze_l2=ft_freqgrandaverage([],allgazetasklate2{:});
gagaze_e2=ft_freqgrandaverage([],allgazetaskearly2{:});

gagaze_b4=ft_freqgrandaverage([],allgazebase4{:});
gagaze_l4=ft_freqgrandaverage([],allgazetasklate4{:});
gagaze_e4=ft_freqgrandaverage([],allgazetaskearly4{:});

gagaze_b6=ft_freqgrandaverage([],allgazebase6{:});
gagaze_l6 =ft_freqgrandaverage([],allgazetasklate6{:});
gagaze_e6=ft_freqgrandaverage([],allgazetaskearly6{:});
%%
for subj=1:length(allgazebase2)
    load2early{subj}=allgazebase2{subj};
    load4early{subj}=allgazebase4{subj};
    load6early{subj}=allgazebase6{subj};
    
    load2late{subj}=allgazebase2{subj};
    load4late{subj}=allgazebase4{subj};
    load6late{subj}=allgazebase6{subj};
    
    load2late{subj}.powspctrm=allgazetasklate2{subj}.powspctrm-allgazebase2{subj}.powspctrm;
    load4late{subj}.powspctrm=allgazetasklate4{subj}.powspctrm-allgazebase4{subj}.powspctrm;
    load6late{subj}.powspctrm=allgazetasklate6{subj}.powspctrm-allgazebase6{subj}.powspctrm;
    
    load2early{subj}.powspctrm=allgazetaskearly2{subj}.powspctrm-allgazebase2{subj}.powspctrm;
    load4early{subj}.powspctrm=allgazetaskearly4{subj}.powspctrm-allgazebase4{subj}.powspctrm;
    load6early{subj}.powspctrm=allgazetaskearly6{subj}.powspctrm-allgazebase6{subj}.powspctrm;
    
end
%%
gagaze_2=ft_freqgrandaverage([],load2early{:});
gagaze_4 =ft_freqgrandaverage([],load4early{:});
gagaze_6=ft_freqgrandaverage([],load6early{:});
%%

% close all
figure;
cfg = [];
% cfg.zlim = [-10 10];
% cfg.ylim = [0 1];
cfg.figure= 'gcf';
subplot(2,2,1);
ft_singleplotTFR(cfg,gagaze_2);
title('baseline');
subplot(2,2,2);
ft_singleplotTFR(cfg,gagaze_4);
title('early');
subplot(2,2,3);
ft_singleplotTFR(cfg,gagaze_6);
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
%%
[stat2early] = ft_freqstatistics(cfg, allgazetaskearly2{:},allgazebase2{:});
cohensd=((stat2early.stat)./sqrt(numel(subjects)));
stat2early.stat=cohensd;
stat2early.stat(stat2early.mask==0)=0;% set everything not relevant to zero

[stat4early] = ft_freqstatistics(cfg, allgazetaskearly4{:},allgazebase4{:});
cohensd=((stat4early.stat)./sqrt(numel(subjects)));
stat4early.stat=cohensd;
stat4early.stat(stat4early.mask==0)=0;% set everything not relevant to zero


[stat6early] = ft_freqstatistics(cfg, allgazetaskearly6{:},allgazebase6{:});
cohensd=((stat6early.stat)./sqrt(numel(subjects)));
stat6early.stat=cohensd;
stat6early.stat(stat6early.mask==0)=0;% set everything not relevant to zero

[stat2late] = ft_freqstatistics(cfg, allgazetasklate2{:},allgazebase2{:});
cohensd=((stat2late.stat)./sqrt(numel(subjects)));
stat2late.stat=cohensd;
stat2late.stat(stat2late.mask==0)=0;

[stat4late] = ft_freqstatistics(cfg, allgazetasklate4{:}, allgazebase4{:});
cohensd=((stat4late.stat)./sqrt(numel(subjects)));
stat4late.stat=cohensd;
stat4late.stat(stat4late.mask==0)=0;

[stat6late] = ft_freqstatistics(cfg, allgazetasklate6{:},allgazebase6{:});
cohensd=((stat6late.stat)./sqrt(numel(subjects)));
stat6late.stat=cohensd;
stat6late.stat(stat6late.mask==0)=0;
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
title('0-1 sec')

subplot(3,2,2);
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
title('1-3sec')

subplot(3,2,3);
ft_singleplotTFR(cfg,stat4early);

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
ft_singleplotTFR(cfg,stat4late);

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
title('1-3sec')

subplot(3,2,5);
ft_singleplotTFR(cfg,stat6early);

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
ft_singleplotTFR(cfg,stat6late);

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
title('1-3sec')
% hgsave(gcf, 'gaze_sternberg.fig', '-v7.3');   % legacy equivalent
%%
save statssternberg stat2early stat4early stat6early stat2late stat4late stat6late -v7.3

save sterngaze allgazebase2 allgazetasklate2 allgazetaskearly2 allgazebase4 allgazetasklate4 allgazetaskearly4 allgazebase6 allgazetasklate6 allgazetaskearly6 -v7.3

save sterngaze_norm allgazebase2_norm allgazetasklate2_norm allgazetaskearly2_norm allgazebase4_norm allgazetasklate4_norm allgazetaskearly4_norm allgazebase6_norm allgazetasklate6_norm allgazetaskearly6_norm -v7.3
%% test main effect
cfg                  = [];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesFunivariate'; % use the dependent samples F-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                                              % permutation distribution.
cfg.neighbours       = [];
cfg.tail             = 1;          % 1 as the F distribution is skewed
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.025;      % alpha level of the permutation test
cfg.numrandomization = 1000;        % number of draws from the permutation distribution

n_U  = numel(subjects);
n_P  = numel(subjects);
n_N = numel(subjects);

clear design
design = zeros(2,3*n_U);
cfg.design(1,:)           = [ones(1,n_U), ones(1,n_P)*2,ones(1,n_N)*3]; % design matrix
cfg.design(2,:)           = [1:n_U,1:n_P, 1:n_N]; 
cfg.ivar                  = 1; % number or list with indices indicating the independent variable(s)
cfg.uvar                  = 2;% units-of-observation (subjects or trials
% [statFsb_early] = ft_freqstatistics(cfg, load2early{:},load4early{:},load6early{:});
[statFsb_late] = ft_freqstatistics(cfg, load2late{:},load4late{:},load6late{:});
statFsb_late.stat(statFsb_late.mask==0)=0;
%%
close all
cfg = [];
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
cfg.colormap = 'YlOrRd';
% cfg.zlim = [-.8 .8];
figure; ft_singleplotTFR(cfg,statFsb_late);
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