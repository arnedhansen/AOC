%% Heatmap for AOC Sternberg COMBINED gaze data SPLIT EARLY AND LATE
% combined data of loads 2, 4, and 6
% Baseline = [-0.5 -0.25]
% Early = 0-1000ms
% Late = 1000-2000ms

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET_sternberg'])
    clc
    disp(upper(['Loading ET data for Subject ' num2str(subj) '/' num2str(length(subjects)) '...']))

    %% Segment data in split latencies and conditions, filter and compute heatmaps
    % Find condition indexes for each trial
    ind2 = find(dataet.trialinfo == 22);
    ind4 = find(dataet.trialinfo == 24);
    ind6 = find(dataet.trialinfo == 26);

    % Common config
    cfg = [];
    cfg.avgovertime  = 'no';
    cfg.keeptrials   = 'yes';
    num_bins = 1000;
    smoothing_factor = 10;

    % BASELINE
    cfg.latency = [-.5 -0.25];
    dataBaseCOMB = ft_selectdata(cfg,dataETlong);
    dataBaseCOMB = horzcat(dataBaseCOMB.trial{:});
    dataBaseCOMBAllsubs{subj} = computeGazeHeatmap(dataBaseCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataBase2 = ft_selectdata(cfg,dataETlong);
    dataBase2 = horzcat(dataBase2.trial{:});
    dataBase2Allsubs{subj} = computeGazeHeatmap(dataBase2, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataBase4 = ft_selectdata(cfg,dataETlong);
    dataBase4 = horzcat(dataBase4.trial{:});
    dataBase4Allsubs{subj} = computeGazeHeatmap(dataBase4, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataBase6 = ft_selectdata(cfg,dataETlong);
    dataBase6 = horzcat(dataBase6.trial{:});
    dataBase6Allsubs{subj} = computeGazeHeatmap(dataBase6, num_bins, smoothing_factor);

    % EARLY
    cfg.latency = [0 1];
    dataEarlyCOMB = ft_selectdata(cfg,dataETlong);
    dataEarlyCOMB = horzcat(dataEarlyCOMB.trial{:});
    dataEarlyCOMBAllsubs{subj} = computeGazeHeatmap(dataEarlyCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataEarly2 = ft_selectdata(cfg,dataETlong);
    dataEarly2 = horzcat(dataEarly2.trial{:});
    dataEarly2Allsubs{subj} = computeGazeHeatmap(dataEarly2, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataEarly4 = ft_selectdata(cfg,dataETlong);
    dataEarly4 = horzcat(dataEarly4.trial{:});
    dataEarly4Allsubs{subj} = computeGazeHeatmap(dataEarly4, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataEarly6 = ft_selectdata(cfg,dataETlong);
    dataEarly6 = horzcat(dataEarly6.trial{:});
    dataEarly6Allsubs{subj} = computeGazeHeatmap(dataEarly6, num_bins, smoothing_factor);

    % LATE
    cfg.latency = [1 3];
    dataLateCOMB = ft_selectdata(cfg,dataETlong);
    dataLateCOMB = horzcat(dataLateCOMB.trial{:});
    dataLateCOMBAllsubs{subj} = computeGazeHeatmap(dataLateCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataLate2 = ft_selectdata(cfg,dataETlong);
    dataLate2 = horzcat(dataLate2.trial{:});
    dataLate2Allsubs{subj} = computeGazeHeatmap(dataLate2, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataLate4 = ft_selectdata(cfg,dataETlong);
    dataLate4 = horzcat(dataLate4.trial{:});
    dataLate4Allsubs{subj} = computeGazeHeatmap(dataLate4, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataLate6 = ft_selectdata(cfg,dataETlong);
    dataLate6 = horzcat(dataLate6.trial{:});
    dataLate6Allsubs{subj} = computeGazeHeatmap(dataLate6, num_bins, smoothing_factor);

end

%% Baseline data
disp('BASELININIG...')
for subj = 1 : length(subjects)
    % EARLY
    dataEarlyCOMBAllsubs_bl{subj}           = dataEarlyCOMBAllsubs{subj};
    dataEarly2Allsubs_bl{subj}              = dataEarly2Allsubs{subj};
    dataEarly4Allsubs_bl{subj}              = dataEarly4Allsubs{subj};
    dataEarly6Allsubs_bl{subj}              = dataEarly6Allsubs{subj};

    dataEarlyCOMBAllsubs_bl{subj}.powspctrm = dataEarlyCOMBAllsubs{subj}.powspctrm - dataBaseCOMBAllsubs{subj}.powspctrm;
    dataEarly2Allsubs_bl{subj}.powspctrm    = dataEarly2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataEarly4Allsubs_bl{subj}.powspctrm    = dataEarly4Allsubs{subj}.powspctrm - dataBase4Allsubs{subj}.powspctrm;
    dataEarly6Allsubs_bl{subj}.powspctrm    = dataEarly6Allsubs{subj}.powspctrm - dataBase6Allsubs{subj}.powspctrm;

    % LATE
    dataLateCOMBAllsubs_bl{subj}            = dataLateCOMBAllsubs{subj};
    dataLate2Allsubs_bl{subj}               = dataLate2Allsubs{subj};
    dataLate4Allsubs_bl{subj}               = dataLate4Allsubs{subj};
    dataLate6Allsubs_bl{subj}               = dataLate6Allsubs{subj};
    dataLateCOMBAllsubs_bl{subj}.powspctrm  = dataLateCOMBAllsubs{subj}.powspctrm - dataBaseCOMBAllsubs{subj}.powspctrm;
    dataLate2Allsubs_bl{subj}.powspctrm     = dataLate2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataLate4Allsubs_bl{subj}.powspctrm     = dataLate4Allsubs{subj}.powspctrm - dataBase4Allsubs{subj}.powspctrm;
    dataLate6Allsubs_bl{subj}.powspctrm     = dataLate6Allsubs{subj}.powspctrm - dataBase6Allsubs{subj}.powspctrm;
end

%% Average across subjects
disp('AVERAGING...')
datBaseGACOMB = ft_freqgrandaverage([],dataBaseCOMBAllsubs{:});
datBaseGA2    = ft_freqgrandaverage([],dataBase2Allsubs{:});
datBaseGA4    = ft_freqgrandaverage([],dataBase4Allsubs{:});
datBaseGA6    = ft_freqgrandaverage([],dataBase6Allsubs{:});

datEarlyGACOMB = ft_freqgrandaverage([],dataEarlyCOMBAllsubs_bl{:});
datEarlyGA2    = ft_freqgrandaverage([],dataEarly2Allsubs_bl{:});
datEarlyGA4    = ft_freqgrandaverage([],dataEarly4Allsubs_bl{:});
datEarlyGA6    = ft_freqgrandaverage([],dataEarly6Allsubs_bl{:});

datLateGACOMB = ft_freqgrandaverage([],dataLateCOMBAllsubs_bl{:});
datLateGA2    = ft_freqgrandaverage([],dataLate2Allsubs_bl{:});
datLateGA4    = ft_freqgrandaverage([],dataLate4Allsubs_bl{:});
datLateGA6    = ft_freqgrandaverage([],dataLate6Allsubs_bl{:});

%% Set up stats
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'montecarlo';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.correctm           = 'cluster';
cfg.clusteralpha       = 0.05;
cfg.clusterstatistic   = 'maxsum';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.025;
cfg.numrandomization   = 100;
cfg.neighbours         = [];

cfg.neighbours=[];
clear design
subj = numel(subjects);
% subj = 20
% subjects = [1:20]
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

%% Compute significant differences
clc; disp('COMPUTING STATS...')

clc; disp('dataEarly2Allsubs...')
[stat2early] = ft_freqstatistics(cfg, dataEarly2Allsubs{:},dataBase2Allsubs{:});
cohensd=((stat2early.stat)./sqrt(numel(subjects)));
stat2early.stat=cohensd;
stat2early.stat(stat2early.mask==0)=0;

clc; disp('dataEarly4Allsubs...')
[stat4early] = ft_freqstatistics(cfg, dataEarly4Allsubs{:},dataBase4Allsubs{:});
cohensd=((stat4early.stat)./sqrt(numel(subjects)));
stat4early.stat=cohensd;
stat4early.stat(stat4early.mask==0)=0;

clc; disp('dataEarly6Allsubs...')
[stat6early] = ft_freqstatistics(cfg, dataEarly6Allsubs{:},dataBase6Allsubs{:});
cohensd=((stat6early.stat)./sqrt(numel(subjects)));
stat6early.stat=cohensd;
stat6early.stat(stat6early.mask==0)=0;

clc; disp('dataLate2Allsubs...')
[stat2late] = ft_freqstatistics(cfg, dataLate2Allsubs{:},dataBase2Allsubs{:});
cohensd=((stat2late.stat)./sqrt(numel(subjects)));
stat2late.stat=cohensd;
stat2late.stat(stat2late.mask==0)=0;

clc; disp('dataLate4Allsubs...')
[stat4late] = ft_freqstatistics(cfg, dataLate4Allsubs{:}, dataBase4Allsubs{:});
cohensd=((stat4late.stat)./sqrt(numel(subjects)));
stat4late.stat=cohensd;
stat4late.stat(stat4late.mask==0)=0;

clc; disp('dataLate6Allsubs...')
[stat6late] = ft_freqstatistics(cfg, dataLate6Allsubs{:},dataBase6Allsubs{:});
cohensd=((stat6late.stat)./sqrt(numel(subjects)));
stat6late.stat=cohensd;
stat6late.stat(stat6late.mask==0)=0;

clc; disp('STATS DONE...')

%% Plot stats OVERVIEW
close all

% Common config
cfg               = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'outline';
cfg.zlim          = 'maxabs';
cfg.zlim          = [-.5 .5];
colMap            = customcolormap_preset('red-white-blue');
cfg.colormap      = colMap;
cfg.figure        = 'gcf';
overallFontSize   = 15;
centerX           = 800 / 2;
centerY           = 600 / 2;

% Open fig
figure;
set(gcf, 'Position', [0 0 2000 1200], 'Color', 'W')

% Early WM load 2
subplot(3,2,1);
ft_singleplotTFR(cfg,stat2early);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 2 EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late WM load 2
subplot(3,2,2);
ft_singleplotTFR(cfg,stat2late);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 2 LATE [1 3] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Early WM load 4
subplot(3,2,3);
ft_singleplotTFR(cfg,stat4early);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 4 EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late WM load 4
subplot(3,2,4);
ft_singleplotTFR(cfg,stat4late);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 4 LATE [1 3] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Early WM load 6
subplot(3,2,5);
ft_singleplotTFR(cfg,stat6early);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 6 EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late WM load 6
subplot(3,2,6);
ft_singleplotTFR(cfg,stat6late);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 6 LATE [1 3] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_stats_OVERVIEW_TZVETAN_CORRECTED_COMPUTATION.png')
