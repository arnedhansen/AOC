%% Heatmap for AOC Nback COMBINED gaze data SPLIT EARLY AND LATE
% combined data of loads 1, 2, and 3
% Baseline = [-0.5 -0.25]
% Early = 0-1000ms
% Late = 1000-2000ms

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET_nback'])
    clc
    disp(upper(['Loading ET data for Subject ' num2str(subj) '/' num2str(length(subjects)) '...']))

    %% Segment data in split latencies and conditions, filter and compute heatmaps
    % Find condition indexes for each trial
    ind2 = find(dataet.trialinfo == 21);
    ind4 = find(dataet.trialinfo == 22);
    ind6 = find(dataet.trialinfo == 23);

    % Common config
    cfg = [];
    cfg.avgovertime  = 'no';
    cfg.keeptrials   = 'yes';
    num_bins = 1000;
    smoothing_factor = 10;

    % BASELINE
    cfg.latency = [-0.5 -0.25];
    dataBaseCOMB = ft_selectdata(cfg,dataETlong);
    dataBaseCOMB = horzcat(dataBaseCOMB.trial{:});
    dataBaseCOMBAllsubs{subj} = computeGazeHeatmap(dataBaseCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataBase1 = ft_selectdata(cfg,dataETlong);
    dataBase1 = horzcat(dataBase1.trial{:});
    dataBase1Allsubs{subj} = computeGazeHeatmap(dataBase1, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataBase2 = ft_selectdata(cfg,dataETlong);
    dataBase2 = horzcat(dataBase2.trial{:});
    dataBase2Allsubs{subj} = computeGazeHeatmap(dataBase2, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataBase3 = ft_selectdata(cfg,dataETlong);
    dataBase3 = horzcat(dataBase3.trial{:});
    dataBase3Allsubs{subj} = computeGazeHeatmap(dataBase3, num_bins, smoothing_factor);

    % EARLY
    cfg.latency = [0 1];
    dataEarlyCOMB = ft_selectdata(cfg,dataETlong);
    dataEarlyCOMB = horzcat(dataEarlyCOMB.trial{:});
    dataEarlyCOMBAllsubs{subj} = computeGazeHeatmap(dataEarlyCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataEarly1 = ft_selectdata(cfg,dataETlong);
    dataEarly1 = horzcat(dataEarly1.trial{:});
    dataEarly1Allsubs{subj} = computeGazeHeatmap(dataEarly1, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataEarly2 = ft_selectdata(cfg,dataETlong);
    dataEarly2 = horzcat(dataEarly2.trial{:});
    dataEarly2Allsubs{subj} = computeGazeHeatmap(dataEarly2, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataEarly3 = ft_selectdata(cfg,dataETlong);
    dataEarly3 = horzcat(dataEarly3.trial{:});
    dataEarly3Allsubs{subj} = computeGazeHeatmap(dataEarly3, num_bins, smoothing_factor);

    % LATE
    cfg.latency = [1 2];
    dataLateCOMB = ft_selectdata(cfg,dataETlong);
    dataLateCOMB = horzcat(dataLateCOMB.trial{:});
    dataLateCOMBAllsubs{subj} = computeGazeHeatmap(dataLateCOMB, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataLate1 = ft_selectdata(cfg,dataETlong);
    dataLate1 = horzcat(dataLate1.trial{:});
    dataLate1Allsubs{subj} = computeGazeHeatmap(dataLate1, num_bins, smoothing_factor);
    cfg.trials = ind4;
    dataLate2 = ft_selectdata(cfg,dataETlong);
    dataLate2 = horzcat(dataLate2.trial{:});
    dataLate2Allsubs{subj} = computeGazeHeatmap(dataLate2, num_bins, smoothing_factor);
    cfg.trials = ind6;
    dataLate3 = ft_selectdata(cfg,dataETlong);
    dataLate3 = horzcat(dataLate3.trial{:});
    dataLate3Allsubs{subj} = computeGazeHeatmap(dataLate3, num_bins, smoothing_factor);

end

%% Baseline data
disp('BASELININIG...')
for subj = 1 : length(subjects)
    % EARLY
    dataEarlyCOMBAllsubs_bl{subj}           = dataEarlyCOMBAllsubs{subj};
    dataEarly1Allsubs_bl{subj}              = dataEarly1Allsubs{subj};
    dataEarly2Allsubs_bl{subj}              = dataEarly2Allsubs{subj};
    dataEarly3Allsubs_bl{subj}              = dataEarly3Allsubs{subj};

    dataEarlyCOMBAllsubs_bl{subj}.powspctrm = dataEarlyCOMBAllsubs{subj}.powspctrm - dataBaseCOMBAllsubs{subj}.powspctrm;
    dataEarly1Allsubs_bl{subj}.powspctrm    = dataEarly1Allsubs{subj}.powspctrm - dataBase1Allsubs{subj}.powspctrm;
    dataEarly2Allsubs_bl{subj}.powspctrm    = dataEarly2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataEarly3Allsubs_bl{subj}.powspctrm    = dataEarly3Allsubs{subj}.powspctrm - dataBase3Allsubs{subj}.powspctrm;

    % LATE
    dataLateCOMBAllsubs_bl{subj}            = dataLateCOMBAllsubs{subj};
    dataLate1Allsubs_bl{subj}               = dataLate1Allsubs{subj};
    dataLate2Allsubs_bl{subj}               = dataLate2Allsubs{subj};
    dataLate3Allsubs_bl{subj}               = dataLate3Allsubs{subj};
    dataLateCOMBAllsubs_bl{subj}.powspctrm  = dataLateCOMBAllsubs{subj}.powspctrm - dataBaseCOMBAllsubs{subj}.powspctrm;
    dataLate1Allsubs_bl{subj}.powspctrm     = dataLate1Allsubs{subj}.powspctrm - dataBase1Allsubs{subj}.powspctrm;
    dataLate2Allsubs_bl{subj}.powspctrm     = dataLate2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataLate3Allsubs_bl{subj}.powspctrm     = dataLate3Allsubs{subj}.powspctrm - dataBase3Allsubs{subj}.powspctrm;
end

%% Average across subjects
disp('AVERAGING...')
datBaseGACOMB = ft_freqgrandaverage([],dataBaseCOMBAllsubs{:});
datBaseGA1    = ft_freqgrandaverage([],dataBase1Allsubs{:});
datBaseGA2    = ft_freqgrandaverage([],dataBase2Allsubs{:});
datBaseGA3    = ft_freqgrandaverage([],dataBase3Allsubs{:});

datEarlyGACOMB = ft_freqgrandaverage([],dataEarlyCOMBAllsubs_bl{:});
datEarlyGA1    = ft_freqgrandaverage([],dataEarly1Allsubs_bl{:});
datEarlyGA2    = ft_freqgrandaverage([],dataEarly2Allsubs_bl{:});
datEarlyGA3    = ft_freqgrandaverage([],dataEarly3Allsubs_bl{:});

datLateGACOMB = ft_freqgrandaverage([],dataLateCOMBAllsubs_bl{:});
datLateGA1    = ft_freqgrandaverage([],dataLate1Allsubs_bl{:});
datLateGA2    = ft_freqgrandaverage([],dataLate2Allsubs_bl{:});
datLateGA3    = ft_freqgrandaverage([],dataLate3Allsubs_bl{:});

%% Plot GRAND AVERAGE HEATMAPS RAW
close all;
overallFontSize = 40;

% Common configuration
centerX = 800 / 2;
centerY = 600 / 2;
colMapBase = customcolormap_preset('white-red');
colMap = customcolormap_preset('red-white-blue');
maxval = max([max(datEarlyGACOMB.powspctrm(:)), max(datLateGACOMB.powspctrm(:))]);

% Plot BASELINE heatmap
figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, datBaseGACOMB);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMapBase);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('Baseline Period [-0.75 -0.25] Gaze Heatmap', 'FontSize', 30)
set(gca, "Clim", [0 max(datBaseGACOMB.powspctrm(:))])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_BASELINE.png');

% Plot EARLY heatmap
figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, datEarlyGACOMB);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('EARLY [0 1] Gaze Heatmap', 'FontSize', 30)
set(gca, "Clim", [-maxval maxval])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_EARLY.png');

% Plot LATE heatmap
figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, datLateGACOMB);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('LATE [1 2] Gaze Heatmap', 'FontSize', 30)
set(gca, "Clim", [-maxval maxval])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_LATE.png');

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
cfg.numrandomization   = 1000;
cfg.neighbours         = [];

clear design
subj   = numel(subjects);
design = [1:subj, 1:subj; ones(1,subj), 2*ones(1,subj)]; % subject IDs; condition labels
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

%% Compute significant differences
clc; disp('COMPUTING STATS...')

clc; disp('dataEarlyCOMBAllsubs...')
[statCOMBearly] = ft_freqstatistics(cfg, dataEarlyCOMBAllsubs{:},dataBaseCOMBAllsubs{:});
cohensd=((statCOMBearly.stat)./sqrt(numel(subjects)));
statCOMBearly.stat=cohensd;
statCOMBearly.stat(statCOMBearly.mask==0)=0;

clc; disp('dataEarly1Allsubs...')
[stat1early] = ft_freqstatistics(cfg, dataEarly1Allsubs{:},dataBase1Allsubs{:});
cohensd=((stat1early.stat)./sqrt(numel(subjects)));
stat1early.stat=cohensd;
stat1early.stat(stat1early.mask==0)=0;

clc; disp('dataEarly2Allsubs...')
[stat2early] = ft_freqstatistics(cfg, dataEarly2Allsubs{:},dataBase2Allsubs{:});
cohensd=((stat2early.stat)./sqrt(numel(subjects)));
stat2early.stat=cohensd;
stat2early.stat(stat2early.mask==0)=0;

clc; disp('dataEarly3Allsubs...')
[stat3early] = ft_freqstatistics(cfg, dataEarly3Allsubs{:},dataBase3Allsubs{:});
cohensd=((stat3early.stat)./sqrt(numel(subjects)));
stat3early.stat=cohensd;
stat3early.stat(stat3early.mask==0)=0;

clc; disp('dataLateCOMBAllsubs...')
[statCOMBlate] = ft_freqstatistics(cfg, dataLateCOMBAllsubs{:},dataBaseCOMBAllsubs{:});
cohensd=((statCOMBlate.stat)./sqrt(numel(subjects)));
statCOMBlate.stat=cohensd;
statCOMBlate.stat(statCOMBlate.mask==0)=0;

clc; disp('dataLate1Allsubs...')
[stat1late] = ft_freqstatistics(cfg, dataLate1Allsubs{:},dataBase1Allsubs{:});
cohensd=((stat1late.stat)./sqrt(numel(subjects)));
stat1late.stat=cohensd;
stat1late.stat(stat1late.mask==0)=0;

clc; disp('dataLate2Allsubs...')
[stat2late] = ft_freqstatistics(cfg, dataLate2Allsubs{:}, dataBase2Allsubs{:});
cohensd=((stat2late.stat)./sqrt(numel(subjects)));
stat2late.stat=cohensd;
stat2late.stat(stat2late.mask==0)=0;

clc; disp('dataLate3Allsubs...')
[stat3late] = ft_freqstatistics(cfg, dataLate3Allsubs{:},dataBase3Allsubs{:});
cohensd=((stat3late.stat)./sqrt(numel(subjects)));
stat3late.stat=cohensd;
stat3late.stat(stat3late.mask==0)=0;

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

% Early 1-back
subplot(3,2,1);
ft_singleplotTFR(cfg,stat1early);
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
title('1-back EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late 1-back
subplot(3,2,2);
ft_singleplotTFR(cfg,stat1late);
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
title('1-back LATE [1 2] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Early 2-back
subplot(3,2,3);
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
title('2-back EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late 2-back
subplot(3,2,4);
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
title('2-back LATE [1 2] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Early 3-back
subplot(3,2,5);
ft_singleplotTFR(cfg,stat3early);
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
title('3-back EARLY [0 1] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late 3-back
subplot(3,2,6);
ft_singleplotTFR(cfg,stat3late);
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
title('3-back LATE [1 2] Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_stats_OVERVIEW.png')

%% Plot stats INDIVIDUAL figures
plotGazeHeatmap(stat1early, '1-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early2');
plotGazeHeatmap(stat1late, '1-back LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late2');
plotGazeHeatmap(stat2early, '2-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early4');
plotGazeHeatmap(stat2late, '2-back LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late4');
plotGazeHeatmap(stat3early, '3-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early6');
plotGazeHeatmap(stat3late, '3-back LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late6');
plotGazeHeatmap(statCOMBearly, 'COMB EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early_COMB');
plotGazeHeatmap(statCOMBlate, 'COMB LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late_COMB');
