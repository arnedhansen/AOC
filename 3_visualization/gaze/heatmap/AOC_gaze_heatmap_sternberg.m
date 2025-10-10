%% Heatmap for AOC Sternberg gaze data

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
    smoothing_factor = 5;

    % BASELINE
    cfg.latency = [-0.5 -0.25];
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

    % LATE
    cfg.latency = [1 2];
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
for subj = 1:length(subjects)
    % LATE
    dataLate2Allsubs_bl{subj}               = dataLate2Allsubs{subj};
    dataLate4Allsubs_bl{subj}               = dataLate4Allsubs{subj};
    dataLate6Allsubs_bl{subj}               = dataLate6Allsubs{subj};
    dataLate2Allsubs_bl{subj}.powspctrm     = dataLate2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataLate4Allsubs_bl{subj}.powspctrm     = dataLate4Allsubs{subj}.powspctrm - dataBase4Allsubs{subj}.powspctrm;
    dataLate6Allsubs_bl{subj}.powspctrm     = dataLate6Allsubs{subj}.powspctrm - dataBase6Allsubs{subj}.powspctrm;
end

%% Average across subjects
disp('AVERAGING...')
datBaseGA2    = ft_freqgrandaverage([],dataBase2Allsubs{:});
datBaseGA4    = ft_freqgrandaverage([],dataBase4Allsubs{:});
datBaseGA6    = ft_freqgrandaverage([],dataBase6Allsubs{:});

datLateGA2    = ft_freqgrandaverage([],dataLate2Allsubs{:});
datLateGA4    = ft_freqgrandaverage([],dataLate4Allsubs{:});
datLateGA6    = ft_freqgrandaverage([],dataLate6Allsubs{:});

datLateGA2_bl    = ft_freqgrandaverage([],dataLate2Allsubs_bl{:});
datLateGA4_bl    = ft_freqgrandaverage([],dataLate4Allsubs_bl{:});
datLateGA6_bl    = ft_freqgrandaverage([],dataLate6Allsubs_bl{:});

%% Plot GRAND AVERAGE HEATMAPS RAW
% close all;
% overallFontSize = 40;
% 
% % Common configuration
% centerX = 800 / 2;
% centerY = 600 / 2;
% colMapRaw = customcolormap_preset('white-red');
% maxval = max([max(datLateGA2.powspctrm(:)), max(datLateGA4.powspctrm(:)) , max(datLateGA6.powspctrm(:))]);
% robustMax = prctile([datLateGA2.powspctrm(:); datLateGA4.powspctrm(:); datLateGA6.powspctrm(:)], 99.9995); 
% 
% % Plot RAW heatmap WM2
% figure;
% set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
% cfg = [];
% cfg.figure = 'gcf';
% cfg.zlim = [0 robustMax];
% ft_singleplotTFR(cfg, datLateGA2);
% xlim([0 800]);
% ylim([0 600]);
% xlabel('Screen Width [px]');
% ylabel('Screen Height [px]');
% colormap(colMapRaw);
% cb = colorbar;
% ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
% hold on
% plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
% set(gca, 'FontSize', overallFontSize);
% title('WM load 2 Gaze Heatmap', 'FontSize', 30)
% 
% % Save
% saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_raw_WM2.png');
% 
% % Plot RAW heatmap WM4
% figure;
% set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
% cfg = [];
% cfg.figure = 'gcf';
% ft_singleplotTFR(cfg, datLateGA4);
% xlim([0 800]);
% ylim([0 600]);
% xlabel('Screen Width [px]');
% ylabel('Screen Height [px]');
% colormap(colMapRaw);
% cb = colorbar;
% ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
% hold on
% plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
% set(gca, 'FontSize', overallFontSize);
% title('WM load 4 Gaze Heatmap', 'FontSize', 30)
% 
% % Save
% saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_raw_WM4.png');
% 
% % Plot RAW heatmap WM6
% figure;
% set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
% cfg = [];
% cfg.figure = 'gcf';
% ft_singleplotTFR(cfg, datLateGA6);
% xlim([0 800]);
% ylim([0 600]);
% xlabel('Screen Width [px]');
% ylabel('Screen Height [px]');
% colormap(colMapRaw);
% cb = colorbar;
% ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
% hold on
% plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
% set(gca, 'FontSize', overallFontSize);
% title('WM load 6 Gaze Heatmap', 'FontSize', 30)
% 
% % Save
% saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_raw_WM6.png');
% 
% % Plot RAW heatmap DIFF (WM6-WM2)
% figure;
% set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
% cfg = [];
% cfg.figure = 'gcf';
% diff = datLateGA6;
% diff.powspctrm = datLateGA6.powspctrm - datLateGA2.powspctrm;
% robustLim = prctile(abs(diff.powspctrm(:)), 99.5);   % robust symmetric limit
% cfg.zlim = [-robustLim robustLim];
% ft_singleplotTFR(cfg, diff);
% xlim([0 800]);
% ylim([0 600]);
% xlabel('Screen Width [px]');
% ylabel('Screen Height [px]');
% colMap = customcolormap_preset('red-white-blue');
% colormap(colMap);
% cb = colorbar;
% ylabel(cb, '\Delta Gaze density [a.u.]', 'FontSize', overallFontSize);
% hold on
% plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
% set(gca, 'FontSize', overallFontSize);
% title('Difference (WM6-WM2) Heatmap', 'FontSize', 30)
% 
% % Save
% saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_raw_diff.png');

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

%% Compute significant differences
clc; disp('COMPUTING STATS...')

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

clc; disp('dataLate6Allsubs...')
[statDIFF] = ft_freqstatistics(cfg, dataLate6Allsubs{:},dataLate2Allsubs{:});
cohensd=((statDIFF.stat)./sqrt(numel(subjects)));
statDIFF.stat=cohensd;
statDIFF.stat(statDIFF.mask==0)=0;

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
set(gcf, 'Position', [0 0 800 1200], 'Color', 'W')

% Late WM load 2
subplot(3,1,1);
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
title('WM load 2 Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late WM load 4
subplot(3,1,2);
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
title('WM load 4 Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late WM load 6
subplot(3,1,3);
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
title('WM load 6 Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_stats_OVERVIEW_onlylate.png')

%% DIFF STATS Heatmap
close all
figure;
set(gcf, 'Position', [0 0 2000 1200], 'Color', 'W');
ft_singleplotTFR(cfg,statDIFF);
overallFontSize = 30;
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
title('Difference (WM6 vs. WM2) Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_stats_DIFF_late6vs2.png')

%% Plot stats INDIVIDUAL figures
% plotGazeHeatmap(stat2early, 'WM load 2 EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_early2');
% plotGazeHeatmap(stat2late, 'WM load 2 LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_late2');
% plotGazeHeatmap(stat4early, 'WM load 4 EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_early4');
% plotGazeHeatmap(stat4late, 'WM load 4 LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_late4');
% plotGazeHeatmap(stat6early, 'WM load 6 EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_early6');
% plotGazeHeatmap(stat6late, 'WM load 6 LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_late6');
% plotGazeHeatmap(statCOMBearly, 'COMB EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_early_COMB');
% plotGazeHeatmap(statCOMBlate, 'COMB LATE [1 2] Gaze Heatmap', 'AOC_gaze_heatmap_sternberg_stats_late_COMB');

