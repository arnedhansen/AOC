%% Heatmap for AOC N-back gaze data

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
    ind1 = find(dataet.trialinfo(:, 1) == 21);
    ind2 = find(dataet.trialinfo(:, 1) == 22);
    ind3 = find(dataet.trialinfo(:, 1) == 23);

    % Common config
    cfg = [];
    cfg.avgovertime  = 'no';
    cfg.keeptrials   = 'yes';
    num_bins = 1000;
    smoothing_factor = 10;

    % BASELINE
    cfg.latency = [-0.5 -0.25];
    cfg.trials = ind1;
    dataBase2 = ft_selectdata(cfg,dataETlong);
    dataBase2 = horzcat(dataBase2.trial{:});
    dataBase1Allsubs{subj} = computeGazeHeatmap(dataBase2, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataBase4 = ft_selectdata(cfg,dataETlong);
    dataBase4 = horzcat(dataBase4.trial{:});
    dataBase2Allsubs{subj} = computeGazeHeatmap(dataBase4, num_bins, smoothing_factor);
    cfg.trials = ind3;
    dataBase6 = ft_selectdata(cfg,dataETlong);
    dataBase6 = horzcat(dataBase6.trial{:});
    dataBase3Allsubs{subj} = computeGazeHeatmap(dataBase6, num_bins, smoothing_factor);

    % Stimulus Presentation Window
    cfg.latency = [0 2];
    cfg.trials = ind1;
    dataStim2 = ft_selectdata(cfg,dataETlong);
    dataStim2 = horzcat(dataStim2.trial{:});
    dataStim1Allsubs{subj} = computeGazeHeatmap(dataStim2, num_bins, smoothing_factor);
    cfg.trials = ind2;
    dataStim4 = ft_selectdata(cfg,dataETlong);
    dataStim4 = horzcat(dataStim4.trial{:});
    dataStim2Allsubs{subj} = computeGazeHeatmap(dataStim4, num_bins, smoothing_factor);
    cfg.trials = ind3;
    dataStim6 = ft_selectdata(cfg,dataETlong);
    dataStim6 = horzcat(dataStim6.trial{:});
    dataStim3Allsubs{subj} = computeGazeHeatmap(dataStim6, num_bins, smoothing_factor);

end

%% Baseline data
disp('BASELININIG...')
for subj = 1:length(subjects)
    % LATE
    dataStim1Allsubs_bl{subj}               = dataStim1Allsubs{subj};
    dataStim2Allsubs_bl{subj}               = dataStim2Allsubs{subj};
    dataStim3Allsubs_bl{subj}               = dataStim3Allsubs{subj};
    dataStim1Allsubs_bl{subj}.powspctrm     = dataStim1Allsubs{subj}.powspctrm - dataBase1Allsubs{subj}.powspctrm;
    dataStim2Allsubs_bl{subj}.powspctrm     = dataStim2Allsubs{subj}.powspctrm - dataBase2Allsubs{subj}.powspctrm;
    dataStim3Allsubs_bl{subj}.powspctrm     = dataStim3Allsubs{subj}.powspctrm - dataBase3Allsubs{subj}.powspctrm;
end

%% Average across subjects
disp('AVERAGING...')
datBaseGA1    = ft_freqgrandaverage([],dataBase1Allsubs{:});
datBaseGA2    = ft_freqgrandaverage([],dataBase2Allsubs{:});
datBaseGA3    = ft_freqgrandaverage([],dataBase3Allsubs{:});

datFullGA1    = ft_freqgrandaverage([],dataStim1Allsubs{:});
datFullGA2    = ft_freqgrandaverage([],dataStim2Allsubs{:});
datFullGA3    = ft_freqgrandaverage([],dataStim3Allsubs{:});

datFullGA1_bl    = ft_freqgrandaverage([],dataStim1Allsubs_bl{:});
datFullGA2_bl    = ft_freqgrandaverage([],dataStim2Allsubs_bl{:});
datFullGA3_bl    = ft_freqgrandaverage([],dataStim3Allsubs_bl{:});

%% Plot GRAND AVERAGE HEATMAPS RAW
close all;
overallFontSize = 40;

% Common configuration
centerX = 800 / 2;
centerY = 600 / 2;
colMapRaw = customcolormap_preset('white-red');
maxval = max([max(datFullGA1.powspctrm(:)), max(datFullGA4.powspctrm(:)) , max(datFullGA6.powspctrm(:))]);
robustMax = prctile([datFullGA1.powspctrm(:); datFullGA4.powspctrm(:); datFullGA6.powspctrm(:)], 99.9995); 

% Plot RAW heatmap 1-back
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
cfg.zlim = [0 robustMax];
ft_singleplotTFR(cfg, datFullGA1);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMapRaw);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 2 Gaze Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_raw_1.png');

% Plot RAW heatmap 2-back
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, datFullGA2);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMapRaw);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 4 Gaze Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_raw_2.png');

% Plot RAW heatmap 3-back
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, datFullGA3);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMapRaw);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('WM load 6 Gaze Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_raw_3.png');

% Plot RAW heatmap DIFF (3-back - 1-back)
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
diff = datFullGA6;
diff.powspctrm = datFullGA6.powspctrm - datFullGA2.powspctrm;
robustLim = prctile(abs(diff.powspctrm(:)), 99.5);   % robust symmetric limit
cfg.zlim = [-robustLim robustLim];
ft_singleplotTFR(cfg, diff);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colMap = customcolormap_preset('red-white-blue');
colormap(colMap);
cb = colorbar;
ylabel(cb, '\Delta Gaze density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('Difference (3-back - 1-back) Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_raw_diff.png');

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

clc; disp('dataStim1Allsubs...')
[stat2late] = ft_freqstatistics(cfg, dataStim1Allsubs{:},dataBase1Allsubs{:});
cohensd=((stat2late.stat)./sqrt(numel(subjects)));
stat2late.stat=cohensd;
stat2late.stat(stat2late.mask==0)=0;

clc; disp('dataStim2Allsubs...')
[stat4late] = ft_freqstatistics(cfg, dataStim2Allsubs{:}, dataBase2Allsubs{:});
cohensd=((stat4late.stat)./sqrt(numel(subjects)));
stat4late.stat=cohensd;
stat4late.stat(stat4late.mask==0)=0;

clc; disp('dataStim3Allsubs...')
[stat6late] = ft_freqstatistics(cfg, dataStim3Allsubs{:},dataBase3Allsubs{:});
cohensd=((stat6late.stat)./sqrt(numel(subjects)));
stat6late.stat=cohensd;
stat6late.stat(stat6late.mask==0)=0;

clc; disp('dataDIFF...')
[statDIFF] = ft_freqstatistics(cfg, dataStim3Allsubs{:},dataStim1Allsubs{:});
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

% Late 1-back
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
title('1-back Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late 2-back
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
title('2-back Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Late 3-back
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
title('3-back Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_stats_OVERVIEW_onlystim.png')

%% DIFF STATS Heatmap
close all
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
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
title('Difference (3-back vs. 1-back) Gaze Heatmap', 'FontSize', overallFontSize*1.25)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_stats_DIFF_stim3vs1.png')

%% Plot stats INDIVIDUAL figures
% plotGazeHeatmap(stat2early, '1-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early2');
% plotGazeHeatmap(stat2late, '1-back LATE [0 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late2');
% plotGazeHeatmap(stat4early, '2-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early4');
% plotGazeHeatmap(stat4late, '2-back LATE [0 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late4');
% plotGazeHeatmap(stat6early, '3-back EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early6');
% plotGazeHeatmap(stat6late, '3-back LATE [0 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late6');
% plotGazeHeatmap(statCOMBearly, 'COMB EARLY [0 1] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_early_COMB');
% plotGazeHeatmap(statCOMBlate, 'COMB LATE [0 2] Gaze Heatmap', 'AOC_gaze_heatmap_nback_stats_late_COMB');
