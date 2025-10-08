%% Heatmap for AOC Sternberg COMBINED gaze data SPLIT EARLY AND LATE
% combined data of loads 2, 4, and 6
% Early = 0-1000ms
% Late = 1000-2000ms

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load data
for subj = 1:10%%%%%length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET_sternberg'])
    clc
    disp(['Loading ET data for Subject ' num2str(subj) '/' num2str(length(subjects)) '...'])

    %% Segment data in split times
    cfg = [];
    cfg.avgovertime  = 'no';
    cfg.keeptrials   = 'yes';
    cfg.latency = [-0.75 -0.25];
    dataBaseline = ft_selectdata(cfg,dataETlong);
    cfg.latency = [0 1];
    dataEarly = ft_selectdata(cfg,dataETlong);
    cfg.latency = [1 2];
    dataLate = ft_selectdata(cfg,dataETlong);

    %% Filter data for out-of-screen data points and zeros from blinks
    for conds = 1:3
        if conds == 1
            data = horzcat(dataBaseline.trial{:});
        elseif conds == 2
            data = horzcat(dataEarly.trial{:});
        elseif conds == 3
            data = horzcat(dataLate.trial{:});
        end

        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices);

        % Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(valid_data, win_size);

        x_positions = data(1, :);
        y_positions = data(2, :);

        %% Create scatterplot for data check
        % figure;
        % scatterhist(x_positions, y_positions, 'Location', 'SouthEast', 'Color', 'k', 'Marker', '.');
        %
        % % Calculate mean values
        % mean_x = mean(x_positions);
        % mean_y = mean(y_positions);
        %
        % % Add mean markers and labels
        % hold on;
        % plot(mean_x, mean_y, 'ro', 'MarkerSize', 10);
        %
        % % Set axis labels
        % xlabel('X Position');
        % ylabel('Y Position');
        % title('Scatterhist of Eye Tracker Data');
        % % Invert y-axis to match the typical eye-tracking coordinate system
        % set(gca, 'YDir','reverse')
        % xlim([0 800]);
        % ylim([0 600]);

        %% Bin and smooth data
        % Create custom grid for heatmap in pixels
        num_bins = 100;
        x_grid_pixels = linspace(0, 800, num_bins);
        y_grid_pixels = linspace(0, 600, num_bins);

        % Bin data
        smoothing_factor = 5;
        binned_data_pixels = histcounts2(x_positions, y_positions, x_grid_pixels, y_grid_pixels);

        % Apply gaussian smoothing
        smoothed_data_pixels(subj,conds, :, :) = imgaussfilt(binned_data_pixels, smoothing_factor);

        % Treat ET data as TFR for stats
        freq = [];
        freq.freq       = linspace(0, 600, 99);
        freq.time       = linspace(0, 800, 99);
        freq.label      = {'et'};
        freq.dimord     = 'chan_freq_time';
        tmp(1,:,:)      = squeeze(smoothed_data_pixels(subj,conds, :, :));
        freq.powspctrm  = tmp;

        if conds == 1
            dataBaselineAll{subj} = freq;
        elseif conds     == 2
            dataEarlyAll{subj} = freq;
        elseif conds == 3
            dataLateAll{subj}  = freq;
        end
    end
end

%% Average across subjects
subject_average = squeeze(mean(smoothed_data_pixels, 1));
datBase  = subject_average(1, :, :);
datEarly = subject_average(2, :, :);
datLate  = subject_average(3, :, :);

%% Plot HEATMAPS
close all;
overallFontSize = 40;

% Common configuration
centerX = 800 / 2;
centerY = 600 / 2;
colMap = customcolormap_preset('white-red');
maxval = max([max(datEarly(:)), max(datLate(:))]);

% Plot BASELINE heatmap
freq.powspctrm(1,:,:) = squeeze(datBase)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, freq);
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
title('Baseline Period [-0.75 -0.25] Gaze Heatmap', 'FontSize', 30)
set(gca, "Clim", [0 max(datBase(:))])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_BASELINE.png');

% Plot EARLY heatmap
freq.powspctrm(1,:,:) = squeeze(datEarly)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, freq);
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
set(gca, "Clim", [0 maxval])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_EARLY.png');

% Plot LATE heatmap
freq.powspctrm(1,:,:) = squeeze(datLate)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, freq);
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
set(gca, "Clim", [0 maxval])

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_LATE.png');

%% Calculate significant differences between baseline and time windows
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'analytic';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 10000;
cfg.neighbours         = [];

clear design
subj = length(subjects);
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

[stat] = ft_freqstatistics(cfg, dataLateAll{:}, dataBaselineAll{:});

% Handle NaNs by replacing them with 0 (or another placeholder value)
stat.stat(isnan(stat.stat)) = 0;  % Replace NaNs with 0
% Check if there are any remaining NaNs
disp(sum(isnan(stat.stat(:))));  % Count NaNs in stat.stat
stat.stat(stat.mask == 0) = 0;  % Mask out all non-significant
statsternberg = stat;
cohensd = 2 * ((statsternberg.stat) ./ sqrt(numel(design)));  % Calculate Cohen's d
statsternberg.stat = cohensd;
% Interpolate NaNs
stat.stat = fillmissing(stat.stat, 'linear', 2);  % Linear interpolation along the 2nd dimension (time/frequency)

%% Plot t-value stats
close all
freq.powspctrm(1,:,:) = squeeze(stat.stat)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

maxval = max(stat.stat(:));
Clim = [-maxval maxval];

figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, freq);
clim(Clim);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
cb = colorbar;
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('Sternberg t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_tvalues_bl_early.png');

%% MONTECARLO
close all

% Calculate significant differences l2 and l6
stat = [];
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'montecarlo';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 10000;
cfg.neighbours         = [];

clear design
subj = length(subjects);
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

[stat] = ft_freqstatistics(cfg, l6g{:}, l2g{:});
stat.stat(stat.mask==0)=0; % mask out all non significant
statsternberg=stat;
cohensd=2*((statsternberg.stat)./sqrt(numel(design)));
statsternberg.stat=cohensd;

% Plot t-value stats
freq.powspctrm(1,:,:)= squeeze(stat.stat)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
cfg = [];
cfg.figure = 'gcf';
ft_singleplotTFR(cfg, freq);
clim(Clim);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
cb = colorbar;
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('Sternberg CBPT t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_cbpt.png');
