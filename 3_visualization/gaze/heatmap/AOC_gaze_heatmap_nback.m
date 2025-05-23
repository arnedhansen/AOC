%% Heatmap for AOC gaze data

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load data
for subj = 1:length(subjects)
    clc
    disp(['Processing Subject ', num2str(subjects{subj})])

    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET_nback'])

    %% Segment data per condition
    ind1 = find(dataet.trialinfo == 21);
    ind2 = find(dataet.trialinfo == 22);
    ind3 = find(dataet.trialinfo == 23);
    cfg = [];
    cfg.latency = [0 2];
    cfg.trials = ind1;
    dataetL1 = ft_selectdata(cfg,dataet);
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind3;
    dataetL3 = ft_selectdata(cfg,dataet);

    %% Filter data for out-of-screen data points and zeros from blinks
    condcounter=0;
    for condition = 1:3
        condcounter=condcounter+1;
        if condition == 1
            data=dataetL1;
            data=horzcat(dataetL1.trial{:});
        elseif condition == 2
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 3
            data=dataetL3;
            data=horzcat(dataetL3.trial{:});
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
        %
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
        smoothed_data_pixels(subj,condcounter, :, :) = imgaussfilt(binned_data_pixels, smoothing_factor);

        % Treat ET data as TFR for stats
        freq = [];
        freq.freq       = linspace(0, 600, 99);
        freq.time       = linspace(0, 800, 99);
        freq.label      = {'et'};
        freq.dimord     = 'chan_freq_time';
        tmp(1,:,:)      = squeeze(smoothed_data_pixels(subj,condcounter, :, :));
        freq.powspctrm  = tmp;

        if condition     == 1
            l1g{subj}    = freq;
        elseif condition == 2
            l2g{subj}    = freq;
        elseif condition == 3
            l3g{subj}    = freq;
        end
    end
end

%% Average across subjects
subject_average = squeeze(mean(smoothed_data_pixels, 1));
l1 = subject_average(1, :, :);
l2 = subject_average(2, :, :);
l3 = subject_average(3, :, :);

%% Calculate significant differences between low and high contrast
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'analytic';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.0005%%%%% 0.05;
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

[stat] = ft_freqstatistics(cfg, l3g{:}, l1g{:});

% Handle NaNs by replacing them with 0 (or another placeholder value)
stat.stat(isnan(stat.stat)) = 0;  % Replace NaNs with 0
% Check if there are any remaining NaNs
disp(sum(isnan(stat.stat(:))));  % Count NaNs in stat.stat
stat.stat(stat.mask == 0) = 0;  % Mask out all non-significant
statnback = stat;
cohensd = 2 * ((statnback.stat) ./ sqrt(numel(design)));  % Calculate Cohen's d
statnback.stat = cohensd;
% Interpolate NaNs
stat.stat = fillmissing(stat.stat, 'linear', 2);  % Linear interpolation along the 2nd dimension (time/frequency)

%% Plot HEATMAPS (1-back, 2-back & 3-back)
close all;
overallFontSize = 40;

% Common configuration
centerX = 800 / 2;
centerY = 600 / 2;
mycolormap = customcolormap_preset('red-white-blue');
maxval = max([max(l1(:)), max(l2(:)), max(l3(:))]);
Clim = [0 maxval];

% Plot 1-back heatmap
freq.powspctrm(1,:,:) = squeeze(l1)';
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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', overallFontSize);
title('1-back Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_1back.png');

% Plot 2-back heatmap
freq.powspctrm(1,:,:) = squeeze(l2)';
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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', overallFontSize);
title('2-back Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_2back.png');

% Plot 3-back heatmap
freq.powspctrm(1,:,:) = squeeze(l3)';
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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', overallFontSize);
title('3-back Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_3back.png');

%% Plot difference heatmap (3-back - 1-back)
close all
diff = l3 - l1;
freq.powspctrm(1,:,:) = squeeze(diff)';
freq.time = x_grid_pixels(1:end-1);
freq.freq = y_grid_pixels(1:end-1);
freq.label = {'et'};
freq.dimord = 'chan_freq_time';

maxval = max(diff(:));
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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('N-back Difference Heatmap (3-back minus 1-back)', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_diff.png');

%% Plot differences heatmap INDIVIDUAL subjects (3-back - 1-back)
% for subj = 1:length(subjects)
%     close all;
%     diff = l3g{subj};
%     diff = l3g{subj}.powspctrm - l1g{subj}.powspctrm;
% 
%     maxval = max(diff(:));
%     Clim = [-maxval maxval];
% 
%     freq.powspctrm(1,:,:) = squeeze(diff)';
%     freq.time = x_grid_pixels(1:end-1);
%     freq.freq = y_grid_pixels(1:end-1);
%     freq.label = {'et'};
%     freq.dimord = 'chan_freq_time';
% 
%     figure;
%     set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
%     cfg = [];
%     cfg.figure = 'gcf';
%     ft_singleplotTFR(cfg, freq);
%     clim(Clim);
%     xlim([0 800]);
%     ylim([0 600]);
%     xlabel('Screen Width [px]');
%     ylabel('Screen Height [px]');
%     colormap(mycolormap);
%     cb = colorbar;
%     ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', overallFontSize);
%     hold on
%     plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
%     set(gca, 'FontSize', overallFontSize);
%     title(['N-back Difference Heatmap Subject ', subjects{subj}, ' (3-back minus 1-back)'], 'FontSize', 30)
% 
%     saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/subjects/AOC_gaze_heatmap_nback_diff_subj', subjects{subj}, '.png']);
% end

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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('N-back t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_tvalues.png');

%% MONTECARLO
close all

% Calculate significant differences l2 and l8
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

[stat] = ft_freqstatistics(cfg, l3g{:}, l1g{:});
stat.stat(stat.mask==0)=0; % mask out all non significant
statnback=stat;
cohensd=2*((statnback.stat)./sqrt(numel(design)));
statnback.stat=cohensd;

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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title('N-back CBPT t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_nback_cbpt.png');
