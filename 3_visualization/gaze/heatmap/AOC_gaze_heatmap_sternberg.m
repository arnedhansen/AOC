%% Heatmap for AOC Sternberg gaze data

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/gaze');
    load([datapath, filesep 'dataET_sternberg'])

    %% Segment data per condition
    ind2 = find(dataet.trialinfo == 52);
    ind4 = find(dataet.trialinfo == 54);
    ind6 = find(dataet.trialinfo == 56);
    cfg = [];
    cfg.latency = [1 2];
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind4;
    dataetL4 = ft_selectdata(cfg,dataet);
    cfg.trials = ind6;
    dataetL6 = ft_selectdata(cfg,dataet);

    %% Filter data for out-of-screen data points and zeros from blinks
    condcounter=0;
    for condition = 1:3
        condcounter=condcounter+1;
        if condition == 1
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 2
            data=dataetL4;
            data=horzcat(dataetL4.trial{:});
        elseif condition == 3
            data=dataetL6;
            data=horzcat(dataetL6.trial{:});
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
            l2g{subj}    = freq;
        elseif condition == 2
            l4g{subj}    = freq;
        elseif condition == 3
            l6g{subj}    = freq;
        end
    end
end

%% Average across subjects
subject_average = squeeze(mean(smoothed_data_pixels, 1));
l2 = subject_average(1, :, :);
l4 = subject_average(2, :, :);
l6 = subject_average(3, :, :);

%% Calculate significant differences between low and high contrast
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

[stat] = ft_freqstatistics(cfg, l6g{:}, l2g{:});

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

%% Plot HEATMAPS (WM load 2, WM load 4 & WM load 6)
close all;

% Common configuration
centerX = 800 / 2;
centerY = 600 / 2;
mycolormap = customcolormap_preset('red-white-blue');
maxval = max([max(l2(:)), max(l4(:)), max(l6(:))]);
Clim = [-maxval maxval];

% Plot WM load 2 heatmap
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
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 25);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', 25);
title('WM load 2 Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_WM2.png');

% Plot WM load 4 heatmap
freq.powspctrm(1,:,:) = squeeze(l4)';
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
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 25);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', 25);
title('WM load 4 Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_WM4.png');

% Plot WM load 6 heatmap
freq.powspctrm(1,:,:) = squeeze(l6)';
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
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 25);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', 25);
title('WM load 6 Heatmap', 'FontSize', 30)

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_WM6.png');

%% Plot difference heatmap (WM load 6 - WM load 2)
close all
diff = l6 - l2;
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
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 25);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', 25);
title('Sternberg Difference Heatmap (WM load 6 minus WM load 2)', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_diff.png');

%% Plot differences heatmap INDIVIDUAL subjects (WM load 6 - WM load 2)
for subj = 1:length(subjects)
    close all;
    diff = l6g{subj};
    diff = l6g{subj}.powspctrm - l2g{subj}.powspctrm;

    maxval = max(diff(:));
    Clim = [-maxval maxval];

    freq.powspctrm(1,:,:) = squeeze(diff)';
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
    ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 25);
    hold on
    plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
    set(gca, 'Fontsize', 25);
    title(['Sternberg Difference Heatmap Subject ', subjects{subj}, ' (WM load 6 minus WM load 2)'], 'FontSize', 30)

    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/subjects/AOC_gaze_heatmap_sternberg_diff_subj', subjects{subj}, '.png']);
end

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
set(gca, 'Fontsize', 25);
title('Sternberg t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_tvalues.png');

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
colormap(mycolormap);
cb = colorbar;
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar
hold on
plot(centerX, centerY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
set(gca, 'Fontsize', 25);
title('Sternberg CBPT t-value Stats', 'FontSize', 30)

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/heatmap/AOC_gaze_heatmap_sternberg_cbpt.png');
