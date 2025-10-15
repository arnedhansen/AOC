%% AOC Gaze Sternberg Scan Path Length

%% Setup
clear
clc
close all

basepath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
dirs = dir(basepath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

fs = 500;  % Hz
t_series = -0.5:1/fs:2;   % full time vector for reference
scanpath_all = nan(length(subjects), length(t_series)-1); % since diff() shortens by 1

%%
for subj = 1:length(subjects)
    datapath = strcat(basepath, subjects{subj}, '/gaze/gaze_series_sternberg_trials.mat');
    if ~isfile(datapath)
        warning('Missing file for subject %s, skipping.', subjects{subj})
        continue
    end

    load(datapath, 'ScanPathSeries', 'ScanPathSeriesT')

    % Concatenate all trials from this subject, resampling to common t-series
    subj_series = nan(numel(ScanPathSeries), length(t_series)-1);

    for trl = 1:numel(ScanPathSeries)
        s = ScanPathSeries{trl};
        t = ScanPathSeriesT{trl};

        if isempty(s) || numel(t) ~= numel(s)
            continue
        end

        % Interpolate to common grid
        try
            s_interp = interp1(t, s, t_series(2:end), 'linear', NaN);
            subj_series(trl,:) = s_interp;
        catch
            continue
        end
    end

    % Mean over trials (within subject)
    subj_mean = nanmean(subj_series,1);

    scanpath_all(subj,:) = subj_mean;
end

% Grand average across subjects
scanpath_grandmean = nanmean(scanpath_all,1);
scanpath_sem = nanstd(scanpath_all,[],1) ./ sqrt(sum(~isnan(scanpath_all),1));

%% Plot
figure
set(gcf, 'Color', 'w', 'Position', [100 100 1000 500])
hold on

fill([t_series(2:end) fliplr(t_series(2:end))], ...
     [scanpath_grandmean+scanpath_sem fliplr(scanpath_grandmean-scanpath_sem)], ...
     [0.85 0.9 1], 'EdgeColor', 'none')

plot(t_series(2:end), scanpath_grandmean, 'k', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Scan Path Length (px)')
title('Grand Average Scan Path Length Time Series')
xlim([t_series(1) t_series(end)])
box on
set(gca, 'FontSize', 14)