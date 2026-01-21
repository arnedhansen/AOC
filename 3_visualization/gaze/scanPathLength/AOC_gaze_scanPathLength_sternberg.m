%% AOC Gaze Scan Path Length — Sternberg
% Loads gaze_series_sternberg_trials (ScanPathSeries), aggregates SPL over time per subject and condition. Plots grand average and per-load with SEM. Saves figures.
%
% Key outputs:
%   Scan-path length time-series figures (grand average, per load)

%% Setup
clear
clc
close all

colors = color_def('AOC');   

basepath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
dirs = dir(basepath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

fs = 500;                                % Hz
t_series = -0.5:1/fs:2;                  % full time vector for reference
T = numel(t_series) - 1;                 % step series align to t(2:end)

scanpath_all          = nan(length(subjects), T);  % subject mean across all trials
scanpath_all_WM2      = nan(length(subjects), T);  % subject mean per load
scanpath_all_WM4      = nan(length(subjects), T);
scanpath_all_WM6      = nan(length(subjects), T);

%% ----------------------
% Aggregate per subject
% ----------------------
for subj = 1:length(subjects)
    datapath = strcat(basepath, subjects{subj}, '/gaze/gaze_series_sternberg_trials.mat');
    if ~isfile(datapath)
        warning('Missing file for subject %s, skipping.', subjects{subj})
        continue
    end

    fprintf('Loading %s ...\n', datapath)
    load(datapath, 'ScanPathSeries', 'ScanPathSeriesT', 'trialinfo')

    % Robust trialinfo parsing -> condition codes (22/24/26)
    if exist('trialinfo','var')
        if size(trialinfo,2) == 2
            conds = trialinfo(:,1);
        elseif size(trialinfo,1) == 2
            conds = trialinfo(1,:)';
        else
            error('Unexpected trialinfo shape for %s', subjects{subj});
        end
    else
        error('trialinfo not found in %s', datapath);
    end

    % Interpolate every trial to common grid
    subj_trials = nan(numel(ScanPathSeries), T);
    for trl = 1:numel(ScanPathSeries)
        s = ScanPathSeries{trl};
        t = ScanPathSeriesT{trl};
        if isempty(s) || numel(t) ~= numel(s)
            continue
        end
        try
            subj_trials(trl,:) = interp1(t, s, t_series(2:end), 'linear', NaN);
        catch
            % leave as NaN
        end
    end

    % Subject means across all trials and per load
    scanpath_all(subj,:) = nanmean(subj_trials, 1);

    idx2 = (conds == 22);
    idx4 = (conds == 24);
    idx6 = (conds == 26);

    if any(idx2), scanpath_all_WM2(subj,:) = nanmean(subj_trials(idx2,:), 1); end
    if any(idx4), scanpath_all_WM4(subj,:) = nanmean(subj_trials(idx4,:), 1); end
    if any(idx6), scanpath_all_WM6(subj,:) = nanmean(subj_trials(idx6,:), 1); end
end

%% ----------------------
% Grand averages + SEM
% ----------------------
grand_all   = nanmean(scanpath_all, 1);
sem_all     = nanstd(scanpath_all, [], 1) ./ sqrt(sum(isfinite(scanpath_all), 1));

grand_WM2   = nanmean(scanpath_all_WM2, 1);
sem_WM2     = nanstd(scanpath_all_WM2, [], 1) ./ sqrt(sum(isfinite(scanpath_all_WM2), 1));

grand_WM4   = nanmean(scanpath_all_WM4, 1);
sem_WM4     = nanstd(scanpath_all_WM4, [], 1) ./ sqrt(sum(isfinite(scanpath_all_WM4), 1));

grand_WM6   = nanmean(scanpath_all_WM6, 1);
sem_WM6     = nanstd(scanpath_all_WM6, [], 1) ./ sqrt(sum(isfinite(scanpath_all_WM6), 1));

t_plot = t_series(2:end);

%% Plot GRAND AVERAGE OVERALL
close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200])
hold on

% Light fill for the overall grand-average SEM (keep your existing vibe)
fill([t_plot fliplr(t_plot)], ...
     [grand_all+sem_all fliplr(grand_all-sem_all)], ...
     [0.85 0.9 1], 'EdgeColor', 'none')

% Per-load shaded error bars (use AOC palette from setup)
% Expecting colors to be an N×3 RGB matrix; adjust if your palette is structured differently.
col2 = colors(1,:); 
col4 = colors(2,:); 
col6 = colors(3,:);

shadedErrorBar(t_plot, grand_WM2, sem_WM2, 'lineProps', {'-','Color',col2,'LineWidth',2}, 'transparent', true);
shadedErrorBar(t_plot, grand_WM4, sem_WM4, 'lineProps', {'-','Color',col4,'LineWidth',2}, 'transparent', true);
shadedErrorBar(t_plot, grand_WM6, sem_WM6, 'lineProps', {'-','Color',col6,'LineWidth',2}, 'transparent', true);

% Black grand-average line on top
plot(t_plot, grand_all, 'k', 'LineWidth', 2.2)

xlabel('Time (s)')
ylabel('Scan Path Length (px)')
title('Sternberg Grand Average Scan Path Length Time Series')
xlim([t_series(1) t_series(end)])
box on
set(gca, 'FontSize', 25)

legend({'All-trials SEM','WM2','WM4','WM6','Grand mean'}, 'Location','northwest')

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/scanPathLength/AOC_gaze_scanPathLength_sternberg_overall.png')

%% Plot GRAND AVERAGE for each condition
close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200])
hold on
col2 = colors(1,:); 
col4 = colors(2,:); 
col6 = colors(3,:);

shadedErrorBar(t_plot, grand_WM2, sem_WM2, 'lineProps', {'-','Color',col2,'LineWidth',2}, 'transparent', true);
shadedErrorBar(t_plot, grand_WM4, sem_WM4, 'lineProps', {'-','Color',col4,'LineWidth',2}, 'transparent', true);
shadedErrorBar(t_plot, grand_WM6, sem_WM6, 'lineProps', {'-','Color',col6,'LineWidth',2}, 'transparent', true);

xlabel('Time [s]')
ylabel('Scan Path Length [px]')
title('Sternberg Grand Average Scan Path Length Time Series')
xlim([t_series(1) t_series(end)])
xline(0)
box on
set(gca, 'FontSize', 25)

legend({'WM2','WM4','WM6'}, 'Location','northeast')

% save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/scanPathLength/AOC_gaze_scanPathLength_sternberg_conds.png')
