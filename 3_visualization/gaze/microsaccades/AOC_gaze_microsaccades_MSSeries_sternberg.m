%% AOC Gaze Microsaccade Series — Sternberg
% Loads gaze_series_sternberg_trials (MSSeries), aggregates MS rate over time across subjects. Plots bar + SEM. Saves figures.
%
% Key outputs:
%   Microsaccade rate time-series figure (grand average ± SEM)

%% Setup
startup
clear
[~, ~, ~ , ~] = setup('AOC');
path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';

dirs     = dir(path);
folders  = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');

%% Load data
clc
for subj = 1:length(subjects)
    disp(['Loading data for Subject ', num2str(subj)])
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'gaze_series_sternberg_trials'])
    meanMSSeriesSubs(subj, :) = nanmean(MSSeries, 1);
end

avgMSSeries = nanmean(meanMSSeriesSubs, 1);
avgSEM = nanstd(meanMSSeriesSubs, 0, 1) ./ sqrt(length(subjects));

%% Plotting Parameters 
colors        = color_def('AOC');
fontSize      = 25;
twin          = [-.5 2];
binWidth      = 0.05; % histogram bin width 50ms
edges         = twin(1):binWidth:twin(2);

%% Visualisation
close all
figure('Color','w','Position',[0 0 2000 1200])

% Per-subject mean ± SEM (counts/trial/bin)
bar(edges(1:end-1), avgMSSeries, 'BarWidth', 1, 'FaceColor', colors(2,:), 'EdgeColor', 'none')
hold on
plot(edges(1:end-1), avgMSSeries + avgSEM, '-', 'LineWidth', 1.5, 'Color', 'k')
plot(edges(1:end-1), avgMSSeries - avgSEM, '-', 'LineWidth', 1.5, 'Color', 'k')
xline(0,'k--','LineWidth',2)
xlim(twin)
xlabel('Time [s]')
ylabel('Averaged Binned Microsaccade Rate [MS/s/bin]')
title('Sternberg Microsaccade Time Series')
set(gca,'FontSize',fontSize)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/microsaccades/AOC_gaze_microsaccades_MSSeries_sternberg.png');