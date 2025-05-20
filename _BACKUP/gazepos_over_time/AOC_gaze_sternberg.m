%% Visualization of gaze for AOC Sternberg
% This script visualizes the average gaze positions across conditions.

%% Setup
clear 
clc 
close all

%% Load gaze data
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_sternberg.mat');
trialinfo = trialinfo - 50; % Adjust trial info

% Identify conditions
conditions = unique(trialinfo);
numConditions = numel(conditions);

%% Process gaze data
% Preallocate cell arrays for storing averages
grandAvgX = cell(1, numConditions);
grandAvgY = cell(1, numConditions);
grandErrX = cell(1, numConditions);
grandErrY = cell(1, numConditions);

% Extract gaze data for each condition
for condIdx = 1:numConditions
    cond = conditions(condIdx);
    condTrialsIdx = find(trialinfo == cond); % Indices of trials in this condition
    
    % Collect all trials' gaze data (handling variable lengths)
    condGazeX = cellfun(@(x) padarray(x(:)', [0, 600 - numel(x)], NaN, 'post'), ...
                         gaze_x(condTrialsIdx), 'UniformOutput', false);
    condGazeY = cellfun(@(y) padarray(y(:)', [0, 600 - numel(y)], NaN, 'post'), ...
                         gaze_y(condTrialsIdx), 'UniformOutput', false);

    % Convert cell arrays to matrices
    condGazeX = vertcat(condGazeX{:});
    condGazeY = vertcat(condGazeY{:});

    % Compute grand averages and standard errors (ignoring NaNs)
    grandAvgX{condIdx} = mean(condGazeX, 1, 'omitnan');
    grandErrX{condIdx} = std(condGazeX, 0, 1, 'omitnan') / sqrt(size(condGazeX, 1));
    grandAvgY{condIdx} = mean(condGazeY, 1, 'omitnan');
    grandErrY{condIdx} = std(condGazeY, 0, 1, 'omitnan') / sqrt(size(condGazeY, 1));
end

%% Convert to Percentage Deviation
gazeXDevs = cellfun(@(x) ((x / 400) - 1) * 100, grandAvgX, 'UniformOutput', false);
gazeYDevs = cellfun(@(y) ((y / 300) - 1) * 100, grandAvgY, 'UniformOutput', false);
gazeXErrs = cellfun(@(x) (x / 400) * 100, grandErrX, 'UniformOutput', false);
gazeYErrs = cellfun(@(y) (y / 300) * 100, grandErrY, 'UniformOutput', false);

%% Plot Gaze Deviation
timeVec = linspace(-0.5, 3, 600); % Adjusted for max trial length
colors = {'b', 'k', 'r'};
transparency = 0.5;

figure('Position', [0, 0, 2000, 1200], 'Color', 'w');

% Plot X-Axis Deviation
subplot(2, 1, 1); hold on;
for condIdx = 1:numConditions
    sbar = shadedErrorBar(timeVec, gazeXDevs{condIdx}, gazeXErrs{condIdx}, ...
        {'Color', colors{condIdx}, 'markerfacecolor', colors{condIdx}});
    set(sbar.patch, 'FaceAlpha', transparency);
end
yline(0, '--', 'Color', [0.5, 0.5, 0.5]);
xline(0, '--', 'Color', [0.5, 0.5, 0.5]);
text(-0.45, -12, 'L', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 12, 'R', 'FontSize', 20, 'FontWeight', 'bold');
xlim([-.5 2]); ylim([-15, 15]);
xlabel('Time [s]', 'FontSize', 20);
ylabel('X Deviation [%]', 'FontSize', 20);
title('X-Axis Deviation from Fixation Cross', 'FontSize', 30);
legend({'WM load 2', 'WM load 4', 'WM load 6'}, 'FontSize', 20);

% Plot Y-Axis Deviation
subplot(2, 1, 2); hold on;
for condIdx = 1:numConditions
    sbar = shadedErrorBar(timeVec, gazeYDevs{condIdx}, gazeYErrs{condIdx}, ...
        {'Color', colors{condIdx}, 'markerfacecolor', colors{condIdx}});
    set(sbar.patch, 'FaceAlpha', transparency);
end
yline(0, '--', 'Color', [0.5, 0.5, 0.5]);
xline(0, '--', 'Color', [0.5, 0.5, 0.5]);
text(-0.45, -12, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 12, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
xlim([-.5 2]); ylim([-15, 15]);
xlabel('Time [s]', 'FontSize', 20);
ylabel('Y Deviation [%]', 'FontSize', 20);
title('Y-Axis Deviation from Fixation Cross', 'FontSize', 30);
legend({'WM load 2', 'WM load 4', 'WM load 6'}, 'FontSize', 20);

% Save Figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_timepoints.png');
