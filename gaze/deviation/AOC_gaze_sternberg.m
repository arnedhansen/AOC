%% Visualization of gaze for AOC Sternberg

% Visualization of:
%   Average gaze positions

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load all eye movements for all subjects
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_sternberg.mat') 
trialinfo = trialinfo-50;
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg.mat')

%% Visualize average X and Y position by conditions for TIMEPOINTS averaged over trials
% Preallocate matrices to hold the average gaze data
conditions = unique(trialinfo); % Updated to use trialinfo
numConditions = length(conditions);
numTrials = length(gaze_x); % Number of trials (400)
timePoints = length(gaze_x{1}); % Assuming all trials have the same length

% Initialize variables to store gaze data
allGazeX = cell(1, numConditions);
allGazeY = cell(1, numConditions);

% Extract and average gaze data for each condition
for condIdx = 1:numConditions
    cond = conditions(condIdx);
    % Get indices of trials corresponding to the current condition
    condTrialsIdx = find(trialinfo == cond); % Updated to use trialinfo

    % Preallocate matrix to hold gaze data for current condition
    condGazeX = zeros(length(condTrialsIdx), timePoints);
    condGazeY = zeros(length(condTrialsIdx), timePoints);

    % Extract gaze data for the current condition
    for j = 1:length(condTrialsIdx)
        condGazeX(j, :) = gaze_x{condTrialsIdx(j)};
        condGazeY(j, :) = gaze_y{condTrialsIdx(j)};
    end

    % Store the condition's gaze data for later averaging
    allGazeX{condIdx} = condGazeX;
    allGazeY{condIdx} = condGazeY;
end

% Compute the grand average and standard error across all trials for each condition
grandAverageGazeX = cell(1, numConditions);
grandErrorGazeX = cell(1, numConditions);
grandAverageGazeY = cell(1, numConditions);
grandErrorGazeY = cell(1, numConditions);

for condIdx = 1:numConditions
    % Get the gaze data for the current condition
    condGazeX = allGazeX{condIdx};
    condGazeY = allGazeY{condIdx};

    % Compute the grand average and standard error for the current condition
    grandAverageGazeX{condIdx} = mean(condGazeX, 1, 'omitnan'); % Mean across trials for each time point
    grandErrorGazeX{condIdx} = std(condGazeX, 0, 1, 'omitnan') / sqrt(size(condGazeX, 1)); % Standard error across trials for each time point
    grandAverageGazeY{condIdx} = mean(condGazeY, 1, 'omitnan'); % Mean across trials for each time point
    grandErrorGazeY{condIdx} = std(condGazeY, 0, 1, 'omitnan') / sqrt(size(condGazeY, 1)); % Standard error across trials for each time point
end

% Convert grand average to percentage deviation
gaze_x_devs = cellfun(@(x) ((x / 400) - 1) * 100, grandAverageGazeX, 'UniformOutput', false);
gaze_y_devs = cellfun(@(y) ((y / 300) - 1) * 100, grandAverageGazeY, 'UniformOutput', false);

% Convert standard error to percentage deviation based on the original scale
gaze_x_errors = cellfun(@(x) (x / 400) * 100, grandErrorGazeX, 'UniformOutput', false);
gaze_y_errors = cellfun(@(y) (y / 300) * 100, grandErrorGazeY, 'UniformOutput', false);

% Plotting
timeVec = linspace(-0.5, 3, timePoints); % Assuming time vector from -0.5 to 3 seconds
close all
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');

% Subplot for X gaze data
colors = {'b', 'k', 'r'};
subplot(2, 1, 1);
hold on;
transparency = 0.5;
for condIdx = 1:length(conditions)
    sbar(condIdx) = shadedErrorBar(timeVec, gaze_x_devs{condIdx}, gaze_x_errors{condIdx}, ...
        {'Color', colors{condIdx}, 'markerfacecolor', colors{condIdx}});
    set(sbar(condIdx).patch, 'FaceAlpha', transparency);
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -12, 'L', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 12, 'R', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ylim([-15 15])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend([sbar.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontName', 'Arial', 'FontSize', 20);
title('Deviation on X-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

% Subplot for Y gaze data
subplot(2, 1, 2);
hold on;
for condIdx = 1:length(conditions)
    sbar(condIdx) = shadedErrorBar(timeVec, gaze_y_devs{condIdx}, gaze_y_errors{condIdx}, ...
        {'Color', colors{condIdx}, 'markerfacecolor', colors{condIdx}});
    set(sbar(condIdx).patch, 'FaceAlpha', transparency);
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -12, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 12, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ylim([-15 15])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend([sbar.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontName', 'Arial', 'FontSize', 20);
title('Deviation on Y-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_timepoints.png');
