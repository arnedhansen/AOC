%% Visualization of gaze for AOC Sternberg

% Visualization of:
%   Average gaze positions
%   etc...

%% Setup
clear
clc
close all
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
close
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load all eye movements for all subjects
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_sternberg.mat')
load('/Volumes/methlab/Students/Arne/AOC/data/features/gaze_matrix_sternberg.mat')

%% Visualize average X and Y position by conditions for TIMEPOINTS averaged over trials

% Preallocate matrices to hold the average gaze data
conditions = unique([gaze_data_sternberg.Condition]);
numConditions = length(conditions);
numSubjects = size(gaze_x, 1); % Number of subjects (100)
numTrials = size(gaze_x, 2); % Number of trials (400)
timePoints = length(gaze_x{1}); % Assuming all trials have the same length

allGazeX = cell(1, numSubjects);
allGazeY = cell(1, numSubjects);

for subj = 1:numSubjects
    subjectAverageGazeX = cell(1, numConditions);
    subjectAverageGazeY = cell(1, numConditions);

    for condIdx = 1:numConditions
        cond = conditions(condIdx);
        % Get indices of trials corresponding to the current condition
        condTrialsIdx = find([gaze_data_sternberg.Condition] == cond);

        % Preallocate matrix to hold gaze data for current condition
        condGazeX = zeros(length(condTrialsIdx), timePoints);
        condGazeY = zeros(length(condTrialsIdx), timePoints);

        % Extract gaze data for the current condition
        for j = 1:length(condTrialsIdx)
            condGazeX(j, :) = gaze_x{subj, condTrialsIdx(j)};
            condGazeY(j, :) = gaze_y{subj, condTrialsIdx(j)};
        end

        % Compute average gaze position for each time point for the current condition
        subjectAverageGazeX{condIdx} = mean(condGazeX, 1);
        subjectAverageGazeY{condIdx} = mean(condGazeY, 1);
    end

    % Store the subject's average gaze data
    allGazeX{subj} = subjectAverageGazeX;
    allGazeY{subj} = subjectAverageGazeY;
    fprintf('Subject %.2d loaded \n', subj)
end

% Compute the grand average and standard error across all subjects for each condition
grandAverageGazeX = cell(1, numConditions);
grandErrorGazeX = cell(1, numConditions);
grandAverageGazeY = cell(1, numConditions);
grandErrorGazeY = cell(1, numConditions);

for condIdx = 1:numConditions
    % Concatenate all subjects' average gaze data for the current condition
    allSubjectGazeX = zeros(numSubjects, timePoints);
    allSubjectGazeY = zeros(numSubjects, timePoints);

    for subj = 1:numSubjects
        allSubjectGazeX(subj, :) = allGazeX{subj}{condIdx};
        allSubjectGazeY(subj, :) = allGazeY{subj}{condIdx};
    end

    % Compute the grand average and standard error for the current condition
    grandAverageGazeX{condIdx} = mean(allSubjectGazeX, 1);
    grandErrorGazeX{condIdx} = std(allSubjectGazeX, 0, 1) / sqrt(numSubjects);
    grandAverageGazeY{condIdx} = mean(allSubjectGazeY, 1);
    grandErrorGazeY{condIdx} = std(allSubjectGazeY, 0, 1) / sqrt(numSubjects);
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
colors = {'b', 'g', 'k', 'r'};
subplot(2, 1, 1);
hold on;
for condIdx = 1:length(conditions)
    shadedErrorBar(timeVec, gaze_x_devs{condIdx}, gaze_x_errors{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -25, 'L', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 5, 'R', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Average Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on X-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

% Subplot for Y gaze data
subplot(2, 1, 2);
hold on;
for condIdx = 1:length(conditions)
    shadedErrorBar(timeVec, gaze_y_devs{condIdx}, gaze_y_errors{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -20, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Average Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on Y-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_timepoints.png');

%% Visualize average X and Y deviation with SLIDING WINDOW by conditions for TIMEPOINTS averaged over trials
clear
% Define the path
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

% Preallocate matrices to hold the average gaze data across all subjects
conditions = [];
numSubjects = length(subjects);
allGazeX = cell(1, numSubjects);
allGazeY = cell(1, numSubjects);
timePoints = [];

% Load data for each subject and accumulate the gaze data
for subj = 1:numSubjects
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load('dataETstern.mat');

    numTrials = length(dataet.trial);
    lGazeX = cell(1, numTrials);
    lGazeY = cell(1, numTrials);
    for i = 1:numTrials
        lGazeX{i} = dataet.trial{1, i}(1, :);
        lGazeY{i} = dataet.trial{1, i}(2, :);
    end

    % Store unique conditions and time points
    if isempty(conditions)
        conditions = unique(dataet.trialinfo);
        timePoints = length(dataet.time{1});
    end

    % Preallocate matrices to hold the average gaze data for the current subject
    subjectAverageGazeX = cell(1, length(conditions));
    subjectAverageGazeY = cell(1, length(conditions));
    subjectErrorGazeX = cell(1, length(conditions));
    subjectErrorGazeY = cell(1, length(conditions));

    for condIdx = 1:length(conditions)
        cond = conditions(condIdx);
        % Get indices of trials corresponding to the current condition
        condTrialsIdx = find(dataet.trialinfo == cond);

        % Preallocate matrix to hold gaze data for current condition
        condGazeX = zeros(length(condTrialsIdx), timePoints);
        condGazeY = zeros(length(condTrialsIdx), timePoints);

        % Extract gaze data for the current condition
        for j = 1:length(condTrialsIdx)
            condGazeX(j, :) = lGazeX{condTrialsIdx(j)};
            condGazeY(j, :) = lGazeY{condTrialsIdx(j)};
        end

        % Compute average and error gaze position for each time point for the current subject
        subjectAverageGazeX{condIdx} = mean(condGazeX, 1);
        subjectErrorGazeX{condIdx} = std(condGazeX, 0, 1) / sqrt(length(condTrialsIdx));
        subjectAverageGazeY{condIdx} = mean(condGazeY, 1);
        subjectErrorGazeY{condIdx} = std(condGazeY, 0, 1) / sqrt(length(condTrialsIdx));
    end

    % Store the subject's average gaze data
    allGazeX{subj} = subjectAverageGazeX;
    allGazeY{subj} = subjectAverageGazeY;
    allErrorGazeX{subj} = subjectErrorGazeX;
    allErrorGazeY{subj} = subjectErrorGazeY;
end

% Compute the grand average and standard error across all subjects for each condition
grandAverageGazeX = cell(1, length(conditions));
grandErrorGazeX = cell(1, length(conditions));
grandAverageGazeY = cell(1, length(conditions));
grandErrorGazeY = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    % Concatenate all subjects' average gaze data for the current condition
    allSubjectGazeX = zeros(numSubjects, timePoints);
    allSubjectGazeY = zeros(numSubjects, timePoints);

    for subj = 1:numSubjects
        allSubjectGazeX(subj, :) = allGazeX{subj}{condIdx};
        allSubjectGazeY(subj, :) = allGazeY{subj}{condIdx};
    end

    % Compute the grand average and standard error for the current condition
    grandAverageGazeX{condIdx} = mean(allSubjectGazeX, 1);
    grandErrorGazeX{condIdx} = std(allSubjectGazeX, 0, 1) / sqrt(numSubjects);
    grandAverageGazeY{condIdx} = mean(allSubjectGazeY, 1);
    grandErrorGazeY{condIdx} = std(allSubjectGazeY, 0, 1) / sqrt(numSubjects);
end

% Convert individual subject values to percentage deviation
for subj = 1:numSubjects
    for condIdx = 1:length(conditions)
        allGazeX{subj}{condIdx} = ((allGazeX{subj}{condIdx} / 400) - 1) * 100;
        allGazeY{subj}{condIdx} = ((allGazeY{subj}{condIdx} / 300) - 1) * 100;
        allErrorGazeX{subj}{condIdx} = (allErrorGazeX{subj}{condIdx} / 400) * 100;
        allErrorGazeY{subj}{condIdx} = (allErrorGazeY{subj}{condIdx} / 300) * 100;
    end
end

% Convert grand average and std to percentage deviation
grandAverageGazeX = cellfun(@(x) ((x / 400) - 1) * 100, grandAverageGazeX, 'UniformOutput', false);
grandAverageGazeY = cellfun(@(y) ((y / 300) - 1) * 100, grandAverageGazeY, 'UniformOutput', false);
grandErrorGazeX = cellfun(@(x) (x / 400) * 100, grandErrorGazeX, 'UniformOutput', false);
grandErrorGazeY = cellfun(@(y) (y / 300) * 100, grandErrorGazeY, 'UniformOutput', false);

save /Volumes/methlab/Students/Arne/AOC/data/features/gaze/gaze_deviation allGazeX allErrorGazeX grandAverageGazeX grandErrorGazeX allGazeY allErrorGazeY grandAverageGazeY grandErrorGazeY

% Define window length
windowLength = 0.1; % 100 ms
windowPoints = round(windowLength / (dataet.time{1}(2) - dataet.time{1}(1))); % Number of points in 100 ms window

% Helper function to compute moving average
moving_average = @(data, win_length) movmean(data, win_length, 2);

% Plotting
timeVec = dataet.time{1};
close all

% Individual participant plots
for subj = 1:numSubjects
    figure;
    set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');

    % Subplot for X gaze data
    subplot(2, 1, 1);
    hold on;
    colors = {'b', 'g', 'k', 'r'};
    for condIdx = 1:length(conditions)
        subjectGazeX = moving_average(allGazeX{subj}{condIdx}, windowPoints);
        subjectErrorX = moving_average(allErrorGazeX{subj}{condIdx}, windowPoints);
        sbar = shadedErrorBar(timeVec, subjectGazeX, subjectErrorX, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
        set(sbar.patch, 'FaceAlpha', 0.1)
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    text(-0.45, max([allGazeX{subj}{:}])*0.75, 'RIGHT', 'FontSize', 20, 'FontWeight', 'bold');
    text(-0.45, -max([allGazeX{subj}{:}])*0.75, 'LEFT', 'FontSize', 20, 'FontWeight', 'bold');
    hold off;
    xlim([-0.5 3])
    ylim([min([allGazeX{subj}{:}])*1.1 max([allGazeX{subj}{:}])*1.1])
    ax = gca;
    ax.XTickLabel = ax.XTickLabel;
    ax.YTickLabel = ax.YTickLabel;
    ax.FontSize = 15;
    xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
    ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
    legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 20, 'Location','best');
    title(['Average Deviation on X-Axis from Fixation Cross - Subject ' num2str(subj)], 'FontName', 'Arial', 'FontSize', 30);

    % Subplot for Y gaze data
    subplot(2, 1, 2);
    hold on;
    for condIdx = 1:length(conditions)
        subjectGazeY = moving_average(allGazeY{subj}{condIdx}, windowPoints);
        subjectErrorY = moving_average(allErrorGazeY{subj}{condIdx}, windowPoints);
        sbar = shadedErrorBar(timeVec, subjectGazeY, subjectErrorY, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
        set(sbar.patch, 'FaceAlpha', 0.1)
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    text(-0.45, max([allGazeY{subj}{:}])*0.75, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
    text(-0.45, -max([allGazeY{subj}{:}])*0.75, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
    hold off;
    xlim([-0.5 3])
    ylim([min([allGazeY{subj}{:}])*1.1 max([allGazeY{subj}{:}])*1.1])
    ax = gca;
    ax.XTickLabel = ax.XTickLabel;
    ax.YTickLabel = ax.YTickLabel;
    ax.FontSize = 15;
    xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
    ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
    title(['Average Deviation on Y-Axis from Fixation Cross - Subject ' num2str(subj)], 'FontName', 'Arial', 'FontSize', 30);

    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_subj' num2str(subj) '.png']);
end

%% Grand average plot
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');

% Subplot for X gaze data
subplot(2, 1, 1);
hold on;
colors = {'b', 'g', 'k', 'r'};
for condIdx = 1:length(conditions)
    grandGazeX{condIdx} = moving_average(grandAverageGazeX{condIdx}, windowPoints);
    grandErrorX{condIdx} = moving_average(grandErrorGazeX{condIdx}, windowPoints);
    sbar = shadedErrorBar(timeVec, grandGazeX{condIdx}, grandErrorX{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
    set(sbar.patch, 'FaceAlpha', 0.1)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, max([grandGazeX{:}])*1.25, 'R', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, -max([grandGazeX{:}])*1.25, 'L', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ylim([min([grandGazeX{:}])-max([grandErrorX{:}])*1.2 max([grandGazeX{:}])+max([grandErrorX{:}])*1.2])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 20, 'Location','best');
title('Average Deviation on X-Axis from Fixation Cross - Grand Average', 'FontName', 'Arial', 'FontSize', 30);

% Subplot for Y gaze data
subplot(2, 1, 2);
hold on;
for condIdx = 1:length(conditions)
    grandGazeY{condIdx} = moving_average(grandAverageGazeY{condIdx}, windowPoints);
    grandErrorY{condIdx} = moving_average(grandErrorGazeY{condIdx}, windowPoints);
    sbar = shadedErrorBar(timeVec, grandGazeY{condIdx}, grandErrorY{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
    set(sbar.patch, 'FaceAlpha', 0.1)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, max([grandGazeY{:}])*1.25,  'UP', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, -max([grandGazeY{:}])*1.25, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 3])
ylim([min([grandGazeX{:}])-max([grandErrorX{:}])*1.2 max([grandGazeX{:}])+max([grandErrorX{:}])*1.2])
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on Y-Axis from Fixation Cross - Grand Average', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_timepoints_grand_average.png');
