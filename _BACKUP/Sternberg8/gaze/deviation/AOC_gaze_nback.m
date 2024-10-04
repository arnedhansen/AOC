%% Visualization of eye events for AOC Nback

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
path = '/Volumes/methlab/Students/Arne/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

% Screen dimensions
screen_width = 800;
screen_height = 600;
center_x = screen_width/2;
center_y = screen_height/2;

%% LOAD GAZE DATA

%% Visualize average X position by conditions for ALL TRIALS
% Get the lengths of all arrays
array_lengths = cellfun(@length, gaze_x);
% Find the maximum length for each condition (column)
max_lengths = max(array_lengths);
% Initialize matrices filled with NaNs for each condition
padded_matrices = cell(1, 3);
for i = 1:3
    padded_matrices{i} = NaN(length(gaze_x), max_lengths(i));
end
% Pad each cell array and place it into the corresponding matrix
for i = 1:size(gaze_x, 1)
    for j = 1:size(gaze_x, 2)
        current_array = gaze_x{i, j};
        padded_matrices{j}(i, 1:length(current_array)) = current_array;
    end
end
gaze_x_matrices = cellfun(@(x) x, padded_matrices, 'UniformOutput', false);
gaze_x_avg = cellfun(@(x) mean(x, 'omitnan'), gaze_x_matrices, 'UniformOutput', false);
gaze_x_dev = cellfun(@(x) ((x/400) - 1) * 100, gaze_x_avg, 'UniformOutput', false);

% Plot the results
close all
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');
hold on;
for i = 1:3
    custom_x = 1:length(gaze_x_dev{i}); % Custom x-axis for each condition
    plot(custom_x, gaze_x_dev{i});
end
plot(custom_x, zeros(size(custom_x)), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 0.5);
hold off;

% Customize x-axis and y-axis ticks and labels
ax = gca;
ax.XTick = get(gca, 'XTick'); % Use default ticks
ax.XTickLabel = ax.XTick / 10000; % Scale down by 10^4 to display scientific notation
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Trials x 10^4', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Average Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on X-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_trials.png');

%% Visualize average X and Y position by conditions for TIMEPOINTS averaged over trials
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
    load('dataET_nback.mat');
    
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
        
        % Compute average gaze position for each time point for the current subject
        subjectAverageGazeX{condIdx} = mean(condGazeX, 1);
        subjectAverageGazeY{condIdx} = mean(condGazeY, 1);
    end
    
    % Store the subject's average gaze data
    allGazeX{subj} = subjectAverageGazeX;
    allGazeY{subj} = subjectAverageGazeY;
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

% Convert grand average to percentage deviation
gaze_x_devs = cellfun(@(x) ((x / 400) - 1) * 100, grandAverageGazeX, 'UniformOutput', false);
gaze_y_devs = cellfun(@(y) ((y / 300) - 1) * 100, grandAverageGazeY, 'UniformOutput', false);

% Convert standard error to percentage deviation based on the original scale
gaze_x_errors = cellfun(@(x) (x / 400) * 100, grandErrorGazeX, 'UniformOutput', false);
gaze_y_errors = cellfun(@(y) (y / 300) * 100, grandErrorGazeY, 'UniformOutput', false);

% Plotting
timeVec = dataet.time{1};
close all
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');

% Subplot for X gaze data
subplot(2, 1, 1);
hold on;
colors = {'b', 'k', 'r'}; 
for condIdx = 1:length(conditions)
    shadedErrorBar(timeVec, gaze_x_devs{condIdx}, gaze_x_errors{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -35, 'L', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 3, 'R', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 2])
ax = gca;
ax.XTickLabel = ax.XTickLabel; 
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20, 'Location','best');
title('Average Deviation on X-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

% Subplot for Y gaze data
subplot(2, 1, 2);
hold on;
for condIdx = 1:length(conditions)
    shadedErrorBar(timeVec, gaze_y_devs{condIdx}, gaze_y_errors{condIdx}, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, 2.5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, -30, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 2])
ax = gca;
ax.XTickLabel = ax.XTickLabel; 
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on Y-Axis from Fixation Cross', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_timepoints.png');

%% Visualize average X and Y deviation with SLIDING WINDOW by conditions for TIMEPOINTS averaged over trials
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
    load('dataET_nback.mat');
    
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
    colors = {'b', 'k', 'r'};
    for condIdx = 1:length(conditions)
        subjectGazeX = moving_average(allGazeX{subj}{condIdx}, windowPoints);
        subjectErrorX = moving_average(allErrorGazeX{subj}{condIdx}, windowPoints);
        sbar = shadedErrorBar(timeVec, subjectGazeX, subjectErrorX, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
        set(sbar.patch, 'FaceAlpha', 0.1)
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
    text(-0.45, -35, 'L', 'FontSize', 20, 'FontWeight', 'bold');
    text(-0.45, 3, 'R', 'FontSize', 20, 'FontWeight', 'bold');
    hold off;
    xlim([-0.5 2])
    ax = gca;
    ax.XTickLabel = ax.XTickLabel; 
    ax.YTickLabel = ax.YTickLabel;
    ax.FontSize = 15;
    xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
    ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
    legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20, 'Location','best');
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
    text(-0.45, 2.5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
    text(-0.45, -30, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
    hold off;
    xlim([-0.5 2])
    ax = gca;
    ax.XTickLabel = ax.XTickLabel; 
    ax.YTickLabel = ax.YTickLabel;
    ax.FontSize = 15;
    xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
    ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
    title(['Average Deviation on Y-Axis from Fixation Cross - Subject ' num2str(subj)], 'FontName', 'Arial', 'FontSize', 30);
    
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_subj' num2str(subj) '.png']);
end

% Grand average plot
figure;
set(gcf, 'Position', [0, 0, 2000, 1200], 'Color', 'w');

% Subplot for X gaze data
subplot(2, 1, 1);
hold on;
colors = {'b', 'k', 'r'};
for condIdx = 1:length(conditions)
    grandGazeX = moving_average(grandAverageGazeX{condIdx}, windowPoints);
    grandErrorX = moving_average(grandErrorGazeX{condIdx}, windowPoints);
    sbar = shadedErrorBar(timeVec, grandGazeX, grandErrorX, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
    set(sbar.patch, 'FaceAlpha', 0.1)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, -35, 'L', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, 3, 'R', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 2])
ax = gca;
ax.XTickLabel = ax.XTickLabel; 
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on X-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
legend({'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20, 'Location','best');
title('Average Deviation on X-Axis from Fixation Cross - Grand Average', 'FontName', 'Arial', 'FontSize', 30);

% Subplot for Y gaze data
subplot(2, 1, 2);
hold on;
for condIdx = 1:length(conditions)
    grandGazeY = moving_average(grandAverageGazeY{condIdx}, windowPoints);
    grandErrorY = moving_average(grandErrorGazeY{condIdx}, windowPoints);
    sbar = shadedErrorBar(timeVec, grandGazeY, grandErrorY, 'lineprops', {'Color', colors{condIdx}, 'DisplayName', ['Condition ' num2str(conditions(condIdx))]});
    set(sbar.patch, 'FaceAlpha', 0.1)
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
text(-0.45, 2.5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.45, -30, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
hold off;
xlim([-0.5 2])
ax = gca;
ax.XTickLabel = ax.XTickLabel; 
ax.YTickLabel = ax.YTickLabel;
ax.FontSize = 15;
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Deviation on Y-Axis [%]', 'FontName', 'Arial', 'FontSize', 20);
title('Average Deviation on Y-Axis from Fixation Cross - Grand Average', 'FontName', 'Arial', 'FontSize', 30);

saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_timepoints_grand_average.png');


%% Plot deviations as DENSITY PLOT
% Compute the deviations for all conditions
deviations = cell(1, length(conditions));
for condIdx = 1:length(conditions)
    % Concatenate all subjects' average gaze data for the current condition
    allSubjectGazeX = zeros(numSubjects, timePoints);
    for subj = 1:numSubjects
        allSubjectGazeX(subj, :) = allGazeX{subj}{condIdx};
    end
    
    % Calculate deviations for the current condition
    deviations{condIdx} = ((allSubjectGazeX / 400) - 1) * 100; % percentage deviation from center (x=400)
end

% Plotting the density of deviations
figure;
set(gcf, 'Position', [0, 0, 800, 600], 'Color', 'w');
hold on;

% Colors for the conditions
colors = {'b', 'k', 'r'};
conditionLabels = {'1-back', '2-back', '3-back'};
legendHandles = [];

% Plot density for each condition
for condIdx = 1:length(conditions)
    % Extract data for current condition
    conditionData = deviations{condIdx}(:); % Flatten the matrix
   
    % Compute density
    [f, xi] = ksdensity(conditionData, 'Bandwidth', 2); 
    
    % Plot density
    h = fill([xi, fliplr(xi)], [f, zeros(size(f))], colors{condIdx}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(xi, f, 'Color', colors{condIdx}, 'LineWidth', 1); 
    h.Annotation.LegendInformation.IconDisplayStyle = 'on';
    legendHandles = [legendHandles, h];
end

legend(legendHandles, conditionLabels, 'Location', 'northeast', 'FontSize', 12);
xlabel('Response Deviation from Target [%]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Density', 'FontName', 'Arial', 'FontSize', 20);
title('Density Plot of Response Deviation from Target', 'FontName', 'Arial', 'FontSize', 25);
xlim([-130 130]);
% ylim([0 0.1]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_timepoints_devplot.png');