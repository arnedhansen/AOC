%% Visualization of gaze deviation for AOC Nback

% Visualizations:
%   Boxplots of euclidean distances for gaze deviation

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

%% Euclidean deviation for N-back task
clc
close all

% Screen dimensions
screen_width = 800;
screen_height = 600;
center_x = screen_width / 2;
center_y = screen_height / 2;

% Compute the deviations for all conditions for X and Y gaze data
deviationsX = cell(1, length(conditions));
deviationsY = cell(1, length(conditions));
euclideanDeviations = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    % Concatenate all subjects' average gaze data for the current condition
    allSubjectGazeX = zeros(numSubjects, timePoints);
    allSubjectGazeY = zeros(numSubjects, timePoints);
    for subj = 1:numSubjects
        allSubjectGazeX(subj, :) = allGazeX{subj}{condIdx};
        allSubjectGazeY(subj, :) = allGazeY{subj}{condIdx};
    end

    % Calculate deviations for the current condition
    deviationsX{condIdx} = allSubjectGazeX;
    deviationsY{condIdx} = allSubjectGazeY;
    
    % Calculate Euclidean deviations for the current condition
    euclideanDeviations{condIdx} = sqrt((allSubjectGazeX - center_x).^2 + (allSubjectGazeY - center_y).^2);
end

% Prepare data for boxplot
dataDeviation = zeros(numSubjects, length(conditions));
for condIdx = 1:length(conditions)
    for subj = 1:numSubjects
        dataDeviation(subj, condIdx) = mean(euclideanDeviations{condIdx}(subj, :), 'omitnan');
    end
end

% Create a vector for conditions
conditionsVec = [];
for i = 1:length(conditions)
    conditionsVec = [conditionsVec; repmat(i, numSubjects, 1)];
end

% Plotting
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'g', 'r'};

% Plot for Euclidean deviations
hold on;
% Boxplots
boxplot(dataDeviation, conditionsVec, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:numSubjects
    plot(1:length(conditions), dataDeviation(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end
% Scatter plot for individual points
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, numSubjects), dataDeviation(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Mean Euclidean Deviation [px]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', {'1-back', '2-back', '3-back'}, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 3.5])
legend(scatterHandles, {'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('N-back Mean Euclidean Gaze Deviation', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_nback_boxplot_euclidean.png');

% Stats
for condIdx = 1:length(conditions)
    euclideanDeviation = [];
    for subj = 1:numSubjects
        euclideanDeviation = [euclideanDeviation; euclideanDeviations{condIdx}(subj, :)];
    end
    means(condIdx) = mean(mean(euclideanDeviation, 2, 'omitnan'));
    Stds(condIdx) = std(mean(euclideanDeviation, 2, 'omitnan'), 'omitnan');
end

% Display the results
disp('Means:');
disp(means)
disp('Stds:');
disp(Stds);
