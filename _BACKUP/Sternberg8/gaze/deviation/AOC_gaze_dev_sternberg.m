%% Visualization of gaze deviation for AOC Sternberg

% Visualizations:
%   Boxplots of euclidean distances for gaze deviation

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

%% Load all eye movements

for subj = 1:length(subjects)
 
end

% Plotting
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'k', 'r'};

% Plot for Euclidean deviations
hold on;
% Boxplots
boxplot(dataDeviation, conditionsVec, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:numSubjects
    plot(1:length(conditions), dataDeviation(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end
% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions)); 
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, numSubjects), dataDeviation(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Mean Euclidean Deviation [px]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 4.5])
set(gca, 'YLim', [0.5 max(dataDeviation(:))*1.15])
legend(scatterHandles, {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Mean Euclidean Gaze Deviation', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_boxplot_euclidean.png');

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
