%% Visualization of gaze deviation for AOC Sternberg

% Visualizations:
%   Boxplots of euclidean distances for gaze deviation

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Load gaze data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    cd(datapath);
    load('gaze_dev_sternberg.mat');
    devs2(subj) = l2gdev;
    devs4(subj) = l4gdev;
    devs6(subj) = l6gdev;
end

%% Plot Euclidean deviations BOXPLOTS
dataDeviation = [devs2', devs4', devs6'];
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'k', 'r'};
hold on;

% Boxplots
boxplot(dataDeviation, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataDeviation(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)), dataDeviation(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Mean Euclidean Deviation [px]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
%set(gca, 'YLim', [min(dataDeviation(:)) * 0.85 max(dataDeviation(:)) * 1.15]);
set(gca, 'YLim', [min(dataDeviation(:)) * 0.85 85])


legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Mean Euclidean Gaze Deviation', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/AOC_dev_sternberg_boxplot_euclidean.png');

% Stats
means = mean(dataDeviation, 'omitnan');
Stds = std(dataDeviation, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);
