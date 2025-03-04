%% Visualization of behavioral data for AOC Sternberg

% Visualizations:
%   Boxplots of reaction times and accuracy

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load RT and Acc data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/behavioral');
    cd(datapath);
    load('acc_sternberg.mat');
    load('rt_sternberg.mat');
    acc2(subj) = l2acc;
    acc4(subj) = l4acc;
    acc6(subj) = l6acc;
    rt2(subj) = l2rt;
    rt4(subj) = l4rt;
    rt6(subj) = l6rt;
end

%% Plot Accuracy BOXPLOTS
dataAcc = [acc2', acc4', acc6'];
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'k', 'r'};
hold on;

% Boxplots
boxplot(dataAcc, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataAcc(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)), dataAcc(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Accuracy [%]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
set(gca, 'YLim', [min(dataAcc(:)) * 0.85 max(dataAcc(:)) * 1.15]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Accuracy', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/behavioral/AOC_acc_sternberg_boxplot.png');

% Stats
means = mean(dataAcc, 'omitnan');
Stds = std(dataAcc, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);

%% Plot Reaction Times BOXPLOTS
dataRT = [rt2', rt4', rt6'];
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
close all
figure;
set(gcf, 'Position', [0, 0, 1000, 600], 'Color', 'w');
colors = {'b', 'k', 'r'};
hold on;

% Boxplots
boxplot(dataRT, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:length(subjects)
    plot(1:length(conditions), dataRT(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, length(subjects)), dataRT(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Reaction Times [ms]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 15);
set(gca, 'LineWidth', 1.5);
set(gca, 'XLim', [0.5 length(conditions) + 0.5]);
set(gca, 'YLim', [min(dataRT(:)) * 0.85 max(dataRT(:)) * 1.15]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
title('Sternberg Reaction Times', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/behavioral/AOC_rt_sternberg_boxplot.png');

% Stats
means = mean(dataRT, 'omitnan');
Stds = std(dataRT, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);

%% Plot Reaction Times vs. Accuracy
dataAcc = [acc2', acc4', acc6'];
dataRT = [rt2', rt4', rt6'];
dataRT = dataRT * 1000;
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
close all
figure;
set(gcf, 'Position', [0, 0, 1200, 800], 'Color', 'w');
colors = {'b', 'k', 'r'};
hold on;

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(dataRT(:, condIdx), dataAcc(:, condIdx), 250, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end
lsline()
xlabel('Reaction Time [ms]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Accuracy [%]', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XLim', [500 1000]);
set(gca, 'YLim', [60 100]);

legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best', 'FontSize', 20);
title('Sternberg Reaction Times vs. Accuracy', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/behavioral/AOC_sternberg_rt_vs_acc.png');
