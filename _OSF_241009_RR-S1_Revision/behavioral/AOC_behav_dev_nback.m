%% Visualization of behavioral data for AOC Nback

% Visualizations:
%   Boxplots of euclidean distances for gaze deviation

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
    load('acc_nback.mat');
    load('rt_nback.mat');
    acc1(subj) = l1acc;
    acc2(subj) = l2acc;
    acc3(subj) = l3acc;
    rt1(subj) = l1rt;
    rt2(subj) = l2rt;
    rt3(subj) = l3rt;
end

%% Plot Accuracy BOXPLOTS
dataAcc = [acc1', acc2', acc3'];
conditions = {'1-back', '2-back', '3-back'};
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
title('N-back Accuracy', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/behavioral/AOC_acc_nback_boxplot.png');

% Stats
means = mean(dataAcc, 'omitnan');
Stds = std(dataAcc, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);

%% Plot Reaction Times BOXPLOTS
dataRT = [rt1', rt2', rt3'];
conditions = {'1-back', '2-back', '3-back'};
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
title('N-back Reaction Times', 'FontName', 'Arial', 'FontSize', 25);
hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/behavioral/AOC_rt_nback_boxplot.png');

% Stats
means = mean(dataRT, 'omitnan');
Stds = std(dataRT, 'omitnan');

% Display the results
disp('Means:');
disp(means);
disp('Stds:');
disp(Stds);
