%% Visualization of gaze deviation for AOC Sternberg

% Visualizations:
%   Time Series of Gaze Position Using a Sliding Window

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

%% Compute gaze metrics
numSubjects = length(subjects);
allGazeX = cell(1, numSubjects);
allGazeY = cell(1, numSubjects);
allErrorGazeX = cell(1, numSubjects);
allErrorGazeY = cell(1, numSubjects);

for subj = 1:numSubjects
    gazePath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(gazePath, 'dataET_sternberg.mat'));

    numTrials = length(dataETlong.trial);
    % Preallocate cell arrays for each trial’s gaze data
    lGazeX = cell(1, numTrials);
    lGazeY = cell(1, numTrials);

    for i = 1:numTrials
        lGazeX{i} = dataETlong.trial{i}(1, :);
        lGazeY{i} = dataETlong.trial{i}(2, :);
    end

    % Get unique conditions (adjust as needed)
    conditions = unique(dataETlong.trialinfo);
    conditions = (conditions-50)/2;
    timePoints = length(dataETlong.time{1});

    subjectAverageGazeX = cell(1, length(conditions));
    subjectAverageGazeY = cell(1, length(conditions));
    subjectErrorGazeX = cell(1, length(conditions));
    subjectErrorGazeY = cell(1, length(conditions));

    for condIdx = 1:length(conditions)
        cond = conditions(condIdx);
        cond = (cond*2)+50;
        condTrialsIdx = find(dataETlong.trialinfo == cond);

        condGazeX = zeros(length(condTrialsIdx), timePoints);
        condGazeY = zeros(length(condTrialsIdx), timePoints);

        for j = 1:length(condTrialsIdx)
            condGazeX(j, :) = lGazeX{condTrialsIdx(j)};
            condGazeY(j, :) = lGazeY{condTrialsIdx(j)};
        end

        subjectAverageGazeX{condIdx} = mean(condGazeX, 1);
        subjectErrorGazeX{condIdx} = std(condGazeX, 0, 1) / sqrt(length(condTrialsIdx));
        subjectAverageGazeY{condIdx} = mean(condGazeY, 1);
        subjectErrorGazeY{condIdx} = std(condGazeY, 0, 1) / sqrt(length(condTrialsIdx));
    end

    allGazeX{subj} = subjectAverageGazeX;
    allGazeY{subj} = subjectAverageGazeY;
    allErrorGazeX{subj} = subjectErrorGazeX;
    allErrorGazeY{subj} = subjectErrorGazeY;
end

% Compute grand average and standard error for each condition
grandAverageGazeX = cell(1, length(conditions));
grandErrorGazeX = cell(1, length(conditions));
grandAverageGazeY = cell(1, length(conditions));
grandErrorGazeY = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    allSubjectGazeX = zeros(numSubjects, timePoints);
    allSubjectGazeY = zeros(numSubjects, timePoints);

    for subj = 1:numSubjects
        allSubjectGazeX(subj, :) = allGazeX{subj}{condIdx};
        allSubjectGazeY(subj, :) = allGazeY{subj}{condIdx};
    end

    grandAverageGazeX{condIdx} = mean(allSubjectGazeX, 1);
    grandErrorGazeX{condIdx} = std(allSubjectGazeX, 0, 1) / sqrt(numSubjects);
    grandAverageGazeY{condIdx} = mean(allSubjectGazeY, 1);
    grandErrorGazeY{condIdx} = std(allSubjectGazeY, 0, 1) / sqrt(numSubjects);
end

% Optionally convert gaze positions to percentage deviation from centre (400,300)
for subj = 1:numSubjects
    for condIdx = 1:length(conditions)
        allGazeX{subj}{condIdx} = ((allGazeX{subj}{condIdx} / 400) - 1) * 100;
        allGazeY{subj}{condIdx} = ((allGazeY{subj}{condIdx} / 300) - 1) * 100;
        allErrorGazeX{subj}{condIdx} = (allErrorGazeX{subj}{condIdx} / 400) * 100;
        allErrorGazeY{subj}{condIdx} = (allErrorGazeY{subj}{condIdx} / 300) * 100;
    end
end

grandAverageGazeX = cellfun(@(x) ((x / 400) - 1) * 100, grandAverageGazeX, 'UniformOutput', false);
grandAverageGazeY = cellfun(@(y) ((y / 300) - 1) * 100, grandAverageGazeY, 'UniformOutput', false);
grandErrorGazeX = cellfun(@(x) (x / 400) * 100, grandErrorGazeX, 'UniformOutput', false);
grandErrorGazeY = cellfun(@(y) (y / 300) * 100, grandErrorGazeY, 'UniformOutput', false);

% Define sliding window parameters
windowLength = 0.05; % 100 ms
dt = dataETlong.time{1}(2) - dataETlong.time{1}(1);
windowPoints = round(windowLength / dt);
moving_average = @(data, win_length) movmean(data, win_length, 2);
timeVec = dataETlong.time{1};

%% Grand Average
close all
figure;
set(gcf, 'Position', [100, 100, 1800, 1000], 'Color', 'W');

% X-Axis
subplot(2,1,1);
hold on;
eb_x = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    grandGazeX_smoothed = moving_average(grandAverageGazeX{condIdx}, windowPoints);
    grandErrorX_smoothed = moving_average(grandErrorGazeX{condIdx}, windowPoints);

    % Create the shaded error bar without preset colour options
    eb_x{condIdx} = shadedErrorBar(timeVec, grandGazeX_smoothed, grandErrorX_smoothed, {'-'}, 0);

    % Set colours for the main line, patch and edges using the corresponding row in colors
    eb_x{condIdx}.mainLine.Color = colors(condIdx, :);
    eb_x{condIdx}.patch.FaceColor = colors(condIdx, :);
    set(eb_x{condIdx}.mainLine, 'LineWidth', 1, 'Color', colors(condIdx, :));
    set(eb_x{condIdx}.edge(1), 'Color', colors(condIdx, :));
    set(eb_x{condIdx}.edge(2), 'Color', colors(condIdx, :));
    set(eb_x{condIdx}.patch, 'FaceAlpha', 0.5);
end

% Adjust plot aesthetics
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
xline(0, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
title('Sternberg Grand Average X Deviation', 'FontSize', 18);
xlabel('Time [s]');
ylabel('X Gaze Deviation [%]');
ylim([-35 35]);
xlim([-0.25 2]);
text(-0.2, -5, 'LEFT', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.2, 5, 'RIGHT', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20)

% Construct legend using only the three main lines
legend([eb_x{1}.mainLine, eb_x{2}.mainLine, eb_x{3}.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'FontName', 'Arial', 'FontSize', 20);
hold off;

% Y-Axis
subplot(2,1,2);
hold on;
eb_y = cell(1, length(conditions));

for condIdx = 1:length(conditions)
    grandGazeY_smoothed = moving_average(grandAverageGazeY{condIdx}, windowPoints);
    grandErrorY_smoothed = moving_average(grandErrorGazeY{condIdx}, windowPoints);

    % Create the shaded error bar without preset colour options
    eb_y{condIdx} = shadedErrorBar(timeVec, grandGazeY_smoothed, grandErrorY_smoothed, {'-'}, 0);

    % Set colours for the main line, patch and edges using the corresponding row in colors
    eb_y{condIdx}.mainLine.Color = colors(condIdx, :);
    eb_y{condIdx}.patch.FaceColor = colors(condIdx, :);
    set(eb_y{condIdx}.mainLine, 'LineWidth', 1, 'Color', colors(condIdx, :));
    set(eb_y{condIdx}.edge(1), 'Color', colors(condIdx, :));
    set(eb_y{condIdx}.edge(2), 'Color', colors(condIdx, :));
    set(eb_y{condIdx}.patch, 'FaceAlpha', 0.5);
end

% Adjust plot aesthetics
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
xline(0, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
title('Sternberg Grand Average Y Deviation', 'FontSize', 18);
xlabel('Time [s]');
ylabel('Y Deviation [%]');
ylim([-25 25]);
xlim([-0.25 2]);
text(-0.2, -5, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
text(-0.2, 5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20)

% Construct legend using only the three main lines
legend([eb_y{1}.mainLine, eb_y{2}.mainLine, eb_y{3}.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'FontName', 'Arial', 'FontSize', 20);

hold off;

% Save
saveas(gcf, fullfile('/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/', ...
    'AOC_dev_sternberg_all.png'));

%% INDIVIDUAL PLOTS
close all
for subj = 1:numSubjects
    figure;
    set(gcf, 'Position', [100, 100, 1800, 1000], 'Color', 'W');

    % X-Axis Plot
    subplot(2,1,1);
    hold on;
    eb = cell(1, length(conditions)); % to store handles for each condition
    for condIdx = 1:length(conditions)
        subjectGazeX = moving_average(allGazeX{subj}{condIdx}, windowPoints);
        subjectErrorX = moving_average(allErrorGazeX{subj}{condIdx}, windowPoints);
        % Create the shaded error bar without predefined colour options
        eb{condIdx} = shadedErrorBar(timeVec, subjectGazeX, subjectErrorX, {'-'}, 0);

        % Manually set the colour for the main line and the patch
        eb{condIdx}.mainLine.Color = colors(condIdx, :);
        eb{condIdx}.patch.FaceColor = colors(condIdx, :);
        set(eb{condIdx}.mainLine, 'LineWidth', 1, 'Color', colors(condIdx, :));
        set(eb{condIdx}.edge(1), 'Color', colors(condIdx, :));
        set(eb{condIdx}.edge(2), 'Color', colors(condIdx, :));
        set(eb{condIdx}.patch, 'FaceAlpha', 0.5);
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
    title(['Average X Deviation - Subject ' subjects{subj}], 'FontSize', 18);
    xlabel('Time [s]');
    ylabel('X Deviation [%]');
    text(0.05, -5, 'LEFT', 'FontSize', 20, 'FontWeight', 'bold');
    text(0.05, 5, 'RIGHT', 'FontSize', 20, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20)
    legend([eb{1}.mainLine, eb{2}.mainLine, eb{3}.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
        'FontName', 'Arial', 'FontSize', 20);
    hold off;

    % Y-Axis Plot
    subplot(2,1,2);
    hold on;
    eb = cell(1, length(conditions));
    for condIdx = 1:length(conditions)
        subjectGazeY = moving_average(allGazeY{subj}{condIdx}, windowPoints);
        subjectErrorY = moving_average(allErrorGazeY{subj}{condIdx}, windowPoints);
        eb{condIdx} = shadedErrorBar(timeVec, subjectGazeY, subjectErrorY, {'-'}, 0);
        eb{condIdx}.mainLine.Color = colors(condIdx, :);
        eb{condIdx}.patch.FaceColor = colors(condIdx, :);
        set(eb{condIdx}.mainLine, 'LineWidth', 1, 'Color', colors(condIdx, :));
        set(eb{condIdx}.edge(1), 'Color', colors(condIdx, :));
        set(eb{condIdx}.edge(2), 'Color', colors(condIdx, :));
        set(eb{condIdx}.patch, 'FaceAlpha', 0.5);
    end
    yline(0, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
    title(['Average Y Deviation - Subject ' subjects{subj}], 'FontSize', 18);
    xlabel('Time [s]');
    ylabel('Y Deviation [%]');
    text(0.05, -5, 'DOWN', 'FontSize', 20, 'FontWeight', 'bold');
    text(0.05, 5, 'UP', 'FontSize', 20, 'FontWeight', 'bold');
    set(gca, 'FontSize', 20)
    legend([eb{1}.mainLine, eb{2}.mainLine, eb{3}.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
        'FontName', 'Arial', 'FontSize', 20);
    hold off;

    % Save
    saveas(gcf, fullfile('/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/', ...
        ['AOC_dev_sternberg_subj' subjects{subj} '.png']));
end
