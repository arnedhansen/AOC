%% AOC N-back Stats Overview
% Boxplots for all variables
% Percentage change barplots
% Alpha Differences vs. Gaze Differences

%% Load data 
clc
clear
close all
data = readtable('/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv');
data.ReactionTime = data.ReactionTime .* 1000;
data.GazeStd = (data.GazeStdX + data.GazeStdY) ./ 2;
variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'GazeStd', 'MSRate', 'Fixations', 'Saccades', 'AlphaPower', 'IAF'};
save_names = {'acc', 'rt', 'gazedev', 'ms', 'blink', 'fix', 'sacc', 'pow', 'iaf'};
colors = color_def('AOC');
subjects = unique(data.ID);

%% BOXPLOTS for each variable
y_axis_labels = {'Accuracy [%]', 'Reaction Time [ms]', 'Gaze Deviation [px]', 'Gaze Std [px]', 'Microsaccade Rate [ms/s]', 'Fixations', 'Saccades', 'Alpha Power [dB]', 'IAF [Hz]'};

font_size = 20;
for i = 1:length(variables)
    close all
    figure;
    set(gcf, 'Position', [100, 200, 800, 800], 'Color', 'w');
    hold on;

    % Set axis limits
    ylim([min(data.(variables{i})) max(data.(variables{i}))])
    xlim([0.5 2.5])

    % Create boxplot 
    boxplot(data.(variables{i}), data.Condition, 'Labels', {'1-back', '2-back', '3-back'}, 'Colors', 'k');
    set(gca, 'FontSize', 30); 

    % Define jitter
    jitterAmount = 0.25; % Adjust jitter if needed
    for cond = 1:3
        cond_data = data(data.Condition == cond, :);
        x_jittered{cond} = cond_data.Condition + (rand(size(cond_data.(variables{i}))) - 0.5) * jitterAmount;
    end

    % Connect points with a line
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 2 % Ensure the subject has data for both conditions
            x = [x_jittered{1}'; x_jittered{2}'];
            x = x(:, subj);
            y = subj_data.(variables{i});
            plot(x, y, '-o', 'Color', [0.5, 0.5, 0.5, 0.5], 'LineWidth', 2);
        end
    end

    % Scatter individual data points
    for cond = 1:3 
        cond_data = data(data.Condition == cond, :);
        scatter(x_jittered{cond}, cond_data.(variables{i}), 300, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(cond, :), 'MarkerFaceAlpha', 0.85, 'jitter', 'off', 'SizeData', 300);
    end
    
    % Add title and labels
    title(['N-back ', variables{i}], "FontSize", 40);
    ylabel(y_axis_labels{i}, "FontSize", 30);
    set(gca, 'Box', 'off');
    ax = gca;
    ax.FontSize = font_size;

    % Save individual subplot
    saveas(gcf, strcat('/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/boxplots/AOC_stats_boxplots_', save_names{i}, '_nback.png'));
end

%% BOXPLOTS OVERVIEW
close all
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

% Unique subject identifiers
subjects = unique(data.ID);

for i = 1:length(variables)
    subplot(3, 3, i);
    hold on;

    % Set axis limits
    ylim([min(data.(variables{i})) max(data.(variables{i}))])
    xlim([0.5 2.5])

    % Create boxplot
    boxplot(data.(variables{i}), data.Condition, 'Labels', {'1-back', '2-back', '3-back'}, 'Colors', 'k');
    set(gca, 'FontSize', 15); 

    % Overlay individual data points and connect them
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 2 % Ensure the subject has data for both conditions
            % X coordinates: condition indices
            x = subj_data.Condition;
            % Y coordinates: variable values
            y = subj_data.(variables{i});

            % Connect points with a line
            plot(x, y, '-o', 'Color', [0.5, 0.5, 0.5, 0.5], 'LineWidth', 0.5);
        end
    end

    % Scatter individual data points with jitter
    jitterAmount = 0.25;  % Amount of jitter in the x-direction
    for cond = 1:3  % Two conditions (Low Contrast = 1, High Contrast = 2)
        % Filter data based on condition
        cond_data = data(data.Condition == cond, :);

        % Add jitter and scatter the points
        x_jittered = cond_data.Condition + (rand(size(cond_data.(variables{i}))) - 0.5) * jitterAmount;  % Add random jitter to x-axis
        
        % Clear any previous scatter points (in case of overlaid dots)
        scatter(x_jittered, cond_data.(variables{i}), 100, 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(cond, :), 'jitter', 'off', 'SizeData', 100);
    end

    % Add title and labels
    title(variables{i}, "FontSize", 20);
    ylabel(y_axis_labels{i}, "FontSize", 15);
    hold off;
end
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/boxplots/AOC_stats_overview_boxplots_nback.png');

%% PERCENTAGE CHANGE BARPLOTS for each variable
close all;
% Preallocate percentage change matrix
percent_change = zeros(length(subjects), length(variables));

% Loop through each variable to calculate percentage change
for i = 1:length(variables)
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 3 % Ensure the subject has data for all conditions
            % Low and high contrast values
            low_value = subj_data{subj_data.Condition == 1, variables{i}};
            high_value = subj_data{subj_data.Condition == 3, variables{i}};

            % Calculate percentage change ((high - low) / low) * 100
            percent_change(subj, i) = ((high_value - low_value) / low_value) * 100;
        else
            percent_change(subj, i) = NaN; % Handle missing data
        end
    end

    % Create a new figure for the individual variable
    figure;
    set(gcf, 'Position', [100, 200, 1600, 1200], 'Color', 'w'); % Adjust size for individual plots
    set(gca, 'FontSize', 15); 
    hold on;

    % Bar plot for each participant
    bar(1:length(subjects), percent_change(:, i), 'FaceColor', 'k', 'EdgeColor', 'none');

    % Formatting
    xlim([0.5, length(subjects) + 0.5]);
    abw = max(abs([min(percent_change(:, i), [], 'omitnan'), max(percent_change(:, i), [], 'omitnan')]));
    ylim([-abw*1.25 abw*1.25]);
    if abw > 100
        ylim([-100 100]);
    end
    xticks(1:length(subjects));
    xticklabels(subjects);
    xlabel('Subjects', 'FontSize', 25); % Set font size for xlabel
    ylabel('% Change', 'FontSize', 25); % Set font size for ylabel
    title(['N-back ', variables{i}], 'FontSize', 30);

    hold off;

    % Save individual bar plot
    saveas(gcf, strcat('/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/barplots/AOC_stats_barplots_', save_names{i}, '_nback.png'));

    % Close the individual figure to free memory
    close(gcf);
end

%% PERCENTAGE CHANGE BARPLOTS OVERVIEW
close all
% Preallocate percentage change matrix
percent_change = zeros(length(subjects), length(variables));

% Loop through each variable to calculate percentage change
for i = 1:length(variables)
    for subj = 1:length(subjects)
        % Extract data for this subject
        subj_data = data(data.ID == subjects(subj), :);

        if height(subj_data) == 3 % Ensure the subject has data for all conditions
            % Low and high contrast values
            low_value = subj_data{subj_data.Condition == 1, variables{i}};
            high_value = subj_data{subj_data.Condition == 3, variables{i}};

            % Calculate percentage change ((high - low) / low) * 100
            percent_change(subj, i) = ((high_value - low_value) / low_value) * 100;
        else
            percent_change(subj, i) = NaN; % Handle missing data
        end
    end
end

% Plot percentage change bar plots
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

for i = 1:length(variables)
    subplot(3, 3, i);
    hold on;

    % Bar plot for each participant
    bar(1:length(subjects), percent_change(:, i), 'FaceColor', 'k', 'EdgeColor', 'none');

    % Formatting
    xlim([0.5, length(subjects) + 0.5]);
    abw = max(abs([min(percent_change(:, i), [], 'omitnan'), max(percent_change(:, i), [], 'omitnan')]));
    ylim([-abw*1.25 abw*1.25]);
    if abw > 100
        ylim([-100 100]);
    end
    xticks(1:length(subjects));
    xticklabels(subjects);
    xlabel('Subjects');
    ylabel('% Change');
    title(variables{i}, 'FontSize', 20);
    hold off;
end
sgtitle('Percentage Change (HC - LC)', 'FontSize', 24);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/barplots/AOC_stats_overview_barplots_percentage_change_nback.png');

%% Plot DIFFERENCES in ALPHA POWER against GAZE METRICS
% MICROSACCADES,GAZE DEVIATION and labels by ET (BLINKS, FIXATIONS, SACCADES)
close all
% Define the gaze metrics and their labels
gaze_metrics = {'GazeDeviation', 'MSRate', 'Fixations', 'Saccades'};
gaze_metric_labels = {'Gaze Deviation [px]', 'Microsaccade Rate [ms/s]', 'Fixations', 'Saccades'};

% Loop through each condition (assuming conditions are 1, 2 and 3)
for cond = 1:3
    % Extract data for the current condition
    cond_data = data(data.Condition == cond, :);
    
    % Create a new figure for the current condition
    figure;
    set(gcf, 'Position', [100, 100, 1200, 800], 'Color', 'w');
    
    % Number of gaze metrics to plot
    n_metrics = numel(gaze_metrics);
    
    % Create a 2x2 subplot for each gaze metric
    for m = 1:n_metrics
        subplot(2, 2, m);
        hold on;
        
        % Scatter plot: x-axis is the gaze metric, y-axis is Alpha Power
        scatter(cond_data.(gaze_metrics{m}), cond_data.AlphaPower, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
        
        % Label axes and title the subplot
        xlabel(gaze_metric_labels{m}, 'FontSize', 15);
        ylabel('Alpha Power [dB]', 'FontSize', 15);
        title([num2str(cond) '-back : Alpha Power vs ' gaze_metric_labels{m}], 'FontSize', 18);
        
        % Add a least-squares fit line
        lsline;
        hold off;
    end
    
    % Add a super-title for the figure
    sgtitle([num2str(cond) '-back: Association between Alpha Power and Gaze Metrics'], 'FontSize', 20);
    
    % Save the figure (adjust the save path as necessary)
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/AOC_stats_overview_associations_pow_' num2str(cond) 'back.png']);
end

%% Plot DIFFERENCES in ALPHA POWER against GAZE METRICS - PERCENTAGE CHANGE
% MICROSACCADES,GAZE DEVIATION and labels by ET (BLINKS, FIXATIONS, SACCADES)
close all
%variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'MSRate', 'Blinks', 'Fixations', 'Saccades', 'AlphaPower', 'IAF'};
alpha_power_diff = percent_change(:, 8);

% Extract the corresponding values for gaze metrics
gaze_deviation_diff = percent_change(:, 3);
ms_rate_diff = percent_change(:, 4);
fix_diff = percent_change(:, 6);
sacc_diff = percent_change(:, 7);

% Define gaze metrics
metrics = {
    'Gaze Deviation Diff. [%]', gaze_deviation_diff;
    'Microsaccade Rate Diff. [%]', ms_rate_diff;
    'Fixation Diff. [%]', fix_diff;
    'Saccade Diff. [%]', sacc_diff
};

% Number of metrics
n_metrics = size(metrics, 1);

% Create figure for Alpha Power
figure;
set(gcf, 'Position', [0, 0, 800, 800], 'Color', 'w');
for i = 1:n_metrics
    subplot(ceil(n_metrics / 2), 2, i);
    hold on;
    scatter(metrics{i, 2}, alpha_power_diff, 36, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
    xlabel(metrics{i, 1}, 'FontSize', 15);
    ylabel('Alpha Power Diff. [%]', 'FontSize', 15);
    xlim([-100, 100]);
    ylim([-100, 100]);
    xline(0, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
    yline(0, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.5, 'Alpha', 0.25);
    title(['Alpha vs ', metrics{i, 1}], 'FontSize', 18);
    lsline()
    hold off;
end
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/AOC_stats_overview_associations_pow_percentage_nback.png');

%% CORRELATION matrix
close all
% Calculate the correlation matrix between variables
corr_matrix = corr(table2array(data(:, variables)), 'Rows', 'pairwise');

% Correlation matrix
% Generate heatmap
figure('Color', 'white');
set(gcf, 'Position', [200, 400, 1200, 1500]);
h = heatmap(corr_matrix);

% Define and add color map
num_colors = 256;
cpoints = [0, 0, 0, 1;      % blue at the lowest value
           0.46, 1, 1, 1;  % transition to white starts
           0.54, 1, 1, 1;  % transition from white ends
           1, 1, 0, 0];    % red at the highest value

% Preallocate the colormap array.
cmap = zeros(num_colors, 3);

% Linearly interpolate to fill in the colormap.
for i=1:size(cmap,1)
    % Normalize the current index to the 0-1 range based on the colormap size.
    val = (i-1)/(num_colors-1);
    % Find the first point where the value is greater than the current normalized value.
    idx = find(cpoints(:,1) >= val, 1);
    if idx == 1
        % Use the first colour for values below the first point.
        cmap(i,:) = cpoints(1,2:4);
    else
        % Linearly interpolate between the two bounding colours.
        range = cpoints(idx,1) - cpoints(idx-1,1);
        frac = (val - cpoints(idx-1,1)) / range;
        cmap(i,:) = (1-frac)*cpoints(idx-1,2:4) + frac*cpoints(idx,2:4);
    end
end
colormap(h, cmap);

% Customization
h.Title = 'Correlation Heatmap';
h.XDisplayLabels = variables; % Assuming varNames contains variable names excluding 'Date', 'Weekday', and 'Overall Score'
h.YDisplayLabels = variables;
h.ColorLimits = [-1, 1]; % Color limits to ensure proper color mapping
h.FontSize = 12; % Increase the size of the x and y axis tick labels
hTitle.FontSize = 25;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/stats/overview/AOC_stats_overview_correlation_matrix_nback.png');
