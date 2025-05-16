%% AOC Split of Nback EEG Data
% Split Nback participants by condition (1-back, 2-back, 3-back) order. Check if gaze metrics differ
% between these groups.

%% Load data
clc
clear
close all
data = readtable('/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv');
data.ReactionTime = data.ReactionTime .* 1000;
data.Condition = data.Condition .* 0.5; % Set condition to 1-3 for loops
data.GazeStd = (data.GazeStdX + data.GazeStdY) ./ 2;
variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'GazeStd', 'MSRate', 'Fixations', 'Saccades', 'AlphaPower', 'IAF'};
save_names = {'acc', 'rt', 'gazedev', 'ms', 'blink', 'fix', 'sacc', 'pow', 'iaf'};
colors = color_def('AOC');
subjects = unique(data.ID);

%% Compute percentage changes
clear percent_struct
% Preallocate structure array
percent_struct(length(subjects)) = struct();

% Loop through each subject
for subj = 1:length(subjects)
    subj_id = subjects(subj);
    subj_data = data(data.ID == subj_id, :);

    % Initialize subject's struct
    percent_struct(subj).ID = subj_id;

    % Only calculate if subject has data for all 3 conditions
    if height(subj_data) == 3
        for i = 1:length(variables)
            low_value = subj_data{subj_data.Condition == 1, variables{i}};
            high_value = subj_data{subj_data.Condition == 3, variables{i}};
            if ~isempty(low_value) && ~isempty(high_value)
                percent_struct(subj).(save_names{i}) = ((high_value - low_value) / low_value) * 100;
            else
                percent_struct(subj).(save_names{i}) = NaN;
            end
            % Add Alpha Power Direction Column
            if strcmp(save_names{i}, 'pow')
                if percent_struct(subj).(save_names{i}) < 0
                    percent_struct(subj).AlphaPowDirection = -1;
                elseif percent_struct(subj).(save_names{i}) == 0 || isnan(percent_struct(subj).(save_names{i}))
                    percent_struct(subj).AlphaPowDirection = 0;
                elseif percent_struct(subj).(save_names{i}) > 0
                    percent_struct(subj).AlphaPowDirection = 1;
                end
            end
        end
    else
        % Fill with NaNs if data is missing
        for i = 1:length(variables)
            percent_struct(subj).(save_names{i}) = NaN;
        end
    end
end
%% PERCENTAGE CHANGE BARPLOT Nback
close all

% Plot percentage change bar plots
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');

for i = 1:length(variables)
    subplot(3, 3, i);
    hold on;

    % Get the relevant percentage change values
    values = [percent_struct.(save_names{i})];

    % Plot each bar with colour based on AlphaPowDirection
    for subj = 1:length(subjects)
        direction = percent_struct(subj).AlphaPowDirection;

        % Determine colour
        if direction == -1
            bar_color = [0 0 1]; % blue
        elseif direction == 1
            bar_color = [1 0 0]; % red
        else
            bar_color = [0 0 0]; % black
        end

        % Plot individual bar
        bar(subj, values(subj), 'FaceColor', bar_color, 'EdgeColor', 'none');
    end

    % Formatting
    xlim([0.5, length(subjects) + 0.5]);
    abw = max(abs([min([percent_struct.(save_names{i})], [], 'omitnan'), max([percent_struct.(save_names{i})], [], 'omitnan')]));
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
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/tests/AOC_split_nback.png');
