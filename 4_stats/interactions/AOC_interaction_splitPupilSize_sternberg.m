%% AOC Interaction — Sternberg: Split by PupilSize
% Loads merged_data_sternberg_trials, plots scatter of each numeric variable vs PupilSizeFull to inspect associations. Saves figures.
%
% Key outputs:
%   Scatter figures (variable vs PupilSizeFull)

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC', 0);
%subjects = subjects(1:10)

% Figure config
color_map = cbrewer('seq', 'Reds', 64, 'spline');
colors = color_def('AOC');
fontSize = 30;

%% Load data
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

%% Visualize weird pupil size distribution
close all

tbl = merged_data_sternberg_trials;
vars = tbl.Properties.VariableNames;

for v = 1:length(vars)
    close all
    varname = vars{v};

    % skip pupil itself and skip non-numeric variables
    if strcmp(varname, 'PupilSizeFull')
        continue
    end
    if ~isnumeric(tbl.(varname))
        continue
    end
    if strcmp(varname, 'Accuracy')
        xlim([-.5 1.5])
    end

    figure();
    set(gcf, 'Position', [0 0 1532 982], 'Color', 'w');

    scatter(tbl.(varname), tbl.PupilSizeFull, 10, colors(3,:), ...
        'filled', 'MarkerFaceAlpha', 0.7);

    xlabel(varname, 'Interpreter', 'none')
    ylabel('PupilSizeFull')
    title(['Pupil size vs. ' varname], 'Interpreter', 'none')

    ylim([0 6500])
    box on
    set(gca, 'FontSize', 30);
    saveas(gcf, ['/Users/Arne/Documents/GitHub/AOC/3_visualization/interactions/pupilSizeGroups/', varname, '_pupilSize.png'])
end