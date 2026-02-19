%% AOC Omnibus — Early vs Late Window Comparison
% Compares SPL between alpha reduction/amplification groups in early (0-0.5s) vs late (1-2s) windows.
% Tests if the effect changes over time using 2-way ANOVA (Time Window × Alpha Group).
%
% Key outputs:
%   Bar plot comparing early vs late windows; 2-way ANOVA results

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/early_late';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Load CSV data
csv_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_omnibus_raincloud_data.csv';
raincloud_data = readtable(csv_path);
sternberg_data = raincloud_data(strcmpi(raincloud_data.Task, 'sternberg'), :);

% Convert ID to string
if isnumeric(sternberg_data.ID)
    sternberg_data.ID = cellstr(num2str(sternberg_data.ID, '%d'));
elseif iscell(sternberg_data.ID)
    sternberg_data.ID = cellfun(@(x) num2str(x), sternberg_data.ID, 'UniformOutput', false);
elseif isstring(sternberg_data.ID)
    sternberg_data.ID = cellstr(sternberg_data.ID);
end

%% Define time windows
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);

early_window = t_plot_full >= 0 & t_plot_full <= 0.5;
late_window = t_plot_full >= 1 & t_plot_full <= 2;

%% Load data: mean alpha power per subject (across all conditions)
subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);
mean_alpha = nan(n_subj, 1);
spl_early = nan(n_subj, 1);
spl_late = nan(n_subj, 1);

for s = 1:n_subj
    subj_id = subj_ids{s};
    
    % Get mean alpha power across all conditions
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    alpha_vals = [];
    for c = 1:height(subj_rows)
        if ~isnan(subj_rows.AlphaPower(c)) && subj_rows.AlphaPower(c) ~= 0
            alpha_vals(end+1) = subj_rows.AlphaPower(c);
        end
    end
    
    if isempty(alpha_vals)
        continue;
    end
    mean_alpha(s) = nanmean(alpha_vals);
    
    % Load gaze data
    datapath_gaze = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features', ...
        subj_id, 'gaze', 'gaze_series_sternberg_trials.mat');
    
    if ~isfile(datapath_gaze)
        continue;
    end
    
    try
        load(datapath_gaze, 'ScanPathSeries', 'ScanPathSeriesT', 'trialinfo');
    catch
        continue;
    end
    
    % Parse trialinfo
    if exist('trialinfo', 'var')
        if size(trialinfo, 1) == 2
            conds = trialinfo(1, :)';
        elseif size(trialinfo, 2) == 2
            conds = trialinfo(:, 1);
        else
            continue;
        end
    else
        continue;
    end
    
    % Compute SPL time series (all conditions)
    subj_trials_full = nan(numel(ScanPathSeries), Tf);
    for trl = 1:numel(ScanPathSeries)
        srl_full = ScanPathSeries{trl};
        tt_full = ScanPathSeriesT{trl};
        if isempty(srl_full) || isempty(tt_full)
            continue;
        end
        try
            subj_trials_full(trl, :) = interp1(tt_full, srl_full, t_plot_full, 'linear', NaN);
        catch
        end
    end
    
    % Average across all trials and compute window means
    subj_mean = nanmean(subj_trials_full, 1);
    spl_early(s) = nanmean(subj_mean(early_window));
    spl_late(s) = nanmean(subj_mean(late_window));
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Prepare data for ANOVA
% Create long format
data_table = [];
for s = 1:n_subj
    if ~isfinite(mean_alpha(s)) || ~isfinite(spl_early(s)) || ~isfinite(spl_late(s))
        continue;
    end
    
    % Determine group
    if mean_alpha(s) < 0
        group = 1; % Reduction
    else
        group = 2; % Amplification
    end
    
    % Early window
    data_table(end+1, :) = [group, 1, spl_early(s)]; % group, window (1=early), SPL
    % Late window
    data_table(end+1, :) = [group, 2, spl_late(s)]; % group, window (2=late), SPL
end

if ~isempty(data_table)
    anova_data = table();
    anova_data.Group = categorical(data_table(:, 1), [1 2], {'Reduction', 'Amplification'});
    anova_data.Window = categorical(data_table(:, 2), [1 2], {'Early (0-0.5s)', 'Late (1-2s)'});
    anova_data.SPL = data_table(:, 3);
    
    %% 2-way ANOVA
    fprintf('\n=== 2-Way ANOVA: Time Window × Alpha Group ===\n');
    try
        [p, tbl, stats] = anovan(anova_data.SPL, ...
            {double(anova_data.Window), double(anova_data.Group)}, ...
            'model', 'interaction', ...
            'varnames', {'Window', 'Group'}, ...
            'display', 'on');
        
        fprintf('Main effect Window: p = %.4f\n', p(1));
        fprintf('Main effect Group: p = %.4f\n', p(2));
        fprintf('Interaction: p = %.4f\n', p(3));
    catch ME
        warning('ANOVA failed: %s', ME.message);
        p = [nan, nan, nan];
    end
    
    %% Compute means and SEMs
    means = nan(2, 2); % windows × groups
    sems = nan(2, 2);
    ns = nan(2, 2);
    
    for w = 1:2
        for g = 1:2
            mask = (data_table(:, 1) == g) & (data_table(:, 2) == w);
            vals = data_table(mask, 3);
            means(w, g) = nanmean(vals);
            sems(w, g) = nanstd(vals) / sqrt(length(vals));
            ns(w, g) = length(vals);
        end
    end
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1000 800])
    
    x_pos = [1, 2];
    x_offset = [-0.15, 0.15];
    width = 0.25;
    
    hold on
    for w = 1:2
        for g = 1:2
            x = x_pos(w) + x_offset(g);
            y = means(w, g);
            err = sems(w, g);
            n = ns(w, g);
            
            % Bar
            if g == 1
                bar(x, y, width, 'FaceColor', colors(1,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
            else
                bar(x, y, width, 'FaceColor', colors(3,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
            end
            
            % Error bar
            plot([x x], [y-err y+err], 'k-', 'LineWidth', 2)
            plot([x-0.05 x+0.05], [y-err y-err], 'k-', 'LineWidth', 2)
            plot([x-0.05 x+0.05], [y+err y+err], 'k-', 'LineWidth', 2)
            
            % Sample size
            text(x, y+err+0.02, sprintf('n=%d', n), ...
                'HorizontalAlignment', 'center', 'FontSize', fontSize-8);
        end
    end
    
    xlim([0.5 2.5])
    set(gca, 'XTick', x_pos, 'XTickLabel', {'Early (0-0.5s)', 'Late (1-2s)'}, ...
        'FontSize', fontSize-2)
    ylabel('Mean Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Early vs Late Window: Alpha Reduction vs Amplification', ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    
    legend({'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
    grid on
    set(gca, 'GridAlpha', 0.2)
    box on
    
    % Add significance annotation
    if all(isfinite(p))
        y_max = max(means(:) + sems(:)) + 0.1;
        text(1.5, y_max, sprintf('Interaction: p=%.4f', p(3)), ...
            'HorizontalAlignment', 'center', 'FontSize', fontSize-6, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    
    saveas(gcf, fullfile(output_dir, 'AOC_omnibus_early_late_window.png'));
    fprintf('\nSaved figure to: %s\n', fullfile(output_dir, 'AOC_omnibus_early_late_window.png'));
end

fprintf('\n=== Done ===\n');
