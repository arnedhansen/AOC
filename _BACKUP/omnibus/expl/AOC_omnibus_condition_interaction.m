%% AOC Omnibus — Condition × Alpha Group Interaction
% Compares the alpha reduction/amplification effect across WM loads (2, 4, 6).
% Tests if the effect size differs across conditions using 2-way ANOVA.
%
% Key outputs:
%   Interaction plot; 2-way ANOVA results; effect sizes per condition

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/condition_interaction';
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

%% Load gaze data and compute mean SPL per condition
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [22, 24, 26];
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);

% Time window for averaging (0-2s task period)
time_window = t_plot_full >= 0 & t_plot_full <= 2;

% Storage: subject × condition × alpha group
subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);
mean_spl_data = nan(n_subj, 3, 2); % subjects × conditions × groups (1=reduction, 2=amplification)
alpha_vals = nan(n_subj, 3); % subjects × conditions (for grouping)

% Process each subject
for s = 1:n_subj
    subj_id = subj_ids{s};
    subj_id_num = str2double(subj_id);
    
    % Get alpha power for this subject per condition
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    for c = 1:3
        cond_rows = subj_rows(strcmp(subj_rows.Condition, conditions{c}), :);
        if ~isempty(cond_rows) && ~isnan(cond_rows.AlphaPower(1))
            alpha_vals(s, c) = cond_rows.AlphaPower(1);
        end
    end
    
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
    
    % Process each condition
    for c = 1:3
        cond_code = cond_codes(c);
        trial_mask = (conds == cond_code);
        
        if ~any(trial_mask)
            continue;
        end
        
        % Compute SPL time series for this condition
        subj_trials_full = nan(numel(ScanPathSeries), Tf);
        for trl = 1:numel(ScanPathSeries)
            if ~trial_mask(trl)
                continue;
            end
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
        
        % Average across trials and time window
        subj_mean = nanmean(subj_trials_full(trial_mask, :), 1);
        mean_spl_window = nanmean(subj_mean(time_window));
        
        if isfinite(mean_spl_window) && isfinite(alpha_vals(s, c)) && alpha_vals(s, c) ~= 0
            % Determine group (1=reduction, 2=amplification)
            if alpha_vals(s, c) < 0
                mean_spl_data(s, c, 1) = mean_spl_window;
            else
                mean_spl_data(s, c, 2) = mean_spl_window;
            end
        end
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Prepare data for ANOVA
% Create long format table
data_table = [];
for c = 1:3
    for g = 1:2
        group_name = {'Reduction', 'Amplification'};
        vals = mean_spl_data(:, c, g);
        vals = vals(isfinite(vals));
        for v = 1:length(vals)
            data_table(end+1, :) = [c, g, vals(v)];
        end
    end
end

% Convert to table format for easier analysis
if ~isempty(data_table)
    anova_data = table();
    anova_data.Condition = categorical(data_table(:, 1), [1 2 3], conditions);
    anova_data.Group = categorical(data_table(:, 2), [1 2], {'Reduction', 'Amplification'});
    anova_data.SPL = data_table(:, 3);
    
    % Compute means and SEMs for plotting
    means = nan(3, 2);
    sems = nan(3, 2);
    ns = nan(3, 2);
    
    for c = 1:3
        for g = 1:2
            vals = mean_spl_data(:, c, g);
            vals = vals(isfinite(vals));
            means(c, g) = nanmean(vals);
            sems(c, g) = nanstd(vals) / sqrt(length(vals));
            ns(c, g) = length(vals);
        end
    end
    
    %% 2-way ANOVA
    fprintf('\n=== 2-Way ANOVA: Condition × Alpha Group ===\n');
    try
        [p, tbl, stats] = anovan(anova_data.SPL, ...
            {double(anova_data.Condition), double(anova_data.Group)}, ...
            'model', 'interaction', ...
            'varnames', {'Condition', 'Group'}, ...
            'display', 'on');
        
        fprintf('Main effect Condition: p = %.4f\n', p(1));
        fprintf('Main effect Group: p = %.4f\n', p(2));
        fprintf('Interaction: p = %.4f\n', p(3));
    catch ME
        warning('ANOVA failed: %s', ME.message);
        p = [nan, nan, nan];
    end
    
    %% Effect sizes per condition
    fprintf('\n=== Effect Sizes (Cohen''s d) per Condition ===\n');
    cohens_d = nan(3, 1);
    for c = 1:3
        red_vals = mean_spl_data(:, c, 1);
        amp_vals = mean_spl_data(:, c, 2);
        red_vals = red_vals(isfinite(red_vals));
        amp_vals = amp_vals(isfinite(amp_vals));
        
        if length(red_vals) >= 2 && length(amp_vals) >= 2
            pooled_std = sqrt(((length(red_vals)-1)*nanvar(red_vals) + ...
                (length(amp_vals)-1)*nanvar(amp_vals)) / ...
                (length(red_vals) + length(amp_vals) - 2));
            if pooled_std > 0
                cohens_d(c) = (nanmean(amp_vals) - nanmean(red_vals)) / pooled_std;
                fprintf('%s: d = %.3f (n_red=%d, n_amp=%d)\n', ...
                    conditions{c}, cohens_d(c), length(red_vals), length(amp_vals));
            end
        end
    end
    
    %% Plot interaction
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1200 800])
    
    x_pos = [1, 2, 3];
    x_offset = [-0.15, 0.15];
    
    hold on
    for c = 1:3
        for g = 1:2
            x = x_pos(c) + x_offset(g);
            y = means(c, g);
            err = sems(c, g);
            n = ns(c, g);
            
            % Error bar
            plot([x x], [y-err y+err], 'k-', 'LineWidth', 2)
            
            % Point
            if g == 1
                scatter(x, y, 150, colors(1,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            else
                scatter(x, y, 150, colors(3,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            end
            
            % Sample size text
            text(x, y+err+0.02, sprintf('n=%d', n), ...
                'HorizontalAlignment', 'center', 'FontSize', fontSize-8);
        end
        
        % Connect points for each group
        plot(x_pos(c) + x_offset, means(c, :), 'k-', 'LineWidth', 1, 'Alpha', 0.3);
    end
    
    % Add effect size annotations
    for c = 1:3
        if isfinite(cohens_d(c))
            text(x_pos(c), max(means(c, :)) + max(sems(c, :)) + 0.05, ...
                sprintf('d=%.2f', cohens_d(c)), ...
                'HorizontalAlignment', 'center', 'FontSize', fontSize-6, 'FontWeight', 'bold');
        end
    end
    
    xlim([0.5 3.5])
    set(gca, 'XTick', x_pos, 'XTickLabel', conditions, 'FontSize', fontSize-2)
    ylabel('Mean Scan Path Length [px] (0-2s)', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Condition × Alpha Group Interaction', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    
    legend({'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
    grid on
    set(gca, 'GridAlpha', 0.2)
    box on
    
    % Add significance annotations
    if all(isfinite(p))
        y_max = max(means(:) + sems(:)) + 0.15;
        text(2, y_max, sprintf('Interaction: p=%.4f', p(3)), ...
            'HorizontalAlignment', 'center', 'FontSize', fontSize-6, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    
    saveas(gcf, fullfile(output_dir, 'AOC_omnibus_condition_interaction.png'));
    fprintf('\nSaved figure to: %s\n', fullfile(output_dir, 'AOC_omnibus_condition_interaction.png'));
end

fprintf('\n=== Done ===\n');
