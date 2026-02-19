%% AOC Omnibus — Continuous Relationship: Alpha Power vs SPL
% Scatter plot showing alpha power change (continuous) vs mean scan path length.
% Tests if the relationship is linear using correlation/regression.
%
% Key outputs:
%   Scatter plots (one per condition + overall); correlation coefficients; regression lines

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/continuous';
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

%% Load gaze data and compute mean SPL
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [22, 24, 26];
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);
time_window = t_plot_full >= 0 & t_plot_full <= 2;

subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);

% Storage: subject × condition
alpha_data = nan(n_subj, 3);
spl_data = nan(n_subj, 3);
mean_alpha_all = nan(n_subj, 1);
mean_spl_all = nan(n_subj, 1);

for s = 1:n_subj
    subj_id = subj_ids{s};
    
    % Get alpha power per condition
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    alpha_vals = [];
    for c = 1:3
        cond_rows = subj_rows(strcmp(subj_rows.Condition, conditions{c}), :);
        if ~isempty(cond_rows) && ~isnan(cond_rows.AlphaPower(1)) && cond_rows.AlphaPower(1) ~= 0
            alpha_data(s, c) = cond_rows.AlphaPower(1);
            alpha_vals(end+1) = cond_rows.AlphaPower(1);
        end
    end
    
    if ~isempty(alpha_vals)
        mean_alpha_all(s) = nanmean(alpha_vals);
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
        
        subj_mean = nanmean(subj_trials_full(trial_mask, :), 1);
        spl_data(s, c) = nanmean(subj_mean(time_window));
    end
    
    % Overall mean SPL (all conditions)
    mean_spl_all(s) = nanmean(spl_data(s, :));
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Plot: one figure per condition + overall
close all

for plot_idx = 1:4
    if plot_idx == 4
        % Overall (mean across conditions)
        x_vals = mean_alpha_all;
        y_vals = mean_spl_all;
        plot_title = 'Mean across all conditions';
        file_name = 'AOC_omnibus_continuous_overall.png';
    else
        % Per condition
        x_vals = alpha_data(:, plot_idx);
        y_vals = spl_data(:, plot_idx);
        plot_title = conditions{plot_idx};
        file_name = sprintf('AOC_omnibus_continuous_%s.png', lower(strrep(conditions{plot_idx}, ' ', '')));
    end
    
    % Remove NaNs
    valid = isfinite(x_vals) & isfinite(y_vals);
    x_plot = x_vals(valid);
    y_plot = y_vals(valid);
    
    if length(x_plot) < 3
        continue;
    end
    
    % Correlation
    [r, p_corr] = corrcoef(x_plot, y_plot);
    r_val = r(1, 2);
    p_val = p_corr(1, 2);
    
    % Linear regression
    p_fit = polyfit(x_plot, y_plot, 1);
    x_fit = linspace(min(x_plot), max(x_plot), 100);
    y_fit = polyval(p_fit, x_fit);
    
    % Plot
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1000 800])
    
    hold on
    
    % Scatter points colored by condition (if per-condition plot)
    if plot_idx < 4
        scatter(x_plot, y_plot, 100, colors(plot_idx,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
    else
        % Color by alpha value (gradient)
        scatter(x_plot, y_plot, 100, x_plot, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        colormap(gca, 'cool');
        cb = colorbar;
        ylabel(cb, 'Alpha Power Change', 'FontSize', fontSize-4);
    end
    
    % Regression line
    plot(x_fit, y_fit, 'k--', 'LineWidth', 2.5);
    
    % Add zero lines
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    xlabel('Alpha Power Change from Baseline', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Mean Scan Path Length [px] (0-2s)', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('%s\nr = %.3f, p = %.4f', plot_title, r_val, p_val), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    % Add equation
    text(0.05, 0.95, sprintf('y = %.3f x + %.3f', p_fit(1), p_fit(2)), ...
        'Units', 'normalized', 'FontSize', fontSize-4, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'top');
    
    saveas(gcf, fullfile(output_dir, file_name));
    fprintf('Saved: %s (r=%.3f, p=%.4f, n=%d)\n', file_name, r_val, p_val, length(x_plot));
end

fprintf('\n=== Done ===\n');
