%% AOC Omnibus — Forest Plot of Effect Sizes
% Creates forest plot showing Cohen's d with confidence intervals for each condition
% and overall. Meta-analytic summary visualization.
%
% Key outputs:
%   Forest plot with effect sizes and confidence intervals per condition

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/forest_plot';
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

%% Time grid
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);
time_window = t_plot_full >= 0 & t_plot_full <= 2;

conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [22, 24, 26];

subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);

% Storage: condition × group
spl_per_cond_group = nan(3, 2, n_subj); % conditions × groups × subjects
alpha_per_cond = nan(n_subj, 3);

for s = 1:n_subj
    subj_id = subj_ids{s};
    
    % Get alpha power per condition
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    for c = 1:3
        cond_rows = subj_rows(strcmp(subj_rows.Condition, conditions{c}), :);
        if ~isempty(cond_rows) && ~isnan(cond_rows.AlphaPower(1)) && cond_rows.AlphaPower(1) ~= 0
            alpha_per_cond(s, c) = cond_rows.AlphaPower(1);
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
        mean_spl = nanmean(subj_mean(time_window));
        
        if isfinite(mean_spl) && isfinite(alpha_per_cond(s, c))
            if alpha_per_cond(s, c) < 0
                spl_per_cond_group(c, 1, s) = mean_spl; % Reduction
            else
                spl_per_cond_group(c, 2, s) = mean_spl; % Amplification
            end
        end
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Compute effect sizes and confidence intervals
cohens_d = nan(4, 1); % 3 conditions + overall
ci_lower = nan(4, 1);
ci_upper = nan(4, 1);
n_red = nan(4, 1);
n_amp = nan(4, 1);

for c = 1:4
    if c == 4
        % Overall: use mean across conditions
        red_vals = [];
        amp_vals = [];
        for cond = 1:3
            red_cond = spl_per_cond_group(cond, 1, :);
            amp_cond = spl_per_cond_group(cond, 2, :);
            % Convert to column vectors and extract finite values
            red_cond = red_cond(:);
            amp_cond = amp_cond(:);
            red_vals = [red_vals; red_cond(isfinite(red_cond))];
            amp_vals = [amp_vals; amp_cond(isfinite(amp_cond))];
        end
    else
        red_vals = spl_per_cond_group(c, 1, :);
        amp_vals = spl_per_cond_group(c, 2, :);
        % Convert to column vectors
        red_vals = red_vals(:);
        amp_vals = amp_vals(:);
        red_vals = red_vals(isfinite(red_vals));
        amp_vals = amp_vals(isfinite(amp_vals));
    end
    
    if length(red_vals) >= 2 && length(amp_vals) >= 2
        n_red(c) = length(red_vals);
        n_amp(c) = length(amp_vals);
        
        % Cohen's d
        pooled_std = sqrt(((length(red_vals)-1)*nanvar(red_vals) + ...
            (length(amp_vals)-1)*nanvar(amp_vals)) / ...
            (length(red_vals) + length(amp_vals) - 2));
        
        if pooled_std > 0
            cohens_d(c) = (nanmean(amp_vals) - nanmean(red_vals)) / pooled_std;
            
            % Confidence interval (using Hedges' g correction for small samples)
            n_eff = (n_red(c) * n_amp(c)) / (n_red(c) + n_amp(c));
            se_d = sqrt((n_red(c) + n_amp(c)) / (n_red(c) * n_amp(c)) + cohens_d(c)^2 / (2 * (n_red(c) + n_amp(c))));
            t_crit = tinv(0.975, n_red(c) + n_amp(c) - 2);
            ci_lower(c) = cohens_d(c) - t_crit * se_d;
            ci_upper(c) = cohens_d(c) + t_crit * se_d;
        end
    end
end

%% Plot forest plot
close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 1200 800])

hold on

% Y positions
y_pos = [4, 3, 2, 1];
labels = [conditions, {'Overall (mean)'}];

% Plot effect sizes and CIs
for i = 1:4
    if isfinite(cohens_d(i))
        % CI line
        plot([ci_lower(i), ci_upper(i)], [y_pos(i), y_pos(i)], 'k-', 'LineWidth', 2)
        
        % Point estimate
        scatter(cohens_d(i), y_pos(i), 200, colors(min(i,3),:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2)
        
        % Label with effect size and CI
        text(ci_upper(i) + 0.1, y_pos(i), ...
            sprintf('d=%.2f [%.2f, %.2f] (n=%d/%d)', ...
            cohens_d(i), ci_lower(i), ci_upper(i), n_red(i), n_amp(i)), ...
            'FontSize', fontSize-4, 'VerticalAlignment', 'middle')
    end
end

% Vertical reference line at d=0
xline(0, '--k', 'LineWidth', 2)

% Effect size reference lines
xline(0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
xline(0.5, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
xline(-0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
xline(-0.5, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1)

ylim([0.5 4.5])
xlim_range = [min([ci_lower; -0.5]) - 0.2, max([ci_upper; 0.5]) + 0.5];
xlim(xlim_range)

% Set YTick in increasing order (MATLAB requirement) with reversed labels to match plot positions
set(gca, 'YTick', sort(y_pos), 'YTickLabel', labels(end:-1:1), 'FontSize', fontSize-2)
xlabel('Cohen''s d (95% CI)', 'FontSize', fontSize, 'FontWeight', 'bold')
title('Forest Plot: Effect Sizes by Condition', 'FontSize', fontSize+2, 'FontWeight', 'bold')
set(gca, 'GridAlpha', 0.2, 'XGrid', 'on', 'YGrid', 'off')
box on

saveas(gcf, fullfile(output_dir, 'AOC_omnibus_forest_plot.png'));
fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_forest_plot.png'));

% Print summary
fprintf('\n=== Effect Size Summary ===\n');
for i = 1:4
    if isfinite(cohens_d(i))
        fprintf('%s: d=%.3f [%.3f, %.3f], n_red=%d, n_amp=%d\n', ...
            labels{i}, cohens_d(i), ci_lower(i), ci_upper(i), n_red(i), n_amp(i));
    end
end

fprintf('\n=== Done ===\n');
