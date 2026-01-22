%% AOC Omnibus — Subject-Level Effect Sizes
% Computes individual Cohen's d values for each subject (comparing their own SPL
% across conditions or time windows). Tests if effect is consistent across subjects.
%
% Key outputs:
%   Distribution of individual effect sizes; histogram; summary statistics

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/subject_effects';
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

% Storage: subject-level SPL per condition
spl_per_cond = nan(n_subj, 3);
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
        spl_per_cond(s, c) = nanmean(subj_mean(time_window));
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Compute individual effect sizes
% For each subject: compare SPL between conditions where alpha < 0 vs > 0
individual_d = nan(n_subj, 1);

for s = 1:n_subj
    % Get conditions with reduction vs amplification
    red_mask = alpha_per_cond(s, :) < 0 & isfinite(spl_per_cond(s, :));
    amp_mask = alpha_per_cond(s, :) > 0 & isfinite(spl_per_cond(s, :));
    
    red_spl = spl_per_cond(s, red_mask);
    amp_spl = spl_per_cond(s, amp_mask);
    
    if length(red_spl) >= 1 && length(amp_spl) >= 1
        % Compute Cohen's d
        pooled_std = sqrt(((length(red_spl)-1)*nanvar(red_spl) + ...
            (length(amp_spl)-1)*nanvar(amp_spl)) / ...
            (length(red_spl) + length(amp_spl) - 2));
        if pooled_std > 0
            individual_d(s) = (nanmean(amp_spl) - nanmean(red_spl)) / pooled_std;
        end
    end
end

% Remove NaNs
valid_d = individual_d(isfinite(individual_d));

if length(valid_d) >= 3
    %% Summary statistics
    fprintf('\n=== Subject-Level Effect Sizes ===\n');
    fprintf('N subjects: %d\n', length(valid_d));
    fprintf('Mean d: %.3f\n', nanmean(valid_d));
    fprintf('SD d: %.3f\n', nanstd(valid_d));
    fprintf('Median d: %.3f\n', nanmedian(valid_d));
    fprintf('Range: [%.3f, %.3f]\n', min(valid_d), max(valid_d));
    
    % Test if mean differs from zero
    [~, p_test] = ttest(valid_d);
    fprintf('One-sample t-test vs 0: p = %.4f\n', p_test);
    
    % Proportion with positive effect
    prop_positive = sum(valid_d > 0) / length(valid_d);
    fprintf('Proportion with positive effect: %.1f%%\n', prop_positive * 100);
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1200 800])
    
    subplot(2,2,1)
    histogram(valid_d, 20, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 1.5)
    xline(0, '--r', 'LineWidth', 2)
    xline(nanmean(valid_d), '--b', 'LineWidth', 2)
    xlabel('Cohen''s d', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Frequency', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Distribution of Individual Effect Sizes', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    legend({'Distribution', 'Zero', sprintf('Mean (%.3f)', nanmean(valid_d))}, ...
        'Location', 'best', 'FontSize', fontSize-6)
    
    subplot(2,2,2)
    boxplot(valid_d, 'Orientation', 'horizontal', 'Colors', 'k', 'Widths', 0.3)
    xline(0, '--r', 'LineWidth', 2)
    xlabel('Cohen''s d', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Box Plot', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    subplot(2,2,3)
    scatter(1:length(valid_d), sort(valid_d), 80, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
    xline(0, '--r', 'LineWidth', 1)
    yline(0, '--r', 'LineWidth', 1)
    xlabel('Subject (sorted)', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Cohen''s d', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Sorted Individual Effect Sizes', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    subplot(2,2,4)
    qqplot(valid_d)
    xlabel('Theoretical Quantiles', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Sample Quantiles (Cohen''s d)', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Q-Q Plot (Normality Check)', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'FontSize', fontSize-2)
    box on
    
    sgtitle(sprintf('Subject-Level Effect Sizes (n=%d)', length(valid_d)), ...
        'FontSize', fontSize+4, 'FontWeight', 'bold')
    
    saveas(gcf, fullfile(output_dir, 'AOC_omnibus_subject_effect_sizes.png'));
    fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_subject_effect_sizes.png'));
end

fprintf('\n=== Done ===\n');
