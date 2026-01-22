%% AOC Omnibus — Behavioral Correlates
% Examines relationship between SPL difference (reduction vs amplification) and
% task performance (accuracy, reaction time). Tests if SPL differences relate to behavior.
%
% Key outputs:
%   Scatter plots: SPL difference vs accuracy/RT; correlation coefficients

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/behavioral';
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

%% Load behavioral data
subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);

mean_alpha = nan(n_subj, 1);
mean_spl = nan(n_subj, 1);
mean_acc = nan(n_subj, 1);
mean_rt = nan(n_subj, 1);

% Time grid for SPL
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);
time_window = t_plot_full >= 0 & t_plot_full <= 2;

for s = 1:n_subj
    subj_id = subj_ids{s};
    subj_id_num = str2double(subj_id);
    
    % Get mean alpha power
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    alpha_vals = [];
    for c = 1:height(subj_rows)
        if ~isnan(subj_rows.AlphaPower(c)) && subj_rows.AlphaPower(c) ~= 0
            alpha_vals(end+1) = subj_rows.AlphaPower(c);
        end
    end
    if ~isempty(alpha_vals)
        mean_alpha(s) = nanmean(alpha_vals);
    end
    
    % Load behavioral data
    behav_path = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features', ...
        subj_id, 'behavioral', 'behavioral_matrix_sternberg_subj.mat');
    
    if isfile(behav_path)
        try
            load(behav_path, 'subj_data_behav');
            if exist('subj_data_behav', 'var')
                acc_vals = [subj_data_behav.Accuracy];
                rt_vals = [subj_data_behav.ReactionTime];
                mean_acc(s) = nanmean(acc_vals);
                mean_rt(s) = nanmean(rt_vals);
            end
        catch
        end
    end
    
    % Load gaze data for SPL
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
    
    % Compute mean SPL
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
    
    subj_mean = nanmean(subj_trials_full, 1);
    mean_spl(s) = nanmean(subj_mean(time_window));
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Compute SPL difference (amplification - reduction) for each subject
% Group subjects
reduction_mask = mean_alpha < 0 & isfinite(mean_alpha) & isfinite(mean_spl);
amplification_mask = mean_alpha > 0 & isfinite(mean_alpha) & isfinite(mean_spl);

red_spl = mean_spl(reduction_mask);
amp_spl = mean_spl(amplification_mask);
red_acc = mean_acc(reduction_mask);
amp_acc = mean_acc(amplification_mask);
red_rt = mean_rt(reduction_mask);
amp_rt = mean_rt(amplification_mask);

% Group means
mean_spl_red = nanmean(red_spl);
mean_spl_amp = nanmean(amp_spl);
spl_diff = mean_spl_amp - mean_spl_red;

%% Plot: SPL vs Accuracy
close all

% Plot 1: SPL vs Accuracy
figure
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');

subplot(1,2,1)
hold on
scatter(red_spl, red_acc, 100, colors(1,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);
scatter(amp_spl, amp_acc, 100, colors(3,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);

% Overall correlation
all_spl = [red_spl; amp_spl];
all_acc = [red_acc; amp_acc];
valid = isfinite(all_spl) & isfinite(all_acc);
if sum(valid) >= 3
    [r, p] = corrcoef(all_spl(valid), all_acc(valid));
    r_val = r(1,2);
    p_val = p(1,2);
    
    % Overall regression line
    p_fit = polyfit(all_spl(valid), all_acc(valid), 1);
    x_fit = linspace(min(all_spl(valid)), max(all_spl(valid)), 100);
    y_fit = polyval(p_fit, x_fit);
    h_overall = plot(x_fit, y_fit, 'k--', 'LineWidth', 2);
    
    title(sprintf('SPL vs Accuracy\nr=%.3f, p=%.4f', r_val, p_val), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold');
end

% Group-specific regression lines
valid_red = isfinite(red_spl) & isfinite(red_acc);
if sum(valid_red) >= 2
    p_fit_red = polyfit(red_spl(valid_red), red_acc(valid_red), 1);
    x_fit_red = linspace(min(red_spl(valid_red)), max(red_spl(valid_red)), 100);
    y_fit_red = polyval(p_fit_red, x_fit_red);
    h_red = plot(x_fit_red, y_fit_red, '-', 'Color', colors(1,:), 'LineWidth', 2);
end

valid_amp = isfinite(amp_spl) & isfinite(amp_acc);
if sum(valid_amp) >= 2
    p_fit_amp = polyfit(amp_spl(valid_amp), amp_acc(valid_amp), 1);
    x_fit_amp = linspace(min(amp_spl(valid_amp)), max(amp_spl(valid_amp)), 100);
    y_fit_amp = polyval(p_fit_amp, x_fit_amp);
    h_amp = plot(x_fit_amp, y_fit_amp, '-', 'Color', colors(3,:), 'LineWidth', 2);
end

xlabel('Mean Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
ylabel('Accuracy [%]', 'FontSize', fontSize, 'FontWeight', 'bold')
if exist('h_overall', 'var') && exist('h_red', 'var') && exist('h_amp', 'var')
    legend([h_red, h_amp, h_overall], {'Reduction', 'Amplification', 'Overall'}, ...
        'Location', 'best', 'FontSize', fontSize-4)
else
    legend({'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
end
set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
box on

% Plot 2: SPL vs RT
subplot(1,2,2)
hold on
scatter(red_spl, red_rt, 100, colors(1,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);
scatter(amp_spl, amp_rt, 100, colors(3,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);

all_rt = [red_rt; amp_rt];
valid = isfinite(all_spl) & isfinite(all_rt);
if sum(valid) >= 3
    [r, p] = corrcoef(all_spl(valid), all_rt(valid));
    r_val = r(1,2);
    p_val = p(1,2);
    
    % Overall regression line
    p_fit = polyfit(all_spl(valid), all_rt(valid), 1);
    x_fit = linspace(min(all_spl(valid)), max(all_spl(valid)), 100);
    y_fit = polyval(p_fit, x_fit);
    h_overall2 = plot(x_fit, y_fit, 'k--', 'LineWidth', 2);
    
    title(sprintf('SPL vs Reaction Time\nr=%.3f, p=%.4f', r_val, p_val), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold');
end

% Group-specific regression lines
valid_red = isfinite(red_spl) & isfinite(red_rt);
if sum(valid_red) >= 2
    p_fit_red = polyfit(red_spl(valid_red), red_rt(valid_red), 1);
    x_fit_red = linspace(min(red_spl(valid_red)), max(red_spl(valid_red)), 100);
    y_fit_red = polyval(p_fit_red, x_fit_red);
    h_red2 = plot(x_fit_red, y_fit_red, '-', 'Color', colors(1,:), 'LineWidth', 2);
end

valid_amp = isfinite(amp_spl) & isfinite(amp_rt);
if sum(valid_amp) >= 2
    p_fit_amp = polyfit(amp_spl(valid_amp), amp_rt(valid_amp), 1);
    x_fit_amp = linspace(min(amp_spl(valid_amp)), max(amp_spl(valid_amp)), 100);
    y_fit_amp = polyval(p_fit_amp, x_fit_amp);
    h_amp2 = plot(x_fit_amp, y_fit_amp, '-', 'Color', colors(3,:), 'LineWidth', 2);
end

xlabel('Mean Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
ylabel('Reaction Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
if exist('h_overall2', 'var') && exist('h_red2', 'var') && exist('h_amp2', 'var')
    legend([h_red2, h_amp2, h_overall2], {'Reduction', 'Amplification', 'Overall'}, ...
        'Location', 'best', 'FontSize', fontSize-4)
else
    legend({'Reduction', 'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
end
set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
box on

sgtitle('Behavioral Correlates of Scan Path Length', ...
    'FontSize', fontSize+4, 'FontWeight', 'bold')

saveas(gcf, fullfile(output_dir, 'AOC_omnibus_behavioral_correlates.png'));
fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_behavioral_correlates.png'));

% Print summary
fprintf('\n=== Behavioral Summary ===\n');
fprintf('Reduction group: n=%d, SPL=%.3f, Acc=%.1f%%, RT=%.3fs\n', ...
    length(red_spl), nanmean(red_spl), nanmean(red_acc), nanmean(red_rt));
fprintf('Amplification group: n=%d, SPL=%.3f, Acc=%.1f%%, RT=%.3fs\n', ...
    length(amp_spl), nanmean(amp_spl), nanmean(amp_acc), nanmean(amp_rt));

fprintf('\n=== Done ===\n');
