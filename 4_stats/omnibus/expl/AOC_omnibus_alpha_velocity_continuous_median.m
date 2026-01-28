B%% AOC Omnibus — Alpha Power × Eye Velocity: Continuous & Median Split Analysis
% Analyzes the relationship between alpha power and eye velocity for both
% n-back and Sternberg tasks using two approaches:
% 1. Continuous correlation: Time-resolved correlation between alpha power and velocity
% 2. Median split: High vs Low alpha reduction groups
%
% Key outputs:
%   Time-resolved correlation plots; Median split comparisons; Statistical tests

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 15;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/continuous_median';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Time grid
fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);
time_series = linspace(-0.5, 2, 51);
T = numel(time_series) - 1;
t_plot = time_series(2:end);

% Velocity computation parameters
polyOrd = 3;  % Savitzky-Golay polynomial order for velocity
velZthr = 4;  % Z-score threshold for velocity outliers

%% Load CSV data
csv_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_omnibus_raincloud_data.csv';
raincloud_data = readtable(csv_path);

% Convert ID to string
if isnumeric(raincloud_data.ID)
    raincloud_data.ID = cellstr(num2str(raincloud_data.ID, '%d'));
elseif iscell(raincloud_data.ID)
    raincloud_data.ID = cellfun(@(x) num2str(x), raincloud_data.ID, 'UniformOutput', false);
elseif isstring(raincloud_data.ID)
    raincloud_data.ID = cellstr(raincloud_data.ID);
end

%% Process each task
tasks = {'nback', 'sternberg'};
task_configs = struct();
task_configs.nback = struct('name', 'nback', 'conditions', {'1-back', '2-back', '3-back'}, ...
    'cond_codes', [21, 22, 23], 'gaze_file', 'gaze_series_nback_trials.mat');
task_configs.sternberg = struct('name', 'sternberg', 'conditions', {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'cond_codes', [22, 24, 26], 'gaze_file', 'gaze_series_sternberg_trials.mat');

for task_idx = 1:length(tasks)
    task_name = tasks{task_idx};
    config = task_configs.(task_name);
    
    fprintf('\n=== Processing %s task ===\n', task_name);
    
    % Filter data for this task
    task_data = raincloud_data(strcmpi(raincloud_data.Task, task_name), :);
    subj_ids = unique(task_data.ID);
    n_subj = length(subj_ids);
    
    fprintf('Found %d subjects for %s\n', n_subj, task_name);
    
    % Storage
    mean_alpha = nan(n_subj, 1);
    vel_trajectories = nan(n_subj, Tf);
    
    %% Load alpha power and gaze data
    valid_subj_count = 0;
    mean_alpha_valid = [];
    vel_trajectories_valid = [];
    
    fprintf('Loading data for %d subjects...\n', n_subj);
    
    % Track statistics
    n_missing_alpha = 0;
    n_missing_gaze = 0;
    n_load_errors = 0;
    n_no_trials = 0;
    n_insufficient_data = 0;
    
    for s = 1:n_subj
        subj_id = subj_ids{s};
        
        % Get mean alpha power across all conditions
        subj_rows = task_data(strcmp(task_data.ID, subj_id), :);
        alpha_vals = [];
        for c = 1:height(subj_rows)
            if ~isnan(subj_rows.AlphaPower(c)) && subj_rows.AlphaPower(c) ~= 0
                alpha_vals(end+1) = subj_rows.AlphaPower(c);
            end
        end
        if isempty(alpha_vals)
            n_missing_alpha = n_missing_alpha + 1;
            continue;
        end
        subj_mean_alpha = nanmean(alpha_vals);
        
        % Load gaze data - construct path explicitly
        base_path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features';
        % Get filename from config (should be just the filename, not a path)
        if isfield(config, 'gaze_file')
            gaze_filename = config.gaze_file;
        else
            % Fallback based on task
            if strcmp(task_name, 'nback')
                gaze_filename = 'gaze_series_nback_trials.mat';
            else
                gaze_filename = 'gaze_series_sternberg_trials.mat';
            end
        end
        % Ensure we only have the filename (no path)
        [~, name, ext] = fileparts(gaze_filename);
        gaze_filename_clean = [name, ext];
        datapath_gaze = fullfile(base_path, subj_id, 'gaze', gaze_filename_clean);
        
        if ~isfile(datapath_gaze)
            n_missing_gaze = n_missing_gaze + 1;
            continue;
        end
        
        % Try to load with better error handling
        try
            % Suppress "Error closing file" warnings - these are usually harmless
            lastwarn(''); % Clear previous warnings
            warning('off', 'MATLAB:load:errorClosingFile');
            load(datapath_gaze, 'gaze_x', 'gaze_y', 'trialinfo');
            warning('on', 'MATLAB:load:errorClosingFile');
            [warnMsg, ~] = lastwarn;
            if contains(warnMsg, 'Error closing file')
                n_load_errors = n_load_errors + 1;
                % Continue anyway - data was loaded successfully
            end
        catch ME
            % Count all errors but don't print warnings for file closing issues
            n_load_errors = n_load_errors + 1;
            continue;
        end
        
        % Check if variables exist
        if ~exist('gaze_x', 'var') || ~exist('gaze_y', 'var')
            n_load_errors = n_load_errors + 1;
            continue;
        end
        
        if isempty(gaze_x) || isempty(gaze_y)
            n_no_trials = n_no_trials + 1;
            continue;
        end
        
        % Compute velocity time series (all conditions)
        n_trials = size(gaze_x, 2);
        subj_trials_full = nan(n_trials, Tf);
        valid_trials = 0;
        
        for trl = 1:n_trials
            try
                X = gaze_x{trl, 1};
                Y = gaze_y{trl, 1};
            catch
                continue;
            end
            
            if isempty(X) || isempty(Y) || length(X) ~= length(Y)
                continue;
            end
            
            % Extract -0.5 to 2s portion from gaze data (full window)
            % Gaze data is typically -0.5 to 2.0s at 500Hz = 1251 samples
            if length(X) == 1251  % Standard case: -0.5 to 2.0s
                % Use all data - no trimming needed
            elseif length(X) > 1250
                % If longer, take first 1251 samples (assuming it starts at -0.5s)
                X = X(1:1251);
                Y = Y(1:1251);
            elseif length(X) < 1251
                % If shorter, pad or interpolate to match expected length
                if length(X) >= 1000
                    % Pad with last value to reach 1251
                    X = [X, repmat(X(end), 1, 1251 - length(X))];
                    Y = [Y, repmat(Y(end), 1, 1251 - length(Y))];
                else
                    % Too short, skip this trial
                    continue;
                end
            end
            
            % Create time vector for this trial
            gaze_time = linspace(-0.5, 2, length(X));
            
            % Compute eye velocity
            [vx, vy] = compute_velocity_sg(X, Y, fs_full, polyOrd);
            
            % Handle outliers
            zvx = (vx - nanmean(vx)) / (nanstd(vx) + eps);
            zvy = (vy - nanmean(vy)) / (nanstd(vy) + eps);
            bad = abs(zvx) > velZthr | abs(zvy) > velZthr;
            if any(bad)
                vx(bad) = NaN;
                vy(bad) = NaN;
                vx = fillmissing(vx, 'linear', 'EndValues', 'nearest');
                vy = fillmissing(vy, 'linear', 'EndValues', 'nearest');
            end
            
            % Total velocity magnitude
            vel = hypot(vx, vy);
            
            % Interpolate to common time grid
            try
                interp_vals = interp1(gaze_time, vel, t_plot_full, 'linear', NaN);
                if sum(isfinite(interp_vals)) > 100  % At least 100 valid points (reasonable threshold)
                    subj_trials_full(trl, :) = interp_vals;
                    valid_trials = valid_trials + 1;
                end
            catch ME
                % Skip this trial
            end
        end
        
        % Average across all trials
        if valid_trials > 0
            subj_mean = nanmean(subj_trials_full, 1);
            % Check if we have sufficient finite values in the mean
            n_valid_points = sum(isfinite(subj_mean));
            if n_valid_points > 100  % At least 100 time points with data
                valid_subj_count = valid_subj_count + 1;
                mean_alpha_valid(valid_subj_count) = subj_mean_alpha;
                vel_trajectories_valid(valid_subj_count, :) = subj_mean;
            else
                n_insufficient_data = n_insufficient_data + 1;
            end
        else
            n_no_trials = n_no_trials + 1;
        end
        
        clear gaze_x gaze_y trialinfo;
    end
    
    % Assign to output variables
    mean_alpha = mean_alpha_valid(:);
    vel_trajectories = vel_trajectories_valid;
    
    fprintf('Valid subjects with both alpha and gaze data: %d\n', valid_subj_count);
    fprintf('  Statistics: Missing alpha=%d, Missing gaze=%d, Load errors=%d, No trials=%d, Insufficient data=%d\n', ...
        n_missing_alpha, n_missing_gaze, n_load_errors, n_no_trials, n_insufficient_data);
    
    if valid_subj_count < 3
        warning('Not enough valid subjects for %s (%d < 3). Skipping.', task_name, valid_subj_count);
        continue;
    end
    
    % Assign to final variables
    mean_alpha = mean_alpha_valid(:);
    vel_trajectories = vel_trajectories_valid;
    
    %% Analysis 1: Time-resolved correlation
    fprintf('\n--- Computing time-resolved correlations ---\n');
    r_vals = nan(1, Tf);
    p_vals = nan(1, Tf);
    
    for t = 1:Tf
        vel_vals = vel_trajectories(:, t);
        valid = isfinite(mean_alpha) & isfinite(vel_vals);
        if sum(valid) >= 3
            [r, p] = corrcoef(mean_alpha(valid), vel_vals(valid));
            r_vals(t) = r(1, 2);
            p_vals(t) = p(1, 2);
        end
    end
    
    % FDR correction
    toi = t_plot_full >= 0 & t_plot_full <= 2;
    p_sub = p_vals(toi & isfinite(p_vals));
    if ~isempty(p_sub)
        q_sub = mafdr(p_sub, 'BHFDR', true);
        sig_mask = false(size(p_vals));
        sig_mask(toi) = q_sub < 0.05;
    else
        sig_mask = false(size(p_vals));
    end
    
    % Interpolate to binned grid for plotting
    r_vals_binned = interp1(t_plot_full, r_vals, t_plot, 'linear', NaN);
    p_vals_binned = interp1(t_plot_full, p_vals, t_plot, 'linear', NaN);
    sig_mask_binned = false(1, T);
    for t = 1:T
        [~, idx] = min(abs(t_plot_full - t_plot(t)));
        if idx <= length(sig_mask)
            sig_mask_binned(t) = sig_mask(idx);
        end
    end
    
    %% Analysis 2: Median split
    fprintf('\n--- Performing median split ---\n');
    
    % Compute median alpha power (more negative = more reduction)
    median_alpha = nanmedian(mean_alpha);
    high_reduction_mask = mean_alpha <= median_alpha;
    low_reduction_mask = mean_alpha > median_alpha;
    
    n_high = sum(high_reduction_mask);
    n_low = sum(low_reduction_mask);
    
    fprintf('High reduction: n=%d (alpha <= %.4f)\n', n_high, median_alpha);
    fprintf('Low reduction: n=%d (alpha > %.4f)\n', n_low, median_alpha);
    
    % Compute group means
    vel_high = vel_trajectories(high_reduction_mask, :);
    vel_low = vel_trajectories(low_reduction_mask, :);
    
    grand_high = nanmean(vel_high, 1);
    grand_low = nanmean(vel_low, 1);
    sem_high = nanstd(vel_high, [], 1) ./ sqrt(sum(isfinite(vel_high), 1));
    sem_low = nanstd(vel_low, [], 1) ./ sqrt(sum(isfinite(vel_low), 1));
    
    % Interpolate to binned grid
    grand_high_binned = interp1(t_plot_full, grand_high, t_plot, 'linear', NaN);
    grand_low_binned = interp1(t_plot_full, grand_low, t_plot, 'linear', NaN);
    sem_high_binned = interp1(t_plot_full, sem_high, t_plot, 'linear', NaN);
    sem_low_binned = interp1(t_plot_full, sem_low, t_plot, 'linear', NaN);
    
    % Statistical comparison at each time point
    cohens_d_binned = nan(1, T);
    p_comp_binned = nan(1, T);
    
    for t = 1:T
        [~, idx] = min(abs(t_plot_full - t_plot(t)));
        high_vals = vel_high(:, idx);
        low_vals = vel_low(:, idx);
        high_vals = high_vals(isfinite(high_vals));
        low_vals = low_vals(isfinite(low_vals));
        
        if length(high_vals) >= 2 && length(low_vals) >= 2
            [~, p, ~, stats] = ttest2(high_vals, low_vals);
            p_comp_binned(t) = p;
            
            pooled_std = sqrt(((length(high_vals)-1)*nanvar(high_vals) + ...
                (length(low_vals)-1)*nanvar(low_vals)) / ...
                (length(high_vals) + length(low_vals) - 2));
            if pooled_std > 0
                cohens_d_binned(t) = (nanmean(high_vals) - nanmean(low_vals)) / pooled_std;
            end
        end
    end
    
    %% Plot: Combined figure
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1600 1200])
    
    % Top row: Time-resolved correlation
    subplot(3,2,1:2)
    hold on
    plot(t_plot, r_vals_binned, 'k-', 'LineWidth', 3)
    xline(0, '--k', 'LineWidth', 1.5)
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Correlation (r)', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('%s: Time-Resolved Correlation (Alpha Power × Eye Velocity)', upper(task_name)), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    % Middle row: Median split comparison
    subplot(3,2,3:4)
    hold on
    
    % Use shadedErrorBar if available, otherwise use fill
    if exist('shadedErrorBar', 'file')
        sebHigh = shadedErrorBar(t_plot, grand_high_binned, sem_high_binned, ...
            'lineProps', {'-','Color',colors(1,:),'LineWidth',3}, 'transparent', true);
        sebLow = shadedErrorBar(t_plot, grand_low_binned, sem_low_binned, ...
            'lineProps', {'-','Color',colors(3,:),'LineWidth',3}, 'transparent', true);
    else
        % Fallback: plot with fill for error bars
        fill([t_plot, fliplr(t_plot)], [grand_high_binned + sem_high_binned, fliplr(grand_high_binned - sem_high_binned)], ...
            colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        sebHigh.mainLine = plot(t_plot, grand_high_binned, '-', 'Color', colors(1,:), 'LineWidth', 3);
        fill([t_plot, fliplr(t_plot)], [grand_low_binned + sem_low_binned, fliplr(grand_low_binned - sem_low_binned)], ...
            colors(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        sebLow.mainLine = plot(t_plot, grand_low_binned, '-', 'Color', colors(3,:), 'LineWidth', 3);
    end
    
    xline(0, '--k', 'LineWidth', 1.5);
    ylim([0 1000])
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Eye Velocity [px/s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('%s: Median Split (High vs Low Alpha Reduction)', upper(task_name)), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    legend([sebHigh.mainLine, sebLow.mainLine], ...
        {sprintf('High Reduction (n=%d) ± SEM', n_high), ...
         sprintf('Low Reduction (n=%d) ± SEM', n_low)}, ...
        'Location', 'best', 'FontSize', fontSize-4)
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    % Bottom row: Effect size (Cohen's d)
    subplot(3,2,5:6)
    hold on
    
    % Shade significant regions
    d_sig = diff([0, p_comp_binned < 0.05, 0]);
    on_sig = find(d_sig == 1);
    off_sig = find(d_sig == -1) - 1;
    
    ylim_range = [-max(abs(cohens_d_binned(isfinite(cohens_d_binned))))*1.2 ...
        max(abs(cohens_d_binned(isfinite(cohens_d_binned))))*1.2];
    if any(isfinite(cohens_d_binned))
        ylim(ylim_range);
    end
    
    for k = 1:numel(on_sig)
        if on_sig(k) <= length(t_plot) && off_sig(k) <= length(t_plot)
            x0 = t_plot(on_sig(k));
            x1 = t_plot(off_sig(k));
            patch([x0 x1 x1 x0], [ylim_range(1) ylim_range(1) ylim_range(2) ylim_range(2)], ...
                [1 1 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        end
    end
    
    plot(t_plot, cohens_d_binned, 'k-', 'LineWidth', 3);
    xline(0, '--k', 'LineWidth', 1.5);
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    yline(-0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Cohen''s d', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('%s: Effect Size (High - Low Reduction)', upper(task_name)), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    % Save figure
    fig_name = sprintf('AOC_omnibus_alpha_velocity_%s.png', task_name);
    fig_path = fullfile(output_dir, fig_name);
    saveas(gcf, fig_path);
    fprintf('\nSaved figure: %s\n', fig_name);
    
    % Print summary statistics
    fprintf('\n=== Summary Statistics (%s) ===\n', task_name);
    fprintf('Mean alpha power: %.4f (SD=%.4f, range=[%.4f, %.4f])\n', ...
        nanmean(mean_alpha), nanstd(mean_alpha), min(mean_alpha), max(mean_alpha));
    fprintf('Median alpha power: %.4f\n', median_alpha);
    fprintf('Max correlation: r=%.3f at t=%.2fs\n', ...
        nanmax(r_vals_binned), t_plot(r_vals_binned == nanmax(r_vals_binned)));
    fprintf('Max effect size: d=%.3f at t=%.2fs\n', ...
        nanmax(abs(cohens_d_binned)), t_plot(abs(cohens_d_binned) == nanmax(abs(cohens_d_binned))));
end

fprintf('\n=== Done ===\n');
fprintf('All figures saved to: %s\n', output_dir);

%% Helper function: Compute velocity using Savitzky-Golay filter
function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
% Convert to double precision (required by sgolayfilt)
X = double(X);
Y = double(Y);

Ts = 1/fs;
L = numel(X);

% Determine frame length (must be odd)
framelen = min(21, L);
if mod(framelen, 2) == 0
    framelen = framelen - 1;
end

minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0
    minLegal = minLegal + 1;
end
if framelen < minLegal
    framelen = minLegal;
end
if framelen > L
    framelen = L - mod(L, 2) + 1;
end

useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / (Ts^1)) * G(:, 2)';  % 1st-derivative kernel
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end
