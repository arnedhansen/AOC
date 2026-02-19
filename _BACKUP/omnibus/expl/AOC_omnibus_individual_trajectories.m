%% AOC Omnibus — Individual Trajectories
% Shows individual subject SPL time series (faded lines) with group means overlaid.
% Visualizes within-group variability and individual differences.
%
% Key outputs:
%   Individual trajectories plot with group means

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/individual_trajectories';
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
time_series = linspace(-0.5, 2, 51);
T = numel(time_series) - 1;
t_plot = time_series(2:end);

fs_full = 500;
t_full = -0.5:1/fs_full:2;
t_plot_full = t_full(2:end);
Tf = numel(t_plot_full);

%% Get mean alpha per subject (across all conditions)
subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);
mean_alpha = nan(n_subj, 1);

for s = 1:n_subj
    subj_id = subj_ids{s};
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
end

%% Load gaze data and compute trajectories
reduction_trajectories = {};
amplification_trajectories = {};

for s = 1:n_subj
    subj_id = subj_ids{s};
    
    if ~isfinite(mean_alpha(s))
        continue;
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
    
    % Average across all trials
    subj_mean = nanmean(subj_trials_full, 1);
    
    % Interpolate to binned grid
    subj_traj = interp1(t_plot_full, subj_mean, t_plot, 'linear', NaN);
    
    % Assign to group
    if mean_alpha(s) < 0
        reduction_trajectories{end+1} = subj_traj;
    else
        amplification_trajectories{end+1} = subj_traj;
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Compute group means
n_red = length(reduction_trajectories);
n_amp = length(amplification_trajectories);

if n_red > 0 && n_amp > 0
    % Stack trajectories
    red_mat = nan(n_red, T);
    amp_mat = nan(n_amp, T);
    
    for i = 1:n_red
        red_mat(i, :) = reduction_trajectories{i};
    end
    for i = 1:n_amp
        amp_mat(i, :) = amplification_trajectories{i};
    end
    
    % Group means and SEMs
    grand_red = nanmean(red_mat, 1);
    sem_red = nanstd(red_mat, [], 1) ./ sqrt(sum(isfinite(red_mat), 1));
    grand_amp = nanmean(amp_mat, 1);
    sem_amp = nanstd(amp_mat, [], 1) ./ sqrt(sum(isfinite(amp_mat), 1));
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1400 900])
    
    hold on
    
    % Plot individual trajectories (faded)
    for i = 1:n_red
        plot(t_plot, red_mat(i, :), '-', 'Color', [colors(1,:) 0.5], 'LineWidth', 0.5);
    end
    for i = 1:n_amp
        plot(t_plot, amp_mat(i, :), '-', 'Color', [colors(3,:) 0.5], 'LineWidth', 0.5);
    end
    
    % Plot group means with SEM shading
    addpath('/Users/Arne/Documents/matlabtools/shadedErrorBar');
    sebRed = shadedErrorBar(t_plot, grand_red, sem_red, ...
        'lineProps', {'-','Color',colors(1,:),'LineWidth',4}, 'transparent', true);
    sebAmp = shadedErrorBar(t_plot, grand_amp, sem_amp, ...
        'lineProps', {'-','Color',colors(3,:),'LineWidth',4}, 'transparent', true);
    
    % Add vertical line at time 0
    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6);
    
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('Individual Trajectories: Alpha Reduction (n=%d) vs Amplification (n=%d)', ...
        n_red, n_amp), 'FontSize', fontSize+2, 'FontWeight', 'bold')
    
    legend([sebRed.mainLine, sebAmp.mainLine], ...
        {sprintf('Reduction (n=%d) ± SEM', n_red), ...
         sprintf('Amplification (n=%d) ± SEM', n_amp)}, ...
        'Location', 'best', 'FontSize', fontSize-4)
        set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)

    % Save
    saveas(gcf, fullfile(output_dir, 'AOC_omnibus_individual_trajectories.png'));
    fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_individual_trajectories.png'));
    fprintf('Reduction: %d subjects, Amplification: %d subjects\n', n_red, n_amp);
end

fprintf('\n=== Done ===\n');
