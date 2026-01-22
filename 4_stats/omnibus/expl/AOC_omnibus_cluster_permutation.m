%% AOC Omnibus — Cluster-Based Permutation Test
% Performs cluster-based permutation test to identify significant time clusters
% where reduction and amplification groups differ in SPL.
%
% Key outputs:
%   Cluster statistics; significant time clusters; visualization

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/cluster_perm';
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

%% Get mean alpha per subject
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

%% Load gaze data
reduction_data = [];
amplification_data = [];

for s = 1:n_subj
    subj_id = subj_ids{s};
    
    if ~isfinite(mean_alpha(s))
        continue;
    end
    
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
    
    % Compute SPL time series
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
    subj_traj = interp1(t_plot_full, subj_mean, t_plot, 'linear', NaN);
    
    if mean_alpha(s) < 0
        reduction_data(end+1, :) = subj_traj;
    else
        amplification_data(end+1, :) = subj_traj;
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Cluster-based permutation test
n_perm = 1000;
alpha_thresh = 0.05; % For cluster formation
cluster_alpha = 0.05; % For cluster significance

fprintf('\n=== Cluster-Based Permutation Test ===\n');
fprintf('Reduction: n=%d, Amplification: n=%d\n', size(reduction_data, 1), size(amplification_data, 1));

if size(reduction_data, 1) >= 2 && size(amplification_data, 1) >= 2
    % Compute t-statistics at each time point
    t_stats = nan(1, T);
    p_vals = nan(1, T);
    
    for t = 1:T
        red_vals = reduction_data(:, t);
        amp_vals = amplification_data(:, t);
        red_vals = red_vals(isfinite(red_vals));
        amp_vals = amp_vals(isfinite(amp_vals));
        
        if length(red_vals) >= 2 && length(amp_vals) >= 2
            [~, p, ~, stats] = ttest2(red_vals, amp_vals);
            t_stats(t) = stats.tstat;
            p_vals(t) = p;
        end
    end
    
    % Find clusters (contiguous time points with p < alpha_thresh)
    sig_mask = p_vals < alpha_thresh;
    clusters = [];
    cluster_stats = [];
    
    in_cluster = false;
    cluster_start = 0;
    for t = 1:T
        if sig_mask(t) && ~in_cluster
            cluster_start = t;
            in_cluster = true;
        elseif ~sig_mask(t) && in_cluster
            cluster_end = t - 1;
            clusters(end+1, :) = [cluster_start, cluster_end];
            in_cluster = false;
        end
    end
    if in_cluster
        clusters(end+1, :) = [cluster_start, T];
    end
    
    % Compute cluster statistics (sum of t-statistics)
    for c = 1:size(clusters, 1)
        cluster_stats(c) = sum(abs(t_stats(clusters(c,1):clusters(c,2))));
    end
    
    % Permutation test
    all_data = [reduction_data; amplification_data];
    n_red = size(reduction_data, 1);
    n_total = size(all_data, 1);
    
    perm_cluster_stats = nan(n_perm, max(size(clusters, 1), 1));
    
    fprintf('Running %d permutations...\n', n_perm);
    for perm = 1:n_perm
        % Random permutation
        perm_idx = randperm(n_total);
        perm_red = all_data(perm_idx(1:n_red), :);
        perm_amp = all_data(perm_idx(n_red+1:end), :);
        
        % Compute t-stats
        perm_t_stats = nan(1, T);
        for t = 1:T
            red_vals = perm_red(:, t);
            amp_vals = perm_amp(:, t);
            red_vals = red_vals(isfinite(red_vals));
            amp_vals = amp_vals(isfinite(amp_vals));
            
            if length(red_vals) >= 2 && length(amp_vals) >= 2
                [~, ~, ~, stats] = ttest2(red_vals, amp_vals);
                perm_t_stats(t) = stats.tstat;
            end
        end
        
        % Find clusters in permuted data
        perm_p_vals = 2 * (1 - tcdf(abs(perm_t_stats), n_total - 2));
        perm_sig_mask = perm_p_vals < alpha_thresh;
        
        perm_clusters = [];
        in_cluster = false;
        cluster_start = 0;
        for t = 1:T
            if perm_sig_mask(t) && ~in_cluster
                cluster_start = t;
                in_cluster = true;
            elseif ~perm_sig_mask(t) && in_cluster
                cluster_end = t - 1;
                perm_clusters(end+1, :) = [cluster_start, cluster_end];
                in_cluster = false;
            end
        end
        if in_cluster
            perm_clusters(end+1, :) = [cluster_start, T];
        end
        
        % Store max cluster statistic
        if ~isempty(perm_clusters)
            perm_cluster_stats(perm) = max(sum(abs(perm_t_stats(perm_clusters(:,1):perm_clusters(:,2))), 2));
        else
            perm_cluster_stats(perm) = 0;
        end
    end
    
    % Compute p-values for clusters
    cluster_p_vals = nan(size(clusters, 1), 1);
    for c = 1:size(clusters, 1)
        cluster_p_vals(c) = sum(perm_cluster_stats >= cluster_stats(c)) / n_perm;
    end
    
    % Report results
    fprintf('\n=== Cluster Results ===\n');
    for c = 1:size(clusters, 1)
        t_start = t_plot(clusters(c, 1));
        t_end = t_plot(clusters(c, 2));
        fprintf('Cluster %d: [%.2f, %.2f] s, statistic=%.2f, p=%.4f', ...
            c, t_start, t_end, cluster_stats(c), cluster_p_vals(c));
        if cluster_p_vals(c) < cluster_alpha
            fprintf(' * SIGNIFICANT\n');
        else
            fprintf('\n');
        end
    end
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1400 900])
    
    subplot(2,1,1)
    hold on
    plot(t_plot, t_stats, 'k-', 'LineWidth', 2.5)
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    xline(0, '--k', 'LineWidth', 1.5)
    
    % Highlight significant clusters
    for c = 1:size(clusters, 1)
        if cluster_p_vals(c) < cluster_alpha
            x0 = t_plot(clusters(c, 1));
            x1 = t_plot(clusters(c, 2));
            y_max = max(abs(t_stats)) * 1.1;
            patch([x0 x1 x1 x0], [-y_max -y_max y_max y_max], ...
                [1 1 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end
    end
    
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('t-statistic', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Cluster-Based Permutation Test', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    subplot(2,1,2)
    hold on
    plot(t_plot, -log10(p_vals), 'k-', 'LineWidth', 2.5)
    yline(-log10(alpha_thresh), '--r', 'LineWidth', 1.5, 'Label', 'Threshold')
    xline(0, '--k', 'LineWidth', 1.5)
    
    % Highlight significant clusters
    for c = 1:size(clusters, 1)
        if cluster_p_vals(c) < cluster_alpha
            x0 = t_plot(clusters(c, 1));
            x1 = t_plot(clusters(c, 2));
            y_max = max(-log10(p_vals(isfinite(p_vals)))) * 1.1;
            patch([x0 x1 x1 x0], [0 0 y_max y_max], ...
                [1 1 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end
    end
    
    xlim([-.5 2])
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel('-log10(p-value)', 'FontSize', fontSize, 'FontWeight', 'bold')
    title('Uncorrected p-values', 'FontSize', fontSize+2, 'FontWeight', 'bold')
    grid on
    set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
    box on
    
    saveas(gcf, fullfile(output_dir, 'AOC_omnibus_cluster_permutation.png'));
    fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_cluster_permutation.png'));
end

fprintf('\n=== Done ===\n');
