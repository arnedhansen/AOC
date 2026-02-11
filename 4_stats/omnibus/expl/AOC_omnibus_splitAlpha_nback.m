%% AOC Omnibus — N-Back: Split by Alpha Power Change
% Loads omnibus_data_rainclouds.csv, splits N-Back participants by alpha power change
% (reduction vs amplification) per condition and overall. Shows time-resolved scan path
% length for each group with statistical comparisons. Saves figures.
%
% Key outputs:
%   SPL time-series figures (reduction vs amplification) per condition and overall

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

% Figure config
fontSize = 25;

% Common reference grid for gaze step-length series
fs_full     = 500;                                  % matches your gaze pipeline
t_full      = -0.5:1/fs_full:2;                     % sample grid incl. endpoints
t_plot_full = t_full(2:end);                        % step series aligns to t_full(2:end)
Tf          = numel(t_plot_full);                   % full-resolution length

% Time series for binned display (matching reference script)
time_series = linspace(-0.5, 2, 51); % 51 points -> 50 steps
T = numel(time_series) - 1; % step series aligns to time_series(2:end)

% Set up paths (cross-platform)
if ispc
    base_data = 'W:\Students\Arne\AOC\data\features';
    output_dir = 'W:\Students\Arne\AOC\figures\interactions\omnibus_alpha_split';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features';
    output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split';
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load CSV data
csv_path = fullfile(base_data, 'AOC_omnibus_raincloud_data.csv');
if ~isfile(csv_path)
    error('CSV file not found: %s', csv_path);
end

fprintf('Loading CSV data from: %s\n', csv_path);
raincloud_data = readtable(csv_path);

% Filter for N-Back task only
nback_data = raincloud_data(strcmpi(raincloud_data.Task, 'nback'), :);

% Convert ID to string for matching
if isnumeric(nback_data.ID)
    nback_data.ID = cellstr(num2str(nback_data.ID, '%d'));
elseif iscell(nback_data.ID)
    nback_data.ID = cellfun(@(x) num2str(x), nback_data.ID, 'UniformOutput', false);
elseif isstring(nback_data.ID)
    nback_data.ID = cellstr(nback_data.ID);
end

fprintf('Loaded %d N-Back rows for %d unique subjects\n', ...
    height(nback_data), length(unique(nback_data.ID)));

%% Define conditions
conditions = {'1-back', '2-back', '3-back'};
cond_codes = [21, 22, 23]; % Condition codes in gaze trialinfo

%% Process each condition separately + overall mean
all_analyses = [conditions, {'Mean across all conditions'}];

for analysis_idx = 1:length(all_analyses)
    analysis_name = all_analyses{analysis_idx};
    fprintf('\n=== Processing: %s ===\n', analysis_name);
    
    % Determine which condition(s) to use
    if strcmp(analysis_name, 'Mean across all conditions')
        use_conditions = conditions;
        cond_label = 'mean';
    else
        use_conditions = {analysis_name};
        cond_label = lower(strrep(analysis_name, ' ', ''));
    end
    
    % Extract alpha power for relevant conditions
    subj_alpha = containers.Map('KeyType', 'char', 'ValueType', 'double');
    subj_ids_in_csv = unique(nback_data.ID);
    
    for s = 1:length(subj_ids_in_csv)
        subj_id = subj_ids_in_csv{s};
        subj_rows = nback_data(strcmp(nback_data.ID, subj_id), :);
        
        % Get alpha power for this subject in relevant conditions
        alpha_vals = [];
        for c = 1:length(use_conditions)
            cond_rows = subj_rows(strcmp(subj_rows.Condition, use_conditions{c}), :);
            if ~isempty(cond_rows) && ~isnan(cond_rows.AlphaPower(1))
                alpha_vals(end+1) = cond_rows.AlphaPower(1);
            end
        end
        
        if ~isempty(alpha_vals)
            mean_alpha = nanmean(alpha_vals);
            if mean_alpha ~= 0
                subj_alpha(subj_id) = mean_alpha;
            else
                fprintf('  Excluding subject %s (alpha power = 0)\n', subj_id);
            end
        end
    end
    
    % Split into reduction (< 0) and amplification (> 0)
    reduction_ids = {};
    amplification_ids = {};
    
    subj_keys = keys(subj_alpha);
    for k = 1:length(subj_keys)
        subj_id = subj_keys{k};
        alpha_val = subj_alpha(subj_id);
        if alpha_val < 0
            reduction_ids{end+1} = subj_id;
        else
            amplification_ids{end+1} = subj_id;
        end
    end
    
    fprintf('  Reduction group: %d subjects\n', length(reduction_ids));
    fprintf('  Amplification group: %d subjects\n', length(amplification_ids));
    
    if length(reduction_ids) < 2 || length(amplification_ids) < 2
        warning('Not enough subjects in one or both groups for %s. Skipping.', analysis_name);
        continue;
    end
    
    %% Load gaze data and compute SPL time series
    scan_reduction = nan(length(subjects), Tf);
    scan_amplification = nan(length(subjects), Tf);
    
    reduction_subj_count = 0;
    amplification_subj_count = 0;
    missing_gaze = {};
    
    for subj = 1:length(subjects)
        subjID_str = subjects{subj};
        
        is_reduction = any(strcmp(reduction_ids, subjID_str));
        is_amplification = any(strcmp(amplification_ids, subjID_str));
        
        if ~is_reduction && ~is_amplification
            continue;
        end
        
        % Load gaze data
        datapath_gaze = fullfile(base_data, subjID_str, 'gaze', 'gaze_series_nback_trials.mat');
        
        if ~isfile(datapath_gaze)
            missing_gaze{end+1} = subjID_str;
            warning('Missing gaze file for subject %s', subjID_str);
            continue;
        end
        
        try
            load(datapath_gaze, 'ScanPathSeries', 'ScanPathSeriesT', 'trialinfo');
        catch ME
            warning('Error loading gaze data for subject %s: %s', subjID_str, ME.message);
            missing_gaze{end+1} = subjID_str;
            continue;
        end
        
        % Parse trialinfo to get condition codes
        if exist('trialinfo', 'var')
            if size(trialinfo, 1) == 2
                conds = trialinfo(1, :)';
            elseif size(trialinfo, 2) == 2
                conds = trialinfo(:, 1);
            else
                warning('Unexpected trialinfo shape for subject %s. Skipping.', subjID_str);
                continue;
            end
        else
            warning('trialinfo not found for subject %s. Skipping.', subjID_str);
            continue;
        end
        
        % Filter trials by condition(s) for this analysis
        if strcmp(analysis_name, 'Mean across all conditions')
            trial_mask = true(size(conds));
        else
            cond_idx = find(strcmp(conditions, analysis_name));
            if isempty(cond_idx)
                continue;
            end
            cond_code = cond_codes(cond_idx);
            trial_mask = (conds == cond_code);
        end
        
        if ~any(trial_mask)
            continue;
        end
        
        % Interpolate each trial to the common full grid
        subj_trials_full = nan(numel(ScanPathSeries), Tf);
        
        for trl = 1:numel(ScanPathSeries)
            if ~trial_mask(trl)
                continue;
            end
            
            srl_full = ScanPathSeries{trl};
            tt_full  = ScanPathSeriesT{trl};
            
            if isempty(srl_full) || isempty(tt_full)
                continue;
            end
            
            try
                subj_trials_full(trl, :) = interp1(tt_full, srl_full, t_plot_full, 'linear', NaN);
            catch
                % leave as NaN
            end
        end
        
        % Average across trials for this subject
        subj_mean = nanmean(subj_trials_full(trial_mask, :), 1);
        
        if is_reduction
            reduction_subj_count = reduction_subj_count + 1;
            scan_reduction(reduction_subj_count, :) = subj_mean;
        end
        
        if is_amplification
            amplification_subj_count = amplification_subj_count + 1;
            scan_amplification(amplification_subj_count, :) = subj_mean;
        end
        
        clear ScanPathSeries ScanPathSeriesT trialinfo;
    end
    
    % Trim to actual number of subjects
    scan_reduction = scan_reduction(1:reduction_subj_count, :);
    scan_amplification = scan_amplification(1:amplification_subj_count, :);
    
    if ~isempty(missing_gaze)
        fprintf('  Warning: Missing gaze data for %d subjects: %s\n', ...
            length(missing_gaze), strjoin(missing_gaze, ', '));
    end
    
    fprintf('  Loaded gaze data: %d reduction, %d amplification subjects\n', ...
        reduction_subj_count, amplification_subj_count);
    
    if reduction_subj_count < 2 || amplification_subj_count < 2
        warning('Not enough subjects with gaze data for %s. Skipping.', analysis_name);
        continue;
    end
    
    %% Statistical comparison: independent t-tests at each time point
    p_vals = nan(1, Tf);
    t_vals = nan(1, Tf);
    cohens_d = nan(1, Tf);
    
    for t = 1:Tf
        red_vals = scan_reduction(:, t);
        amp_vals = scan_amplification(:, t);
        
        red_vals = red_vals(isfinite(red_vals));
        amp_vals = amp_vals(isfinite(amp_vals));
        
        if length(red_vals) >= 2 && length(amp_vals) >= 2
            [~, p, ~, stats] = ttest2(red_vals, amp_vals);
            p_vals(t) = p;
            t_vals(t) = stats.tstat;
            
            pooled_std = sqrt(((length(red_vals)-1)*nanvar(red_vals) + ...
                (length(amp_vals)-1)*nanvar(amp_vals)) / ...
                (length(red_vals) + length(amp_vals) - 2));
            if pooled_std > 0
                cohens_d(t) = (nanmean(amp_vals) - nanmean(red_vals)) / pooled_std;
            end
        end
    end
    
    % FDR correction across time points
    toi = t_plot_full >= -0.5 & t_plot_full <= 2;
    valid_mask = toi & isfinite(p_vals);
    p_sub = p_vals(valid_mask);
    if ~isempty(p_sub)
        q_sub = mafdr(p_sub, 'BHFDR', true);
        sig_mask = false(size(p_vals));
        sig_mask(valid_mask) = q_sub < 0.05;
    else
        sig_mask = false(size(p_vals));
    end
    
    %% Prepare data for plotting (binned grid)
    scan_reduction_binned = nan(size(scan_reduction, 1), T);
    scan_amplification_binned = nan(size(scan_amplification, 1), T);
    
    for s = 1:size(scan_reduction, 1)
        scan_reduction_binned(s, :) = interp1(t_plot_full, scan_reduction(s, :), ...
            time_series(2:end), 'linear', NaN);
    end
    
    for s = 1:size(scan_amplification, 1)
        scan_amplification_binned(s, :) = interp1(t_plot_full, scan_amplification(s, :), ...
            time_series(2:end), 'linear', NaN);
    end
    
    % Grand averages and SEM (binned)
    grand_reduction = nanmean(scan_reduction_binned, 1);
    grand_amplification = nanmean(scan_amplification_binned, 1);
    sem_reduction = nanstd(scan_reduction_binned, [], 1) ./ sqrt(sum(isfinite(scan_reduction_binned), 1));
    sem_amplification = nanstd(scan_amplification_binned, [], 1) ./ sqrt(sum(isfinite(scan_amplification_binned), 1));
    t_plot = time_series(2:end);
    
    % Cohen's d for binned data (for plotting)
    d_vals_binned = nan(1, T);
    for t = 1:T
        red_vals = scan_reduction_binned(:, t);
        amp_vals = scan_amplification_binned(:, t);
        red_vals = red_vals(isfinite(red_vals));
        amp_vals = amp_vals(isfinite(amp_vals));
        
        if length(red_vals) >= 2 && length(amp_vals) >= 2
            pooled_std = sqrt(((length(red_vals)-1)*nanvar(red_vals) + ...
                (length(amp_vals)-1)*nanvar(amp_vals)) / ...
                (length(red_vals) + length(amp_vals) - 2));
            if pooled_std > 0
                d_vals_binned(t) = (nanmean(amp_vals) - nanmean(red_vals)) / pooled_std;
            end
        end
    end
    
    % Interpolate significance mask to binned grid
    sig_mask_binned = false(1, T);
    for t = 1:T
        t_val = time_series(t+1);
        [~, idx] = min(abs(t_plot_full - t_val));
        if idx <= length(sig_mask)
            sig_mask_binned(t) = sig_mask(idx);
        end
    end
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1512 900])
    
    % Top subplot: Scan Path Length
    subplot(4,1,1:3)
    hold on
    
    sebRed = shadedErrorBar(t_plot, grand_reduction, sem_reduction, ...
        'lineProps', {'-','Color',colors(1,:),'LineWidth',3.5}, 'transparent', true);
    sebAmp = shadedErrorBar(t_plot, grand_amplification, sem_amplification, ...
        'lineProps', {'-','Color',colors(3,:),'LineWidth',3.5}, 'transparent', true);
    
    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6)
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Alpha', 0.5)
    
    ylabel('Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('N-Back SPL: Alpha Reduction vs Amplification (%s)', analysis_name), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    xlim([-.5 2])
    set(gca, 'FontSize', fontSize-2)
    box on
    
    legend([sebRed.mainLine, sebAmp.mainLine], ...
        {sprintf('Reduction (n=%d) ± SEM', reduction_subj_count), ...
         sprintf('Amplification (n=%d) ± SEM', amplification_subj_count)}, ...
        'Location','northeast', 'FontSize', fontSize*0.7, 'FontWeight', 'normal', ...
        'Box', 'on', 'EdgeColor', [0.7 0.7 0.7])
    
    set(gca, 'XTick', -0.5:0.5:2, 'XColor', 'none')
    set(gca, 'YTickLabelMode', 'auto')
    
    % Bottom subplot: Cohen's d with significance
    subplot(4,1,4)
    hold on
    
    y = d_vals_binned;
    ylab = 'Cohen''s d';
    
    ylim_range = [-max(abs(y(isfinite(y))))*1.15 max(abs(y(isfinite(y))))*1.15];
    if any(isfinite(y))
        ylim(ylim_range);
    else
        ylim_range = [-1 1];
        ylim(ylim_range);
    end
    
    d_sig = diff([0, sig_mask_binned, 0]);
    on_sig = find(d_sig == 1);
    off_sig = find(d_sig == -1) - 1;
    
    for k = 1:numel(on_sig)
        if on_sig(k) <= length(t_plot) && off_sig(k) <= length(t_plot)
            x0 = t_plot(on_sig(k));
            x1 = t_plot(off_sig(k));
            patch([x0 x1 x1 x0], [ylim_range(1) ylim_range(1) ylim_range(2) ylim_range(2)], ...
                [1 1 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'LineStyle', 'none')
        end
    end
    
    h = plot(t_plot, y, 'k-', 'LineWidth', 3.5);
    
    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6)
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'Alpha', 0.7)
    yline(0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'Alpha', 0.5)
    yline(-0.2, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'Alpha', 0.5)
    yline(0.5, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'Alpha', 0.5)
    yline(-0.5, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'Alpha', 0.5)
    
    xlim([-.5 2])
    set(gca, 'FontSize', fontSize-2)
    xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    ylabel(ylab, 'FontSize', fontSize, 'FontWeight', 'bold')
    
    if max(abs(ylim_range)) <= 1
        yticks([-1 -0.5 -0.2 0 0.2 0.5 1]);
    else
        yticks([-.5 -.25 0 .25 .5]);
    end
    
    if any(sig_mask_binned)
        sig_periods = [];
        for k = 1:numel(on_sig)
            if on_sig(k) <= length(t_plot) && off_sig(k) <= length(t_plot)
                sig_periods(end+1, :) = [t_plot(on_sig(k)), t_plot(off_sig(k))];
            end
        end
        if ~isempty(sig_periods)
            mid_time = mean(sig_periods, 2);
            for k = 1:size(sig_periods, 1)
                text(mid_time(k), ylim_range(2)*0.85, '*', ...
                    'FontSize', fontSize, 'HorizontalAlignment', 'center', ...
                    'Color', [0.8 0.6 0], 'FontWeight', 'bold')
            end
        end
    end
    
    box on
    
    % Save figure
    fig_name = sprintf('AOC_omnibus_splitAlpha_nback_%s.png', cond_label);
    fig_path = fullfile(output_dir, fig_name);
    saveas(gcf, fig_path);
    fprintf('  Saved figure: %s\n', fig_name);
    
end

fprintf('\n=== Done ===\n');
fprintf('All figures saved to: %s\n', output_dir);
