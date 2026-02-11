%% AOC Omnibus — Sternberg: Split by Alpha Power Change (Microsaccade Rate)
% Loads omnibus_data_rainclouds.csv, splits Sternberg participants by alpha power change
% (reduction vs amplification) per condition and overall. Shows time-resolved microsaccade
% rate for each group with statistical comparisons. Saves figures.
%
% Key outputs:
%   Microsaccade rate time-series figures (reduction vs amplification) per condition and overall

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

% Figure config
fontSize = 25;

% MSSeries time grid: 50 bins of 50ms from -0.5 to 2.0 s
binWidth = 0.05; % 50 ms
twin = [-0.5 2];
edges = twin(1):binWidth:twin(2);
T = numel(edges) - 1; % 50 bins
t_plot = edges(1:end-1) + binWidth/2; % bin centers

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

% Filter for Sternberg task only
sternberg_data = raincloud_data(strcmpi(raincloud_data.Task, 'sternberg'), :);

% Convert ID to string for matching
if isnumeric(sternberg_data.ID)
    sternberg_data.ID = cellstr(num2str(sternberg_data.ID, '%d'));
elseif iscell(sternberg_data.ID)
    sternberg_data.ID = cellfun(@(x) num2str(x), sternberg_data.ID, 'UniformOutput', false);
elseif isstring(sternberg_data.ID)
    sternberg_data.ID = cellstr(sternberg_data.ID);
end

fprintf('Loaded %d Sternberg rows for %d unique subjects\n', ...
    height(sternberg_data), length(unique(sternberg_data.ID)));

%% Define conditions
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [22, 24, 26]; % Condition codes in gaze data

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
    subj_ids_in_csv = unique(sternberg_data.ID);
    
    for s = 1:length(subj_ids_in_csv)
        subj_id = subj_ids_in_csv{s};
        subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
        
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
    
    %% Load gaze data and extract MSSeries
    ms_reduction = nan(length(subjects), T);
    ms_amplification = nan(length(subjects), T);
    
    reduction_subj_count = 0;
    amplification_subj_count = 0;
    missing_gaze = {};
    
    for subj = 1:length(subjects)
        subjID_str = subjects{subj};
        
        % Check if this subject is in either group
        is_reduction = any(strcmp(reduction_ids, subjID_str));
        is_amplification = any(strcmp(amplification_ids, subjID_str));
        
        if ~is_reduction && ~is_amplification
            continue;
        end
        
        % Load gaze data
        datapath_gaze = fullfile(base_data, subjID_str, 'gaze', 'gaze_series_sternberg_trials.mat');
        
        if ~isfile(datapath_gaze)
            missing_gaze{end+1} = subjID_str;
            warning('Missing gaze file for subject %s', subjID_str);
            continue;
        end
        
        try
            load(datapath_gaze, 'MSSeries', 'trialinfo');
        catch ME
            warning('Error loading gaze data for subject %s: %s', subjID_str, ME.message);
            missing_gaze{end+1} = subjID_str;
            continue;
        end
        
        if ~exist('MSSeries', 'var') || isempty(MSSeries)
            warning('MSSeries not found for subject %s. Skipping.', subjID_str);
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
        
        % Filter trials by condition(s)
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
        
        % Ensure trial_mask matches MSSeries dimensions
        n_trials_ms = size(MSSeries, 1);
        if length(trial_mask) > n_trials_ms
            trial_mask = trial_mask(1:n_trials_ms);
        elseif length(trial_mask) < n_trials_ms
            trial_mask(end+1:n_trials_ms) = false;
        end
        
        % Average MSSeries across selected trials
        subj_mean = nanmean(MSSeries(trial_mask, :), 1);
        
        % Ensure correct number of bins
        if length(subj_mean) ~= T
            % Interpolate to match expected bins if needed
            orig_t = linspace(twin(1), twin(2), length(subj_mean)+1);
            orig_centers = orig_t(1:end-1) + diff(orig_t)/2;
            subj_mean = interp1(orig_centers, subj_mean, t_plot, 'linear', NaN);
        end
        
        % Assign to appropriate group
        if is_reduction
            reduction_subj_count = reduction_subj_count + 1;
            ms_reduction(reduction_subj_count, :) = subj_mean;
        end
        
        if is_amplification
            amplification_subj_count = amplification_subj_count + 1;
            ms_amplification(amplification_subj_count, :) = subj_mean;
        end
        
        clear MSSeries trialinfo;
    end
    
    % Trim to actual number of subjects
    ms_reduction = ms_reduction(1:reduction_subj_count, :);
    ms_amplification = ms_amplification(1:amplification_subj_count, :);
    
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
    
    %% Statistical comparison: independent t-tests at each time bin
    p_vals = nan(1, T);
    t_vals = nan(1, T);
    cohens_d = nan(1, T);
    
    for t = 1:T
        red_vals = ms_reduction(:, t);
        amp_vals = ms_amplification(:, t);
        
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
    
    % FDR correction
    valid_mask = isfinite(p_vals);
    p_sub = p_vals(valid_mask);
    if ~isempty(p_sub)
        q_sub = mafdr(p_sub, 'BHFDR', true);
        sig_mask = false(size(p_vals));
        sig_mask(valid_mask) = q_sub < 0.05;
    else
        sig_mask = false(size(p_vals));
    end
    
    %% Grand averages and SEM
    grand_reduction = nanmean(ms_reduction, 1);
    grand_amplification = nanmean(ms_amplification, 1);
    sem_reduction = nanstd(ms_reduction, [], 1) ./ sqrt(sum(isfinite(ms_reduction), 1));
    sem_amplification = nanstd(ms_amplification, [], 1) ./ sqrt(sum(isfinite(ms_amplification), 1));
    
    %% Plot
    close all
    figure
    set(gcf, 'Color', 'w', 'Position', [0 0 1512 900])
    
    % Top subplot: Microsaccade Rate
    subplot(4,1,1:3)
    hold on
    
    grid on
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', '--', 'MinorGridAlpha', 0.05)
    
    sebRed = shadedErrorBar(t_plot, grand_reduction, sem_reduction, ...
        'lineProps', {'-','Color',colors(1,:),'LineWidth',3.5}, 'transparent', true);
    sebAmp = shadedErrorBar(t_plot, grand_amplification, sem_amplification, ...
        'lineProps', {'-','Color',colors(3,:),'LineWidth',3.5}, 'transparent', true);
    
    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6)
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Alpha', 0.5)
    
    ylabel('Microsaccade Rate [MS/s]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('Sternberg MS Rate: Alpha Reduction vs Amplification (%s)', analysis_name), ...
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
    
    grid on
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', '--', 'MinorGridAlpha', 0.05)
    
    y = cohens_d;
    ylab = 'Cohen''s d';
    
    ylim_range = [-max(abs(y(isfinite(y))))*1.15 max(abs(y(isfinite(y))))*1.15];
    if any(isfinite(y))
        ylim(ylim_range);
    else
        ylim_range = [-1 1];
        ylim(ylim_range);
    end
    
    % Shade significant regions
    d_sig = diff([0, sig_mask, 0]);
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
    
    if any(sig_mask)
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
    fig_name = sprintf('AOC_omnibus_splitAlpha_msrate_%s.png', cond_label);
    fig_path = fullfile(output_dir, fig_name);
    saveas(gcf, fig_path);
    fprintf('  Saved figure: %s\n', fig_name);
    
end

fprintf('\n=== Done ===\n');
fprintf('All figures saved to: %s\n', output_dir);
