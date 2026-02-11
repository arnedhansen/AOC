%% AOC Omnibus — Sternberg: Within-Subject Trial Split by Alpha Power Change
% Splits trials within each subject by baselined alpha power (reduction vs amplification),
% averages SPL time-course per group per subject, then performs paired statistics across
% subjects. Uses AlphaPowerLateBL from merged_data_sternberg_trials.mat.
%
% Key difference from AOC_omnibus_splitAlpha_sternberg.m:
%   - Between-subject version splits SUBJECTS by mean alpha change
%   - This script splits TRIALS within each subject by trial-level alpha change
%
% Key outputs:
%   SPL time-series figures (reduction vs amplification trials) per condition and overall
%   Paired Cohen's dz over time with FDR-corrected significance

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

% Time series for binned display
time_series = linspace(-0.5, 2, 51); % 51 points -> 50 steps
T = numel(time_series) - 1;

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

%% Load trial-level merged data
merged_path = fullfile(base_data, 'merged_data_sternberg_trials.mat');
if ~isfile(merged_path)
    error('Merged data file not found: %s', merged_path);
end

fprintf('Loading merged trial data from: %s\n', merged_path);
load(merged_path, 'merged_data_sternberg_trials');

fprintf('Loaded %d trials for %d unique subjects\n', ...
    height(merged_data_sternberg_trials), length(unique(merged_data_sternberg_trials.ID)));

%% Define conditions
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [2, 4, 6]; % Condition codes in merged data (raw 22/24/26 minus 20)

%% Process each condition separately + overall mean
all_analyses = [conditions, {'Mean across all conditions'}];

for analysis_idx = 1:length(all_analyses)
    analysis_name = all_analyses{analysis_idx};
    fprintf('\n=== Processing: %s ===\n', analysis_name);
    
    if strcmp(analysis_name, 'Mean across all conditions')
        cond_label = 'mean';
    else
        cond_label = lower(strrep(analysis_name, ' ', ''));
    end
    
    %% Per-subject: split trials by alpha and compute SPL
    scan_reduction = nan(length(subjects), Tf);   % per-subject mean SPL for reduction trials
    scan_amplification = nan(length(subjects), Tf); % per-subject mean SPL for amplification trials
    valid_subj_count = 0;
    n_red_trials_total = 0;
    n_amp_trials_total = 0;
    
    for subj = 1:length(subjects)
        subjID_str = subjects{subj};
        subjID_num = str2double(subjID_str);
        
        fprintf('  Subject %s (%d/%d)... ', subjID_str, subj, length(subjects));
        
        %% Get trial-level alpha power for this subject
        rows = merged_data_sternberg_trials(merged_data_sternberg_trials.ID == subjID_num, :);
        
        if isempty(rows)
            fprintf('no merged data. Skipping.\n');
            continue;
        end
        
        % Filter by condition if needed
        if ~strcmp(analysis_name, 'Mean across all conditions')
            cond_idx = find(strcmp(conditions, analysis_name));
            cond_code = cond_codes(cond_idx);
            rows = rows(rows.Condition == cond_code, :);
        end
        
        if isempty(rows)
            fprintf('no trials for this condition. Skipping.\n');
            continue;
        end
        
        % Get alpha power (baselined) and trial numbers
        alpha_bl = rows.AlphaPowerLateBL;
        trial_nums = rows.Trial;
        
        % Split by alpha change: < 0 = reduction, > 0 = amplification
        valid = isfinite(alpha_bl) & alpha_bl ~= 0;
        red_trials = trial_nums(valid & alpha_bl < 0);
        amp_trials = trial_nums(valid & alpha_bl > 0);
        
        if length(red_trials) < 2 || length(amp_trials) < 2
            fprintf('too few trials (red=%d, amp=%d). Skipping.\n', ...
                length(red_trials), length(amp_trials));
            continue;
        end
        
        %% Load gaze data
        datapath_gaze = fullfile(base_data, subjID_str, 'gaze', 'gaze_series_sternberg_trials.mat');
        
        if ~isfile(datapath_gaze)
            fprintf('missing gaze file. Skipping.\n');
            continue;
        end
        
        try
            load(datapath_gaze, 'ScanPathSeries', 'ScanPathSeriesT', 'trialinfo');
        catch ME
            fprintf('error loading gaze: %s. Skipping.\n', ME.message);
            continue;
        end
        
        if isempty(ScanPathSeries)
            fprintf('empty ScanPathSeries. Skipping.\n');
            continue;
        end
        
        % Get trial numbers from gaze data
        % Note: trialinfo is saved as dataETlong.trialinfo' (transposed),
        % so it is [2 x nTrials] with row 1 = condition, row 2 = trial number
        if exist('trialinfo', 'var')
            if size(trialinfo, 1) == 2
                gaze_trial_nums = trialinfo(2, :)'; % row 2 = Trial number (transposed format)
            elseif size(trialinfo, 2) >= 2
                gaze_trial_nums = trialinfo(:, 2); % column 2 = Trial number (original format)
            else
                fprintf('unexpected trialinfo shape. Skipping.\n');
                continue;
            end
        else
            fprintf('trialinfo not found. Skipping.\n');
            continue;
        end
        
        %% Interpolate each trial to common grid
        n_gaze_trials = numel(ScanPathSeries);
        subj_trials_full = nan(n_gaze_trials, Tf);
        
        for trl = 1:n_gaze_trials
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
        
        %% Match trials and compute group means
        red_mask = ismember(gaze_trial_nums, red_trials);
        amp_mask = ismember(gaze_trial_nums, amp_trials);
        
        if sum(red_mask) < 2 || sum(amp_mask) < 2
            fprintf('too few matched gaze trials (red=%d, amp=%d). Skipping.\n', ...
                sum(red_mask), sum(amp_mask));
            continue;
        end
        
        % Average SPL across trials within each group
        valid_subj_count = valid_subj_count + 1;
        scan_reduction(valid_subj_count, :) = nanmean(subj_trials_full(red_mask, :), 1);
        scan_amplification(valid_subj_count, :) = nanmean(subj_trials_full(amp_mask, :), 1);
        
        n_red_trials_total = n_red_trials_total + sum(red_mask);
        n_amp_trials_total = n_amp_trials_total + sum(amp_mask);
        
        fprintf('OK (red=%d, amp=%d trials)\n', sum(red_mask), sum(amp_mask));
        
        clear ScanPathSeries ScanPathSeriesT trialinfo;
    end
    
    % Trim to actual number of valid subjects
    scan_reduction = scan_reduction(1:valid_subj_count, :);
    scan_amplification = scan_amplification(1:valid_subj_count, :);
    
    fprintf('\n  Valid subjects: %d\n', valid_subj_count);
    fprintf('  Total trials: %d reduction, %d amplification\n', ...
        n_red_trials_total, n_amp_trials_total);
    
    if valid_subj_count < 3
        warning('Not enough valid subjects for %s (%d < 3). Skipping.', analysis_name, valid_subj_count);
        continue;
    end
    
    %% Statistical comparison: PAIRED t-tests at each time point
    p_vals = nan(1, Tf);
    t_vals = nan(1, Tf);
    cohens_dz = nan(1, Tf);
    
    for t = 1:Tf
        red_vals = scan_reduction(:, t);
        amp_vals = scan_amplification(:, t);
        
        % Only use subjects with valid data in both groups
        valid = isfinite(red_vals) & isfinite(amp_vals);
        
        if sum(valid) >= 3
            % Paired t-test (within-subject)
            [~, p, ~, stats] = ttest(amp_vals(valid), red_vals(valid));
            p_vals(t) = p;
            t_vals(t) = stats.tstat;
            
            % Cohen's dz for paired data: mean(diff) / std(diff)
            diffs = amp_vals(valid) - red_vals(valid);
            if nanstd(diffs) > 0
                cohens_dz(t) = nanmean(diffs) / nanstd(diffs);
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
    scan_reduction_binned = nan(valid_subj_count, T);
    scan_amplification_binned = nan(valid_subj_count, T);
    
    for s = 1:valid_subj_count
        scan_reduction_binned(s, :) = interp1(t_plot_full, scan_reduction(s, :), ...
            time_series(2:end), 'linear', NaN);
        scan_amplification_binned(s, :) = interp1(t_plot_full, scan_amplification(s, :), ...
            time_series(2:end), 'linear', NaN);
    end
    
    % Grand averages and SEM (binned)
    grand_reduction = nanmean(scan_reduction_binned, 1);
    grand_amplification = nanmean(scan_amplification_binned, 1);
    sem_reduction = nanstd(scan_reduction_binned, [], 1) ./ sqrt(sum(isfinite(scan_reduction_binned), 1));
    sem_amplification = nanstd(scan_amplification_binned, [], 1) ./ sqrt(sum(isfinite(scan_amplification_binned), 1));
    t_plot = time_series(2:end);
    
    % Cohen's dz for binned data
    dz_vals_binned = nan(1, T);
    for t = 1:T
        red_vals = scan_reduction_binned(:, t);
        amp_vals = scan_amplification_binned(:, t);
        valid = isfinite(red_vals) & isfinite(amp_vals);
        
        if sum(valid) >= 3
            diffs = amp_vals(valid) - red_vals(valid);
            if nanstd(diffs) > 0
                dz_vals_binned(t) = nanmean(diffs) / nanstd(diffs);
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
    
    % Add subtle grid
    grid on
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', '--', 'MinorGridAlpha', 0.05)
    
    % Plot with thicker lines
    sebRed = shadedErrorBar(t_plot, grand_reduction, sem_reduction, ...
        'lineProps', {'-','Color',colors(1,:),'LineWidth',3.5}, 'transparent', true);
    sebAmp = shadedErrorBar(t_plot, grand_amplification, sem_amplification, ...
        'lineProps', {'-','Color',colors(3,:),'LineWidth',3.5}, 'transparent', true);
    
    % Add vertical line at time 0 (stimulus onset)
    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6)
    
    % Add horizontal reference line
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Alpha', 0.5)
    
    ylabel('Scan Path Length [px]', 'FontSize', fontSize, 'FontWeight', 'bold')
    title(sprintf('Sternberg SPL (Within-Subject): Reduction vs Amplification Trials (%s)', analysis_name), ...
        'FontSize', fontSize+2, 'FontWeight', 'bold')
    xlim([-.5 2])
    set(gca, 'FontSize', fontSize-2)
    box on
    
    % Legend with trial counts
    legend([sebRed.mainLine, sebAmp.mainLine], ...
        {sprintf('Reduction trials ± SEM (n=%d subj)', valid_subj_count), ...
         sprintf('Amplification trials ± SEM (n=%d subj)', valid_subj_count)}, ...
        'Location','northeast', 'FontSize', fontSize*0.7, 'FontWeight', 'normal', ...
        'Box', 'on', 'EdgeColor', [0.7 0.7 0.7])
    
    set(gca, 'XTick', -0.5:0.5:2, 'XColor', 'none')
    set(gca, 'YTickLabelMode', 'auto')
    
    % Bottom subplot: Cohen's dz with significance
    subplot(4,1,4)
    hold on
    
    % Add subtle grid
    grid on
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', '--', 'MinorGridAlpha', 0.05)
    
    y = dz_vals_binned;
    ylab = 'Cohen''s d_z';
    
    % Calculate ylim
    ylim_range = [-max(abs(y(isfinite(y))))*1.15 max(abs(y(isfinite(y))))*1.15];
    if any(isfinite(y))
        ylim(ylim_range);
    else
        ylim_range = [-1 1];
        ylim(ylim_range);
    end
    
    % Shade significant regions
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
    
    % Plot line on top
    h = plot(t_plot, y, 'k-', 'LineWidth', 3.5);
    
    % Reference lines
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
    
    % Annotate significant periods
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
    fig_name = sprintf('AOC_omnibus_splitAlpha_withinsubject_%s.png', cond_label);
    fig_path = fullfile(output_dir, fig_name);
    saveas(gcf, fig_path);
    fprintf('  Saved figure: %s\n', fig_name);
    
    %% Print summary statistics
    fprintf('\n--- Summary Statistics (%s) ---\n', analysis_name);
    fprintf('Valid subjects: %d\n', valid_subj_count);
    fprintf('Total trials: %d reduction, %d amplification\n', ...
        n_red_trials_total, n_amp_trials_total);
    
    % Overall effect
    overall_diffs = nanmean(scan_amplification_binned, 2) - nanmean(scan_reduction_binned, 2);
    valid_diffs = overall_diffs(isfinite(overall_diffs));
    if ~isempty(valid_diffs)
        [~, p_overall] = ttest(valid_diffs);
        dz_overall = mean(valid_diffs) / std(valid_diffs);
        fprintf('Overall paired t-test: t(%d)=%.3f, p=%.4f, dz=%.3f\n', ...
            length(valid_diffs)-1, mean(valid_diffs)/(std(valid_diffs)/sqrt(length(valid_diffs))), ...
            p_overall, dz_overall);
    end
    
    % Peak effect
    if any(isfinite(dz_vals_binned))
        [max_dz, max_idx] = max(abs(dz_vals_binned));
        fprintf('Peak effect size: dz=%.3f at t=%.2fs\n', ...
            dz_vals_binned(max_idx), t_plot(max_idx));
    end
    
    % Significant time windows
    if any(sig_mask_binned)
        fprintf('Significant time windows (FDR q<0.05):\n');
        for k = 1:numel(on_sig)
            if on_sig(k) <= length(t_plot) && off_sig(k) <= length(t_plot)
                fprintf('  %.2f - %.2f s\n', t_plot(on_sig(k)), t_plot(off_sig(k)));
            end
        end
    else
        fprintf('No significant time windows after FDR correction.\n');
    end
    
end

fprintf('\n=== Done ===\n');
fprintf('All figures saved to: %s\n', output_dir);
