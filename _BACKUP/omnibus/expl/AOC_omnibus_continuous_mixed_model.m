%% AOC Omnibus — Continuous Mixed Model (Alpha Power as Continuous Predictor)
% Uses alpha power as continuous predictor in a mixed model framework.
% More nuanced than binary split. Saves data for R/Python analysis.
%
% Key outputs:
%   CSV file with data for mixed model analysis; summary statistics

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 25;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/stats/omnibus';
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

% Create time bins for mixed model
time_bins = 0:0.1:2; % 100ms bins
n_bins = length(time_bins) - 1;

conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
cond_codes = [22, 24, 26];

subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);

%% Prepare long-format data for mixed model
mixed_model_data = table();

row_idx = 1;
for s = 1:n_subj
    subj_id = subj_ids{s};
    
    % Get alpha power per condition
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    alpha_per_cond = nan(1, 3);
    for c = 1:3
        cond_rows = subj_rows(strcmp(subj_rows.Condition, conditions{c}), :);
        if ~isempty(cond_rows) && ~isnan(cond_rows.AlphaPower(1)) && cond_rows.AlphaPower(1) ~= 0
            alpha_per_cond(c) = cond_rows.AlphaPower(1);
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
        if ~isfinite(alpha_per_cond(c))
            continue;
        end
        
        cond_code = cond_codes(c);
        trial_mask = (conds == cond_code);
        
        if ~any(trial_mask)
            continue;
        end
        
        % Compute SPL time series
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
        
        % Bin into time windows
        for b = 1:n_bins
            bin_mask = t_plot_full >= time_bins(b) & t_plot_full < time_bins(b+1);
            spl_bin = nanmean(subj_mean(bin_mask));
            
            if isfinite(spl_bin)
                mixed_model_data.ID(row_idx) = {subj_id};
                mixed_model_data.Condition(row_idx) = {conditions{c}};
                mixed_model_data.AlphaPower(row_idx) = alpha_per_cond(c);
                mixed_model_data.TimeBin(row_idx) = (time_bins(b) + time_bins(b+1)) / 2;
                mixed_model_data.SPL(row_idx) = spl_bin;
                row_idx = row_idx + 1;
            end
        end
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Save data for R/Python analysis
if height(mixed_model_data) > 0
    output_file = fullfile(output_dir, 'AOC_omnibus_mixed_model_data.csv');
    writetable(mixed_model_data, output_file);
    fprintf('\nSaved mixed model data: %s\n', output_file);
    fprintf('Rows: %d, Subjects: %d, Conditions: %d, Time bins: %d\n', ...
        height(mixed_model_data), length(unique(mixed_model_data.ID)), ...
        length(unique(mixed_model_data.Condition)), length(unique(mixed_model_data.TimeBin)));
    
    % Summary statistics
    fprintf('\n=== Summary Statistics ===\n');
    fprintf('Alpha Power: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', ...
        nanmean(mixed_model_data.AlphaPower), nanstd(mixed_model_data.AlphaPower), ...
        min(mixed_model_data.AlphaPower), max(mixed_model_data.AlphaPower));
    fprintf('SPL: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', ...
        nanmean(mixed_model_data.SPL), nanstd(mixed_model_data.SPL), ...
        min(mixed_model_data.SPL), max(mixed_model_data.SPL));
    
    % Correlation
    [r, p] = corrcoef(mixed_model_data.AlphaPower, mixed_model_data.SPL);
    fprintf('Overall correlation (AlphaPower × SPL): r=%.3f, p=%.4f\n', r(1,2), p(1,2));
    
    fprintf('\n=== Model Formula Suggestion (R) ===\n');
    fprintf('library(lme4)\n');
    fprintf('model <- lmer(SPL ~ AlphaPower * Condition * TimeBin + (1|ID), data=dat)\n');
    fprintf('summary(model)\n');
end

fprintf('\n=== Done ===\n');
