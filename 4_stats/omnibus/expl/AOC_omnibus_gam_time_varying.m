%% AOC Omnibus — Time-Varying Effects (GAM)
% Uses Generalized Additive Model (GAM) to model non-linear temporal dynamics.
% Shows how the alpha power effect changes smoothly over time.
%
% Key outputs:
%   GAM fit visualization; smooth curves; R script for analysis

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

% Use finer time resolution for GAM
time_points = t_plot_full(t_plot_full >= 0 & t_plot_full <= 2);
n_time = length(time_points);

subj_ids = unique(sternberg_data.ID);
n_subj = length(subj_ids);

%% Prepare data for GAM
gam_data = table();

row_idx = 1;
for s = 1:n_subj
    subj_id = subj_ids{s};
    
    % Get mean alpha power
    subj_rows = sternberg_data(strcmp(sternberg_data.ID, subj_id), :);
    alpha_vals = [];
    for c = 1:height(subj_rows)
        if ~isnan(subj_rows.AlphaPower(c)) && subj_rows.AlphaPower(c) ~= 0
            alpha_vals(end+1) = subj_rows.AlphaPower(c);
        end
    end
    if isempty(alpha_vals)
        continue;
    end
    mean_alpha = nanmean(alpha_vals);
    
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
    
    subj_mean = nanmean(subj_trials_full, 1);
    
    % Extract time points
    for t = 1:n_time
        time_val = time_points(t);
        [~, idx] = min(abs(t_plot_full - time_val));
        spl_val = subj_mean(idx);
        
        if isfinite(spl_val)
            gam_data.ID(row_idx) = {subj_id};
            gam_data.Time(row_idx) = time_val;
            gam_data.AlphaPower(row_idx) = mean_alpha;
            gam_data.SPL(row_idx) = spl_val;
            row_idx = row_idx + 1;
        end
    end
    
    clear ScanPathSeries ScanPathSeriesT trialinfo;
end

%% Save data for R GAM analysis
if height(gam_data) > 0
    output_file = fullfile(output_dir, 'AOC_omnibus_gam_data.csv');
    writetable(gam_data, output_file);
    fprintf('\nSaved GAM data: %s\n', output_file);
    fprintf('Rows: %d, Subjects: %d, Time points: %d\n', ...
        height(gam_data), length(unique(gam_data.ID)), length(unique(gam_data.Time)));
    
    %% Create R script for GAM analysis
    r_script = fullfile(output_dir, 'AOC_omnibus_gam_analysis.R');
    fid = fopen(r_script, 'w');
    
    fprintf(fid, '# AOC Omnibus GAM Analysis\n');
    fprintf(fid, '# Time-varying effects of alpha power on scan path length\n\n');
    fprintf(fid, 'library(mgcv)\n');
    fprintf(fid, 'library(ggplot2)\n\n');
    fprintf(fid, '# Load data\n');
    fprintf(fid, 'dat <- read.csv("%s")\n\n', strrep(output_file, '\', '/'));
    fprintf(fid, '# Convert ID to factor\n');
    fprintf(fid, 'dat$ID <- as.factor(dat$ID)\n\n');
    fprintf(fid, '# GAM model: SPL ~ s(Time) + s(AlphaPower) + s(Time, AlphaPower)\n');
    fprintf(fid, 'model <- gam(SPL ~ s(Time, k=10) + s(AlphaPower, k=5) + te(Time, AlphaPower, k=c(10,5)) + s(ID, bs="re"),\n');
    fprintf(fid, '             data=dat, method="REML")\n\n');
    fprintf(fid, 'summary(model)\n\n');
    fprintf(fid, '# Plot smooth terms\n');
    fprintf(fid, 'pdf("%s/gam_smooth_terms.pdf", width=12, height=8)\n', strrep(output_dir, '\', '/'));
    fprintf(fid, 'par(mfrow=c(2,2))\n');
    fprintf(fid, 'plot(model, pages=1, residuals=TRUE, pch=1, cex=0.5)\n');
    fprintf(fid, 'dev.off()\n\n');
    fprintf(fid, '# Predictions over time for different alpha power values\n');
    fprintf(fid, 'alpha_vals <- seq(min(dat$AlphaPower), max(dat$AlphaPower), length=50)\n');
    fprintf(fid, 'time_vals <- seq(0, 2, length=100)\n');
    fprintf(fid, 'pred_grid <- expand.grid(Time=time_vals, AlphaPower=alpha_vals, ID=levels(dat$ID)[1])\n');
    fprintf(fid, 'pred_grid$SPL <- predict(model, newdata=pred_grid, type="response")\n\n');
    fprintf(fid, '# Plot predictions\n');
    fprintf(fid, 'p <- ggplot(pred_grid, aes(x=Time, y=SPL, color=AlphaPower)) +\n');
    fprintf(fid, '  geom_line(size=1) +\n');
    fprintf(fid, '  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +\n');
    fprintf(fid, '  labs(x="Time [s]", y="Scan Path Length [px]", color="Alpha Power") +\n');
    fprintf(fid, '  theme_bw()\n');
    fprintf(fid, 'ggsave("%s/gam_predictions.png", p, width=10, height=6, dpi=300)\n', strrep(output_dir, '\', '/'));
    
    fclose(fid);
    
    fprintf('Created R script: %s\n', r_script);
    fprintf('\n=== R GAM Model Formula ===\n');
    fprintf('SPL ~ s(Time) + s(AlphaPower) + te(Time, AlphaPower) + s(ID, bs="re")\n');
    fprintf('This models:\n');
    fprintf('  - Smooth effect of time\n');
    fprintf('  - Smooth effect of alpha power\n');
    fprintf('  - Tensor product (interaction) of time × alpha power\n');
    fprintf('  - Random intercepts per subject\n');
end

fprintf('\n=== Done ===\n');
