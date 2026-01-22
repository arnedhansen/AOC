%% AOC Omnibus — Behavioral Correlates (Alpha Power)
% Examines relationship between alpha power change and task performance 
% (accuracy, reaction time). Tests if alpha power differences relate to behavior.
%
% Key outputs:
%   Scatter plots: Alpha power vs accuracy/RT; correlation coefficients

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 20;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/behavioral';
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end

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
mean_acc = nan(n_subj, 1);
mean_rt = nan(n_subj, 1);

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
end

%% Filter valid data (exclude zeros and NaNs)
valid_mask = isfinite(mean_alpha) & mean_alpha ~= 0 & isfinite(mean_acc) & isfinite(mean_rt);

alpha_valid = mean_alpha(valid_mask);
acc_valid = mean_acc(valid_mask);
rt_valid = mean_rt(valid_mask);

% Group by reduction vs amplification
reduction_mask = alpha_valid < 0;
amplification_mask = alpha_valid > 0;

red_alpha = alpha_valid(reduction_mask);
amp_alpha = alpha_valid(amplification_mask);
red_acc = acc_valid(reduction_mask);
amp_acc = acc_valid(amplification_mask);
red_rt = rt_valid(reduction_mask);
amp_rt = rt_valid(amplification_mask);

%% Plot: Alpha Power vs Accuracy
close all

% Plot 1: Alpha Power vs Accuracy
figure
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');

subplot(1,2,1)
hold on
scatter(red_alpha, red_acc, 100, colors(1,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);
scatter(amp_alpha, amp_acc, 100, colors(3,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);

% Group-specific regression lines and correlations
valid_red = isfinite(red_alpha) & isfinite(red_acc);
r_red = nan;
p_red = nan;
if sum(valid_red) >= 2
    [r, p] = corrcoef(red_alpha(valid_red), red_acc(valid_red));
    r_red = r(1,2);
    p_red = p(1,2);
    
    p_fit_red = polyfit(red_alpha(valid_red), red_acc(valid_red), 1);
    x_fit_red = linspace(min(red_alpha(valid_red)), max(red_alpha(valid_red)), 100);
    y_fit_red = polyval(p_fit_red, x_fit_red);
    h_red = plot(x_fit_red, y_fit_red, '-', 'Color', colors(1,:), 'LineWidth', 2);
end

valid_amp = isfinite(amp_alpha) & isfinite(amp_acc);
r_amp = nan;
p_amp = nan;
if sum(valid_amp) >= 2
    [r, p] = corrcoef(amp_alpha(valid_amp), amp_acc(valid_amp));
    r_amp = r(1,2);
    p_amp = p(1,2);
    
    p_fit_amp = polyfit(amp_alpha(valid_amp), amp_acc(valid_amp), 1);
    x_fit_amp = linspace(min(amp_alpha(valid_amp)), max(amp_alpha(valid_amp)), 100);
    y_fit_amp = polyval(p_fit_amp, x_fit_amp);
    h_amp = plot(x_fit_amp, y_fit_amp, '-', 'Color', colors(3,:), 'LineWidth', 2);
end

% Title with group-specific correlations
title_str = 'Alpha Power vs Accuracy';
if isfinite(r_red) && isfinite(r_amp)
    title_str = sprintf('%s\nR: r=%.3f, p=%.4f | A: r=%.3f, p=%.4f', ...
        title_str, r_red, p_red, r_amp, p_amp);
elseif isfinite(r_red)
    title_str = sprintf('%s\nR: r=%.3f, p=%.4f', title_str, r_red, p_red);
elseif isfinite(r_amp)
    title_str = sprintf('%s\nA: r=%.3f, p=%.4f', title_str, r_amp, p_amp);
end
title(title_str, 'FontSize', fontSize, 'FontWeight', 'bold');

% Add vertical line at zero
xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

xlabel('Alpha Power Change from Baseline', 'FontSize', fontSize, 'FontWeight', 'bold')
ylabel('Accuracy [%]', 'FontSize', fontSize, 'FontWeight', 'bold')
if exist('h_red', 'var') && exist('h_amp', 'var')
    legend([h_red, h_amp], {'Reduction', 'Amplification'}, ...
        'Location', 'best', 'FontSize', fontSize-4)
elseif exist('h_red', 'var')
    legend(h_red, {'Reduction'}, 'Location', 'best', 'FontSize', fontSize-4)
elseif exist('h_amp', 'var')
    legend(h_amp, {'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
end
grid on
set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
box on

% Plot 2: Alpha Power vs RT
subplot(1,2,2)
hold on
scatter(red_alpha, red_rt, 100, colors(1,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);
scatter(amp_alpha, amp_rt, 100, colors(3,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.995);

% Group-specific regression lines and correlations
valid_red = isfinite(red_alpha) & isfinite(red_rt);
r_red = nan;
p_red = nan;
if sum(valid_red) >= 2
    [r, p] = corrcoef(red_alpha(valid_red), red_rt(valid_red));
    r_red = r(1,2);
    p_red = p(1,2);
    
    p_fit_red = polyfit(red_alpha(valid_red), red_rt(valid_red), 1);
    x_fit_red = linspace(min(red_alpha(valid_red)), max(red_alpha(valid_red)), 100);
    y_fit_red = polyval(p_fit_red, x_fit_red);
    h_red2 = plot(x_fit_red, y_fit_red, '-', 'Color', colors(1,:), 'LineWidth', 2);
end

valid_amp = isfinite(amp_alpha) & isfinite(amp_rt);
r_amp = nan;
p_amp = nan;
if sum(valid_amp) >= 2
    [r, p] = corrcoef(amp_alpha(valid_amp), amp_rt(valid_amp));
    r_amp = r(1,2);
    p_amp = p(1,2);
    
    p_fit_amp = polyfit(amp_alpha(valid_amp), amp_rt(valid_amp), 1);
    x_fit_amp = linspace(min(amp_alpha(valid_amp)), max(amp_alpha(valid_amp)), 100);
    y_fit_amp = polyval(p_fit_amp, x_fit_amp);
    h_amp2 = plot(x_fit_amp, y_fit_amp, '-', 'Color', colors(3,:), 'LineWidth', 2);
end

% Title with group-specific correlations
title_str = 'Alpha Power vs Reaction Time';
if isfinite(r_red) && isfinite(r_amp)
    title_str = sprintf('%s\nR: r=%.3f, p=%.4f | A: r=%.3f, p=%.4f', ...
        title_str, r_red, p_red, r_amp, p_amp);
elseif isfinite(r_red)
    title_str = sprintf('%s\nR: r=%.3f, p=%.4f', title_str, r_red, p_red);
elseif isfinite(r_amp)
    title_str = sprintf('%s\nA: r=%.3f, p=%.4f', title_str, r_amp, p_amp);
end
title(title_str, 'FontSize', fontSize, 'FontWeight', 'bold');

% Add vertical line at zero
xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

xlabel('Alpha Power Change from Baseline', 'FontSize', fontSize, 'FontWeight', 'bold')
ylabel('Reaction Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold')
if exist('h_red2', 'var') && exist('h_amp2', 'var')
    legend([h_red2, h_amp2], {'Reduction', 'Amplification'}, ...
        'Location', 'best', 'FontSize', fontSize-4)
elseif exist('h_red2', 'var')
    legend(h_red2, {'Reduction'}, 'Location', 'best', 'FontSize', fontSize-4)
elseif exist('h_amp2', 'var')
    legend(h_amp2, {'Amplification'}, 'Location', 'best', 'FontSize', fontSize-4)
end
grid on
set(gca, 'GridAlpha', 0.2, 'FontSize', fontSize-2)
box on

sgtitle('Behavioral Correlates of Alpha Power', ...
    'FontSize', fontSize+4, 'FontWeight', 'bold')

saveas(gcf, fullfile(output_dir, 'AOC_omnibus_behavioral_correlates_Alpha.png'));
fprintf('\nSaved figure: %s\n', fullfile(output_dir, 'AOC_omnibus_behavioral_correlates_Alpha.png'));

% Print summary
fprintf('\n=== Behavioral Summary ===\n');
fprintf('Reduction group: n=%d, Alpha=%.4f, Acc=%.1f%%, RT=%.3fs\n', ...
    length(red_alpha), nanmean(red_alpha), nanmean(red_acc), nanmean(red_rt));
fprintf('Amplification group: n=%d, Alpha=%.4f, Acc=%.1f%%, RT=%.3fs\n', ...
    length(amp_alpha), nanmean(amp_alpha), nanmean(amp_acc), nanmean(amp_rt));

% Group-specific correlations summary
fprintf('\n=== Group-Specific Correlations ===\n');
fprintf('Reduction group:\n');
if sum(isfinite(red_alpha) & isfinite(red_acc)) >= 2
    [r, p] = corrcoef(red_alpha(isfinite(red_alpha) & isfinite(red_acc)), ...
        red_acc(isfinite(red_alpha) & isfinite(red_acc)));
    fprintf('  Alpha Power × Accuracy: r=%.3f, p=%.4f\n', r(1,2), p(1,2));
end
if sum(isfinite(red_alpha) & isfinite(red_rt)) >= 2
    [r, p] = corrcoef(red_alpha(isfinite(red_alpha) & isfinite(red_rt)), ...
        red_rt(isfinite(red_alpha) & isfinite(red_rt)));
    fprintf('  Alpha Power × Reaction Time: r=%.3f, p=%.4f\n', r(1,2), p(1,2));
end
fprintf('Amplification group:\n');
if sum(isfinite(amp_alpha) & isfinite(amp_acc)) >= 2
    [r, p] = corrcoef(amp_alpha(isfinite(amp_alpha) & isfinite(amp_acc)), ...
        amp_acc(isfinite(amp_alpha) & isfinite(amp_acc)));
    fprintf('  Alpha Power × Accuracy: r=%.3f, p=%.4f\n', r(1,2), p(1,2));
end
if sum(isfinite(amp_alpha) & isfinite(amp_rt)) >= 2
    [r, p] = corrcoef(amp_alpha(isfinite(amp_alpha) & isfinite(amp_rt)), ...
        amp_rt(isfinite(amp_alpha) & isfinite(amp_rt)));
    fprintf('  Alpha Power × Reaction Time: r=%.3f, p=%.4f\n', r(1,2), p(1,2));
end

fprintf('\n=== Done ===\n');
