%% AOC BCEA Visualization — Sternberg & N-Back
% Loads dataET for both tasks, computes Bivariate Contour Ellipse Area
% (BCEA) per subject x condition in the full retention window [0 2] s.
% BCEA = 2 * k * pi * sigma_x * sigma_y * sqrt(1 - rho^2)
%   where k = -log(1 - P) for proportion P (68% and 95%).
% Reference: Crossland & Rubin (2002); Castet & Crossland (2012).
%
% Produces:
%   1) 2D screen plots with overlaid 68% & 95% confidence ellipses per condition
%   2) Bar plots of mean BCEA ± SEM across conditions (with subject dots)
%
% Key outputs:
%   AOC_gaze_BCEA_sternberg.png       — ellipse overlay
%   AOC_gaze_BCEA_nback.png           — ellipse overlay
%   AOC_gaze_BCEA_bar_sternberg.png   — bar plot BCEA by condition
%   AOC_gaze_BCEA_bar_nback.png       — bar plot BCEA by condition

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

screenW  = 800;
screenH  = 600;
centreX  = 400;
centreY  = 300;
blink_win = 50;
bounds_x = [0 screenW];
bounds_y = [0 screenH];
t_full   = [0 2];

P_inner = 0.68;  % inner ellipse proportion
P_outer = 0.95;  % outer ellipse proportion
k_inner = -log(1 - P_inner);
k_outer = -log(1 - P_outer);

figDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/convexHull/';
mkdir(figDir);

numSubjects = length(subjects);
overallFontSize = 20;
theta = linspace(0, 2*pi, 200);  % for drawing ellipses

%% Helper: compute BCEA and ellipse parameters from gaze x,y
%  Returns bcea_68, bcea_95, ellipse struct (for plotting)
compute_bcea = @(gx, gy) local_bcea(double(gx), double(gy), k_inner, k_outer, theta);

%% ====================================================================
%  STERNBERG
%  ====================================================================
condLabels = {'WM load 2', 'WM load 4', 'WM load 6'};
condCodes  = [22, 24, 26];
nConds     = length(condCodes);

bcea68_subj_stern = nan(numSubjects, nConds);
bcea95_subj_stern = nan(numSubjects, nConds);

% Store grand-pooled gaze per condition
gaze_pool_x = cell(1, nConds);
gaze_pool_y = cell(1, nConds);
for c = 1:nConds; gaze_pool_x{c} = []; gaze_pool_y{c} = []; end

for subj = 1:numSubjects
    gazePath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(gazePath, 'dataET_sternberg.mat'));

    for c = 1:nConds
        condIdx = find(dataETlong.trialinfo(:,1) == condCodes(c));
        all_x = []; all_y = [];

        for ti = 1:length(condIdx)
            trl = condIdx(ti);
            raw_dat = dataETlong.trial{trl}(1:3, :);
            t = dataETlong.time{trl};
            raw_dat(2,:) = screenH - raw_dat(2,:);
            inb = raw_dat(1,:) >= bounds_x(1) & raw_dat(1,:) <= bounds_x(2) & ...
                  raw_dat(2,:) >= bounds_y(1) & raw_dat(2,:) <= bounds_y(2);
            raw_dat(:, ~inb) = NaN;
            raw_dat = remove_blinks(raw_dat, blink_win);
            idx_full = t >= t_full(1) & t <= t_full(2);
            x = raw_dat(1, idx_full); y = raw_dat(2, idx_full);
            valid = isfinite(x) & isfinite(y);
            all_x = [all_x, x(valid)]; all_y = [all_y, y(valid)];
        end

        if numel(all_x) >= 10
            [bcea68_subj_stern(subj,c), bcea95_subj_stern(subj,c)] = ...
                compute_bcea(all_x, all_y);
        end

        gaze_pool_x{c} = [gaze_pool_x{c}, all_x];
        gaze_pool_y{c} = [gaze_pool_y{c}, all_y];
    end
    clc; disp(['BCEA Sternberg — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 1: Sternberg — Overlaid 68% & 95% BCEA ellipses
close all
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
hold on;

legendHandles = gobjects(1, nConds);

for c = 1:nConds
    gx = double(gaze_pool_x{c}); gy = double(gaze_pool_y{c});
    if numel(gx) < 10; continue; end

    mx = mean(gx); my = mean(gy);
    sx = std(gx);  sy = std(gy);
    rho = corr(gx(:), gy(:));
    cov_mat = [sx^2, rho*sx*sy; rho*sx*sy, sy^2];
    [V, D] = eig(cov_mat);

    for ki = [k_outer, k_inner]  % draw outer first (behind)
        r = sqrt(2 * ki * diag(D));  % radii
        ell = V * [r(1)*cos(theta); r(2)*sin(theta)];
        ex = ell(1,:) + mx;
        ey = ell(2,:) + my;

        if ki == k_outer
            fill(ex, ey, colors(c,:), 'FaceAlpha', 0.08, ...
                'EdgeColor', colors(c,:), 'LineWidth', 1.5, 'LineStyle', '--');
        else
            fill(ex, ey, colors(c,:), 'FaceAlpha', 0.20, ...
                'EdgeColor', colors(c,:), 'LineWidth', 2.5, 'LineStyle', '-');
            legendHandles(c) = plot(ex, ey, '-', 'Color', colors(c,:), 'LineWidth', 2.5);
        end
    end
end

plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');
xlim([-20 screenW+20]); ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('Sternberg — BCEA (68% solid, 95% dashed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(legendHandles, condLabels, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal; xlim([-20 screenW+20]); ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_sternberg.png'));

%% Plot 2: Sternberg — Violin + box + dot plot of BCEA (68%) per condition
close all
plot_violin_box(bcea68_subj_stern, condLabels, colors, ...
    'BCEA 68% [px^2]', 'Sternberg — BCEA (68%) by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_bar68_sternberg.png'));

%% Plot 2b: Sternberg — Violin + box + dot plot of BCEA (95%) per condition
close all
plot_violin_box(bcea95_subj_stern, condLabels, colors, ...
    'BCEA 95% [px^2]', 'Sternberg — BCEA (95%) by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_bar95_sternberg.png'));

%% ====================================================================
%  N-BACK
%  ====================================================================
condLabels_nb = {'1-back', '2-back', '3-back'};
condCodes_nb  = [21, 22, 23];
nConds_nb     = length(condCodes_nb);

bcea68_subj_nb = nan(numSubjects, nConds_nb);
bcea95_subj_nb = nan(numSubjects, nConds_nb);
gaze_pool_x_nb = cell(1, nConds_nb);
gaze_pool_y_nb = cell(1, nConds_nb);
for c = 1:nConds_nb; gaze_pool_x_nb{c} = []; gaze_pool_y_nb{c} = []; end

for subj = 1:numSubjects
    gazePath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(gazePath, 'dataET_nback.mat'));

    for c = 1:nConds_nb
        condIdx = find(dataETlong.trialinfo(:,1) == condCodes_nb(c));
        all_x = []; all_y = [];

        for ti = 1:length(condIdx)
            trl = condIdx(ti);
            raw_dat = dataETlong.trial{trl}(1:3, :);
            t = dataETlong.time{trl};
            raw_dat(2,:) = screenH - raw_dat(2,:);
            inb = raw_dat(1,:) >= bounds_x(1) & raw_dat(1,:) <= bounds_x(2) & ...
                  raw_dat(2,:) >= bounds_y(1) & raw_dat(2,:) <= bounds_y(2);
            raw_dat(:, ~inb) = NaN;
            raw_dat = remove_blinks(raw_dat, blink_win);
            idx_full = t >= t_full(1) & t <= t_full(2);
            x = raw_dat(1, idx_full); y = raw_dat(2, idx_full);
            valid = isfinite(x) & isfinite(y);
            all_x = [all_x, x(valid)]; all_y = [all_y, y(valid)];
        end

        if numel(all_x) >= 10
            [bcea68_subj_nb(subj,c), bcea95_subj_nb(subj,c)] = ...
                compute_bcea(all_x, all_y);
        end

        gaze_pool_x_nb{c} = [gaze_pool_x_nb{c}, all_x];
        gaze_pool_y_nb{c} = [gaze_pool_y_nb{c}, all_y];
    end
    clc; disp(['BCEA N-Back — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 3: N-Back — Overlaid BCEA ellipses
close all
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
hold on;

legendHandles_nb = gobjects(1, nConds_nb);

for c = 1:nConds_nb
    gx = double(gaze_pool_x_nb{c}); gy = double(gaze_pool_y_nb{c});
    if numel(gx) < 10; continue; end

    mx = mean(gx); my = mean(gy);
    sx = std(gx);  sy = std(gy);
    rho = corr(gx(:), gy(:));
    cov_mat = [sx^2, rho*sx*sy; rho*sx*sy, sy^2];
    [V, D] = eig(cov_mat);

    for ki = [k_outer, k_inner]
        r = sqrt(2 * ki * diag(D));
        ell = V * [r(1)*cos(theta); r(2)*sin(theta)];
        ex = ell(1,:) + mx; ey = ell(2,:) + my;

        if ki == k_outer
            fill(ex, ey, colors(c,:), 'FaceAlpha', 0.08, ...
                'EdgeColor', colors(c,:), 'LineWidth', 1.5, 'LineStyle', '--');
        else
            fill(ex, ey, colors(c,:), 'FaceAlpha', 0.20, ...
                'EdgeColor', colors(c,:), 'LineWidth', 2.5, 'LineStyle', '-');
            legendHandles_nb(c) = plot(ex, ey, '-', 'Color', colors(c,:), 'LineWidth', 2.5);
        end
    end
end

plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');
xlim([-20 screenW+20]); ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('N-Back — BCEA (68% solid, 95% dashed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(legendHandles_nb, condLabels_nb, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal; xlim([-20 screenW+20]); ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_nback.png'));

%% Plot 4: N-Back — Violin + box + dot plot of BCEA (68%) per condition
close all
plot_violin_box(bcea68_subj_nb, condLabels_nb, colors, ...
    'BCEA 68% [px^2]', 'N-Back — BCEA (68%) by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_bar68_nback.png'));

%% Plot 4b: N-Back — Violin + box + dot plot of BCEA (95%) per condition
close all
plot_violin_box(bcea95_subj_nb, condLabels_nb, colors, ...
    'BCEA 95% [px^2]', 'N-Back — BCEA (95%) by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_BCEA_bar95_nback.png'));

%% Print summary statistics
fprintf('\n=== Sternberg: BCEA 68%% (px²) ===\n');
for c = 1:nConds
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels{c}, ...
        nanmean(bcea68_subj_stern(:,c)), ...
        nanstd(bcea68_subj_stern(:,c)) / sqrt(sum(~isnan(bcea68_subj_stern(:,c)))));
end
fprintf('\n=== N-Back: BCEA 68%% (px²) ===\n');
for c = 1:nConds_nb
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels_nb{c}, ...
        nanmean(bcea68_subj_nb(:,c)), ...
        nanstd(bcea68_subj_nb(:,c)) / sqrt(sum(~isnan(bcea68_subj_nb(:,c)))));
end

%% ====================================================================
%  LOCAL FUNCTION
%  ====================================================================
function [bcea68, bcea95, ell_struct] = local_bcea(gx, gy, k68, k95, ~)
    sx  = std(gx, 'omitnan');
    sy  = std(gy, 'omitnan');
    rho = corr(gx(:), gy(:), 'rows', 'complete');
    denom = sqrt(1 - rho^2);
    bcea68 = 2 * k68 * pi * sx * sy * denom;
    bcea95 = 2 * k95 * pi * sx * sy * denom;
    ell_struct = struct('mx', mean(gx), 'my', mean(gy), 'sx', sx, 'sy', sy, 'rho', rho);
end

function plot_violin_box(data, condLabels, colors, yLabel, titleStr, fs)
    nConds = size(data, 2);
    all_vals = data(:); all_vals = all_vals(~isnan(all_vals));
    Q1g = prctile(all_vals, 25); Q3g = prctile(all_vals, 75);
    IQRg = Q3g - Q1g;
    lo_fence = Q1g - 1.5 * IQRg; hi_fence = Q3g + 1.5 * IQRg;

    figure; set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w'); hold on;
    vw = 0.30; bxw = 0.18; dot_off = 0.08;

    for c = 1:nConds
        v = data(:,c); v = v(~isnan(v)); if isempty(v); continue; end
        [f, xi] = ksdensity(v, 'NumPoints', 100);
        f = f / max(f) * vw;
        fill([c - f, fliplr(repmat(c, 1, length(f)))], [xi, fliplr(xi)], ...
            colors(c,:), 'FaceAlpha', 0.30, 'EdgeColor', colors(c,:), 'LineWidth', 1);
        q1 = prctile(v,25); q3 = prctile(v,75); med = median(v);
        p5 = prctile(v,5); p95 = prctile(v,95);
        plot([c c], [p5 q1], '-k', 'LineWidth', 1.2);
        plot([c c], [q3 p95], '-k', 'LineWidth', 1.2);
        plot(c+[-bxw/2,bxw/2], [p5 p5], '-k', 'LineWidth', 1.2);
        plot(c+[-bxw/2,bxw/2], [p95 p95], '-k', 'LineWidth', 1.2);
        rectangle('Position', [c-bxw/2, q1, bxw, q3-q1], ...
            'FaceColor', [colors(c,:), 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
        plot(c+[-bxw/2,bxw/2], [med med], '-k', 'LineWidth', 2);
        rng(42); n = numel(v); jit = dot_off + 0.12*(rand(n,1)-0.5);
        is_out = v < lo_fence | v > hi_fence;
        scatter(c+jit(~is_out), v(~is_out), 25, colors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
        if any(is_out)
            scatter(c+jit(is_out), v(is_out), 50, 'MarkerEdgeColor', colors(c,:), ...
                'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        end
    end
    inliers = all_vals(all_vals >= lo_fence & all_vals <= hi_fence);
    if ~isempty(inliers)
        rng_y = max(inliers) - min(inliers);
        ylim([max(0, min(inliers)-0.15*rng_y), max(inliers)+0.15*rng_y]);
    end
    set(gca, 'XTick', 1:nConds, 'XTickLabel', condLabels, 'FontSize', fs);
    ylabel(yLabel, 'FontSize', fs); title(titleStr, 'FontSize', fs+2);
    box off; set(gca, 'FontSize', fs); hold off;
end
