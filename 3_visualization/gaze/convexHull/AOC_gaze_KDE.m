%% AOC KDE Iso-Contour Visualization — Sternberg & N-Back
% Loads dataET for both tasks, computes 2D kernel density estimates of gaze
% positions per condition in the full retention window [0 2] s. Draws
% iso-contours enclosing 50% and 90% of the density mass. Computes the
% area within the 90% contour per subject x condition.
%
% Produces:
%   1) 2D screen plots with overlaid KDE contours per condition
%   2) Bar plots of mean 90%-contour area ± SEM across conditions
%
% Key outputs:
%   AOC_gaze_KDE_sternberg.png       — contour overlay
%   AOC_gaze_KDE_nback.png           — contour overlay
%   AOC_gaze_KDE_bar_sternberg.png   — bar plot area by condition
%   AOC_gaze_KDE_bar_nback.png       — bar plot area by condition

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

% KDE grid resolution
nGrid    = 200;
xEdges   = linspace(0, screenW, nGrid+1);
yEdges   = linspace(0, screenH, nGrid+1);
xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
[Xg, Yg] = meshgrid(xCenters, yCenters);
dx = xCenters(2) - xCenters(1);
dy = yCenters(2) - yCenters(1);
cellArea = dx * dy;  % area of each grid cell in px²

% Contour proportions
prop_inner = 0.50;   % 50% iso-contour
prop_outer = 0.90;   % 90% iso-contour

figDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/convexHull/';
mkdir(figDir);

numSubjects = length(subjects);
overallFontSize = 20;

%% Helper: find density threshold enclosing a given proportion of mass
find_threshold = @(density, prop) local_find_threshold(density, prop);

%% ====================================================================
%  STERNBERG
%  ====================================================================
condLabels = {'WM load 2', 'WM load 4', 'WM load 6'};
condCodes  = [22, 24, 26];
nConds     = length(condCodes);

kde90_area_stern = nan(numSubjects, nConds);
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

        % Per-subject KDE to get 90% contour area
        if numel(all_x) >= 50
            try
                dens = ksdensity([double(all_x(:)), double(all_y(:))], ...
                    [Xg(:), Yg(:)], 'Function', 'pdf');
                dens = reshape(dens, size(Xg));
                thresh90 = find_threshold(dens, prop_outer);
                kde90_area_stern(subj, c) = sum(dens(:) >= thresh90) * cellArea;
            catch
                kde90_area_stern(subj, c) = NaN;
            end
        end

        gaze_pool_x{c} = [gaze_pool_x{c}, all_x];
        gaze_pool_y{c} = [gaze_pool_y{c}, all_y];
    end
    clc; disp(['KDE Sternberg — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 1: Sternberg — Overlaid KDE contours (grand-pooled)
close all
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
hold on;

legendHandles = gobjects(1, nConds);

for c = 1:nConds
    gx = double(gaze_pool_x{c}(:));
    gy = double(gaze_pool_y{c}(:));
    if numel(gx) < 50; continue; end

    dens = ksdensity([gx, gy], [Xg(:), Yg(:)], 'Function', 'pdf');
    dens = reshape(dens, size(Xg));

    thresh50 = find_threshold(dens, prop_inner);
    thresh90 = find_threshold(dens, prop_outer);

    % 90% contour (outer, dashed)
    contour(Xg, Yg, dens, [thresh90 thresh90], '--', ...
        'LineColor', colors(c,:), 'LineWidth', 1.5);
    % 50% contour (inner, solid)
    [~, h] = contour(Xg, Yg, dens, [thresh50 thresh50], '-', ...
        'LineColor', colors(c,:), 'LineWidth', 2.5);
    legendHandles(c) = h;
end

plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');
xlim([-20 screenW+20]); ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('Sternberg — KDE Contours (50% solid, 90% dashed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(legendHandles, condLabels, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal; xlim([-20 screenW+20]); ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_KDE_sternberg.png'));

%% Plot 2: Sternberg — Violin + box + dot plot of 90%-contour area
close all
plot_violin_box(kde90_area_stern, condLabels, colors, ...
    'KDE 90% Contour Area [px^2]', 'Sternberg — KDE 90% Area by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_KDE_bar_sternberg.png'));

%% ====================================================================
%  N-BACK
%  ====================================================================
condLabels_nb = {'1-back', '2-back', '3-back'};
condCodes_nb  = [21, 22, 23];
nConds_nb     = length(condCodes_nb);

kde90_area_nb = nan(numSubjects, nConds_nb);
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

        if numel(all_x) >= 50
            try
                dens = ksdensity([double(all_x(:)), double(all_y(:))], ...
                    [Xg(:), Yg(:)], 'Function', 'pdf');
                dens = reshape(dens, size(Xg));
                thresh90 = find_threshold(dens, prop_outer);
                kde90_area_nb(subj, c) = sum(dens(:) >= thresh90) * cellArea;
            catch
                kde90_area_nb(subj, c) = NaN;
            end
        end

        gaze_pool_x_nb{c} = [gaze_pool_x_nb{c}, all_x];
        gaze_pool_y_nb{c} = [gaze_pool_y_nb{c}, all_y];
    end
    clc; disp(['KDE N-Back — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 3: N-Back — Overlaid KDE contours
close all
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
hold on;

legendHandles_nb = gobjects(1, nConds_nb);

for c = 1:nConds_nb
    gx = double(gaze_pool_x_nb{c}(:));
    gy = double(gaze_pool_y_nb{c}(:));
    if numel(gx) < 50; continue; end

    dens = ksdensity([gx, gy], [Xg(:), Yg(:)], 'Function', 'pdf');
    dens = reshape(dens, size(Xg));

    thresh50 = find_threshold(dens, prop_inner);
    thresh90 = find_threshold(dens, prop_outer);

    contour(Xg, Yg, dens, [thresh90 thresh90], '--', ...
        'LineColor', colors(c,:), 'LineWidth', 1.5);
    [~, h] = contour(Xg, Yg, dens, [thresh50 thresh50], '-', ...
        'LineColor', colors(c,:), 'LineWidth', 2.5);
    legendHandles_nb(c) = h;
end

plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');
xlim([-20 screenW+20]); ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('N-Back — KDE Contours (50% solid, 90% dashed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(legendHandles_nb, condLabels_nb, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal; xlim([-20 screenW+20]); ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_KDE_nback.png'));

%% Plot 4: N-Back — Violin + box + dot plot of 90%-contour area
close all
plot_violin_box(kde90_area_nb, condLabels_nb, colors, ...
    'KDE 90% Contour Area [px^2]', 'N-Back — KDE 90% Area by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_KDE_bar_nback.png'));

%% Print summary statistics
fprintf('\n=== Sternberg: KDE 90%% Contour Area (px²) ===\n');
for c = 1:nConds
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels{c}, ...
        nanmean(kde90_area_stern(:,c)), ...
        nanstd(kde90_area_stern(:,c)) / sqrt(sum(~isnan(kde90_area_stern(:,c)))));
end
fprintf('\n=== N-Back: KDE 90%% Contour Area (px²) ===\n');
for c = 1:nConds_nb
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels_nb{c}, ...
        nanmean(kde90_area_nb(:,c)), ...
        nanstd(kde90_area_nb(:,c)) / sqrt(sum(~isnan(kde90_area_nb(:,c)))));
end

%% ====================================================================
%  LOCAL FUNCTION
%  ====================================================================
function thresh = local_find_threshold(dens, prop)
%LOCAL_FIND_THRESHOLD  Find density threshold enclosing a given proportion.
    vals = sort(dens(:), 'descend');
    cumMass = cumsum(vals) / sum(vals);
    idx = find(cumMass >= prop, 1, 'first');
    if isempty(idx)
        thresh = min(vals);
    else
        thresh = vals(idx);
    end
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
