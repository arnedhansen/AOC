%% AOC Convex Hull Visualization (95th-Percentile Trimmed) — Sternberg & N-Back
% Loads dataET for both tasks, trims gaze to 95th percentile of distance
% from centroid, then computes per-subject convex hulls in the full
% retention window [0 2] s per condition. Produces:
%   1) 2D screen plots with overlaid grand-average convex hulls per condition
%   2) Bar plots of mean CHA ± SEM across conditions (with subject dots)
%
% Key outputs:
%   AOC_gaze_convexHull_sternberg.png  — spatial hull overlay
%   AOC_gaze_convexHull_nback.png      — spatial hull overlay
%   AOC_gaze_convexHull_bar_sternberg.png — bar plot CHA by condition
%   AOC_gaze_convexHull_bar_nback.png     — bar plot CHA by condition

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
pctile_keep = 95;  % keep gaze samples within this percentile of centroid distance

figDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/convexHull/';
mkdir(figDir);

numSubjects = length(subjects);
overallFontSize = 20;

%% ====================================================================
%  STERNBERG
%  ====================================================================
%% Load & compute per-subject convex hulls — Sternberg
condLabels = {'WM load 2', 'WM load 4', 'WM load 6'};
condCodes  = [22, 24, 26];
nConds     = length(condCodes);

% Store per-subject CHA and hull vertices for each condition
cha_subj_stern = nan(numSubjects, nConds);            % area per subj x cond
hull_x_stern   = cell(numSubjects, nConds);           % hull vertices
hull_y_stern   = cell(numSubjects, nConds);
gaze_pool_x    = cell(1, nConds);                     % pooled gaze for grand hull
gaze_pool_y    = cell(1, nConds);
for c = 1:nConds; gaze_pool_x{c} = []; gaze_pool_y{c} = []; end

for subj = 1:numSubjects
    gazePath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(gazePath, 'dataET_sternberg.mat'));

    nTrials = size(dataETlong.trialinfo, 1);

    for c = 1:nConds
        % Collect all valid gaze samples for this subject x condition
        condIdx = find(dataETlong.trialinfo(:,1) == condCodes(c));
        all_x = [];
        all_y = [];

        for ti = 1:length(condIdx)
            trl = condIdx(ti);
            raw_dat = dataETlong.trial{trl}(1:3, :);
            t = dataETlong.time{trl};

            % Invert Y
            raw_dat(2,:) = screenH - raw_dat(2,:);

            % Filter out-of-bounds
            inb = raw_dat(1,:) >= bounds_x(1) & raw_dat(1,:) <= bounds_x(2) & ...
                  raw_dat(2,:) >= bounds_y(1) & raw_dat(2,:) <= bounds_y(2);
            raw_dat(:, ~inb) = NaN;

            % Remove blinks
            raw_dat = remove_blinks(raw_dat, blink_win);

            % Full window
            idx_full = t >= t_full(1) & t <= t_full(2);
            x = raw_dat(1, idx_full);
            y = raw_dat(2, idx_full);

            % Keep valid
            valid = isfinite(x) & isfinite(y);
            all_x = [all_x, x(valid)];
            all_y = [all_y, y(valid)];
        end

        % 95th-percentile trimming from gaze centroid
        if ~isempty(all_x)
            cx = mean(double(all_x)); cy = mean(double(all_y));
            d = sqrt((double(all_x) - cx).^2 + (double(all_y) - cy).^2);
            keep = d <= prctile(d, pctile_keep);
            all_x = all_x(keep);
            all_y = all_y(keep);
        end

        % Compute convex hull for this subject x condition (on trimmed data)
        upts = unique([all_x(:), all_y(:)], 'rows');
        if size(upts, 1) >= 3
            try
                [k, cha_subj_stern(subj, c)] = convhull(double(all_x(:)), double(all_y(:)));
                hull_x_stern{subj, c} = all_x(k);
                hull_y_stern{subj, c} = all_y(k);
            catch
                cha_subj_stern(subj, c) = NaN;
            end
        end

        % Pool trimmed data for grand-average hull
        gaze_pool_x{c} = [gaze_pool_x{c}, all_x];
        gaze_pool_y{c} = [gaze_pool_y{c}, all_y];
    end

    clc
    disp(['Sternberg — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 1: Sternberg — Overlaid grand-average convex hulls on screen
close all
figure;
set(gcf, 'position', [0, 0, 1512, 982], 'color', 'w');
hold on;

hullHandles = gobjects(1, nConds);

for c = 1:nConds
    xp = gaze_pool_x{c};
    yp = gaze_pool_y{c};

    upts = unique([xp(:), yp(:)], 'rows');
    if size(upts, 1) >= 3
        try
            k = convhull(double(xp(:)), double(yp(:)));

            % Fill hull
            fill(xp(k), yp(k), colors(c,:), ...
                'FaceAlpha', 0.15, 'EdgeColor', colors(c,:), ...
                'LineWidth', 2.5, 'LineStyle', '-');
            hullHandles(c) = plot(xp(k), yp(k), '-', ...
                'Color', colors(c,:), 'LineWidth', 2.5);
        catch
            hullHandles(c) = plot(NaN, NaN, '-', 'Color', colors(c,:));
        end
    end
end

% Screen centre marker
plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');

% Screen boundary
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');

% Cosmetics
xlim([-20 screenW+20]);
ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('Sternberg — Convex Hull (95th Percentile Trimmed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(hullHandles, condLabels, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal
xlim([-20 screenW+20]);
ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_sternberg.png'));

%% Plot 2: Sternberg — Violin + box + dot plot of CHA per condition
close all
plot_violin_box(cha_subj_stern, condLabels, colors, ...
    'Convex Hull Area [px^2]', ...
    'Sternberg — Convex Hull Area by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_bar_sternberg.png'));

%% Plot 3: Sternberg — Individual subject hulls overlaid per condition
close all
figure;
set(gcf, 'Position', [100, 100, 1800, 500], 'Color', 'w');

for c = 1:nConds
    subplot(1, nConds, c);
    hold on;

    for subj = 1:numSubjects
        if ~isempty(hull_x_stern{subj, c})
            fill(hull_x_stern{subj, c}, hull_y_stern{subj, c}, colors(c,:), ...
                'FaceAlpha', 0.08, 'EdgeColor', colors(c,:), ...
                'EdgeAlpha', 0.3, 'LineWidth', 0.8);
        end
    end

    % Screen centre
    plot(centreX, centreY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
    rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1, 'LineStyle', '--');

    xlim([-20 screenW+20]);
    ylim([-20 screenH+20]);
    set(gca, 'YDir', 'normal');
    xlabel('X [px]');
    ylabel('Y [px]');
    title(condLabels{c}, 'FontSize', overallFontSize);
    set(gca, 'FontSize', overallFontSize - 4);
    axis equal
    xlim([-20 screenW+20]);
    ylim([-20 screenH+20]);
    hold off;
end

sgtitle('Sternberg — Individual Subject Convex Hulls', 'FontSize', overallFontSize + 2);
saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_individual_sternberg.png'));

%% ====================================================================
%  N-BACK
%  ====================================================================
%% Load & compute per-subject convex hulls — N-Back
condLabels_nb = {'1-back', '2-back', '3-back'};
condCodes_nb  = [21, 22, 23];
nConds_nb     = length(condCodes_nb);

cha_subj_nb = nan(numSubjects, nConds_nb);
hull_x_nb   = cell(numSubjects, nConds_nb);
hull_y_nb   = cell(numSubjects, nConds_nb);
gaze_pool_x_nb = cell(1, nConds_nb);
gaze_pool_y_nb = cell(1, nConds_nb);
for c = 1:nConds_nb; gaze_pool_x_nb{c} = []; gaze_pool_y_nb{c} = []; end

for subj = 1:numSubjects
    gazePath = fullfile(path, subjects{subj}, 'gaze');
    load(fullfile(gazePath, 'dataET_nback.mat'));

    nTrials = size(dataETlong.trialinfo, 1);

    for c = 1:nConds_nb
        condIdx = find(dataETlong.trialinfo(:,1) == condCodes_nb(c));
        all_x = [];
        all_y = [];

        for ti = 1:length(condIdx)
            trl = condIdx(ti);
            raw_dat = dataETlong.trial{trl}(1:3, :);
            t = dataETlong.time{trl};

            % Invert Y
            raw_dat(2,:) = screenH - raw_dat(2,:);

            % Filter out-of-bounds
            inb = raw_dat(1,:) >= bounds_x(1) & raw_dat(1,:) <= bounds_x(2) & ...
                  raw_dat(2,:) >= bounds_y(1) & raw_dat(2,:) <= bounds_y(2);
            raw_dat(:, ~inb) = NaN;

            % Remove blinks
            raw_dat = remove_blinks(raw_dat, blink_win);

            % Full window
            idx_full = t >= t_full(1) & t <= t_full(2);
            x = raw_dat(1, idx_full);
            y = raw_dat(2, idx_full);

            valid = isfinite(x) & isfinite(y);
            all_x = [all_x, x(valid)];
            all_y = [all_y, y(valid)];
        end

        % 95th-percentile trimming from gaze centroid
        if ~isempty(all_x)
            cx = mean(double(all_x)); cy = mean(double(all_y));
            d = sqrt((double(all_x) - cx).^2 + (double(all_y) - cy).^2);
            keep = d <= prctile(d, pctile_keep);
            all_x = all_x(keep);
            all_y = all_y(keep);
        end

        % Compute convex hull for this subject x condition (on trimmed data)
        upts = unique([all_x(:), all_y(:)], 'rows');
        if size(upts, 1) >= 3
            try
                [k, cha_subj_nb(subj, c)] = convhull(double(all_x(:)), double(all_y(:)));
                hull_x_nb{subj, c} = all_x(k);
                hull_y_nb{subj, c} = all_y(k);
            catch
                cha_subj_nb(subj, c) = NaN;
            end
        end

        % Pool trimmed data for grand-average hull
        gaze_pool_x_nb{c} = [gaze_pool_x_nb{c}, all_x];
        gaze_pool_y_nb{c} = [gaze_pool_y_nb{c}, all_y];
    end

    clc
    disp(['N-Back — Subject ' num2str(subj) '/' num2str(numSubjects) ' done.'])
end

%% Plot 4: N-Back — Overlaid grand-average convex hulls on screen
close all
figure;
set(gcf, 'position', [0, 0, 1512, 982], 'color', 'w');
hold on;

hullHandles_nb = gobjects(1, nConds_nb);

for c = 1:nConds_nb
    xp = gaze_pool_x_nb{c};
    yp = gaze_pool_y_nb{c};

    upts = unique([xp(:), yp(:)], 'rows');
    if size(upts, 1) >= 3
        try
            k = convhull(double(xp(:)), double(yp(:)));
            fill(xp(k), yp(k), colors(c,:), ...
                'FaceAlpha', 0.15, 'EdgeColor', colors(c,:), ...
                'LineWidth', 2.5, 'LineStyle', '-');
            hullHandles_nb(c) = plot(xp(k), yp(k), '-', ...
                'Color', colors(c,:), 'LineWidth', 2.5);
        catch
            hullHandles_nb(c) = plot(NaN, NaN, '-', 'Color', colors(c,:));
        end
    end
end

plot(centreX, centreY, '+', 'MarkerSize', 20, 'LineWidth', 2.5, 'Color', 'k');
rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 1, 'LineStyle', '--');

xlim([-20 screenW+20]);
ylim([-20 screenH+20]);
set(gca, 'YDir', 'normal');
xlabel('Screen X [px]', 'FontSize', overallFontSize);
ylabel('Screen Y [px]', 'FontSize', overallFontSize);
title('N-Back — Convex Hull (95th Percentile Trimmed)', 'FontSize', overallFontSize + 2);
set(gca, 'FontSize', overallFontSize);
legend(hullHandles_nb, condLabels_nb, 'FontSize', overallFontSize - 2, 'Location', 'northeast');
axis equal
xlim([-20 screenW+20]);
ylim([-20 screenH+20]);
hold off;

saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_nback.png'));

%% Plot 5: N-Back — Violin + box + dot plot of CHA per condition
close all
plot_violin_box(cha_subj_nb, condLabels_nb, colors, ...
    'Convex Hull Area [px^2]', ...
    'N-Back — Convex Hull Area by Condition', overallFontSize);
saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_bar_nback.png'));

%% Plot 6: N-Back — Individual subject hulls overlaid per condition
close all
figure;
set(gcf, 'Position', [100, 100, 1800, 500], 'Color', 'w');

for c = 1:nConds_nb
    subplot(1, nConds_nb, c);
    hold on;

    for subj = 1:numSubjects
        if ~isempty(hull_x_nb{subj, c})
            fill(hull_x_nb{subj, c}, hull_y_nb{subj, c}, colors(c,:), ...
                'FaceAlpha', 0.08, 'EdgeColor', colors(c,:), ...
                'EdgeAlpha', 0.3, 'LineWidth', 0.8);
        end
    end

    plot(centreX, centreY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
    rectangle('Position', [0, 0, screenW, screenH], 'EdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1, 'LineStyle', '--');

    xlim([-20 screenW+20]);
    ylim([-20 screenH+20]);
    set(gca, 'YDir', 'normal');
    xlabel('X [px]');
    ylabel('Y [px]');
    title(condLabels_nb{c}, 'FontSize', overallFontSize);
    set(gca, 'FontSize', overallFontSize - 4);
    axis equal
    xlim([-20 screenW+20]);
    ylim([-20 screenH+20]);
    hold off;
end

sgtitle('N-Back — Individual Subject Convex Hulls', 'FontSize', overallFontSize + 2);
saveas(gcf, fullfile(figDir, 'AOC_gaze_convexHull_individual_nback.png'));

%% Print summary statistics
fprintf('\n=== Sternberg: Convex Hull Area (px²) ===\n');
for c = 1:nConds
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels{c}, ...
        nanmean(cha_subj_stern(:,c)), ...
        nanstd(cha_subj_stern(:,c)) / sqrt(sum(~isnan(cha_subj_stern(:,c)))));
end

fprintf('\n=== N-Back: Convex Hull Area (px²) ===\n');
for c = 1:nConds_nb
    fprintf('  %s:  Mean = %.1f,  SEM = %.1f\n', condLabels_nb{c}, ...
        nanmean(cha_subj_nb(:,c)), ...
        nanstd(cha_subj_nb(:,c)) / sqrt(sum(~isnan(cha_subj_nb(:,c)))));
end

%% ====================================================================
%  LOCAL FUNCTION — Violin + Box + Dot plot with outlier handling
%  ====================================================================
function plot_violin_box(data, condLabels, colors, yLabel, titleStr, fs)
    nConds = size(data, 2);
    all_vals = data(:); all_vals = all_vals(~isnan(all_vals));
    Q1g = prctile(all_vals, 25); Q3g = prctile(all_vals, 75);
    IQRg = Q3g - Q1g;
    lo_fence = Q1g - 1.5 * IQRg; hi_fence = Q3g + 1.5 * IQRg;

    figure; set(gcf, 'Position', [100, 100, 800, 600], 'Color', 'w'); hold on;
    vw = 0.30; bw = 0.18; dot_off = 0.08;

    for c = 1:nConds
        v = data(:,c); v = v(~isnan(v)); if isempty(v); continue; end

        % Half-violin (left side)
        [f, xi] = ksdensity(v, 'NumPoints', 100);
        f = f / max(f) * vw;
        fill([c - f, fliplr(repmat(c, 1, length(f)))], [xi, fliplr(xi)], ...
            colors(c,:), 'FaceAlpha', 0.30, 'EdgeColor', colors(c,:), 'LineWidth', 1);

        % Box + whiskers (5th–95th)
        q1 = prctile(v,25); q3 = prctile(v,75); med = median(v);
        p5 = prctile(v,5); p95 = prctile(v,95);
        plot([c c], [p5 q1],  '-k', 'LineWidth', 1.2);
        plot([c c], [q3 p95], '-k', 'LineWidth', 1.2);
        plot(c + [-bw/2, bw/2], [p5 p5],   '-k', 'LineWidth', 1.2);
        plot(c + [-bw/2, bw/2], [p95 p95], '-k', 'LineWidth', 1.2);
        rectangle('Position', [c-bw/2, q1, bw, q3-q1], ...
            'FaceColor', [colors(c,:), 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
        plot(c + [-bw/2, bw/2], [med med], '-k', 'LineWidth', 2);

        % Dots — inliers filled, outliers open circles
        rng(42); n = numel(v);
        jit = dot_off + 0.12 * (rand(n,1) - 0.5);
        is_out = v < lo_fence | v > hi_fence;
        scatter(c + jit(~is_out), v(~is_out), 25, colors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
        if any(is_out)
            scatter(c + jit(is_out), v(is_out), 50, 'MarkerEdgeColor', colors(c,:), ...
                'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        end
    end

    % Y-axis clipped to inlier range
    inliers = all_vals(all_vals >= lo_fence & all_vals <= hi_fence);
    if ~isempty(inliers)
        rng_y = max(inliers) - min(inliers);
        yl = max(0, min(inliers) - 0.15*rng_y);
        yh = max(inliers) + 0.15*rng_y;
        ylim([yl, yh]);
    end
    set(gca, 'XTick', 1:nConds, 'XTickLabel', condLabels, 'FontSize', fs);
    ylabel(yLabel, 'FontSize', fs); title(titleStr, 'FontSize', fs+2);
    box off; set(gca, 'FontSize', fs); hold off;
end
