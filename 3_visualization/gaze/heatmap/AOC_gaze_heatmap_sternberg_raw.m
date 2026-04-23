%% AOC Gaze Heatmap - Sternberg (Raw, Non-Baselined)
% Raw heatmaps without baseline correction.
% Window: 1-2 s. Difference map: WM6 minus WM2.

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'gaze', 'heatmap');
if ~isfolder(figDir), mkdir(figDir); end

%% Parameters
latencyWindow = [1 2];
num_bins = 100;
smoothing_factor = 5;
x_edges = linspace(0, 800, num_bins);
y_edges = linspace(0, 600, num_bins);
x_centers = x_edges(2:end);
y_centers = y_edges(2:end);

%% Subject-level heatmaps using direct 2D binning
nSub = length(subjects);
heatmaps = nan(nSub, 3, num_bins - 1, num_bins - 1);

for subj = 1:nSub
    datapath = fullfile(path, subjects{subj}, 'gaze');
    load([datapath, filesep 'dataET_sternberg']);
    clc
    disp(upper(['Loading ET data for subject ' num2str(subj) '/' num2str(nSub) '...']))

    trialinfo = dataETlong.trialinfo;
    if size(trialinfo, 2) > 1
        trialinfo = trialinfo(:, 1);
    end

    condCodes = [22 24 26];
    cfg = [];
    cfg.avgovertime = 'no';
    cfg.keeptrials = 'yes';
    cfg.latency = latencyWindow;

    for condIdx = 1:3
        cfg.trials = find(trialinfo == condCodes(condIdx));
        dataSel = ft_selectdata(cfg, dataETlong);
        dataSel = horzcat(dataSel.trial{:});

        x_positions = dataSel(1, :);
        y_positions = dataSel(2, :);

        valid_idx = x_positions >= 0 & x_positions <= 800 & ...
                    y_positions >= 0 & y_positions <= 600 & ...
                    ~(x_positions == 0 & y_positions == 0);
        x_positions = x_positions(valid_idx);
        y_positions = y_positions(valid_idx);

        binned_data = histcounts2(x_positions, y_positions, x_edges, y_edges);
        heatmaps(subj, condIdx, :, :) = imgaussfilt(binned_data, smoothing_factor);
    end
end

%% Grand averages
grandMaps = squeeze(mean(heatmaps, 1, 'omitnan'));
map2 = squeeze(grandMaps(1, :, :));
map4 = squeeze(grandMaps(2, :, :));
map6 = squeeze(grandMaps(3, :, :));
diffMap = map6 - map2;

%% Plot settings
close all
overallFontSize = 35;
centerX = 400;
centerY = 300;
rawColMap = customcolormap_preset('white-red');
diffColMap = customcolormap_preset('red-white-blue');
robustMax = prctile([map2(:); map4(:); map6(:)], 99.5);
robustLim = prctile(abs(diffMap(:)), 99.9);

% plotRawMap(map2, x_centers, y_centers, rawColMap, [0 robustMax], centerX, centerY, overallFontSize, ...
%     'WM load 2 Gaze Heatmap (Raw, 1-2 s)', fullfile(figDir, 'AOC_gaze_heatmap_sternberg_raw_WM2.png'));
% plotRawMap(map4, x_centers, y_centers, rawColMap, [0 robustMax], centerX, centerY, overallFontSize, ...
%     'WM load 4 Gaze Heatmap (Raw, 1-2 s)', fullfile(figDir, 'AOC_gaze_heatmap_sternberg_raw_WM4.png'));
% plotRawMap(map6, x_centers, y_centers, rawColMap, [0 robustMax], centerX, centerY, overallFontSize, ...
%     'WM load 6 Gaze Heatmap (Raw, 1-2 s)', fullfile(figDir, 'AOC_gaze_heatmap_sternberg_raw_WM6.png'));
plotRawMap(diffMap, x_centers, y_centers, diffColMap, [-robustLim robustLim], centerX, centerY, overallFontSize, ...
    'Difference Heatmap: WM6 - WM2 (Raw, 1-2 s)', fullfile(figDir, 'AOC_gaze_heatmap_sternberg_raw_diff.png'));

function plotRawMap(mapData, x_centers, y_centers, cmap, climVals, centerX, centerY, fontSize, panelTitle, outPath)
figure('Position', [0 0 1512 982], 'Color', 'w');
imagesc(x_centers, y_centers, mapData');
axis tight
xlim([0 800]);
ylim([0 600]);
yticks([0 150 300 450 600]);
yticklabels({'600','450','300','150','0'});
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(cmap);
clim(climVals);
cb = colorbar;
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', fontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', fontSize);
title(panelTitle, 'FontSize', 30);
saveas(gcf, outPath);
end
