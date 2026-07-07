function plotGazeHeatmap(data, titleString, savePath)
% Common config
cfg               = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'outline';
cfg.zlim          = 'maxabs';
cfg.zlim          = [-.5 .5];
colMap            = customcolormap_preset('red-white-blue');
cfg.colormap      = colMap;
cfg.figure        = 'gcf';
overallFontSize   = 30;

% Set up figure
close all
figure;
set(gcf, 'Position', [0, 0, 1600, 1000], 'Color', 'W');
ft_singleplotTFR(cfg,data);
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
colormap(colMap);
colB = colorbar;
colB.LineWidth = 1;
colB.Ticks = [-.5 0 .5];
title(colB,'Effect size \itd')
hold on
centerX = 800 / 2;
centerY = 600 / 2;
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', overallFontSize);
title(titleString)

% Save
[~, paths] = setup('AOC', 0);
saveBasepath = fullfile(paths.figures, 'gaze', 'heatmap');
if ~isfolder(saveBasepath)
    mkdir(saveBasepath);
end
saveas(gcf, fullfile(saveBasepath, strcat(savePath, '.png')));