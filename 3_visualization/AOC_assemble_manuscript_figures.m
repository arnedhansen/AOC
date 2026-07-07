%% AOC Assemble Manuscript Results Figures (RR-S2)
% Combines exported panel PNGs into publication-ready composite figures
% matching the Results section of the Stage-2 manuscript (Figures 3 to 8).
%
% Prerequisite outputs (run these first if files are missing):
%   4_stats/AOC_stats_glmms_rainclouds.py
%   3_visualization/gaze/deviation/AOC_gaze_dev_*.m
%   3_visualization/gaze/microsaccades/AOC_gaze_microsaccades_*.m
%   3_visualization/gaze/heatmap/AOC_gaze_heatmap_*_CBPT_outline.m
%   3_visualization/eeg/ersd/AOC_eeg_ersd_*.m
%   splits/AOC_split_AlphaAmpRed_GazeDev.m
%   splits/AOC_split_AlphaAmpRed_MS.m
%
% Output: figures/manuscript/X.png

%% Setup
startup
[~, paths, ~, ~] = setup('AOC', 0);

exportDpi = 600;
outDir = fullfile(paths.figures, 'manuscript');
if ~isfolder(outDir)
    mkdir(outDir);
end

rainDir = fullfile(paths.figures, 'stats', 'rainclouds');
gazeDevDir = fullfile(paths.figures, 'gaze', 'deviation');
msDir = fullfile(paths.figures, 'gaze', 'microsaccades');
heatmapDir = fullfile(paths.figures, 'gaze', 'heatmap');
ersdDir = fullfile(paths.figures, 'eeg', 'ersd');
splitDir = fullfile(paths.figures, 'splits', 'SplitERSERD');

fprintf('Assembling manuscript figures...\n');
fprintf('Output directory: %s\n', outDir);

%% Figure 3: Behavioral rainclouds (2 x 2)
fig3 = struct();
fig3.name = 'Figure3';
fig3.nrow = 2;
fig3.ncol = 2;
fig3.figSize = [0 0 1512 982];
fig3.panels = {
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_rt_nback.png'), ...
        'A', 'N-back Reaction Time');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_rt_sternberg.png'), ...
        'B', 'Sternberg Reaction Time');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_acc_nback.png'), ...
        'C', 'N-back Accuracy');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_acc_sternberg.png'), ...
        'D', 'Sternberg Accuracy');
    };

%% Figure 4: Gaze deviation and microsaccade rate rainclouds (2 x 2)
fig4 = struct();
fig4.name = 'Figure4';
fig4.nrow = 2;
fig4.ncol = 2;
fig4.figSize = [0 0 1512 982];
fig4.panels = {
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_gazedev_nback.png'), ...
        'A', 'N-back Gaze Deviation');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_gazedev_sternberg.png'), ...
        'B', 'Sternberg Gaze Deviation');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ms_nback.png'), ...
        'C', 'N-back Microsaccades');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ms_sternberg.png'), ...
        'D', 'Sternberg Microsaccades');
    };

%% Figure 5: Gaze time courses and density maps (3 x 2)
fig5 = struct();
fig5.name = 'Figure5';
fig5.nrow = 3;
fig5.ncol = 2;
fig5.figSize = [0 0 1512 round(982 * 1.5)];
fig5.panels = {
    panelSpec(fullfile(gazeDevDir, 'AOC_gaze_dev_timecourse_nback.png'), ...
        'A', 'N-back Gaze Deviation Time Course');
    panelSpec(fullfile(gazeDevDir, 'AOC_gaze_dev_timecourse_sternberg.png'), ...
        'B', 'Sternberg Gaze Deviation Time Course');
    panelSpec(fullfile(msDir, 'AOC_gaze_microsaccades_timecourse_nback.png'), ...
        'C', 'N-back Microsaccades Time Course');
    panelSpec(fullfile(msDir, 'AOC_gaze_microsaccades_timecourse_sternberg.png'), ...
        'D', 'Sternberg Microsaccades Time Course');
    panelSpec(fullfile(heatmapDir, 'AOC_gaze_heatmap_nback_CBPT_outline_stat_full_d.png'), ...
        'E', 'N-back Gaze Density Heatmap');
    panelSpec(fullfile(heatmapDir, 'AOC_gaze_heatmap_sternberg_CBPT_outline_stat_full_d.png'), ...
        'F', 'Sternberg Gaze Density Heatmap');
    };

%% Figure 6: ERSD time courses, topographies, rainclouds (3 x 2)
fig6 = struct();
fig6.name = 'Figure6';
fig6.nrow = 3;
fig6.ncol = 2;
fig6.figSize = [0 0 1512 round(982 * 1.5)];
% Row 2 (topoplots) is set to NaN so its height is derived automatically
% from the actual (white-border-trimmed) image content at column width,
% instead of a manually guessed fraction. This removes the dead space
% above/below the topoplots without shrinking them, while rows 1 and 3
% share whatever vertical space remains.
fig6.rowHeights = [1, NaN, 1];
fig6.panels = {
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_nback_timecourse.png'), ...
        'A', 'N-back ERS/ERD Time Course');
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_sternberg_timecourse.png'), ...
        'B', 'Sternberg ERS/ERD Time Course');
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_nback_topos.png'), ...
        'C', 'N-back ERS/ERD Topographies', 1.00, true);
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_sternberg_topos.png'), ...
        'D', 'Sternberg ERS/ERD Topographies', 1.00, true);
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ersd_nback.png'), ...
        'E', 'N-back ERS/ERD');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ersd_sternberg.png'), ...
        'F', 'Sternberg ERS/ERD');
    };

%% Figure 7: Baselined gaze metrics rainclouds (2 x 2)
fig7 = struct();
fig7.name = 'Figure7';
fig7.nrow = 2;
fig7.ncol = 2;
fig7.figSize = [0 0 1512 982];
fig7.panels = {
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_gazedev_bl_nback.png'), ...
        'A', 'N-back Baselined Gaze Deviation');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_gazedev_bl_sternberg.png'), ...
        'B', 'Sternberg Baselined Gaze Deviation');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ms_bl_nback.png'), ...
        'C', 'N-back Baselined Microsaccades');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_ms_bl_sternberg.png'), ...
        'D', 'Sternberg Baselined Microsaccades');
    };

%% Figure 8: Exploratory alpha split gaze coupling (2 x 2)
fig8 = struct();
fig8.name = 'Figure8';
fig8.nrow = 2;
fig8.ncol = 2;
fig8.figSize = [0 0 1512 982];
fig8.panels = {
    panelSpec(fullfile(splitDir, ...
        'AOC_splitERSERD_GazeDev_timecourse_nback_splitMedian_GazeDev_pct_CBPT.png'), ...
        'A', 'N-back Gaze Deviation Split');
    panelSpec(fullfile(splitDir, ...
        'AOC_splitERSERD_GazeDev_timecourse_sternberg_split0_GazeDev_pct_CBPT.png'), ...
        'B', 'Sternberg Gaze Deviation Split');
    panelSpec(fullfile(splitDir, ...
        'AOC_splitERSERD_MS_timecourse_nback_splitMedian_MS_pct_CBPT.png'), ...
        'C', 'N-back Microsaccades Split');
    panelSpec(fullfile(splitDir, ...
        'AOC_splitERSERD_MS_timecourse_sternberg_split0_MS_pct_CBPT.png'), ...
        'D', 'Sternberg Microsaccades Split');
    };

%% Assemble all figures
figSpecs = {fig3, fig4, fig5, fig6, fig7, fig8};
for iFig = 1:numel(figSpecs)
    spec = figSpecs{iFig};
    outPng = fullfile(outDir, ['AOC_manuscript_' spec.name '.png']);
    assembleManuscriptFigure(spec, outPng, exportDpi);
    fprintf('Saved %s\n', outPng);
end

fprintf('\nDone. Manuscript composites saved to:\n  %s\n', outDir);

%% Local functions
function p = panelSpec(imagePath, letter, titleText, imageScale, trimWhite)
p = struct('file', imagePath, 'letter', letter, 'title', titleText);
if nargin >= 4 && ~isempty(imageScale)
    p.imageScale = imageScale;
end
if nargin >= 5 && ~isempty(trimWhite)
    p.trimWhite = trimWhite;
end
end

function assembleManuscriptFigure(spec, outPng, exportDpi)
missing = {};
for i = 1:numel(spec.panels)
    if ~isfile(spec.panels{i}.file)
        missing{end + 1} = spec.panels{i}.file; %#ok<AGROW>
    end
end
if ~isempty(missing)
    fprintf('\nMissing panel files for %s:\n', spec.name);
    for i = 1:numel(missing)
        fprintf('  %s\n', missing{i});
    end
    error('AOC_assemble_manuscript_figures:MissingPanels', ...
        'Cannot assemble %s until all panel PNGs exist.', spec.name);
end

fig = figure('Position', spec.figSize, 'Color', 'w');

if isfield(spec, 'rowHeights') && numel(spec.rowHeights) == spec.nrow
    layout = manuscriptPanelLayout(spec, spec.figSize);
    for i = 1:numel(spec.panels)
        p = spec.panels{i};
        slot = layout.slots(i);
        renderManuscriptPanel(fig, p, slot);
    end
else
    % Tighten vertical spacing for multi-row composites (e.g. Figure 5)
    if spec.nrow >= 3
        tileSpacing = 'tight';
    else
        tileSpacing = 'compact';
    end
    tl = tiledlayout(spec.nrow, spec.ncol, 'TileSpacing', tileSpacing, 'Padding', 'compact');

    for i = 1:numel(spec.panels)
        p = spec.panels{i};
        ax = nexttile(tl);
        renderPanelImage(ax, p);
    end
end

drawnow;
exportgraphics(fig, outPng, 'Resolution', exportDpi, 'BackgroundColor', 'white');
close(fig);
end

function layout = manuscriptPanelLayout(spec, figSize)
pad = struct('left', 0.04, 'right', 0.04, 'top', 0.02, 'bottom', 0.02);
gapX = 0.04;
gapY = 0.012;
titleH = 0.035;

usableW = 1 - pad.left - pad.right;
usableH = 1 - pad.top - pad.bottom;
colW = (usableW - (spec.ncol - 1) * gapX) / spec.ncol;

titleTotal = spec.nrow * titleH;
gapTotal = (spec.nrow - 1) * gapY;
imgAreaH = usableH - titleTotal - gapTotal;

rowWeights = spec.rowHeights(:)';
autoRows = find(isnan(rowWeights));
manualRows = find(~isnan(rowWeights));

rowImgH = zeros(1, spec.nrow);
figAspect = figSize(3) / figSize(4);
for r = autoRows
    panelIdx = (r - 1) * spec.ncol + 1;
    p = spec.panels{panelIdx};
    img = loadPanelImage(p);
    aspect = size(img, 2) / size(img, 1);
    imageScale = 1;
    if isfield(p, 'imageScale') && ~isempty(p.imageScale)
        imageScale = p.imageScale;
    end
    rowImgH(r) = colW * imageScale * figAspect / aspect;
end

remainingH = imgAreaH - sum(rowImgH(autoRows));
manualWeights = rowWeights(manualRows);
manualWeights = manualWeights / sum(manualWeights);
rowImgH(manualRows) = manualWeights * remainingH;

rowTops = zeros(1, spec.nrow);
yCursor = 1 - pad.top;
for r = 1:spec.nrow
    rowTops(r) = yCursor - titleH - rowImgH(r);
    yCursor = rowTops(r) - gapY;
end

emptySlot = struct('colLeft', 0, 'colW', 0, 'titlePos', zeros(1, 4), 'imgPos', zeros(1, 4));
layout = struct('slots', repmat(emptySlot, 1, numel(spec.panels)));
for i = 1:numel(spec.panels)
    row = ceil(i / spec.ncol);
    col = mod(i - 1, spec.ncol) + 1;

    layout.slots(i).colLeft = pad.left + (col - 1) * (colW + gapX);
    layout.slots(i).colW = colW;
    layout.slots(i).titlePos = [layout.slots(i).colLeft, rowTops(row) + rowImgH(row), colW, titleH];
    layout.slots(i).imgPos = [layout.slots(i).colLeft, rowTops(row), colW, rowImgH(row)];
end
end

function renderManuscriptPanel(fig, p, slot)
imageScale = 1;
if isfield(p, 'imageScale') && ~isempty(p.imageScale)
    imageScale = p.imageScale;
end

img = loadPanelImage(p);
imgPos = slot.imgPos;

if imageScale < 1
    % Row height already matches this content's aspect ratio (set in
    % manuscriptPanelLayout), so only the width needs horizontal
    % centering here; no vertical centering/gap is introduced.
    scaledW = slot.colW * imageScale;
    imgPos = [ ...
        slot.colLeft + (slot.colW - scaledW) / 2, ...
        slot.imgPos(2), ...
        scaledW, slot.imgPos(4)];
end

axTitle = axes(fig, 'Position', slot.titlePos, 'Visible', 'off', 'Color', 'w'); %#ok<LAXES>
labelStr = sprintf('\\bf{%s} | %s', p.letter, p.title);
title(axTitle, labelStr, 'Interpreter', 'tex', ...
    'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

ax = axes(fig, 'Position', imgPos, 'Color', 'w'); %#ok<LAXES>
renderPanelImage(ax, p, false, img);
end

function renderPanelImage(ax, p, includeTitle, img)
if nargin < 3
    includeTitle = true;
end
if nargin < 4
    img = loadPanelImage(p);
end
image(ax, img);
axis(ax, 'image');
axis(ax, 'off');
set(ax, 'Color', 'w');

if includeTitle
    labelStr = sprintf('\\bf{%s} | %s', p.letter, p.title);
    title(ax, labelStr, 'Interpreter', 'tex', ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
end

function img = loadPanelImage(p)
img = imread(p.file);
if isfield(p, 'trimWhite') && p.trimWhite
    img = trimWhiteBorders(img);
end
end

function imgOut = trimWhiteBorders(imgIn)
imgOut = imgIn;
if isempty(imgIn)
    return;
end

if size(imgIn, 3) == 1
    gray = imgIn;
else
    gray = rgb2gray(imgIn(:, :, 1:3));
end

% Keep non-background content, then crop with a small safety padding.
mask = gray < 248;
if ~any(mask(:))
    return;
end

[rows, cols] = find(mask);
r1 = min(rows);
r2 = max(rows);
c1 = min(cols);
c2 = max(cols);

pad = 6;
r1 = max(1, r1 - pad);
r2 = min(size(imgIn, 1), r2 + pad);
c1 = max(1, c1 - pad);
c2 = min(size(imgIn, 2), c2 + pad);

imgOut = imgIn(r1:r2, c1:c2, :);
end
