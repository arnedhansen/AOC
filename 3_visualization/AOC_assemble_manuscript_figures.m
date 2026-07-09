%% AOC Assemble Manuscript Results Figures (AOC RR-S2)
% Combines exported panel PNGs into publication-ready composite figures
% matching the Stage-2 manuscript (Figures 1 to 9)

%% Setup
startup
[~, paths, ~, ~] = setup('AOC', 0);

exportDpi = 600;
outDir = fullfile(paths.figures, 'manuscript');
if ~isfolder(outDir)
    mkdir(outDir);
end

paradigmsDir = fullfile(paths.figures, 'paradigms');
powerDir = fullfile(paths.figures, 'power_analysis');
rainDir = fullfile(paths.figures, 'stats', 'rainclouds');
gazeDevDir = fullfile(paths.figures, 'gaze', 'deviation');
msDir = fullfile(paths.figures, 'gaze', 'microsaccades');
heatmapDir = fullfile(paths.figures, 'gaze', 'heatmap');
ersdDir = fullfile(paths.figures, 'eeg', 'ersd');
splitDir = fullfile(paths.figures, 'splits', 'SplitERSERD');

fprintf('Assembling manuscript figures...\n');
fprintf('Output directory: %s\n', outDir);

%% Figure 1: Task paradigms (2 x 1)
fig1 = struct();
fig1.name = 'Figure1';
fig1.nrow = 2;
fig1.ncol = 1;
fig1.figSize = [0 0 1512 982];
fig1.rowHeights = [1, 1];
fig1.panels = {
    panelSpec(fullfile(paradigmsDir, 'AOC_paradigms_nback.png'), ...
        'A', 'N-back Task Paradigm');
    panelSpec(fullfile(paradigmsDir, 'AOC_paradigms_sternberg.png'), ...
        'B', 'Sternberg Task Paradigm');
    };

%% Figure 2: Power analysis (1 x 3)
fig2 = struct();
fig2.name = 'Figure2';
fig2.nrow = 1;
fig2.ncol = 3;
fig2.figSize = [0 0 1512 982];
fig2GazePanel = panelSpec(fullfile(powerDir, 'AOC_power_analysis_gaze.png'), ...
    'C', 'Alpha ~ Gaze');
fig2GazePanel.titleInterpreter = 'none';
fig2.panels = {
    panelSpec(fullfile(powerDir, 'AOC_power_analysis_sternberg.png'), ...
        'A', 'Sternberg Power Analysis');
    panelSpec(fullfile(powerDir, 'AOC_power_analysis_nback.png'), ...
        'B', 'N-back Power Analysis');
    fig2GazePanel;
    };

%% Figure 3: Behavioral rainclouds (2 x 2)
fig3 = struct();
fig3.name = 'Figure3';
fig3.nrow = 2;
fig3.ncol = 2;
fig3.figSize = [0 0 1512 982];
fig3.rowHeights = [1, 1];
fig3.colGap = 0.012;
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
fig4.rowHeights = [1, 1];
fig4.colGap = 0.012;
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
fig5.rowHeights = [1, 1, NaN];
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
        'E', 'N-back Gaze Density Heatmap', 0.9);
    panelSpec(fullfile(heatmapDir, 'AOC_gaze_heatmap_sternberg_CBPT_outline_stat_full_d.png'), ...
        'F', 'Sternberg Gaze Density Heatmap', 0.9);
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
% Tighter gap below time courses; extra padding below topo row titles
fig6.rowGaps = [0.004, 0.012];
fig6.titleImageGap = [0, 0.012, 0];
fig6.titleShiftY = [-0.03, 0, 0];
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
fig7.rowHeights = [1, 1];
fig7.colGap = 0.012;
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

%% Figure 8: Baseline window power spectrum rainclouds (1 x 2)
fig8 = struct();
fig8.name = 'Figure8';
fig8.nrow = 1;
fig8.ncol = 2;
fig8.figSize = [0 0 1512 982];
fig8.panels = {
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_pow_baselineWindow_nback.png'), ...
        'A', 'N-back Baseline Window Alpha Power');
    panelSpec(fullfile(rainDir, 'AOC_stats_rainclouds_pow_baselineWindow_sternberg.png'), ...
        'B', 'Sternberg Baseline Window Alpha Power');
    };

%% Figure 9: Exploratory alpha split gaze coupling (2 x 2)
fig9 = struct();
fig9.name = 'Figure9';
fig9.nrow = 2;
fig9.ncol = 2;
fig9.figSize = [0 0 1512 982];
fig9.panels = {
    panelSpec(fullfile(splitDir, 'AOC_splitERSERD_timecourse_nback_GazeDev.png'), ...
        'A', 'N-back Gaze Deviation Split');
    panelSpec(fullfile(splitDir, 'AOC_splitERSERD_timecourse_sternberg_GazeDev.png'), ...
        'B', 'Sternberg Gaze Deviation Split');
    panelSpec(fullfile(splitDir, 'AOC_splitERSERD_timecourse_nback_MS.png'), ...
        'C', 'N-back Microsaccades Split');
    panelSpec(fullfile(splitDir, 'AOC_splitERSERD_timecourse_sternberg_MS.png'), ...
        'D', 'Sternberg Microsaccades Split');
    };

%% Assemble all figures
figSpecs = {fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9};
for iFig = 1:numel(figSpecs)
    spec = figSpecs{iFig};
    outPng = fullfile(outDir, ['AOC_manuscript_' spec.name '.png']);
    assembleManuscriptFigure(spec, outPng, exportDpi);
    fprintf('Saved %s\n', outPng);
end

fprintf('\nDone. Manuscript composites saved to:\n  %s\n', outDir);

%% Local functions
function p = panelSpec(imagePath, letter, titleText, imageScale, trimWhite, whiteBackground, hideTitle, titleInterpreter)
p = struct('file', imagePath, 'letter', letter, 'title', titleText);
if nargin >= 4 && ~isempty(imageScale)
    p.imageScale = imageScale;
end
if nargin >= 5 && ~isempty(trimWhite)
    p.trimWhite = trimWhite;
end
if nargin >= 6 && ~isempty(whiteBackground)
    p.whiteBackground = whiteBackground;
end
if nargin >= 7 && ~isempty(hideTitle)
    p.hideTitle = hideTitle;
end
if nargin >= 8 && ~isempty(titleInterpreter)
    p.titleInterpreter = titleInterpreter;
end
end

function assembleManuscriptFigure(spec, outPng, exportDpi)
missing = {};
for i = 1:numel(spec.panels)
    if ~isfile(spec.panels{i}.file)
        missing{end + 1} = spec.panels{i}.file;
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
    if isfield(spec, 'tileSpacing') && ~isempty(spec.tileSpacing)
        tileSpacing = spec.tileSpacing;
    elseif spec.nrow >= 3
        tileSpacing = 'tight';
    else
        tileSpacing = 'compact';
    end
    if ~ischar(tileSpacing) && ~isstring(tileSpacing)
        error('AOC_assemble_manuscript_figures:InvalidTileSpacing', ...
            'tileSpacing must be ''loose'', ''compact'', ''tight'', or ''none''.');
    end
    tileSpacing = char(tileSpacing);
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
if isfield(spec, 'colGap') && ~isempty(spec.colGap)
    gapX = spec.colGap;
else
    gapX = 0.04;
end
gapY = 0.012;
titleH = 0.035;

usableW = 1 - pad.left - pad.right;
usableH = 1 - pad.top - pad.bottom;
colW = (usableW - (spec.ncol - 1) * gapX) / spec.ncol;

if isfield(spec, 'rowGaps') && numel(spec.rowGaps) == spec.nrow - 1
    rowGaps = spec.rowGaps(:)';
else
    rowGaps = repmat(gapY, 1, spec.nrow - 1);
end

if isfield(spec, 'titleImageGap') && numel(spec.titleImageGap) == spec.nrow
    titleImageGap = spec.titleImageGap(:)';
else
    titleImageGap = zeros(1, spec.nrow);
end

if isfield(spec, 'titlePadTop') && numel(spec.titlePadTop) == spec.nrow
    titlePadTop = spec.titlePadTop(:)';
else
    titlePadTop = zeros(1, spec.nrow);
end

if isfield(spec, 'titleShiftY') && numel(spec.titleShiftY) == spec.nrow
    titleShiftY = spec.titleShiftY(:)';
else
    titleShiftY = zeros(1, spec.nrow);
end

titleTotal = spec.nrow * titleH;
gapTotal = sum(rowGaps);
titleImageGapTotal = sum(titleImageGap);
titlePadTopTotal = sum(titlePadTop);
imgAreaH = usableH - titleTotal - gapTotal - titleImageGapTotal - titlePadTopTotal;

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
    rowTops(r) = yCursor - titlePadTop(r) - titleH - titleImageGap(r) - rowImgH(r);
    if r < spec.nrow
        yCursor = rowTops(r) - rowGaps(r);
    end
end

emptySlot = struct('colLeft', 0, 'colW', 0, 'titlePos', zeros(1, 4), 'imgPos', zeros(1, 4));
layout = struct('slots', repmat(emptySlot, 1, numel(spec.panels)));
for i = 1:numel(spec.panels)
    row = ceil(i / spec.ncol);
    col = mod(i - 1, spec.ncol) + 1;

    layout.slots(i).colLeft = pad.left + (col - 1) * (colW + gapX);
    layout.slots(i).colW = colW;
    layout.slots(i).titlePos = [layout.slots(i).colLeft, ...
        rowTops(row) + rowImgH(row) + titleImageGap(row), colW, titleH];
    layout.slots(i).titleShiftY = titleShiftY(row);
    layout.slots(i).imgPos = [layout.slots(i).colLeft, rowTops(row), colW, rowImgH(row)];
    p = spec.panels{i};
    hasScale = isfield(p, 'imageScale') && ~isempty(p.imageScale) && p.imageScale ~= 1;
    layout.slots(i).layoutIncludesScale = ismember(row, autoRows) && hasScale;
end
end

function renderManuscriptPanel(fig, p, slot)
imageScale = 1;
if isfield(p, 'imageScale') && ~isempty(p.imageScale)
    imageScale = p.imageScale;
end

img = loadPanelImage(p);
imgPos = slot.imgPos;

if isfield(slot, 'layoutIncludesScale') && slot.layoutIncludesScale
    scaledW = slot.colW * imageScale;
    imgPos = [ ...
        slot.colLeft + (slot.colW - scaledW) / 2, ...
        slot.imgPos(2), ...
        scaledW, slot.imgPos(4)];
elseif imageScale ~= 1
    scaledW = slot.colW * imageScale;
    scaledH = slot.imgPos(4) * imageScale;
    imgPos = [ ...
        slot.colLeft + (slot.colW - scaledW) / 2, ...
        slot.imgPos(2) + (slot.imgPos(4) - scaledH) / 2, ...
        scaledW, scaledH];
end

labelStr = panelLabelStr(p);
titleInterpreter = panelTitleInterpreter(p);
titlePos = slot.titlePos;
if isfield(slot, 'titleShiftY') && ~isempty(slot.titleShiftY)
    titlePos(2) = titlePos(2) + slot.titleShiftY;
end
if ~(isfield(p, 'hideTitle') && p.hideTitle)
    annotation(fig, 'textbox', titlePos, ...
        'String', labelStr, ...
        'Interpreter', titleInterpreter, ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', ...
        'FitBoxToText', 'off', ...
        'Margin', 0);
end

ax = axes(fig, 'Position', imgPos, 'Color', 'w');
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

if includeTitle && ~(isfield(p, 'hideTitle') && p.hideTitle)
    labelStr = panelLabelStr(p);
    title(ax, labelStr, 'Interpreter', panelTitleInterpreter(p), ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
end

function labelStr = panelLabelStr(p)
if isfield(p, 'titleInterpreter') && strcmpi(p.titleInterpreter, 'none')
    labelStr = sprintf('%s | %s', p.letter, p.title);
else
    labelStr = sprintf('\\bf{%s} | %s', p.letter, p.title);
end
end

function interpreter = panelTitleInterpreter(p)
if isfield(p, 'titleInterpreter') && ~isempty(p.titleInterpreter)
    interpreter = p.titleInterpreter;
else
    interpreter = 'tex';
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
