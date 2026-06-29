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
%   splits/AOC_split_AlphaAmpRed.m
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
fig6.panels = {
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_nback_timecourse.png'), ...
        'A', 'N-back ERS/ERD Time Course');
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_sternberg_timecourse.png'), ...
        'B', 'Sternberg ERS/ERD Time Course');
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_nback_topos.png'), ...
        'C', 'N-back ERS/ERD Topographies');
    panelSpec(fullfile(ersdDir, 'AOC_eeg_ersd_sternberg_topos.png'), ...
        'D', 'Sternberg ERS/ERD Topographies');
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

%% Figure 8: Exploratory alpha split gaze coupling (single panel)
fig8 = struct();
fig8.name = 'Figure8';
fig8.nrow = 1;
fig8.ncol = 1;
fig8.figSize = [0 0 1512 982];
fig8.panels = {
    panelSpec(fullfile(splitDir, ...
        'AOC_splitERSERD_timecourse_sternberg_split0_gaze_deviation_pct_CBPT.png'), ...
        'A', 'Sternberg Gaze Deviation by Alpha Response Direction');
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
function p = panelSpec(imagePath, letter, titleText)
p = struct('file', imagePath, 'letter', letter, 'title', titleText);
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
% Tighten vertical spacing for multi-row composites (e.g. Figures 5 and 6)
if spec.nrow >= 3
    tileSpacing = 'tight';
else
    tileSpacing = 'compact';
end
tl = tiledlayout(spec.nrow, spec.ncol, 'TileSpacing', tileSpacing, 'Padding', 'compact');

for i = 1:numel(spec.panels)
    p = spec.panels{i};
    ax = nexttile(tl);
    img = imread(p.file);
    image(ax, img);
    axis(ax, 'image');
    axis(ax, 'off');
    set(ax, 'Color', 'w');

    labelStr = sprintf('\\bf{%s} | %s', p.letter, p.title);
    title(ax, labelStr, 'Interpreter', 'tex', ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

drawnow;
exportgraphics(fig, outPng, 'Resolution', exportDpi, 'BackgroundColor', 'white');
close(fig);
end
