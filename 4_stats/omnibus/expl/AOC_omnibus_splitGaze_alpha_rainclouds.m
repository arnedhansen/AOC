%% AOC — Split by Gaze Deviation / Scan Path Length → Alpha Power Rainclouds
% Splits subjects (between-subject) and trials (within-subject) by gaze
% deviation or scan path length using a median split, then visualises
% alpha power distributions as raincloud plots.
%
% PART 1 — Between-Subject Median Split (LONG data)
%   Loads merged_data_*_LONG.mat, computes mean gaze per subject (collapsed
%   across conditions), median-splits into low / high groups, and shows
%   alpha power rainclouds with independent-samples t-test + Cohen's d.
%
% PART 2 — Within-Subject Trial-Level Median Split (trials data)
%   Loads merged_data_*_trials.mat, median-splits trials within each
%   subject by trial-level gaze, computes subject-level means for low/high
%   groups, runs paired t-test + Cohen's dz. Produces two figures per
%   analysis: (A) subject-level means, (B) all-trials distribution.
%
% Split variables:  GazeDeviation / ScanPathLength  (full epoch, [0-2]s)
% Alpha DVs:        AlphaPower (raw, full epoch) + AlphaPower_FOOOF_bl (FOOOF, baselined)
%                   Trial-level: AlphaPowerFull (raw) + AlphaPower_FOOOF_bl (FOOOF)
% Tasks:            N-back, Sternberg
%
% Note on FOOOF alpha in the trials data:
%   FOOOF alpha is computed at the condition level (not trial level).
%   In the trials table it is repeated for every trial of a given
%   subject×condition. Within-subject trial splits using FOOOF alpha
%   therefore reflect condition-composition differences rather than
%   true trial-by-trial alpha variation.
%
% Key outputs:
%   24 raincloud PNG figures in output_dir
%   Console statistics for every analysis

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');
rng(12345);  % reproducible jitter across all figures

% Paths (cross-platform)
if ispc
    base_data  = 'W:\Students\Arne\AOC\data\features';
    output_dir = 'W:\Students\Arne\AOC\figures\interactions\omnibus_alpha_split\gazedev';
else
    base_data  = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features';
    output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/gazedev';
end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Figure config
fontSize = 20;
pal = [colors(1,:); colors(3,:)];  % blue → low gaze, red → high gaze

% ====================================================================
%% PART 1: BETWEEN-SUBJECT MEDIAN SPLIT (LONG DATA)
% ====================================================================
fprintf('\n%s\n  PART 1: Between-Subject Median Split (LONG data)\n%s\n', ...
    repmat('=',1,64), repmat('=',1,64));

taskCfgs = struct( ...
    'name',    {'nback',                      'sternberg'}, ...
    'file',    {'merged_data_nback_LONG.mat', 'merged_data_sternberg_LONG.mat'}, ...
    'varname', {'merged_data_nback_LONG',     'merged_data_sternberg_LONG'});

splitVars   = {'GazeDeviation', 'ScanPathLength'};
splitLabels = {'Gaze Deviation', 'Scan Path Length'};

% Alpha DVs (same for both tasks, full epoch):
alphaVars   = {'AlphaPower', 'AlphaPower_FOOOF_bl'};
alphaLabels = {'Alpha Power [\muV^2/Hz]', 'Alpha Power (FOOOF, BL)'};
alphaSaves  = {'alpha', 'alpha_fooof_bl'};

for ti = 1:numel(taskCfgs)
    fprintf('\n--- Task: %s ---\n', upper(taskCfgs(ti).name));

    S   = load(fullfile(base_data, taskCfgs(ti).file));
    dat = S.(taskCfgs(ti).varname);

    allIDs = [dat.ID];
    uIDs   = unique(allIDs);
    nSubj  = numel(uIDs);

    for si = 1:numel(splitVars)
        splitVar = splitVars{si};
        splitLab = splitLabels{si};

        % --- Collapse conditions: mean split variable per subject ---
        splitPerSubj = nan(nSubj, 1);
        for s = 1:nSubj
            mask = allIDs == uIDs(s);
            splitPerSubj(s) = mean([dat(mask).(splitVar)], 'omitnan');
        end

        % --- Median split ---
        validSplit = isfinite(splitPerSubj);
        medVal  = median(splitPerSubj(validSplit));
        lowIdx  = splitPerSubj <= medVal & validSplit;
        highIdx = splitPerSubj >  medVal & validSplit;

        fprintf('  Split by %s: median = %.2f, low n=%d, high n=%d\n', ...
            splitVar, medVal, sum(lowIdx), sum(highIdx));

        for ai = 1:numel(alphaVars)
            alphaVar  = alphaVars{ai};
            alphaLab  = alphaLabels{ai};
            alphaSave = alphaSaves{ai};

            % --- Collapse conditions: mean alpha per subject ---
            alphaPerSubj = nan(nSubj, 1);
            for s = 1:nSubj
                mask = allIDs == uIDs(s);
                alphaPerSubj(s) = mean([dat(mask).(alphaVar)], 'omitnan');
            end

            % Keep subjects with valid data in both split + alpha
            validBoth = lowIdx  & isfinite(alphaPerSubj) | ...
                        highIdx & isfinite(alphaPerSubj);
            alphaLow  = alphaPerSubj(lowIdx  & isfinite(alphaPerSubj));
            alphaHigh = alphaPerSubj(highIdx & isfinite(alphaPerSubj));

            if numel(alphaLow) < 3 || numel(alphaHigh) < 3
                warning('Too few subjects for %s × %s (%s). Skipping.', ...
                    splitVar, alphaVar, taskCfgs(ti).name);
                continue
            end

            % --- Statistics: two-sample t-test + Cohen's d ---
            [~, pVal, ~, stats] = ttest2(alphaLow, alphaHigh);
            pooled_sd = sqrt(((numel(alphaLow)-1)*var(alphaLow) + ...
                (numel(alphaHigh)-1)*var(alphaHigh)) / ...
                (numel(alphaLow) + numel(alphaHigh) - 2));
            cohens_d = (mean(alphaHigh) - mean(alphaLow)) / max(pooled_sd, eps);

            fprintf('    %s: t(%d)=%.3f, p=%.4f, d=%.3f  (low n=%d, high n=%d)\n', ...
                alphaVar, stats.df, stats.tstat, pVal, cohens_d, ...
                numel(alphaLow), numel(alphaHigh));

            % --- Raincloud figure ---
            close all
            fig = figure; set(fig, 'Color', 'w', 'Position', [0 0 1400 1200]);
            ax = gca; hold(ax, 'on');

            drawOneCloud(ax, alphaLow,  1, pal(1,:), -0.20, 0.40, 0.20, 30, 0.50);
            drawOneCloud(ax, alphaHigh, 2, pal(2,:), -0.20, 0.40, 0.20, 30, 0.50);

            % Y range + bracket
            allV = [alphaLow(:); alphaHigh(:)];
            yR   = max(allV) - min(allV);
            if yR == 0, yR = 1; end
            brkY = max(allV) + 0.08 * yR;
            ylim(ax, [min(allV) - 0.05*yR, brkY + 0.22*yR]);
            addSigBracket(ax, 1, 2, brkY, pVal, cohens_d, 'd', fontSize);

            % Labels
            set(ax, 'XTick', [1 2], 'XTickLabel', ...
                {sprintf('Low %s\n(n=%d)', splitLab, numel(alphaLow)), ...
                 sprintf('High %s\n(n=%d)', splitLab, numel(alphaHigh))});
            ylabel(ax, alphaLab, 'FontSize', fontSize);
            title(ax, sprintf('%s: %s split by %s (Between-Subject)', ...
                upper(taskCfgs(ti).name), strrep(alphaVar,'_',' '), splitLab), ...
                'FontSize', fontSize, 'FontWeight', 'bold');
            ax.YGrid = 'on'; ax.GridAlpha = 0.35; ax.Box = 'off';
            set(ax, 'FontSize', fontSize-2);
            xlim(ax, [0 3]);

            % Save
            splitSave = lower(strrep(splitVar, ' ', ''));
            fname = sprintf('AOC_splitGaze_%s_%s_%s_between.png', ...
                splitSave, alphaSave, taskCfgs(ti).name);
            saveas(fig, fullfile(output_dir, fname));
            fprintf('    Saved: %s\n', fname);
            close(fig);
        end
    end
end

% ====================================================================
%% PART 2: WITHIN-SUBJECT TRIAL-LEVEL MEDIAN SPLIT (TRIALS DATA)
% ====================================================================
fprintf('\n%s\n  PART 2: Within-Subject Trial-Level Median Split (trials data)\n%s\n', ...
    repmat('=',1,64), repmat('=',1,64));

taskCfgs_trl = struct( ...
    'name',    {'nback',                        'sternberg'}, ...
    'file',    {'merged_data_nback_trials.mat', 'merged_data_sternberg_trials.mat'}, ...
    'varname', {'merged_data_nback_trials',     'merged_data_sternberg_trials'});

splitVars_trl   = {'GazeDeviationFull', 'ScanPathLengthFull'};
splitLabels_trl = {'Gaze Deviation', 'Scan Path Length'};

% Alpha DVs (same for both tasks, full epoch):
alphaVars_trl   = {'AlphaPowerFull', 'AlphaPower_FOOOF_bl'};
alphaLabels_trl = {'Alpha Power [\muV^2/Hz]', 'Alpha Power (FOOOF, BL)'};
alphaSaves_trl  = {'alphafull', 'alpha_fooof_bl'};

for ti = 1:numel(taskCfgs_trl)
    fprintf('\n--- Task: %s ---\n', upper(taskCfgs_trl(ti).name));

    S = load(fullfile(base_data, taskCfgs_trl(ti).file));
    T = S.(taskCfgs_trl(ti).varname);                 % MATLAB table

    uIDs  = unique(T.ID);
    nSubj = numel(uIDs);

    for si = 1:numel(splitVars_trl)
        splitVar = splitVars_trl{si};
        splitLab = splitLabels_trl{si};

        for ai = 1:numel(alphaVars_trl)
            alphaVar  = alphaVars_trl{ai};
            alphaLab  = alphaLabels_trl{ai};
            alphaSave = alphaSaves_trl{ai};

            % --- Per-subject: median split trials, compute means ---
            subjMeanLow  = nan(nSubj, 1);
            subjMeanHigh = nan(nSubj, 1);
            allTrialLow  = [];
            allTrialHigh = [];
            validSubj    = 0;
            nLowTotal    = 0;
            nHighTotal   = 0;

            for s = 1:nSubj
                rows      = T(T.ID == uIDs(s), :);
                gazeVals  = rows.(splitVar);
                alphaVals = rows.(alphaVar);

                valid     = isfinite(gazeVals) & isfinite(alphaVals);
                gazeVals  = gazeVals(valid);
                alphaVals = alphaVals(valid);

                if numel(gazeVals) < 4, continue; end

                med      = median(gazeVals);
                lowMask  = gazeVals <= med;
                highMask = gazeVals >  med;

                if sum(lowMask) < 2 || sum(highMask) < 2, continue; end

                validSubj = validSubj + 1;
                subjMeanLow(validSubj)  = mean(alphaVals(lowMask));
                subjMeanHigh(validSubj) = mean(alphaVals(highMask));

                allTrialLow  = [allTrialLow;  alphaVals(lowMask)];   %#ok<AGROW>
                allTrialHigh = [allTrialHigh; alphaVals(highMask)];   %#ok<AGROW>
                nLowTotal    = nLowTotal  + sum(lowMask);
                nHighTotal   = nHighTotal + sum(highMask);
            end

            subjMeanLow  = subjMeanLow(1:validSubj);
            subjMeanHigh = subjMeanHigh(1:validSubj);

            if validSubj < 3
                warning('Too few subjects for %s × %s (%s). Skipping.', ...
                    splitVar, alphaVar, taskCfgs_trl(ti).name);
                continue
            end

            % --- Statistics: paired t-test + Cohen's dz ---
            [~, pVal, ~, stats] = ttest(subjMeanHigh, subjMeanLow);
            diffs = subjMeanHigh - subjMeanLow;
            dz    = mean(diffs) / max(std(diffs), eps);

            fprintf('  %s × %s: n=%d subj (%d + %d trials), t(%d)=%.3f, p=%.4f, dz=%.3f\n', ...
                splitVar, alphaVar, validSubj, nLowTotal, nHighTotal, ...
                stats.df, stats.tstat, pVal, dz);

            splitClean = strrep(splitVar, 'Full', '');
            splitSave  = lower(strrep(splitVar, ' ', ''));

            % ============================================================
            % Figure A: Subject-level means raincloud
            % ============================================================
            close all
            fig = figure; set(fig, 'Color', 'w', 'Position', [0 0 1500 1200]);
            ax = gca; hold(ax, 'on');

            drawOneCloud(ax, subjMeanLow,  1, pal(1,:), -0.20, 0.40, 0.20, 40, 0.50);
            drawOneCloud(ax, subjMeanHigh, 2, pal(2,:), -0.20, 0.40, 0.20, 40, 0.50);

            allV = [subjMeanLow; subjMeanHigh];
            yR   = max(allV) - min(allV);
            if yR == 0, yR = 1; end
            brkY = max(allV) + 0.08 * yR;
            ylim(ax, [min(allV) - 0.05*yR, brkY + 0.22*yR]);
            addSigBracket(ax, 1, 2, brkY, pVal, dz, 'd_z', fontSize);

            set(ax, 'XTick', [1 2], 'XTickLabel', ...
                {sprintf('Low %s', splitClean), ...
                 sprintf('High %s', splitClean)});
            ylabel(ax, alphaLab, 'FontSize', fontSize);
            title(ax, sprintf('%s: %s split by %s\n(Within-Subject, n=%d, Subject Means)', ...
                upper(taskCfgs_trl(ti).name), strrep(alphaVar,'_',' '), ...
                splitClean, validSubj), ...
                'FontSize', fontSize-2, 'FontWeight', 'bold');
            ax.YGrid = 'on'; ax.GridAlpha = 0.35; ax.Box = 'off';
            set(ax, 'FontSize', fontSize-2);
            xlim(ax, [0 3]);

            fname = sprintf('AOC_splitGaze_%s_%s_%s_within_subjMeans.png', ...
                splitSave, alphaSave, taskCfgs_trl(ti).name);
            saveas(fig, fullfile(output_dir, fname));
            fprintf('    Saved: %s\n', fname);
            close(fig);

            % ============================================================
            % Figure B: All-trials raincloud
            % ============================================================
            fig = figure; set(fig, 'Color', 'w', 'Position', [0 0 1500 1200]);
            ax = gca; hold(ax, 'on');

            drawOneCloud(ax, allTrialLow,  1, pal(1,:), -0.20, 0.40, 0.20, 8, 0.12);
            drawOneCloud(ax, allTrialHigh, 2, pal(2,:), -0.20, 0.40, 0.20, 8, 0.12);

            % Use robust range (1st–99th pctl) for y limits
            allV  = [allTrialLow; allTrialHigh];
            q_lo  = prctile(allV, 1);
            q_hi  = prctile(allV, 99);
            yR    = q_hi - q_lo;
            if yR == 0, yR = 1; end
            brkY  = q_hi + 0.08 * yR;
            ylim(ax, [q_lo - 0.05*yR, brkY + 0.22*yR]);
            addSigBracket(ax, 1, 2, brkY, pVal, dz, 'd_z', fontSize);

            set(ax, 'XTick', [1 2], 'XTickLabel', ...
                {sprintf('Low %s\n(%d trials)', splitClean, nLowTotal), ...
                 sprintf('High %s\n(%d trials)', splitClean, nHighTotal)});
            ylabel(ax, alphaLab, 'FontSize', fontSize);
            title(ax, sprintf('%s: %s split by %s\n(Within-Subject, n=%d, All Trials)', ...
                upper(taskCfgs_trl(ti).name), strrep(alphaVar,'_',' '), ...
                splitClean, validSubj), ...
                'FontSize', fontSize-2, 'FontWeight', 'bold');
            ax.YGrid = 'on'; ax.GridAlpha = 0.35; ax.Box = 'off';
            set(ax, 'FontSize', fontSize-2);
            xlim(ax, [0 3]);

            fname = sprintf('AOC_splitGaze_%s_%s_%s_within_trials.png', ...
                splitSave, alphaSave, taskCfgs_trl(ti).name);
            saveas(fig, fullfile(output_dir, fname));
            fprintf('    Saved: %s\n', fname);
            close(fig);
        end
    end
end

%% Summary
fprintf('\n%s\n  Done — all figures saved to:\n  %s\n%s\n', ...
    repmat('=',1,64), output_dir, repmat('=',1,64));

% ====================================================================
%% LOCAL FUNCTIONS
% ====================================================================

function drawOneCloud(ax, yvals, xpos, col, cloudOffset, maxViolW, boxW, dotSize, dotAlpha)
%DRAWONECLOUD  Half-violin + boxplot + jittered dots for one group.
%   Matches the Python raincloud style: KDE cloud to the LEFT of the
%   group position, centred boxplot, jittered scatter dots.

    yvals = yvals(isfinite(yvals));
    if numel(yvals) < 3, return; end

    % --- Half-violin (KDE) extending to the left ---
    [f, xi] = ksdensity(yvals, 'NumPoints', 200);
    if max(f) == 0, return; end
    f_scaled = f / max(f) * maxViolW;

    x_right = repmat(xpos + cloudOffset, size(xi));
    x_left  = xpos + cloudOffset - f_scaled;

    fill(ax, [x_right, fliplr(x_left)], [xi, fliplr(xi)], col, ...
        'FaceAlpha', 0.60, 'EdgeColor', 'none');

    % --- Boxplot (manual: 5th/95th whiskers) ---
    q = prctile(yvals, [5 25 50 75 95]);

    % Box body (Q1–Q3)
    patch(ax, xpos + boxW/2*[-1 1 1 -1], [q(2) q(2) q(4) q(4)], ...
        col, 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineWidth', 1);
    % Median line
    line(ax, xpos + boxW/2*[-1 1], [q(3) q(3)], ...
        'Color', 'k', 'LineWidth', 2);
    % Whiskers
    line(ax, [xpos xpos], [q(1) q(2)], 'Color', 'k', 'LineWidth', 1);
    line(ax, [xpos xpos], [q(4) q(5)], 'Color', 'k', 'LineWidth', 1);
    % Caps
    capW = boxW / 3;
    line(ax, xpos + capW*[-1 1]/2, [q(1) q(1)], 'Color', 'k', 'LineWidth', 1);
    line(ax, xpos + capW*[-1 1]/2, [q(5) q(5)], 'Color', 'k', 'LineWidth', 1);

    % --- Jittered dots ---
    jit = (rand(size(yvals)) - 0.5) * boxW;
    scatter(ax, xpos + jit, yvals, dotSize, col, 'filled', ...
        'MarkerFaceAlpha', dotAlpha, 'MarkerEdgeColor', 'none');
end


function addSigBracket(ax, x1, x2, yBrk, pVal, effSize, effLabel, fontSize)
%ADDSIGBRACKET  Draw a significance bracket with stars + effect size.

    yR   = diff(ylim(ax));
    head = 0.015 * yR;

    % Bracket lines
    line(ax, [x1 x1], [yBrk yBrk+head], 'Color', 'k', 'LineWidth', 1.5, 'Clipping', 'off');
    line(ax, [x1 x2], [yBrk+head yBrk+head], 'Color', 'k', 'LineWidth', 1.5, 'Clipping', 'off');
    line(ax, [x2 x2], [yBrk yBrk+head], 'Color', 'k', 'LineWidth', 1.5, 'Clipping', 'off');

    % Text label
    sigStr = sigLabel(pVal);
    txt = sprintf('%s  %s = %.2f', sigStr, effLabel, effSize);
    text(ax, (x1+x2)/2, yBrk + head*3, txt, ...
        'HorizontalAlignment', 'center', 'FontSize', fontSize*0.60, ...
        'FontWeight', 'bold', 'Clipping', 'off');
end


function s = sigLabel(p)
%SIGLABEL  Convert p-value to significance stars or 'n.s.'.
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'n.s.';
    end
end
