%% AOC Gaze Microsaccade Suppression — Sternberg & N-Back
% Detects microsaccades per trial using detect_microsaccades.m (Engbert &
% Kliegl, 2003). The velocity threshold is GLOBAL — computed over the full
% computation window passed to detect_microsaccades.
%
% Builds binary spike trains from detected onsets, convolves with a
% Gaussian kernel to produce smoothed microsaccade rate time courses, and
% plots per-condition rate curves with SEM shading and raster plots.
%
% Edge effects from the Gaussian smoothing kernel are handled by computing
% over a wider window (t_comp) and cropping to the display window (t_win).
%
% Key outputs (per task):
%   1. Smoothed MS rate time courses per WM load condition (± SEM)
%   1b. Percentage-change baseline-corrected rate time courses
%   2. Raster plots of microsaccade onsets per condition
%   3. Combined figure (rate + raster)
%   3b. Combined figure with percentage-change rate

%% Setup
startup
setup('AOC');
clear
close all
clc

datapath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
figpath  = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/gaze/microsaccades/';
mkdir(figpath);

dirs     = dir(datapath);
folders  = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');
nSubj    = length(subjects);

%% Parameters
fsample   = 500;        % Hz
screenW   = 800;
screenH   = 600;
blink_win = 50;         % blink removal window (samples, ±100 ms)

% Gaussian smoothing kernel for continuous rate estimation
sigma_ms   = 50;                                     % kernel SD in ms
sigma_samp = round(sigma_ms / (1000 / fsample));     % SD in samples
kHalf      = 3 * sigma_samp;
x_kern     = -kHalf : kHalf;
gKernel    = exp(-x_kern.^2 / (2 * sigma_samp^2));
gKernel    = gKernel / sum(gKernel);                 % normalise

% Computation window (wider than display to avoid convolution edge effects)
% Kernel half-width = 3*sigma = 150 ms; 500 ms padding per side is safe.
t_comp     = [-1.0 2.5];                                % seconds
n_comp     = round(diff(t_comp) * fsample) + 1;         % samples in computation window
t_comp_vec = linspace(t_comp(1), t_comp(2), n_comp);    % computation time axis

% Display time window (what is shown in figures)
t_win = [-0.5 2];                                        % seconds
[~, crop_start] = min(abs(t_comp_vec - t_win(1)));
[~, crop_end]   = min(abs(t_comp_vec - t_win(2)));
crop_idx = crop_start : crop_end;
n_samp   = length(crop_idx);
t_vec    = t_comp_vec(crop_idx);                         % display time axis

% Baseline period for percentage-change normalisation
bl_win = [-0.5 0];                                       % seconds (pre-stimulus)
bl_idx = t_vec >= bl_win(1) & t_vec <= bl_win(2);       % logical index into t_vec

% Plotting
colors   = color_def('AOC');
fontSize = 20;

%% Task definitions
taskNames  = {'sternberg', 'nback'};
taskFiles  = {'dataET_sternberg', 'dataET_nback'};
condCodes  = {[2 4 6], [1 2 3]};
condLabels = {{'WM 2', 'WM 4', 'WM 6'}, {'1-back', '2-back', '3-back'}};

%% ====================================================================
%  Process both tasks
%  ====================================================================
for iTask = 1:2
    tName   = taskNames{iTask};
    tFile   = taskFiles{iTask};
    cCodes  = condCodes{iTask};
    cLabels = condLabels{iTask};
    nConds  = length(cCodes);

    fprintf('\n=== Processing %s ===\n', tName);

    % Preallocate: subjects x timepoints x conditions
    subjRates = nan(nSubj, n_samp, nConds);
    rasterAll = cell(nConds, 1);

    for subj = 1:nSubj
        fprintf('  Subject %d/%d\n', subj, nSubj);
        spath = fullfile(datapath, subjects{subj}, 'gaze');

        % Load ET data
        fpath = fullfile(spath, [tFile, '.mat']);
        if ~isfile(fpath)
            warning('File not found for subject %s', subjects{subj});
            continue
        end
        S = load(fpath);

        % Prefer dataETlong (extended time range incl. baseline)
        if isfield(S, 'dataETlong')
            et = S.dataETlong;
        elseif isfield(S, 'dataet')
            et = S.dataet;
        else
            warning('No ET struct for subject %s, task %s', subjects{subj}, tName);
            continue
        end

        nTrials = length(et.trial);
        condSpikes = cell(nConds, 1);

        for trl = 1:nTrials
            raw = et.trial{trl};
            t   = et.time{trl};

            % Keep x, y, pupil; invert Y
            raw = raw(1:min(3, size(raw, 1)), :);
            raw(2, :) = screenH - raw(2, :);

            % Out-of-bounds → NaN
            oob = raw(1,:) < 0 | raw(1,:) > screenW | ...
                  raw(2,:) < 0 | raw(2,:) > screenH;
            raw(:, oob) = NaN;

            % Remove blinks
            raw = remove_blinks(raw, blink_win);

            % Extract COMPUTATION window (wider than display for edge-free smoothing)
            idx_win = t >= t_comp(1) & t <= t_comp(2);
            gx = raw(1, idx_win);
            gy = raw(2, idx_win);

            if sum(isfinite(gx) & isfinite(gy)) < 50
                continue
            end

            % Detect microsaccades using detect_microsaccades
            % Global threshold: computed over the full computation window
            velData = [gx; gy];
            [~, msDetails] = detect_microsaccades(fsample, velData, length(gx));

            % Build binary spike train from onset indices
            spikeVec = zeros(1, length(gx));
            if ~isempty(msDetails.Onset)
                onsets = msDetails.Onset(msDetails.Onset >= 1 & msDetails.Onset <= length(gx));
                spikeVec(onsets) = 1;
            end

            % Pad/trim to fixed computation-window length
            if length(spikeVec) >= n_comp
                spikeVec = spikeVec(1:n_comp);
            else
                spikeVec(end+1 : n_comp) = 0;
            end

            % Determine condition
            if size(et.trialinfo, 2) > 1
                cond = et.trialinfo(trl, 1) - 20;
            else
                cond = et.trialinfo(trl) - 20;
            end

            cIdx = find(cCodes == cond);
            if isempty(cIdx); continue; end

            condSpikes{cIdx}(end+1, :) = spikeVec;
        end

        % Per-subject: average spike trains → rate → smooth → crop
        for c = 1:nConds
            if size(condSpikes{c}, 1) >= 3
                % Mean spike probability × fsample = rate (Hz)
                rate = mean(condSpikes{c}, 1) * fsample;

                % Smooth with Gaussian kernel (over full computation window)
                smoothed = conv(rate, gKernel, 'same');

                % Crop to display window (edge-effect-free)
                subjRates(subj, :, c) = smoothed(crop_idx);

                % Accumulate rasters (cropped to display window)
                rasterAll{c} = [rasterAll{c}; condSpikes{c}(:, crop_idx)];
            end
        end

        clear S et
    end

    %% Grand averages across subjects (raw Hz)
    grandMean = squeeze(nanmean(subjRates, 1));                        % n_samp x nConds
    nValid    = squeeze(sum(~isnan(subjRates(:, 1, :)), 1));           % nConds x 1
    grandSEM  = squeeze(nanstd(subjRates, 0, 1)) ./ sqrt(nValid');    % n_samp x nConds
    grandAll  = nanmean(grandMean, 2);                                 % grand avg across conditions

    %% Percentage-change baseline normalisation (per subject, per condition)
    subjRates_pct = nan(size(subjRates));
    for subj = 1:nSubj
        for c = 1:nConds
            rate_ts = subjRates(subj, :, c);
            if all(isnan(rate_ts)); continue; end
            bl_mean = nanmean(rate_ts(bl_idx));
            if bl_mean == 0 || isnan(bl_mean); continue; end
            subjRates_pct(subj, :, c) = (rate_ts - bl_mean) ./ bl_mean * 100;
        end
    end

    % Grand averages (percentage change)
    grandMean_pct = squeeze(nanmean(subjRates_pct, 1));                          % n_samp x nConds
    nValid_pct    = squeeze(sum(~isnan(subjRates_pct(:, 1, :)), 1));             % nConds x 1
    grandSEM_pct  = squeeze(nanstd(subjRates_pct, 0, 1)) ./ sqrt(nValid_pct');  % n_samp x nConds
    grandAll_pct  = nanmean(grandMean_pct, 2);                                   % grand avg across conditions

    %% ================================================================
    %  FIGURE 1: Smoothed MS Rate Time Courses per Condition (raw Hz)
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);
    hold on

    % Per-condition lines with SEM shading
    h = gobjects(nConds, 1);
    for c = 1:nConds
        mu  = grandMean(:, c);
        sem = grandSEM(:, c);

        % SEM ribbon (excluded from legend)
        fill([t_vec, fliplr(t_vec)], ...
             [(mu + sem)', fliplr((mu - sem)')], ...
             colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        h(c) = plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]');
    ylabel('Microsaccade Rate [Hz]');
    title([upper(tName(1)), tName(2:end), ' — Microsaccade Rate Time Course']);
    legend(h, cLabels, 'Location', 'best', 'FontSize', fontSize - 4);
    set(gca, 'FontSize', fontSize);
    hold off

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_rate_', tName, '.png']));

    %% ================================================================
    %  FIGURE 1b: Percentage-Change MS Rate Time Courses per Condition
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);
    hold on

    % Per-condition lines with SEM shading
    h_pct = gobjects(nConds, 1);
    for c = 1:nConds
        mu  = grandMean_pct(:, c);
        sem = grandSEM_pct(:, c);

        % SEM ribbon
        fill([t_vec, fliplr(t_vec)], ...
             [(mu + sem)', fliplr((mu - sem)')], ...
             colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        h_pct(c) = plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(t_win);
    xlabel('Time [s]');
    ylabel('Microsaccade Rate [% change]');
    title([upper(tName(1)), tName(2:end), ' — Microsaccade Rate (Percentage Change)']);
    legend(h_pct, cLabels, 'Location', 'best', 'FontSize', fontSize - 4);
    set(gca, 'FontSize', fontSize);
    hold off

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_rate_pct_', tName, '.png']));

    %% ================================================================
    %  FIGURE 2: Raster Plots per Condition
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);

    for c = 1:nConds
        subplot(1, nConds, c);
        hold on

        raster = rasterAll{c};
        if isempty(raster)
            title(cLabels{c}, 'FontSize', fontSize);
            continue
        end

        nR = size(raster, 1);

        % Subsample if too many trials for readability
        maxShow = 300;
        if nR > maxShow
            rIdx = sort(randperm(nR, maxShow));
            raster = raster(rIdx, :);
            nR = maxShow;
        end

        % Plot each microsaccade onset as a dot
        nSampPlot = min(size(raster, 2), n_samp);
        [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
        if ~isempty(trialIdx)
            plot(t_vec(sampleIdx), trialIdx, '.', ...
                'Color', colors(c, :), 'MarkerSize', 10);
        end

        xline(0, 'k--', 'LineWidth', 1.5);
        xlim(t_win);
        ylim([0 nR + 1]);
        xlabel('Time [s]');
        if c == 1; ylabel('Trial'); end
        title(cLabels{c}, 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize - 4, 'YDir', 'reverse');
        hold off
    end

    sgtitle([upper(tName(1)), tName(2:end), ' — Microsaccade Raster'], ...
        'FontSize', fontSize + 2);

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_raster_', tName, '.png']));

    %% ================================================================
    %  FIGURE 2b: Collapsed Raster (all conditions pooled)
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);
    hold on

    % Pool all conditions
    rasterPooled = [];
    for c = 1:nConds
        if ~isempty(rasterAll{c})
            rasterPooled = [rasterPooled; rasterAll{c}];
        end
    end

    if ~isempty(rasterPooled)
        nR = size(rasterPooled, 1);

        % Subsample if too many trials
        maxShow = 500;
        if nR > maxShow
            rIdx = sort(randperm(nR, maxShow));
            rasterPooled = rasterPooled(rIdx, :);
            nR = maxShow;
        end

        nSampPlot = min(size(rasterPooled, 2), n_samp);
        [trialIdx, sampleIdx] = find(rasterPooled(:, 1:nSampPlot));
        if ~isempty(trialIdx)
            plot(t_vec(sampleIdx), trialIdx, '.', ...
                'Color', [0.3 0.3 0.3], 'MarkerSize', 10);
        end
    end

    xline(0, 'k--', 'LineWidth', 1.5);
    xlim(t_win);
    ylim([0 nR + 1]);
    xlabel('Time [s]');
    ylabel('Trial');
    title([upper(tName(1)), tName(2:end), ' — Microsaccade Raster (All Conditions)']);
    set(gca, 'FontSize', fontSize, 'YDir', 'reverse');
    hold off

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_raster_collapsed_', tName, '.png']));

    %% ================================================================
    %  FIGURE 3: Combined (Rate on top + Raster on bottom)
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);

    % --- Top panel: Smoothed rate ---
    ax1 = subplot(2, 1, 1);
    hold on
    for c = 1:nConds
        mu  = grandMean(:, c);
        sem = grandSEM(:, c);
        fill([t_vec, fliplr(t_vec)], ...
             [(mu + sem)', fliplr((mu - sem)')], ...
             colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
    end
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlim(t_win);
    ylabel('MS Rate [Hz]');
    title([upper(tName(1)), tName(2:end), ' — Microsaccade Dynamics']);
    legend(cLabels, 'Location', 'best', 'FontSize', fontSize - 4);
    set(gca, 'FontSize', fontSize, 'XTickLabel', []);
    hold off

    % --- Bottom panel: Raster grouped by condition ---
    ax2 = subplot(2, 1, 2);
    hold on
    yOff      = 0;
    yTicks    = [];
    yTickLbls = {};

    for c = 1:nConds
        raster = rasterAll{c};
        if isempty(raster); continue; end

        nR = min(size(raster, 1), 100);
        raster = raster(1:nR, :);

        nSampPlot = min(size(raster, 2), n_samp);
        [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
        if ~isempty(trialIdx)
            plot(t_vec(sampleIdx), trialIdx + yOff, '.', ...
                'Color', colors(c, :), 'MarkerSize', 10);
        end

        yTicks(end+1)    = yOff + nR / 2;
        yTickLbls{end+1} = cLabels{c};
        yOff = yOff + nR + 5;  % gap between conditions
    end

    xline(0, 'k--', 'LineWidth', 1.5);
    xlim(t_win);
    ylim([0 yOff]);
    xlabel('Time [s]');
    ylabel('Trials');
    set(gca, 'FontSize', fontSize, 'YDir', 'reverse', ...
        'YTick', yTicks, 'YTickLabel', yTickLbls);
    hold off

    linkaxes([ax1, ax2], 'x');

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_combined_', tName, '.png']));

    %% ================================================================
    %  FIGURE 3b: Combined Percentage-Change (Rate on top + Raster on bottom)
    %  ================================================================
    figure('Color', 'w', 'Position', [0 0 1512 982]);

    % --- Top panel: Percentage-change rate ---
    ax1p = subplot(2, 1, 1);
    hold on
    for c = 1:nConds
        mu  = grandMean_pct(:, c);
        sem = grandSEM_pct(:, c);
        fill([t_vec, fliplr(t_vec)], ...
             [(mu + sem)', fliplr((mu - sem)')], ...
             colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        plot(t_vec, mu, '-', 'Color', colors(c, :), 'LineWidth', 2.5);
    end
    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(t_win);
    ylabel('MS Rate [% change]');
    title([upper(tName(1)), tName(2:end), ' — Microsaccade Dynamics (Percentage Change)']);
    legend(cLabels, 'Location', 'best', 'FontSize', fontSize - 4);
    set(gca, 'FontSize', fontSize, 'XTickLabel', []);
    hold off

    % --- Bottom panel: Raster grouped by condition ---
    ax2p = subplot(2, 1, 2);
    hold on
    yOff      = 0;
    yTicks    = [];
    yTickLbls = {};

    for c = 1:nConds
        raster = rasterAll{c};
        if isempty(raster); continue; end

        nR = min(size(raster, 1), 100);
        raster = raster(1:nR, :);

        nSampPlot = min(size(raster, 2), n_samp);
        [trialIdx, sampleIdx] = find(raster(:, 1:nSampPlot));
        if ~isempty(trialIdx)
            plot(t_vec(sampleIdx), trialIdx + yOff, '.', ...
                'Color', colors(c, :), 'MarkerSize', 10);
        end

        yTicks(end+1)    = yOff + nR / 2;
        yTickLbls{end+1} = cLabels{c};
        yOff = yOff + nR + 5;
    end

    xline(0, 'k--', 'LineWidth', 1.5);
    xlim(t_win);
    ylim([0 yOff]);
    xlabel('Time [s]');
    ylabel('Trials');
    set(gca, 'FontSize', fontSize, 'YDir', 'reverse', ...
        'YTick', yTicks, 'YTickLabel', yTickLbls);
    hold off

    linkaxes([ax1p, ax2p], 'x');

    saveas(gcf, fullfile(figpath, ['AOC_gaze_microsaccades_combined_pct_', tName, '.png']));
end

fprintf('\n=== All figures saved to %s ===\n', figpath);
