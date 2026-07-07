%% AOC Gaze Microsaccade Suppression — N-Back
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
% Key outputs:
%   1. Percentage-change MS rate time courses per WM load condition
%   2. Raw MS rate time courses (Hz)
%   3. Raster plots of microsaccade onsets per condition
%   4. Combined figure (rate + raster)

%% Setup
startup
[~, paths, colors, ~] = setup('AOC');
addpath(paths.seb_path);

featPath = paths.features;
datapath = featPath;
figpath  = fullfile(paths.figures, 'gaze', 'microsaccades');
if ~isfolder(figpath), mkdir(figpath); end

dirs     = dir(datapath);
folders  = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subj_dirs = {folders.name};
subj_dirs = exclude_subjects(subj_dirs, 'AOC');
nSubj    = length(subj_dirs);

tName   = 'nback';
tFile   = 'dataET_nback';
cCodes  = [1 2 3];
cLabels = {'1-back', '2-back', '3-back'};
cLabels_legend = strcat(' ', cLabels);
nConds  = length(cCodes);

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
% Kernel half-width = 3*sigma = 150 ms; extend left to include baseline period.
t_comp     = [-1.5 2.5];                                % seconds
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
bl_win = [-1.5 -0.5];                                    % seconds (pre-stimulus)
bl_idx_comp = t_comp_vec >= bl_win(1) & t_comp_vec <= bl_win(2);

% Plotting
fontSize = 50;
smooth_sec = 0.05;
idx_viable = (t_vec >= 0) & (t_vec <= 2);
outlier_k_iqr = 1.5;
max_interp_gap_sec = 0.20;
min_subject_coverage = 0.85;
winsor_pct = 10;
win_sm = max(1, round(smooth_sec * fsample));

%% Process N-back
fprintf('\n[GAZE MS NBACK] Processing %s\n', tName);

% Preallocate: subjects x timepoints x conditions
subjRates = nan(nSubj, n_samp, nConds);
subjRates_pct = nan(nSubj, n_samp, nConds);
rasterAll = cell(nConds, 1);

for subj = 1:nSubj
    fprintf('  Subject %d/%d\n', subj, nSubj);
    spath = fullfile(datapath, subj_dirs{subj}, 'gaze');

    % Load ET data
    fpath = fullfile(spath, [tFile, '.mat']);
    if ~isfile(fpath)
        warning('File not found for subject %s', subj_dirs{subj});
        continue
    end
    S = load(fpath);

    % Prefer dataETlong (extended time range incl. baseline)
    if isfield(S, 'dataETlong')
        et = S.dataETlong;
    elseif isfield(S, 'dataet')
        et = S.dataet;
    else
        warning('No ET struct for subject %s, task %s', subj_dirs{subj}, tName);
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

            % Percentage-change baseline on full computation axis, then crop
            bl_mean = nanmean(smoothed(bl_idx_comp));
            if isfinite(bl_mean) && bl_mean > 0
                smoothed_pct = (smoothed - bl_mean) ./ bl_mean * 100;
            else
                smoothed_pct = nan(size(smoothed));
            end

            % Crop to display window (edge-effect-free)
            subjRates(subj, :, c) = smoothed(crop_idx);
            subjRates_pct(subj, :, c) = smoothed_pct(crop_idx);

            % Accumulate rasters (cropped to display window)
            rasterAll{c} = [rasterAll{c}; condSpikes{c}(:, crop_idx)];
        end
    end

    clear S et
end

%% Plot baselined microsaccade rate
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:nConds
    X = squeeze(subjRates_pct(:, :, c));
    [m, se] = summarize_subject_tc(X, idx_viable, outlier_k_iqr, ...
        max_interp_gap_sec, fsample, min_subject_coverage, win_sm, winsor_pct);
    eb = shadedErrorBar(t_vec, m, se, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end
yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
xline(0, '--k');
xlabel('Time [s]');
ylabel('Microsaccade Rate [%]');
xlim(t_win);
set(gca, 'FontSize', fontSize - 4);
box off
legend_handles = gobjects(1, nConds);
for c = 1:nConds
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cLabels_legend, 'Location', 'northeast', 'FontSize', fontSize * 0.666, 'Box', 'off');
saveas(gcf, fullfile(figpath, 'AOC_gaze_microsaccades_timecourse_nback.png'));

%% Plot raw microsaccade rate
close all
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:nConds
    X = squeeze(subjRates(:, :, c));
    [m, se] = summarize_subject_tc(X, idx_viable, outlier_k_iqr, ...
        max_interp_gap_sec, fsample, min_subject_coverage, win_sm, winsor_pct);
    eb = shadedErrorBar(t_vec, m, se, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end
yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Microsaccade Rate [Hz]');
xlim(t_win);
ylim([0 4]);
set(gca, 'FontSize', fontSize - 4);
box off
legend_handles = gobjects(1, nConds);
for c = 1:nConds
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cLabels_legend, 'Location', 'northeast', 'FontSize', fontSize * 0.666, 'Box', 'off');
saveas(gcf, fullfile(figpath, 'AOC_gaze_microsaccades_timecourse_nback_raw.png'));

%% Raster plots per condition
close all
figure('Position', [0 0 1512 982], 'Color', 'w');

for c = 1:nConds
    subplot(1, nConds, c);
    hold on

    raster = rasterAll{c};
    if isempty(raster)
        title(cLabels{c}, 'FontSize', fontSize - 4);
        set(gca, 'FontSize', fontSize - 4);
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

    xline(0, '--k');
    xlim(t_win);
    ylim([0 nR + 1]);
    xlabel('Time [s]');
    if c == 1
        ylabel('Trial');
    end
    title(cLabels{c}, 'FontSize', fontSize - 4);
    set(gca, 'FontSize', fontSize - 4, 'YDir', 'reverse');
    hold off
end

saveas(gcf, fullfile(figpath, 'AOC_gaze_microsaccades_raster_nback.png'));

%% Combined figure (rate + raster)
close all
figure('Position', [0 0 1512 982], 'Color', 'w');

ax1 = subplot(2, 1, 1);
hold on
for c = 1:nConds
    X = squeeze(subjRates(:, :, c));
    [m, se] = summarize_subject_tc(X, idx_viable, outlier_k_iqr, ...
        max_interp_gap_sec, fsample, min_subject_coverage, win_sm, winsor_pct);
    eb = shadedErrorBar(t_vec, m, se, 'lineProps', {'-'}, 'transparent', true);
    set(eb.mainLine, 'Color', colors(c, :), 'LineWidth', 2.5);
    set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.20);
    set(eb.edge(1), 'Color', 'none');
    set(eb.edge(2), 'Color', 'none');
end
xline(0, '--k');
xlim(t_win);
ylabel('Microsaccade Rate [Hz]');
legend_handles = gobjects(1, nConds);
for c = 1:nConds
    legend_handles(c) = patch(nan, nan, colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.60);
end
legend(legend_handles, cLabels_legend, 'Location', 'northeast', 'FontSize', fontSize * 0.666, 'Box', 'off');
set(gca, 'FontSize', fontSize - 4, 'XTickLabel', []);
box off
hold off

ax2 = subplot(2, 1, 2);
hold on
yOff      = 0;
yTicks    = [];
yTickLbls = {};

for c = 1:nConds
    raster = rasterAll{c};
    if isempty(raster)
        continue
    end

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

xline(0, '--k');
xlim(t_win);
ylim([0 yOff]);
xlabel('Time [s]');
ylabel('Trials');
set(gca, 'FontSize', fontSize - 4, 'YDir', 'reverse', ...
    'YTick', yTicks, 'YTickLabel', yTickLbls);
hold off

linkaxes([ax1, ax2], 'x');

saveas(gcf, fullfile(figpath, 'AOC_gaze_microsaccades_combined_nback.png'));

fprintf('\n[GAZE MS NBACK] Figures saved to %s\n', figpath);

%%
function [m, se] = summarize_subject_tc(X, idx_viable, outlier_k_iqr, ...
    max_interp_gap_sec, fs, min_subject_coverage, win_sm, winsor_pct)

X(~isfinite(X)) = NaN;

med_metric = median(X(:, idx_viable), 2, 'omitnan');
[X, ~] = exclude_outlier_trajectories(X, med_metric, outlier_k_iqr);

max_interp_gap_smp = max(1, round(max_interp_gap_sec * fs));
for s = 1:size(X, 1)
    X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
end

subj_cov = mean(isfinite(X(:, idx_viable)), 2);
keep_cov = subj_cov >= min_subject_coverage;
X = X(keep_cov, :);

if win_sm > 1
    X = movmean(X, win_sm, 2, 'omitnan');
end

[m, se] = winsorized_nanmean_se(X, winsor_pct);
end

function [X_keep, keep_subj] = exclude_outlier_trajectories(X, metric, k_iqr)
med_m = median(metric, 'omitnan');
iqr_m = iqr(metric);
if ~isfinite(iqr_m) || iqr_m == 0
    low_m = -inf;
    high_m = inf;
else
    low_m = med_m - k_iqr * iqr_m;
    high_m = med_m + k_iqr * iqr_m;
end
keep_subj = isfinite(metric) & (metric >= low_m) & (metric <= high_m);
X_keep = X(keep_subj, :);
end

function x = fill_short_nan_gaps(x, max_gap_smp)
if isempty(x)
    return
end
valid = isfinite(x);
if all(~valid) || all(valid)
    return
end

n = numel(x);
i = 1;
while i <= n
    if valid(i)
        i = i + 1;
        continue
    end
    j = i;
    while j <= n && ~valid(j)
        j = j + 1;
    end
    gap_start = i;
    gap_end = j - 1;
    gap_len = gap_end - gap_start + 1;

    left = gap_start - 1;
    right = gap_end + 1;
    if left >= 1 && right <= n && valid(left) && valid(right) && gap_len <= max_gap_smp
        x(gap_start:gap_end) = interp1([left right], [x(left) x(right)], gap_start:gap_end);
    end
    i = j;
end
end

function [m, se] = winsorized_nanmean_se(X, pct)
nT = size(X, 2);
m = nan(1, nT);
se = nan(1, nT);
for t = 1:nT
    xt = X(:, t);
    xt = xt(isfinite(xt));
    if isempty(xt)
        continue
    end
    xt = winsorize_values(xt, pct);
    m(t) = mean(xt);
    n = numel(xt);
    if n > 1
        se(t) = std(xt, 0) / sqrt(n);
    end
end
end

function xw = winsorize_values(x, pct)
xw = x;
if numel(xw) < 5 || pct <= 0
    return
end
lo = prctile(xw, pct);
hi = prctile(xw, 100 - pct);
xw(xw < lo) = lo;
xw(xw > hi) = hi;
end
