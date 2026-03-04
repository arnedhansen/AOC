%% AOC Gaze Time Courses — Supplements
% Time-course figures for gaze deviation, microsaccades, BCEA, gaze velocity, SPL.
% X-axis: time; Y-axis: metric unit. Three lines with shaded error bars per condition.
% Both Sternberg and N-back tasks.
%
% Output: 10 figures (5 metrics x 2 tasks)

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');
colors = color_def('AOC');
subjects = exclude_subjects(subjects, 'AOC');

basepath = path;
if ~contains(path, 'features')
    basepath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
end

fs = 500;
binDur = 0.05;
t_series = -0.5:1/fs:2;
n_bins = length(t_series(1):binDur:t_series(end)) - 1;
t_plot = linspace(-0.5 + binDur/2, 2 - binDur/2, n_bins);
centreX = 400;
centreY = 300;
k95 = -log(1 - 0.95);

% Figure output (assumes basepath = .../AOC/data/features)
proj_root = fileparts(fileparts(fileparts(basepath)));
figdir = fullfile(proj_root, 'figures', 'gaze', 'timeCourses');
if ~exist(figdir, 'dir'), mkdir(figdir); end

%% Task config
tasks = {'sternberg', 'nback'};
cond_labels_sternberg = {'WM2', 'WM4', 'WM6'};
cond_labels_nback = {'1-back', '2-back', '3-back'};

%% Aggregate per task
for itask = 1:2
    task = tasks{itask};
    if itask == 1
        cond_labels = cond_labels_sternberg;
    else
        cond_labels = cond_labels_nback;
    end

    % Preallocate: [nSubj x nConds x nBins]
    nSubj = length(subjects);
    dev_subj    = nan(nSubj, 3, n_bins);
    ms_subj     = nan(nSubj, 3, n_bins);
    bcea_subj   = nan(nSubj, 3, n_bins);
    vel_subj    = nan(nSubj, 3, n_bins);
    spl_subj    = nan(nSubj, 3, n_bins);

    datapath_base = fullfile(basepath, '%s', 'gaze', ['gaze_series_' task '_trials.mat']);
    has_MSSeries = strcmp(task, 'sternberg');

    for subj = 1:nSubj
        datapath = sprintf(datapath_base, subjects{subj});
        if ~isfile(datapath)
            warning('Missing %s for subject %s, skipping.', task, subjects{subj});
            continue
        end

        load(datapath, 'gaze_x', 'gaze_y', 'trialinfo', 'ScanPathSeries', 'ScanPathSeriesT');
        if has_MSSeries
            load(datapath, 'MSSeries');
        end

        % Condition indices (robust parsing): map to 1=low, 2=med, 3=high
        if size(trialinfo, 2) == 2
            conds_raw = trialinfo(:, 1);
        elseif size(trialinfo, 1) == 2
            conds_raw = trialinfo(1, :)';
        else
            conds_raw = trialinfo(:);
        end
        u = unique(conds_raw);
        u = sort(u(isfinite(u)));
        if numel(u) < 3
            u = [u; repmat(u(end), 3 - numel(u), 1)];
        end
        conds = nan(size(conds_raw));
        conds(conds_raw == u(1)) = 1;
        conds(conds_raw == u(2)) = 2;
        conds(conds_raw == u(3)) = 3;
        idx1 = (conds == 1);
        idx2 = (conds == 2);
        idx3 = (conds == 3);
        idxc = {idx1, idx2, idx3};

        nTrials = numel(ScanPathSeries);
        binSamp = round(binDur * fs);

        % Per-trial arrays for metrics requiring computation
        dev_trials = nan(nTrials, n_bins);
        ms_trials  = nan(nTrials, n_bins);
        bcea_trials = nan(nTrials, n_bins);
        vel_trials = nan(nTrials, n_bins);
        spl_trials = nan(nTrials, n_bins);

        for trl = 1:nTrials
            gx = gaze_x{trl};
            gy = gaze_y{trl};
            if isempty(gx) || isempty(gy)
                continue
            end
            gx = gx(:)';
            gy = gy(:)';
            L = min(numel(gx), numel(gy));
            if L < 50
                continue
            end
            gx = gx(1:L);
            gy = gy(1:L);

            % Gaze deviation (Euclidean from center) -> bin by mean
            dev_raw = sqrt((gx - centreX).^2 + (gy - centreY).^2);
            dev_binned = bin_series(dev_raw, fs, binDur, 'mean');
            nvb = min(n_bins, numel(dev_binned));
            dev_trials(trl, 1:nvb) = dev_binned(1:nvb);

            % Microsaccades (MSSeries: count per bin -> rate = count/binDur)
            if has_MSSeries
                n_ms = min(n_bins, size(MSSeries, 2));
                ms_trials(trl, 1:n_ms) = MSSeries(trl, 1:n_ms) / binDur;
            else
                vel_full = [gx; gy];
                valid_xy = isfinite(gx) & isfinite(gy);
                T_full = sum(valid_xy) / fs;
                if T_full > 0
                    [~, msf] = detect_microsaccades(fs, vel_full, L);
                    ms_onsets = zeros(1, L);
                    if ~isempty(msf) && isfield(msf, 'Onset') && ~isempty(msf.Onset)
                        onsets = msf.Onset(msf.Onset >= 1 & msf.Onset <= L);
                        ms_onsets(onsets) = 1;
                    end
                    for b = 1:n_bins
                        i1 = (b-1)*binSamp + 1;
                        i2 = min(b*binSamp, L);
                        n_ev = sum(ms_onsets(i1:i2));
                        n_val = sum(valid_xy(i1:i2));
                        eff_sec = n_val / fs;
                        if eff_sec > 0
                            ms_trials(trl, b) = n_ev / eff_sec;
                        end
                    end
                end
            end

            % Sliding-window BCEA (200 ms, step 50 ms)
            win_samp = round(0.2 * fs);
            for b = 1:n_bins
                t_center = t_plot(b);
                i_center = round((t_center - t_series(1)) * fs) + 1;
                i1 = max(1, i_center - win_samp/2);
                i2 = min(L, i_center + win_samp/2 - 1);
                xw = gx(i1:i2);
                yw = gy(i1:i2);
                ok = isfinite(xw) & isfinite(yw);
                xw = double(xw(ok));
                yw = double(yw(ok));
                if numel(xw) >= 10
                    sx = std(xw);
                    sy = std(yw);
                    rho = corr(xw(:), yw(:));
                    if isnan(rho), rho = 0; end
                    bcea_trials(trl, b) = 2 * k95 * pi * sx * sy * sqrt(max(0, 1 - rho^2));
                end
            end

            % Gaze velocity (Savitzky-Golay) -> bin by mean
            [vx, vy] = compute_velocity_sg(gx, gy, fs, 3);
            spd = hypot(vx, vy);
            % Optional: attenuate saccade spikes (|z|>4) via interpolation
            velZthr = 4;
            zspd = (spd - nanmean(spd)) / max(nanstd(spd), eps);
            bad = abs(zspd) > velZthr;
            if any(bad)
                spd(bad) = NaN;
                spd = fillmissing(spd, 'linear', 'EndValues', 'nearest');
            end
            vel_binned = bin_series(spd, fs, binDur, 'mean');
            nvb = min(n_bins, numel(vel_binned));
            vel_trials(trl, 1:nvb) = vel_binned(1:nvb);

            % SPL (step length per sample, interpolate to common grid, then bin)
            st = ScanPathSeriesT{trl};
            ss = ScanPathSeries{trl};
            if ~isempty(st) && numel(st) == numel(ss)
                try
                    spl_interp = interp1(st, ss, t_series(2:end), 'linear', NaN);
                    spl_binned = bin_series(spl_interp, fs, binDur, 'sum');
                    nvb = min(n_bins, numel(spl_binned));
                    spl_trials(trl, 1:nvb) = spl_binned(1:nvb);
                catch
                    % leave NaN
                end
            end
        end

        % Subject mean per condition
        for c = 1:3
            idx = idxc{c};
            if any(idx)
                dev_subj(subj, c, :) = nanmean(dev_trials(idx, :), 1);
                ms_subj(subj, c, :)  = nanmean(ms_trials(idx, :), 1);
                bcea_subj(subj, c, :) = nanmean(bcea_trials(idx, :), 1);
                vel_subj(subj, c, :)  = nanmean(vel_trials(idx, :), 1);
                spl_subj(subj, c, :)  = nanmean(spl_trials(idx, :), 1);
            end
        end
    end

    %% Grand average and SEM per condition
    grand_dev  = squeeze(nanmean(dev_subj, 1));
    sem_dev    = squeeze(nanstd(dev_subj, [], 1)) ./ sqrt(sum(isfinite(dev_subj), 1));
    grand_ms   = squeeze(nanmean(ms_subj, 1));
    sem_ms     = squeeze(nanstd(ms_subj, [], 1)) ./ sqrt(sum(isfinite(ms_subj), 1));
    grand_bcea = squeeze(nanmean(bcea_subj, 1));
    sem_bcea   = squeeze(nanstd(bcea_subj, [], 1)) ./ sqrt(sum(isfinite(bcea_subj), 1));
    grand_vel  = squeeze(nanmean(vel_subj, 1));
    sem_vel    = squeeze(nanstd(vel_subj, [], 1)) ./ sqrt(sum(isfinite(vel_subj), 1));
    grand_spl  = squeeze(nanmean(spl_subj, 1));
    sem_spl    = squeeze(nanstd(spl_subj, [], 1)) ./ sqrt(sum(isfinite(spl_subj), 1));

    %% Replace SEM NaN/Inf with 0 for plotting
    sem_dev(~isfinite(sem_dev)) = 0;
    sem_ms(~isfinite(sem_ms)) = 0;
    sem_bcea(~isfinite(sem_bcea)) = 0;
    sem_vel(~isfinite(sem_vel)) = 0;
    sem_spl(~isfinite(sem_spl)) = 0;

    %% Plot each metric
    metrics = {'deviation', 'microsaccades', 'BCEA', 'velocity', 'SPL'};
    grand_all = {grand_dev, grand_ms, grand_bcea, grand_vel, grand_spl};
    sem_all  = {sem_dev, sem_ms, sem_bcea, sem_vel, sem_spl};
    ylabels = {'Gaze deviation [px]', 'Microsaccade rate [MS/s]', 'BCEA [px^2]', 'Gaze velocity [px/s]', 'Scan path length [px]'};

    for im = 1:5
        grand = grand_all{im};
        sem = sem_all{im};
        figure('Position', [0 0 1512 982], 'Color', 'w');
        hold on
        for c = 1:3
            shadedErrorBar(t_plot, grand(c, :), sem(c, :), ...
                'lineProps', {'-', 'Color', colors(c,:), 'LineWidth', 2}, 'transparent', true);
        end
        xline(0, 'k--', 'LineWidth', 1);
        xlabel('Time [s]');
        ylabel(ylabels{im});
        title([upper(task(1)) task(2:end) ' — ' strrep(metrics{im}, '_', ' ')]);
        xlim([t_series(1) t_series(end)]);
        box on
        set(gca, 'FontSize', 20);
        legend(cond_labels, 'Location', 'northeast');
        saveas(gcf, fullfile(figdir, ['AOC_gaze_timeCourse_' metrics{im} '_' task '.png']));
    end
end

%% Helper functions
function out = bin_series(x, fs, binDur, how)
    if iscolumn(x), x = x.'; end
    binSamp = round(binDur * fs);
    nb = floor(numel(x) / binSamp);
    out = nan(1, nb);
    for b = 1:nb
        i1 = (b-1)*binSamp + 1;
        i2 = b*binSamp;
        seg = x(i1:i2);
        switch lower(how)
            case 'mean', out(b) = mean(seg, 'omitnan');
            case 'sum',  out(b) = nansum(seg);
            otherwise, error('how must be ''mean'' or ''sum''.');
        end
    end
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
    X = double(X(:)');
    Y = double(Y(:)');
    Ts = 1/fs;
    L = numel(X);
    framelen = min(21, L);
    if mod(framelen, 2) == 0, framelen = framelen - 1; end
    minLegal = polyOrd + 3;
    if mod(minLegal, 2) == 0, minLegal = minLegal + 1; end
    if framelen < minLegal, framelen = minLegal; end
    if framelen > L, framelen = L - mod(L, 2) + 1; end
    useFallback = framelen < 5;

    if ~useFallback
        Xs = sgolayfilt(X, polyOrd, framelen);
        Ys = sgolayfilt(Y, polyOrd, framelen);
        [~, G] = sgolay(polyOrd, framelen);
        d1 = (factorial(1) / (Ts^1)) * G(:, 2)';
        vx = conv(Xs, d1, 'same');
        vy = conv(Ys, d1, 'same');
    else
        vx = gradient(X) * fs;
        vy = gradient(Y) * fs;
    end
end
