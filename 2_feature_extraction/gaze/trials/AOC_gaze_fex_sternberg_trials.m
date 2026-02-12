%% AOC Gaze Feature Extraction — Sternberg (Trial-Level)
% Trial-wise gaze deviation, SPL, pupil, MS rate in early/late/full windows with trial-specific baselines; MSSeries, ScanPathSeries. Loads dataET_sternberg. Saves gaze_matrix_sternberg_trials.mat.
%
% Extracted features:
%   GazeDeviation, ScanPathLength, PupilSize, MSRate (Early/Late/Full, baselined)
%   MSSeries, ScanPathSeries, ScanPathSeriesT, ScanPathSeriesBins. Windows: BL [-0.5 -0.25], early [0 1], late [1 2], full [0 2] s

%% Setup
startup
setup('AOC');

path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

gaze_data_sternberg_trials = struct('Trial', {}, 'ID', {}, 'Condition', {}, ...
    'GazeDeviationEarly', {}, 'GazeDeviationEarlyBL', {}, ...
    'GazeDeviationLate', {}, 'GazeDeviationLateBL', {}, ...
    'GazeDeviationFull', {}, 'GazeDeviationFullBL', {}, ...
    'ScanPathLengthEarly', {}, 'ScanPathLengthEarlyBL', {}, ...
    'ScanPathLengthLate', {}, 'ScanPathLengthLateBL', {}, ...
    'ScanPathLengthFull', {}, 'ScanPathLengthFullBL', {}, ...
    'PupilSizeEarly', {}, 'PupilSizeEarlyBL', {}, ...
    'PupilSizeLate', {}, 'PupilSizeLateBL', {}, ...
    'PupilSizeFull', {}, 'PupilSizeFullBL', {}, ...
    'MSRateEarly', {}, 'MSRateEarlyBL', {}, ...
    'MSRateLate', {}, 'MSRateLateBL', {}, ...
    'MSRateFull', {}, 'MSRateFullBL', {}, ...
    'ConvexHullAreaEarly', {}, 'ConvexHullAreaEarlyBL', {}, ...
    'ConvexHullAreaLate', {}, 'ConvexHullAreaLateBL', {}, ...
    'ConvexHullAreaFull', {}, 'ConvexHullAreaFullBL', {});

%% Parameters
fsample = 500; % Hz
screenW = 800; screenH = 600;
centreX = 400; centreY = 300;
blink_win = 50; % samples (+/-) for removal in your helper
min_valid_samples = 100; % per window after cleaning
bounds_x = [0 screenW];
bounds_y = [0 screenH];

t_base = [-0.5 -0.25];
t_early = [0 1];
t_late = [1 2];
t_full = [0 2];
t_series = [-0.5 2]; % for scan-path time series

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET_sternberg'])
    nTrials = size(dataETlong.trialinfo,1);

    %% Preallocate per-subject arrays
    ID              = nan(nTrials,1);
    Trial           = nan(nTrials,1);
    Condition       = nan(nTrials,1);

    GazeDeviationEarly    = nan(nTrials,1);
    GazeDeviationEarlyBL  = nan(nTrials,1);
    GazeDeviationLate     = nan(nTrials,1);
    GazeDeviationLateBL   = nan(nTrials,1);
    GazeDeviationFull     = nan(nTrials,1);
    GazeDeviationFullBL   = nan(nTrials,1);

    ScanPathLengthEarly   = nan(nTrials,1);
    ScanPathLengthEarlyBL = nan(nTrials,1);
    ScanPathLengthLate    = nan(nTrials,1);
    ScanPathLengthLateBL  = nan(nTrials,1);
    ScanPathLengthFull    = nan(nTrials,1);
    ScanPathLengthFullBL  = nan(nTrials,1);

    PupilSizeEarly        = nan(nTrials,1);
    PupilSizeEarlyBL      = nan(nTrials,1);   % percent
    PupilSizeLate         = nan(nTrials,1);
    PupilSizeLateBL       = nan(nTrials,1);   % percent
    PupilSizeFull         = nan(nTrials,1);
    PupilSizeFullBL       = nan(nTrials,1);   % percent

    MSRateEarly           = nan(nTrials,1);   % events/s
    MSRateEarlyBL         = nan(nTrials,1);   % dB
    MSRateLate            = nan(nTrials,1);   % events/s
    MSRateLateBL          = nan(nTrials,1);   % dB
    MSRateFull            = nan(nTrials,1);   % events/s
    MSRateFullBL          = nan(nTrials,1);   % dB
    MSSeries              = nan(nTrials,50);   % events/s/bin

    ConvexHullAreaEarly   = nan(nTrials,1);   % px²
    ConvexHullAreaEarlyBL = nan(nTrials,1);   % dB
    ConvexHullAreaLate    = nan(nTrials,1);   % px²
    ConvexHullAreaLateBL  = nan(nTrials,1);   % dB
    ConvexHullAreaFull    = nan(nTrials,1);   % px²
    ConvexHullAreaFullBL  = nan(nTrials,1);   % dB

    ScanPathSeriesT       = cell(nTrials,1);
    ScanPathSeries        = cell(nTrials,1);

    gaze_x = cell(1,nTrials);
    gaze_y = cell(1,nTrials);

    trial_num = dataETlong.trialinfo(:,2);        % global trial ID
    cond_code = dataETlong.trialinfo(:,1) - 20;   % 22/24/26

    %% Get trial-by-trial gaze data
    for trl = 1:nTrials
        close all
        raw_dat = dataETlong.trial{trl};           % rows: [x;y;pupil]
        t   = dataETlong.time{trl};

        % Invert Y and keep only first 3 rows
        raw_dat = raw_dat(1:3,:);
        raw_dat(2,:) = screenH - raw_dat(2,:);

        % Filter out-of-bounds
        inb = raw_dat(1,:) >= bounds_x(1) & raw_dat(1,:) <= bounds_x(2) & ...
            raw_dat(2,:) >= bounds_y(1) & raw_dat(2,:) <= bounds_y(2);
        raw_dat(:,~inb) = NaN;

        % Remove blinks
        raw_dat = remove_blinks(raw_dat, blink_win);

        % indices
        idx_base  = t >= t_base(1)  & t <= t_base(2);
        idx_early = t >= t_early(1) & t <= t_early(2);
        idx_late  = t >= t_late(1)  & t <= t_late(2);
        idx_full  = t >= t_full(1)  & t <= t_full(2);
        idx_series= t >= t_series(1)& t <= t_series(2);

        % Extract windows
        dat_base  = raw_dat(:, idx_base);
        dat_early = raw_dat(:, idx_early);
        dat_late  = raw_dat(:, idx_late);
        dat_full  = raw_dat(:, idx_full);
        dat_series= raw_dat(:, idx_series);

        % Minimum valid-sample threshold per window (post-cleaning)
        ok_base  = sum(all(isfinite(dat_base(1:2,:)),1))  >= min_valid_samples;
        ok_early = sum(all(isfinite(dat_early(1:2,:)),1)) >= min_valid_samples;
        ok_late  = sum(all(isfinite(dat_late(1:2,:)),1))  >= min_valid_samples;
        ok_full  = sum(all(isfinite(dat_full(1:2,:)),1))  >= min_valid_samples;

        %% Baseline metrics (per trial)
        % Gaze deviation baseline (mean Euclidean distance from centre)
        if ok_base
            dx_b = dat_base(1,:) - centreX;
            dy_b = dat_base(2,:) - centreY;
            gaze_dev_base = nanmean( sqrt(dx_b.^2 + dy_b.^2) );
        else
            gaze_dev_base = NaN;
        end

        % Scan path length baseline (sum of step lengths within window)
        if ok_base
            xb = dat_base(1,:); yb = dat_base(2,:);
            dxb = diff(xb); dyb = diff(yb);
            spl_base = nansum( sqrt(dxb.^2 + dyb.^2) );
        else
            spl_base = NaN;
        end

        % Pupil size baseline (mean, raw_dat units)
        if ok_base
            pupil_base = nanmean(dat_base(3,:));
        else
            pupil_base = NaN;
        end

        % MS rate baseline (events/s) within baseline window
        if ok_base
            vel_base = [dat_base(1,:); dat_base(2,:)];
            T_base = sum(isfinite(vel_base(1,:)) & isfinite(vel_base(2,:))) / fsample;
            if T_base > 0
                [~, msb] = detect_microsaccades(fsample, vel_base, size(vel_base,2));
                ms_count_b = numel(msb.Onset); % assume your struct has onsets
                ms_rate_base = ms_count_b / T_base;
                if ~isfinite(ms_rate_base) || ms_rate_base <= 0
                    ms_rate_base = NaN;
                end
            else
                ms_rate_base = NaN;
            end
        else
            ms_rate_base = NaN;
        end

        % Convex hull area baseline (px²)
        if ok_base
            valid_b = isfinite(dat_base(1,:)) & isfinite(dat_base(2,:));
            xb_v = dat_base(1,valid_b); yb_v = dat_base(2,valid_b);
            upts_b = unique([xb_v(:), yb_v(:)], 'rows');
            if size(upts_b, 1) >= 3
                try
                    [~, cha_base] = convhull(xb_v(:), yb_v(:));
                catch
                    cha_base = NaN;
                end
            else
                cha_base = NaN;
            end
        else
            cha_base = NaN;
        end

        %% Early window
        if ok_early
            xe = dat_early(1,:); ye = dat_early(2,:);
            dxe = xe - centreX; dye = ye - centreY;
            gd_early= nanmean( sqrt(dxe.^2 + dye.^2) );

            dxe_s = diff(xe); dye_s = diff(ye);
            spl_early= nansum( sqrt(dxe_s.^2 + dye_s.^2) );

            pup_early= nanmean(dat_early(3,:));

            vel_early= [xe; ye];
            T_early= sum(isfinite(xe) & isfinite(ye)) / fsample;
            if T_early> 0
                [~, msei] = detect_microsaccades(fsample, vel_early, size(vel_early,2));
                ms_early= numel(msei.Onset) / T_early;  % events/s
                if ~isfinite(ms_early); ms_early= NaN; end
            else
                ms_early= NaN;
            end

            % Convex hull area (px²)
            valid_e = isfinite(xe) & isfinite(ye);
            xe_v = xe(valid_e); ye_v = ye(valid_e);
            upts_e = unique([xe_v(:), ye_v(:)], 'rows');
            if size(upts_e, 1) >= 3
                try
                    [~, cha_early] = convhull(xe_v(:), ye_v(:));
                catch
                    cha_early = NaN;
                end
            else
                cha_early = NaN;
            end
        else
            gd_early= NaN; spl_early= NaN; pup_early= NaN; ms_early= NaN; cha_early= NaN;
        end

        % Baseline-correct early (dB for GD/SPL/MS/CHA; % for pupil)
        if isfinite(gd_early) && isfinite(gaze_dev_base) && gaze_dev_base > 0
            gd_early_bl = 10*log10(gd_early/ gaze_dev_base);
        else
            gd_early_bl = NaN;
        end

        if isfinite(spl_early) && isfinite(spl_base) && spl_base > 0
            spl_early_bl = 10*log10(spl_early/ spl_base);
        else
            spl_early_bl = NaN;
        end

        if isfinite(ms_early) && isfinite(ms_rate_base) && ms_rate_base > 0
            ms_early_bl = 10*log10(ms_early/ ms_rate_base);
        else
            ms_early_bl = NaN;
        end

        if isfinite(pup_early) && isfinite(pupil_base) && pupil_base ~= 0
            pup_early_bl = 100 * (pup_early- pupil_base) / pupil_base;
        else
            pup_early_bl = NaN;
        end

        if isfinite(cha_early) && isfinite(cha_base) && cha_base > 0
            cha_early_bl = 10*log10(cha_early / cha_base);
        else
            cha_early_bl = NaN;
        end

        %% Late window
        if ok_late
            xl = dat_late(1,:); yl = dat_late(2,:);
            dxl = xl - centreX; dyl = yl - centreY;
            gd_late = nanmean( sqrt(dxl.^2 + dyl.^2) );

            dxl_s = diff(xl); dyl_s = diff(yl);
            spl_late = nansum( sqrt(dxl_s.^2 + dyl_s.^2) );

            pup_late = nanmean(dat_late(3,:));

            vel_late = [xl; yl];
            T_late = sum(isfinite(xl) & isfinite(yl)) / fsample;
            if T_late> 0
                [~, msl] = detect_microsaccades(fsample, vel_late, size(vel_late,2));
                ms_late = numel(msl.Onset) / T_late;  % events/s
                if ~isfinite(ms_late); ms_late = NaN; end
            else
                ms_late = NaN;
            end

            % Convex hull area (px²)
            valid_l = isfinite(xl) & isfinite(yl);
            xl_v = xl(valid_l); yl_v = yl(valid_l);
            upts_l = unique([xl_v(:), yl_v(:)], 'rows');
            if size(upts_l, 1) >= 3
                try
                    [~, cha_late] = convhull(xl_v(:), yl_v(:));
                catch
                    cha_late = NaN;
                end
            else
                cha_late = NaN;
            end
        else
            gd_late = NaN; spl_late = NaN; pup_late = NaN; ms_late = NaN; cha_late = NaN;
        end

        % Baseline-correct late (dB for GD/SPL/MS/CHA; % for pupil)
        if isfinite(gd_late) && isfinite(gaze_dev_base) && gaze_dev_base > 0
            gd_late_bl = 10*log10(gd_late/ gaze_dev_base);
        else
            gd_late_bl = NaN;
        end

        if isfinite(spl_late) && isfinite(spl_base) && spl_base > 0
            spl_late_bl = 10*log10(spl_late/ spl_base);
        else
            spl_late_bl = NaN;
        end

        if isfinite(ms_late) && isfinite(ms_rate_base) && ms_rate_base > 0
            ms_late_bl = 10*log10(ms_late/ ms_rate_base);
        else
            ms_late_bl = NaN;
        end

        if isfinite(pup_late) && isfinite(pupil_base) && pupil_base ~= 0
            pup_late_bl = 100 * (pup_late- pupil_base) / pupil_base;
        else
            pup_late_bl = NaN;
        end

        if isfinite(cha_late) && isfinite(cha_base) && cha_base > 0
            cha_late_bl = 10*log10(cha_late / cha_base);
        else
            cha_late_bl = NaN;
        end

        %% Full window
        if ok_full
            xf = dat_full(1,:); yf = dat_full(2,:);
            dxf = xf - centreX; dyf = yf - centreY;
            gd_full = nanmean( sqrt(dxf.^2 + dyf.^2) );

            dxf_s = diff(xf); dyf_s = diff(yf);
            spl_full = nansum( sqrt(dxf_s.^2 + dyf_s.^2) );

            pup_full = nanmean(dat_full(3,:));

            vel_full = [xf; yf];
            T_full = sum(isfinite(xf) & isfinite(yf)) / fsample;
            if T_full> 0
                [~, msf] = detect_microsaccades(fsample, vel_full, size(vel_full,2));
                ms_full = numel(msf.Onset) / T_full;  % events/s
                if ~isfinite(ms_full); ms_full = NaN; end
            else
                ms_full = NaN;
            end

            % Convex hull area (px²)
            valid_f = isfinite(xf) & isfinite(yf);
            xf_v = xf(valid_f); yf_v = yf(valid_f);
            upts_f = unique([xf_v(:), yf_v(:)], 'rows');
            if size(upts_f, 1) >= 3
                try
                    [~, cha_full] = convhull(xf_v(:), yf_v(:));
                catch
                    cha_full = NaN;
                end
            else
                cha_full = NaN;
            end
        else
            gd_full = NaN; spl_full = NaN; pup_full = NaN; ms_full = NaN; cha_full = NaN;
        end

        % Baseline-correct full (dB for GD/SPL/MS/CHA; % for pupil)
        if isfinite(gd_full) && isfinite(gaze_dev_base) && gaze_dev_base > 0
            gd_full_bl = 10*log10(gd_full/ gaze_dev_base);
        else
            gd_full_bl = NaN;
        end

        if isfinite(spl_full) && isfinite(spl_base) && spl_base > 0
            spl_full_bl = 10*log10(spl_full/ spl_base);
        else
            spl_full_bl = NaN;
        end

        if isfinite(ms_full) && isfinite(ms_rate_base) && ms_rate_base > 0
            ms_full_bl = 10*log10(ms_full/ ms_rate_base);
        else
            ms_full_bl = NaN;
        end

        if isfinite(pup_full) && isfinite(pupil_base) && pupil_base ~= 0
            pup_full_bl = 100 * (pup_full- pupil_base) / pupil_base;
        else
            pup_full_bl = NaN;
        end

        if isfinite(cha_full) && isfinite(cha_base) && cha_base > 0
            cha_full_bl = 10*log10(cha_full / cha_base);
        else
            cha_full_bl = NaN;
        end

        %% Scan-path time series [-0.5, 2]
        xs = dat_series(1,:); ys = dat_series(2,:);
        ts = t(idx_series);

        % Keep indices that are finite for both x and y
        valid_series = isfinite(xs) & isfinite(ys);
        xs(~valid_series) = NaN; ys(~valid_series) = NaN;

        % Step-wise path length (px) aligned to ts(2:end)
        dxs = diff(xs); dys = diff(ys);
        step_series = sqrt(dxs.^2 + dys.^2);    % one value per sample step
        time_steps  = ts(2:end);

        % Bin parameters
        win_samp = round(0.05 * fsample);       % 50 ms → 25 samples
        n_bins   = floor(numel(step_series) / win_samp);

        % Preallocate binned series and times
        step_series_bin = nan(1, n_bins);
        time_bins       = nan(1, n_bins);

        for b = 1:n_bins
            idx_start = (b-1)*win_samp + 1;
            idx_end   = b*win_samp;
            segment   = step_series(idx_start:idx_end);
            step_series_bin(b) = nansum(segment);               % summed step length per 50 ms
            time_bins(b) = mean(time_steps(idx_start:idx_end));  % centre time of that bin
        end

        % Store binned data
        ScanPathSeriesT{trl}    = time_steps;
        ScanPathSeriesBins{trl} = step_series_bin;   % averaged step lengths per 50 ms bin
        ScanPathSeries{trl}     = step_series;         % their corresponding time centres

        %% Microsaccade Time Series [-0.5, 2]
        % Parameters for MS time course
        bin_size    = 0.05; % 50 ms
        bin_samples = round(bin_size * fsample);
        n_bins      = length(t_series(1):bin_size:t_series(2))-1;
        MSTime      = linspace(-0.5 + bin_size/2, 2 - bin_size/2, n_bins);

        % Per-trial microsaccade series
        % (inside your trial loop; 'trl' is the trial index within participant)
        xpos_series = dat_series(1,:);
        ypos_series = dat_series(2,:);
        vel_full    = [xpos_series; ypos_series];

        valid_xy = isfinite(xpos_series) & isfinite(ypos_series);
        T_full   = sum(valid_xy) / fsample;

        if T_full > 0
            [~, msf] = detect_microsaccades(fsample, vel_full, size(vel_full,2));

            % Binary onset vector
            ms_onsets = zeros(1, length(xpos_series));
            if ~isempty(msf) && isfield(msf, 'Onset') && ~isempty(msf.Onset)
                % Clip to bounds just in case
                onsets = msf.Onset(msf.Onset >= 1 & msf.Onset <= length(ms_onsets));
                ms_onsets(onsets) = 1;
            end

            % Compute events/s per bin, robust to missing data
            MSSeries_trl = nan(1, n_bins);
            for b = 1:n_bins
                idx_start = (b-1)*bin_samples + 1;
                idx_end   = b*bin_samples;

                % Guard against any rounding edge cases
                if idx_end > length(ms_onsets)
                    idx_end = length(ms_onsets);
                end

                n_events = sum(ms_onsets(idx_start:idx_end));
                n_valid  = sum(valid_xy(idx_start:idx_end));
                eff_sec  = n_valid / fsample;        % effective bin duration

                if eff_sec > 0
                    MSSeries_trl(b) = n_events;   % events per second
                else
                    MSSeries_trl(b) = NaN;
                end
            end

            MSSeries(trl, :) = MSSeries_trl;
            % clc
            % disp(['MS Series computed for trial ', num2str(trl)])
            % disp(nansum(MSSeries_trl))
            % pause(0.5)
        else
            % clc
            % disp('T_full < 0!')
            % pause(0.5)
            MSSeries(trl, :) = nan(1, n_bins);
        end

        %% Save trial-wise values
        ID(trl)        = str2double(subjects{subj});
        Trial(trl)     = trial_num(trl);
        Condition(trl) = cond_code(trl);

        GazeDeviationEarly(trl)    = gd_early;
        GazeDeviationEarlyBL(trl)  = gd_early_bl;
        GazeDeviationLate(trl)     = gd_late;
        GazeDeviationLateBL(trl)   = gd_late_bl;
        GazeDeviationFull(trl)     = gd_full;
        GazeDeviationFullBL(trl)   = gd_full_bl;

        ScanPathLengthEarly(trl)   = spl_early;
        ScanPathLengthEarlyBL(trl) = spl_early_bl; 
        ScanPathLengthLate(trl)    = spl_late;
        ScanPathLengthLateBL(trl)  = spl_late_bl;
        ScanPathLengthFull(trl)    = spl_full;
        ScanPathLengthFullBL(trl)  = spl_full_bl;

        PupilSizeEarly(trl)        = pup_early;     % raw_dat units
        PupilSizeEarlyBL(trl)      = pup_early_bl;  % percent
        PupilSizeLate(trl)         = pup_late;
        PupilSizeLateBL(trl)       = pup_late_bl;
        PupilSizeFull(trl)         = pup_full;
        PupilSizeFullBL(trl)       = pup_full_bl;

        MSRateEarly(trl)           = ms_early;      % events/s
        MSRateEarlyBL(trl)         = ms_early_bl;   % dB
        MSRateLate(trl)            = ms_late;       % events/s
        MSRateLateBL(trl)          = ms_late_bl;    % dB
        MSRateFull(trl)            = ms_full;       % events/s
        MSRateFullBL(trl)          = ms_full_bl;    % dB

        ConvexHullAreaEarly(trl)   = cha_early;     % px²
        ConvexHullAreaEarlyBL(trl) = cha_early_bl;  % dB
        ConvexHullAreaLate(trl)    = cha_late;      % px²
        ConvexHullAreaLateBL(trl)  = cha_late_bl;   % dB
        ConvexHullAreaFull(trl)    = cha_full;      % px²
        ConvexHullAreaFullBL(trl)  = cha_full_bl;   % dB

        % Also store cleaned gaze x/y for the series window
        gaze_x{trl} = xs;
        gaze_y{trl} = ys;
    end

    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trials = struct( ...
        'ID', num2cell(ID), ...
        'Trial', num2cell(Trial), ...
        'Condition', num2cell(Condition), ...
        'GazeDeviationEarly', num2cell(GazeDeviationEarly), ...
        'GazeDeviationEarlyBL', num2cell(GazeDeviationEarlyBL), ...
        'GazeDeviationLate', num2cell(GazeDeviationLate), ...
        'GazeDeviationLateBL', num2cell(GazeDeviationLateBL), ...
        'GazeDeviationFull', num2cell(GazeDeviationFull), ...
        'GazeDeviationFullBL', num2cell(GazeDeviationFullBL), ...
        'ScanPathLengthEarly', num2cell(ScanPathLengthEarly), ...
        'ScanPathLengthEarlyBL', num2cell(ScanPathLengthEarlyBL), ...
        'ScanPathLengthLate', num2cell(ScanPathLengthLate), ...
        'ScanPathLengthLateBL', num2cell(ScanPathLengthLateBL), ...
        'ScanPathLengthFull', num2cell(ScanPathLengthFull), ...
        'ScanPathLengthFullBL', num2cell(ScanPathLengthFullBL), ...
        'PupilSizeEarly', num2cell(PupilSizeEarly), ...
        'PupilSizeEarlyBL', num2cell(PupilSizeEarlyBL), ...
        'PupilSizeLate', num2cell(PupilSizeLate), ...
        'PupilSizeLateBL', num2cell(PupilSizeLateBL), ...
        'PupilSizeFull', num2cell(PupilSizeFull), ...
        'PupilSizeFullBL', num2cell(PupilSizeFullBL), ...
        'MSRateEarly', num2cell(MSRateEarly), ...
        'MSRateEarlyBL', num2cell(MSRateEarlyBL), ...
        'MSRateLate', num2cell(MSRateLate), ...
        'MSRateLateBL', num2cell(MSRateLateBL), ...
        'MSRateFull', num2cell(MSRateFull), ...
        'MSRateFullBL', num2cell(MSRateFullBL), ...
        'ConvexHullAreaEarly', num2cell(ConvexHullAreaEarly), ...
        'ConvexHullAreaEarlyBL', num2cell(ConvexHullAreaEarlyBL), ...
        'ConvexHullAreaLate', num2cell(ConvexHullAreaLate), ...
        'ConvexHullAreaLateBL', num2cell(ConvexHullAreaLateBL), ...
        'ConvexHullAreaFull', num2cell(ConvexHullAreaFull), ...
        'ConvexHullAreaFullBL', num2cell(ConvexHullAreaFullBL));

    %% Save data
    savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_sternberg_trials subj_data_gaze_trials

    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_sternberg_trials = [gaze_data_sternberg_trials; subj_data_gaze_trials];

    % Also save per-trial gaze series and trialinfo (for convenience)
    trialinfo = dataETlong.trialinfo';
    save([savepath 'gaze_series_sternberg_trials.mat'], 'gaze_x', 'gaze_y', 'trialinfo', 'ScanPathSeriesT', 'ScanPathSeries', 'ScanPathSeriesBins', 'MSSeries');
end

% Grand save across subjects
save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/gaze_matrix_sternberg_trials gaze_data_sternberg_trials
