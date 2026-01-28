%% AOC Omnibus — Alpha Power × Eye Velocity: Coherence Analysis
% Computes coherence between alpha power and eye velocity over time for
% Sternberg task. All conditions are combined (no condition differences).
% Uses entire stimulus presentation period (0-2 seconds).
%
% Key outputs:
%   Time-resolved coherence plot; Group statistics

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

fontSize = 20;
output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/omnibus_alpha_split/coherence';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Time-frequency parameters
fs = 500;  % Hz (EEG & ET)
t_win = [-.5 2];  % Analysis window: entire stimulus presentation period
t_vec = t_win(1):1/fs:t_win(2);
t_vec = t_vec(2:end);  % Remove first sample to match typical alignment
nTime = length(t_vec);

% Coherence parameters
coh_win_sec = 0.25;  % 250 ms window for coherence estimation
coh_win_samp = round(coh_win_sec * fs);  % 250ms * 500Hz = 125 samples
coh_overlap = 0.5;  % 50% overlap
coh_step = round(coh_win_samp * (1 - coh_overlap));  % Step = 50% of window = 62.5 samples ≈ 62 samples = 0.125s

% Alpha band
alpha_band = [8 14];  % Hz

% Velocity computation parameters
polyOrd = 3;  % Savitzky-Golay polynomial order for velocity
velZthr = 4;  % Z-score threshold for velocity outliers

%% Process Sternberg task only
task_name = 'sternberg';
config = struct('name', 'sternberg', 'gaze_file', 'gaze_series_sternberg_trials.mat', ...
    'tfr_file', 'tfr_stern_trials.mat');

% Occipital ROI labels
roi_labels = {};
try
    elecPath = fullfile(path, subjects{1}, 'eeg');
    cd(elecPath);
    load('power_stern_late_trials.mat', 'powload2_late');
    allLabs = powload2_late.label;
    for i = 1:numel(allLabs)
        L = allLabs{i};
        if contains(L, {'O'}) || contains(L, {'I'})
            roi_labels{end+1} = L;
        end
    end
catch
    warning('Could not load electrode labels, using default occipital channels');
    roi_labels = {'O1', 'Oz', 'O2', 'PO3', 'POz', 'PO4'};
end

fprintf('\n=== Processing %s task ===\n', task_name);

% Storage for coherence across subjects
coh_all_subj = {};  % Cell array: one per subject
coh_time_centers = [];  % Time centers of coherence windows

%% Loop subjects
n_valid = 0;
for s = 1:length(subjects)
    subj_id = subjects{s};
    fprintf('  Processing subject %s (%d/%d)...\n', subj_id, s, length(subjects));

    % Load EEG TFR data
    datapath_eeg = fullfile(path, subj_id, 'eeg');
    if ~exist(datapath_eeg, 'dir')
        continue;
    end

    try
        cd(datapath_eeg);
        load(config.tfr_file, 'tfr_all');
    catch
        warning('  Missing TFR file for subject %s', subj_id);
        continue;
    end

    % Check data structure
    if ~isfield(tfr_all, 'powspctrm')
        warning('  TFR data does not have powspctrm field for subject %s', subj_id);
        continue;
    end

    % Find occipital channel indices
    chan_idx = [];
    for i = 1:length(roi_labels)
        idx = find(strcmp(tfr_all.label, roi_labels{i}));
        if ~isempty(idx)
            chan_idx = [chan_idx, idx];
        end
    end

    if isempty(chan_idx)
        warning('  No occipital channels found for subject %s', subj_id);
        continue;
    end

    % Find frequency indices for alpha band (8-14 Hz)
    freq_idx = find(tfr_all.freq >= alpha_band(1) & tfr_all.freq <= alpha_band(2));
    if isempty(freq_idx)
        warning('  Alpha frequency band not found for subject %s', subj_id);
        continue;
    end

    % Find time indices for the desired window (-0.5 to 2 seconds)
    time_idx = find(tfr_all.time >= t_win(1) & tfr_all.time <= t_win(2));
    if isempty(time_idx)
        warning('  Time window %.1f-%.1fs not found for subject %s (available: %.1f to %.1fs)', ...
            t_win(1), t_win(2), subj_id, min(tfr_all.time), max(tfr_all.time));
        continue;
    end

    % Extract alpha power: average over channels and frequencies
    % powspctrm dimensions: [rpt_chan_freq_time] = [n_trials, n_channels, n_freqs, n_times]
    alpha_power = tfr_all.powspctrm(:, chan_idx, freq_idx, time_idx);
    % Size after selection: [n_trials, n_channels, n_freqs, n_times]
    
    % Average over channels (dimension 2) and frequencies (dimension 3) simultaneously
    % Reshape to combine channel and frequency dimensions, then average
    [n_trials, n_chans, n_freqs, n_times] = size(alpha_power);
    alpha_power = reshape(alpha_power, [n_trials, n_chans * n_freqs, n_times]);
    alpha_power = nanmean(alpha_power, 2);  % Average over combined chan*freq: [n_trials, 1, n_times]
    alpha_power = squeeze(alpha_power);     % Remove singleton: [n_trials, n_times]
    
    % Ensure correct orientation: trials in rows, time in columns
    if size(alpha_power, 1) ~= n_trials || size(alpha_power, 2) ~= length(time_idx)
        % Check if transposed
        if size(alpha_power, 2) == n_trials && size(alpha_power, 1) == length(time_idx)
            alpha_power = alpha_power';
        else
            warning('  Dimension mismatch for subject %s: alpha_power is %s, expected [%d %d]', ...
                subj_id, mat2str(size(alpha_power)), n_trials, length(time_idx));
            continue;
        end
    end

    % Get time vector for this window (actual TFR time points)
    alpha_time = tfr_all.time(time_idx);
    
    % Number of trials (should be first dimension)
    n_trials_eeg = size(alpha_power, 1);
    
    % Final verification: alpha_power should be [n_trials, n_times]
    if size(alpha_power, 2) ~= length(alpha_time)
        warning('  Dimension check failed for subject %s: alpha_power=[%d %d], alpha_time length=%d', ...
            subj_id, size(alpha_power, 1), size(alpha_power, 2), length(alpha_time));
        continue;
    end

    % Load gaze data
    datapath_gaze = fullfile(path, subj_id, 'gaze');
    if ~exist(datapath_gaze, 'dir')
        continue;
    end

    try
        cd(datapath_gaze);
        warning('off', 'MATLAB:load:errorClosingFile');
        load(config.gaze_file, 'gaze_x', 'gaze_y', 'trialinfo');
        warning('on', 'MATLAB:load:errorClosingFile');
    catch
        warning('  Missing gaze file for subject %s', subj_id);
        continue;
    end

    % Check if variables exist
    if ~exist('gaze_x', 'var') || ~exist('gaze_y', 'var')
        continue;
    end

    % Match trials between EEG and gaze
    n_trials_gaze = size(gaze_x, 2);
    n_trials = min(n_trials_eeg, n_trials_gaze);

    if n_trials == 0
        continue;
    end

    % Storage for this subject's coherence - use maximum possible windows
    % Calculate based on expected data length
    % Alpha data should be -0.5 to 2s, which at TFR resolution might be ~50 time points (0.05s steps)
    % After interpolation to match gaze (1251 samples at 500Hz), we'll have 1251 samples
    % Use a reasonable maximum: 1251 samples allows for many windows
    max_expected_length = 1251;  % Expected after alignment (-0.5 to 2s * 500Hz = 1251 samples)
    n_windows_max = floor((max_expected_length - coh_win_samp) / coh_step) + 1;
    if n_windows_max < 1
        n_windows_max = 1;  % At least 1 window
    end
    coh_subj = nan(n_trials, n_windows_max);

    %% Process each trial
    valid_trials = 0;
    for tr = 1:n_trials
        % Extract alpha power time series for this trial
        try
            alpha_ts = alpha_power(tr, :);  % [1 x time]
        catch
            continue;
        end

        % Extract gaze data for this trial
        try
            X = gaze_x{tr, 1};
            Y = gaze_y{tr, 1};
        catch
            continue;
        end

        if isempty(X) || isempty(Y) || length(X) ~= length(Y)
            continue;
        end

        % Extract -0.5 to 2s portion from gaze data (full window)
        % Gaze data is typically -0.5 to 2.0s at 500Hz = 1251 samples
        % We need all samples (1 to 1251) for -0.5 to 2s = 1251 samples
        if length(X) == 1251  % Standard case: -0.5 to 2.0s
            % Use all data - no trimming needed
            % X and Y already contain -0.5 to 2s
        elseif length(X) > 1250
            % If longer, take first 1251 samples (assuming it starts at -0.5s)
            X = X(1:1251);
            Y = Y(1:1251);
        elseif length(X) < 1251
            % If shorter, pad or interpolate to match expected length
            % This shouldn't happen normally, but handle it gracefully
            if length(X) >= 1000
                % Pad with last value to reach 1251
                X = [X, repmat(X(end), 1, 1251 - length(X))];
                Y = [Y, repmat(Y(end), 1, 1251 - length(Y))];
            else
                % Too short, skip this trial
                continue;
            end
        end

        % Create time vectors for interpolation
        % Gaze data: -0.5 to 2.0s at 500Hz (1251 samples, 0.002s spacing)
        % Alpha data: TFR time points (typically 0.05s steps = ~50 points)
        % Strategy: Interpolate ALPHA to match GAZE resolution (500Hz) to preserve temporal detail
        gaze_time = linspace(t_win(1), t_win(2), length(X));  % 1251 points at 500Hz
        
        % Interpolate alpha to match gaze time points (upsample alpha to 500Hz)
        if length(X) ~= length(alpha_ts)
            if length(alpha_ts) > 1 && length(X) > 1 && length(alpha_time) == length(alpha_ts)
                % Interpolate alpha to match gaze time points (upsample to 500Hz)
                alpha_ts = interp1(alpha_time, alpha_ts, gaze_time, 'linear', 'extrap');
                % Now alpha_ts has same length as X (1251 samples at 500Hz)
            else
                continue;
            end
        end

        % Remove NaNs
        %valid_idx = isfinite(X) & isfinite(Y) & isfinite(alpha_ts);
        %if sum(valid_idx) < coh_win_samp
        %    continue;
        %end

        % Compute eye velocity
        [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd);

        % Handle outliers
        zvx = (vx - nanmean(vx)) / (nanstd(vx) + eps);
        zvy = (vy - nanmean(vy)) / (nanstd(vy) + eps);
        bad = abs(zvx) > velZthr | abs(zvy) > velZthr;
        if any(bad)
            vx(bad) = NaN;
            vy(bad) = NaN;
            vx = fillmissing(vx, 'linear', 'EndValues', 'nearest');
            vy = fillmissing(vy, 'linear', 'EndValues', 'nearest');
        end

        % Total velocity magnitude
        vel = hypot(vx, vy);

        % Ensure same length as alpha (both should be 1251 samples at 500Hz after interpolation)
        if length(vel) ~= length(alpha_ts)
            if length(vel) > 1 && length(alpha_ts) > 1
                % Interpolate velocity to match alpha time points (should match gaze_time now)
                vel_time = linspace(t_win(1), t_win(2), length(vel));
                gaze_time_aligned = linspace(t_win(1), t_win(2), length(alpha_ts));
                vel = interp1(vel_time, vel, gaze_time_aligned, 'linear', 'extrap');
            else
                continue;
            end
        end

        % Apply valid mask
        %alpha_clean = alpha_ts(valid_idx);
        %vel_clean = vel(valid_idx);
        alpha_clean = alpha_ts;
        vel_clean = vel;

        if length(alpha_clean) < coh_win_samp
            continue;
        end

        % Compute time-resolved coherence using sliding windows
        % Ensure we can fit at least one full window
        if length(alpha_clean) < coh_win_samp
            continue;
        end
        
        % Calculate number of windows that fit
        % Last window should end before or at the last sample
        n_windows_trial = floor((length(alpha_clean) - coh_win_samp) / coh_step) + 1;
        coh_trial = nan(1, n_windows_trial);

        for w = 1:n_windows_trial
            win_start = (w - 1) * coh_step + 1;
            win_end = win_start + coh_win_samp - 1;

            % Ensure window doesn't exceed data length
            if win_end > length(alpha_clean)
                % Try to use the last possible window ending at the last sample
                if w == 1
                    % Can't even fit one window
                    break;
                end
                % Use the last valid window
                break;
            end

            alpha_win = alpha_clean(win_start:win_end);
            vel_win = vel_clean(win_start:win_end);

            % Remove any remaining NaNs
            valid_win = isfinite(alpha_win) & isfinite(vel_win);
            if sum(valid_win) < coh_win_samp * 0.8  % At least 80% valid
                continue;
            end

            alpha_win = alpha_win(valid_win);
            vel_win = vel_win(valid_win);

            % Validate signals have sufficient length and variance
            if length(alpha_win) < 64  % Need at least 64 samples for reasonable coherence estimation
                continue;
            end
            
            % Validate signals have sufficient variance (avoid constant signals)
            var_alpha = nanvar(alpha_win);
            var_vel = nanvar(vel_win);
            if var_alpha < eps || var_vel < eps
                % One or both signals are constant - skip this window
                continue;
            end
            
            % Check if signals are identical (would give coherence = 1 artificially)
            if max(abs(alpha_win - vel_win)) < eps * max(abs(alpha_win))
                % Signals are essentially identical - skip to avoid artificial coherence = 1
                continue;
            end

            % Compute coherence using mscohere
            % After interpolation, both alpha and gaze are at 500Hz (fs)
            % Use explicit window parameters for better control
            try
                % For mscohere, use a shorter window to allow multiple segments
                % If data length is 125 samples, use ~64 sample window for Welch's method
                mscohere_win_len = min(64, floor(length(alpha_win) / 2));  % Use half the data length or 64, whichever is smaller
                if mod(mscohere_win_len, 2) == 0
                    mscohere_win_len = mscohere_win_len - 1;  % Make odd for better spectral estimation
                end
                mscohere_win_len = max(mscohere_win_len, 33);  % Minimum reasonable window size
                
                % Window parameters: use Hamming window with 50% overlap
                noverlap = round(mscohere_win_len * 0.5);  % 50% overlap
                nfft = 2^nextpow2(mscohere_win_len);  % Next power of 2 for FFT efficiency
                window = hamming(mscohere_win_len);  % Hamming window
                
                [Cxy, f] = mscohere(alpha_win, vel_win, window, noverlap, nfft, fs);
                
                % Filter to alpha band frequencies only (8-14 Hz)
                alpha_freq_idx = f >= alpha_band(1) & f <= alpha_band(2);
                if any(alpha_freq_idx)
                    coh_val = mean(Cxy(alpha_freq_idx));  % Average coherence in alpha band
                else
                    % Fallback: use all frequencies if alpha band not in range
                    coh_val = mean(Cxy);
                end
                
                % Clamp coherence to [0, 0.999] to avoid atanh(1) = Inf in Fisher z-transform
                coh_val = max(0, min(0.999, coh_val));
                coh_trial(w) = coh_val;
            catch
                % Fallback: simple correlation
                if length(alpha_win) > 10
                    r = corrcoef(alpha_win, vel_win);
                    if ~isnan(r(1,2))
                        coh_val = abs(r(1,2));
                        % Clamp to [0, 0.999] to avoid atanh(1) = Inf
                        coh_val = max(0, min(0.999, coh_val));
                        coh_trial(w) = coh_val;
                    end
                end
            end
        end

        if any(isfinite(coh_trial))
            % Store coherence - use actual length of coh_trial
            n_coh = length(coh_trial);
            if n_coh > 0
                % Store as many as we can fit
                n_store = min(n_coh, size(coh_subj, 2));
                if n_store > 0
                    coh_subj(tr, 1:n_store) = coh_trial(1:n_store);
                    valid_trials = valid_trials + 1;
                end
            end
        end
    end

    if valid_trials > 0
        % Find maximum number of windows actually used across all trials
        max_windows_used = 0;
        for tr = 1:n_trials
            if any(isfinite(coh_subj(tr, :)))
                last_valid = find(isfinite(coh_subj(tr, :)), 1, 'last');
                if ~isempty(last_valid) && last_valid > max_windows_used
                    max_windows_used = last_valid;
                end
            end
        end

        if max_windows_used > 0
            % Trim to actual used windows
            coh_subj = coh_subj(:, 1:max_windows_used);

            n_valid = n_valid + 1;
            coh_all_subj{end+1} = coh_subj;
        end
    end

    clear gaze_x gaze_y trialinfo tfr_all
end

fprintf('Valid subjects: %d\n', n_valid);

if n_valid == 0
    warning('No valid subjects for %s task. Skipping.', task_name);
    return;
end

%% Aggregate across subjects

% Find maximum number of windows across all subjects
max_windows_all = 0;
for s = 1:n_valid
    n_w = size(coh_all_subj{s}, 2);
    if n_w > max_windows_all
        max_windows_all = n_w;
    end
end

% Compute time centers for all windows
% After interpolation, both alpha and gaze are at 500Hz (1251 samples from -0.5 to 2.0s)
% Windows are calculated at 500Hz, so time centers are straightforward
% With 200ms windows and 50% overlap (100ms step), centers are at:
% -0.4, -0.3, -0.2, ..., 1.9 (starting from -0.5s, first center is -0.5 + 0.1 = -0.4)
coh_time_centers = zeros(1, max_windows_all);
for w = 1:max_windows_all
    win_start_sample = (w - 1) * coh_step + 1;  % Starting sample (1-based, at 500Hz)
    win_center_sample = win_start_sample + coh_win_samp / 2;  % Center sample
    
    % Convert sample index to time: data starts at t_win(1) = -0.5s, sampled at fs (500Hz)
    % Sample 1 corresponds to t = -0.5s, so sample n corresponds to t = -0.5 + (n-1)/fs
    coh_time_centers(w) = (win_center_sample - 1) / fs + t_win(1);
end

% Average coherence across trials per subject, then across subjects
n_windows = max_windows_all;
coh_subj_means = nan(n_valid, n_windows);

for s = 1:n_valid
    coh_subj_trials = coh_all_subj{s};
    % Align to common time grid (handle different lengths)
    n_windows_subj = size(coh_subj_trials, 2);
    n_use = min(n_windows_subj, n_windows);
    if n_use > 0
        coh_subj_means(s, 1:n_use) = nanmean(coh_subj_trials(:, 1:n_use), 1);
    end
end

% Group statistics
% Use Fisher z-transform for proper averaging of coherence values
% Coherence is bounded [0,1], so Fisher z-transform makes distribution more normal
% IMPORTANT: Clamp coherence to [0, 0.999] to avoid atanh(1) = Inf
coh_subj_means_clamped = coh_subj_means;
coh_subj_means_clamped(coh_subj_means_clamped >= 1.0) = 0.999;
coh_subj_means_clamped(coh_subj_means_clamped < 0) = 0;
coh_subj_means_z = atanh(coh_subj_means_clamped);  % Fisher z-transform
coh_mean_z = nanmean(coh_subj_means_z, 1);
coh_sem_z = nanstd(coh_subj_means_z, [], 1) / sqrt(n_valid);
coh_mean = tanh(coh_mean_z);  % Transform back to coherence scale
coh_sem = tanh(coh_sem_z);  % Approximate SEM on coherence scale (not exact but reasonable)

% Statistical testing: one-sample t-test against 0 at each time point
p_vals = nan(1, n_windows);
t_vals = nan(1, n_windows);
for w = 1:n_windows
    valid_data = coh_subj_means_z(:, w);
    valid_data = valid_data(isfinite(valid_data));
    if length(valid_data) > 2
        [~, p_vals(w), ~, stats] = ttest(valid_data);
        t_vals(w) = stats.tstat;
    end
end

% FDR correction for multiple comparisons
q_vals = nan(size(p_vals));
valid_p = p_vals(isfinite(p_vals));
if ~isempty(valid_p)
    q_vals(isfinite(p_vals)) = mafdr(valid_p, 'BHFDR', true);
end
sig_mask = q_vals < 0.05;

%% Plot
close all
figure('Position', [0 0 1400 800], 'Color', 'w');

% Main plot
hold on

% Shaded error bar
fill([coh_time_centers, fliplr(coh_time_centers)], ...
    [coh_mean - coh_sem, fliplr(coh_mean + coh_sem)], ...
    colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Shade significant regions
if any(sig_mask)
    d_sig = diff([0, sig_mask, 0]);
    on_sig = find(d_sig == 1);
    off_sig = find(d_sig == -1) - 1;
    
    y_max_plot = max(coh_mean + coh_sem);
    if isfinite(y_max_plot) && y_max_plot > 0
        ylim_sig = [0, y_max_plot * 1.1];
    else
        ylim_sig = [0, 1];
    end
    
    for k = 1:numel(on_sig)
        if on_sig(k) <= length(coh_time_centers) && off_sig(k) <= length(coh_time_centers)
            x0 = coh_time_centers(on_sig(k));
            x1 = coh_time_centers(off_sig(k));
            patch([x0 x1 x1 x0], [ylim_sig(1) ylim_sig(1) ylim_sig(2) ylim_sig(2)], ...
                [1 1 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        end
    end
end

% Mean line
h_mean = plot(coh_time_centers, coh_mean, '-', 'Color', colors(1,:), 'LineWidth', 3);
uistack(h_mean, 'top');  % Ensure line is on top of significance shading

xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.6);
xlim([-.25 1.8]);
xticks([-.25 0 .25 .5 .75 1 1.25 1.5 1.75])

% Set ylim only if we have valid data
if any(isfinite(coh_mean)) && any(isfinite(coh_sem))
    y_max = max(coh_mean + coh_sem);
    if isfinite(y_max) && y_max > 0
        ylim([0, y_max * 1.1]);
    end
end

xlabel('Time [s]', 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel('Coherence (Alpha Power × Eye Velocity)', 'FontSize', fontSize, 'FontWeight', 'bold');
title(sprintf('%s: Alpha-Velocity Coherence (n=%d subjects)', upper(task_name), n_valid), ...
    'FontSize', fontSize+2, 'FontWeight', 'bold');
set(gca, 'GridAlpha', 0.15, 'FontSize', fontSize-2);

% Save figure
fig_name = sprintf('AOC_omnibus_alpha_velocity_coherence_%s.png', task_name);
saveas(gcf, fullfile(output_dir, fig_name));
fprintf('Saved: %s\n', fig_name);

%% Summary statistics
fprintf('\n=== Summary Statistics (%s) ===\n', task_name);
fprintf('Mean coherence: %.4f (SD=%.4f, range=[%.4f, %.4f])\n', ...
    nanmean(coh_mean), nanstd(coh_mean), min(coh_mean), max(coh_mean));
fprintf('Max coherence: %.4f at t=%.2fs\n', max(coh_mean), coh_time_centers(coh_mean == max(coh_mean)));

% Report significant time points
n_sig = sum(sig_mask);
if n_sig > 0
    fprintf('Significant time points (FDR-corrected p<0.05): %d/%d\n', n_sig, n_windows);
    sig_times = coh_time_centers(sig_mask);
    fprintf('Significant time range: %.3f to %.3f s\n', min(sig_times), max(sig_times));
else
    fprintf('No significant time points found (FDR-corrected p<0.05)\n');
end

fprintf('\n=== Done ===\n');

%% Helper function: Compute velocity using Savitzky-Golay filter
function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
% Convert to double precision (required by sgolayfilt)
X = double(X);
Y = double(Y);

Ts = 1/fs;
L = numel(X);

% Determine frame length (must be odd)
framelen = min(21, L);
if mod(framelen, 2) == 0
    framelen = framelen - 1;
end

minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0
    minLegal = minLegal + 1;
end
if framelen < minLegal
    framelen = minLegal;
end
if framelen > L
    framelen = L - mod(L, 2) + 1;
end

useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / (Ts^1)) * G(:, 2)';  % 1st-derivative kernel
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end
