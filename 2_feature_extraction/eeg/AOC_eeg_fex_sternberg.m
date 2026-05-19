%% AOC EEG Feature Extraction — Sternberg
% Computes subject-level NON-FOOOF EEG features from preprocessed Sternberg EEG.
% This script is the non-FOOOF branch in the split pipeline and saves
% `eeg_data_sternberg` (MAT + CSV). FOOOF features are produced in
% AOC_eeg_fex_sternberg_TFR.m.
%
% Extracted features:
%   Power Spectrum (Early [0 1], Late [1 2], Full [0 2]) + Baseline [-0.5 -0.25]
%   Baselined spectra (dB) for each window
%   IAF (condition-wise from full window)
%   Alpha power in IAF band (early/late/full; raw + dB)
%   Lateralization index (late baselined)

%% POWSPCTRM (Baseline + Early/Late/Full)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

logDir = paths.logs;
scriptName = 'AOC_eeg_fex_sternberg';

for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG FEX - STERNBERG] Windowed power spectrum extraction for Subject %d / %d \n', subj, length(subjects))
        % Load data
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg   % time-domain FieldTrip data struct (used for tf transforms)

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo(:, 1) == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo(:, 1) == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo(:, 1) == 26); % WM load 6

        % Retention windows [s] — suffix mapping for saved variables
        window.base  = [-0.5 -0.25];
        window.early = [0 1];
        window.late  = [1 2];
        window.full  = [0 2];

        % ----------------------
        % Time-frequency transform (raw) per condition (keeps time for window selection)
        % mtmconvol: Hanning taper, fixed 500 ms window length at each foi
        % ----------------------
        cfg            = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 2:2:40;
        cfg.t_ftimwin  = ones(size(cfg.foi)) * 0.5;  % 500 ms window for all frequencies
        cfg.toi        = -1.5:0.05:3;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR); % chan_freq_time
        cfg.trials = ind4; tfr4 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind6; tfr6 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = window.base;
        cfgb.baselinetype = 'db';
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr4_bl = ft_freqbaseline(cfgb, tfr4);
        tfr6_bl = ft_freqbaseline(cfgb, tfr6);

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [2 40];
        pow2_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr2));
        pow4_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr4));
        pow6_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr6));
        pow2_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr2_bl));
        pow4_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr4_bl));
        pow6_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr6_bl));

        pow2_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr2));
        pow4_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr4));
        pow6_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr6));
        pow2_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr2_bl));
        pow4_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr4_bl));
        pow6_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr6_bl));

        pow2_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr2));
        pow4_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr4));
        pow6_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr6));
        pow2_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr2_bl));
        pow4_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr4_bl));
        pow6_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr6_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_stern_windows.mat', ...
            'pow2_raw_early','pow4_raw_early','pow6_raw_early', ...
            'pow2_raw_late','pow4_raw_late','pow6_raw_late', ...
            'pow2_raw_full','pow4_raw_full','pow6_raw_full', ...
            'pow2_bl_early','pow4_bl_early','pow6_bl_early', ...
            'pow2_bl_late','pow4_bl_late','pow6_bl_late', ...
            'pow2_bl_full','pow4_bl_full','pow6_bl_full')

    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Define channels
subj = 1;
datapath = fullfile(path, subjects{subj}, 'eeg');
cd(datapath);
load('power_stern_windows.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(pow2_raw_full.label)
    label = pow2_raw_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;
% Left and right channels
left_channels = {};
right_channels = {};
for i = 1:length(channels)
    try
        ch = channels{i};
        % Find the first numeric part in the channel name
        numStr = regexp(ch, '\d+', 'match');
        % Convert the first numerical token to a number
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch;
        else
            right_channels{end+1} = ch;
        end
    catch ME
        ME.message
        disp(['Midline channel: ', ch])
    end
end

% Load data and calculate alpha power, IAF and lateralization index
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'IAF', {}, 'Lateralization', {}, ...
    'AlphaPower_raw_early', {}, 'AlphaPower_raw_late', {}, 'AlphaPower_raw_full', {}, ...
    'AlphaPower_bl_early', {}, 'AlphaPower_bl_late', {}, 'AlphaPower_bl_full', {});

for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG FEX - STERNBERG] Alpha power, IAF and lateralization for Subject %d / %d \n', subj, length(subjects))
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath);
        load('power_stern_windows.mat');

        % Channel selection
        channelIdx = find(ismember(pow2_raw_full.label, channels));

        % Extract FULL-window power spectra for selected channels (used for IAF)
        powspctrm2 = mean(pow2_raw_full.powspctrm(channelIdx, :), 1);
        powspctrm4 = mean(pow4_raw_full.powspctrm(channelIdx, :), 1);
        powspctrm6 = mean(pow6_raw_full.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(pow2_raw_full.freq >= alphaRange(1) & pow2_raw_full.freq <= alphaRange(2));

        % Calculate IAF for WM load 2
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = pow2_raw_full.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(pow2_raw_full.freq > (IAF2-4) & pow2_raw_full.freq < (IAF2+2));
            powerIAF2 = mean(powspctrm2(IAF_range2));
        end

        % Calculate IAF for WM load 4
        alphaPower4 = powspctrm4(alphaIndices);
        [pks4,locs] = findpeaks(alphaPower4);
        if isempty(pks4)
            IAF4 = NaN;
            IAF_range4 = NaN;
            powerIAF4 = NaN;
        else
            [~, ind] = max(pks4);
            IAF4 = pow4_raw_full.freq(alphaIndices(locs(ind)));
            IAF_range4 = find(pow4_raw_full.freq > (IAF4-4) & pow4_raw_full.freq < (IAF4+2));
            powerIAF4 = mean(powspctrm4(IAF_range4));
        end

        % Calculate IAF for WM load 6
        alphaPower6 = powspctrm6(alphaIndices);
        [pks6,locs] = findpeaks(alphaPower6);
        if isempty(pks6)
            IAF6 = NaN;
            IAF_range6 = NaN;
            powerIAF6 = NaN;
        else
            [~, ind] = max(pks6);
            IAF6 = pow6_raw_full.freq(alphaIndices(locs(ind)));
            IAF_range6 = find(pow6_raw_full.freq > (IAF6-4) & pow6_raw_full.freq < (IAF6+2));
            powerIAF6 = mean(powspctrm6(IAF_range6));
        end

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF4 == alphaRange(1) || IAF4 == alphaRange(2)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if IAF6 == alphaRange(1) || IAF6 == alphaRange(2)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF4 > max(pks4)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if powerIAF6 > max(pks6)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Compute lateralization index on LATE BASELINED spectra (dB)
        powloads = {pow2_bl_late, pow4_bl_late, pow6_bl_late};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx2 = LatIdx(1);
        LatIdx4 = LatIdx(2);
        LatIdx6 = LatIdx(3);

        % Alpha power in IAF band (raw + baselined) for early/late/full
        IAF_band2 = [IAF2-4 IAF2+2]; if any(isnan(IAF_band2)); IAF_band2 = alphaRange; end
        IAF_band4 = [IAF4-4 IAF4+2]; if any(isnan(IAF_band4)); IAF_band4 = alphaRange; end
        IAF_band6 = [IAF6-4 IAF6+2]; if any(isnan(IAF_band6)); IAF_band6 = alphaRange; end

        AlphaPower_raw_early = [robust_roi_pow(pow2_raw_early, channelIdx, IAF_band2); robust_roi_pow(pow4_raw_early, channelIdx, IAF_band4); robust_roi_pow(pow6_raw_early, channelIdx, IAF_band6)];
        AlphaPower_raw_late  = [robust_roi_pow(pow2_raw_late,  channelIdx, IAF_band2); robust_roi_pow(pow4_raw_late,  channelIdx, IAF_band4); robust_roi_pow(pow6_raw_late,  channelIdx, IAF_band6)];
        AlphaPower_raw_full  = [robust_roi_pow(pow2_raw_full,  channelIdx, IAF_band2); robust_roi_pow(pow4_raw_full,  channelIdx, IAF_band4); robust_roi_pow(pow6_raw_full,  channelIdx, IAF_band6)];
        AlphaPower_bl_early  = [robust_roi_pow(pow2_bl_early, channelIdx, IAF_band2); robust_roi_pow(pow4_bl_early, channelIdx, IAF_band4); robust_roi_pow(pow6_bl_early, channelIdx, IAF_band6)];
        AlphaPower_bl_late   = [robust_roi_pow(pow2_bl_late,  channelIdx, IAF_band2); robust_roi_pow(pow4_bl_late,  channelIdx, IAF_band4); robust_roi_pow(pow6_bl_late,  channelIdx, IAF_band6)];
        AlphaPower_bl_full   = [robust_roi_pow(pow2_bl_full,  channelIdx, IAF_band2); robust_roi_pow(pow4_bl_full,  channelIdx, IAF_band4); robust_roi_pow(pow6_bl_full,  channelIdx, IAF_band6)];

        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'IAF', num2cell([IAF2; IAF4; IAF6]), 'Lateralization', num2cell([LatIdx2; LatIdx4; LatIdx6]), ...
            'AlphaPower_raw_early', num2cell(AlphaPower_raw_early), ...
            'AlphaPower_raw_late', num2cell(AlphaPower_raw_late), ...
            'AlphaPower_raw_full', num2cell(AlphaPower_raw_full), ...
            'AlphaPower_bl_early', num2cell(AlphaPower_bl_early), ...
            'AlphaPower_bl_late', num2cell(AlphaPower_bl_late), ...
            'AlphaPower_bl_full', num2cell(AlphaPower_bl_full));
        % Save
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_sternberg_subj subj_data_eeg
        save alpha_power_sternberg powerIAF2 powerIAF4 powerIAF6
        save IAF_sternberg IAF2 IAF4 IAF6
        save lateralization_sternberg LatIdx2 LatIdx4 LatIdx6
        eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: WM2: %f Hz (Power: %f), WM4: %f Hz (Power: %f), ' ...
            'WM6: %f Hz (Power: %f) | Lateralization: %f %f %f \n'], subjects{subj}, ...
            IAF2, powerIAF2, IAF4, powerIAF4, IAF6, powerIAF6, LatIdx2, LatIdx4, LatIdx6);
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end
save(fullfile(paths.features, 'AOC_eeg_matrix_sternberg.mat'), 'eeg_data_sternberg')
writetable(struct2table(eeg_data_sternberg), fullfile(paths.features, 'AOC_eeg_matrix_sternberg.csv'))

function v = robust_roi_pow(S, channelIdx, band)
fmask = S.freq >= band(1) & S.freq <= band(2);
if ~any(fmask)
    v = NaN;
    return
end

x = S.powspctrm(channelIdx, fmask);
x = x(:);
x = x(isfinite(x));
% Hard plausibility guard to suppress catastrophic numeric explosions.
x = x(abs(x) <= 1e4);
if numel(x) >= 8
    q1 = prctile(x, 25);
    q3 = prctile(x, 75);
    iqr_v = q3 - q1;
    if isfinite(iqr_v) && iqr_v > 0
        lo = q1 - 3 * iqr_v;
        hi = q3 + 3 * iqr_v;
        x = x(x >= lo & x <= hi);
    end
end
if isempty(x)
    v = NaN;
else
    v = mean(x, 'omitnan');
end
end