%% AOC EEG Feature Extraction — Sternberg
% Computes subject-level power spectra in early/late/full windows (raw + baselined),
% IAF, lateralization, and FOOOF-based alpha (raw + baselined) from preprocessed
% Sternberg EEG. Saves per-subject and grand results.
%
% Extracted features:
%   Power Spectrum (Early [0 1], Late [1 2], Full [0 2]) + Baseline [-0.5 -0.25]
%   Baselined spectra (dB) for each window
%   IAF (subject-level, pooled across conditions, full window)
%   Alpha power in IAF band (early/late/full; raw + dB)
%   Lateralization index (late baselined)
%   FOOOF-based spectra (model fit - aperiodic) for each window + baselined (absolute, log space)
%   FOOOF alpha power (8-14 Hz) full + baselined (full/early/late)

%% POWSPCTRM (Baseline + Early/Late/Full) + FOOOF window spectra
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Setup logging
if ispc == 1
    logDir = 'W:\Students\Arne\AOC\data\controls\logs';
else
    logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
end
scriptName = 'AOC_eeg_fex_sternberg';

for subj = 1:length(subjects)
    try
        clc
        disp(['Processing windowed POWSPCTRM + FOOOF (subject-level) for Subject AOC ', num2str(subjects{subj})])
        % Load data
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg   % time-domain FieldTrip data struct (used for tf transforms)

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo(:, 1) == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo(:, 1) == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo(:, 1) == 26); % WM load 6

        % Windows
        t_base = [-0.5 -0.25];
        t_early = [0 1];
        t_late  = [1 2];
        t_full  = [0 2];

        % ----------------------
        % Time-frequency transform (raw) per condition (keeps time for window selection)
        % ----------------------
        cfg            = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 3:1:30;
        cfg.t_ftimwin  = ones(size(cfg.foi))*1;     % 1 s windows
        cfg.toi        = -1.5:0.05:3;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR); % chan_freq_time
        cfg.trials = ind4; tfr4 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind6; tfr6 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = t_base;
        cfgb.baselinetype = 'db';
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr4_bl = ft_freqbaseline(cfgb, tfr4);
        tfr6_bl = ft_freqbaseline(cfgb, tfr6);

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [3 30];
        pow2_early     = remove_time_dimension(select_data(t_early, freq_range, tfr2));
        pow4_early     = remove_time_dimension(select_data(t_early, freq_range, tfr4));
        pow6_early     = remove_time_dimension(select_data(t_early, freq_range, tfr6));
        pow2_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr2_bl));
        pow4_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr4_bl));
        pow6_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr6_bl));

        pow2_late      = remove_time_dimension(select_data(t_late, freq_range, tfr2));
        pow4_late      = remove_time_dimension(select_data(t_late, freq_range, tfr4));
        pow6_late      = remove_time_dimension(select_data(t_late, freq_range, tfr6));
        pow2_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr2_bl));
        pow4_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr4_bl));
        pow6_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr6_bl));

        pow2_full      = remove_time_dimension(select_data(t_full, freq_range, tfr2));
        pow4_full      = remove_time_dimension(select_data(t_full, freq_range, tfr4));
        pow6_full      = remove_time_dimension(select_data(t_full, freq_range, tfr6));
        pow2_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr2_bl));
        pow4_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr4_bl));
        pow6_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr6_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_stern_windows.mat', ...
            'pow2_early','pow4_early','pow6_early', ...
            'pow2_late','pow4_late','pow6_late', ...
            'pow2_full','pow4_full','pow6_full', ...
            'pow2_early_bl','pow4_early_bl','pow6_early_bl', ...
            'pow2_late_bl','pow4_late_bl','pow6_late_bl', ...
            'pow2_full_bl','pow4_full_bl','pow6_full_bl')

        % ----------------------
        % FOOOF spectra per window (subject-level, no sliding)
        % Output convention: pow*_fooof = (model fit - aperiodic) in log space
        % Baselining: absolute difference (window - baseline) in log space
        %
        % NOTE:
        % This path uses one FFT per broad latency window and is more sensitive
        % to occasional catastrophic fits than the sliding-window TFR pipeline.
        % Hard guards are applied below before saving.
        % ----------------------
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = freq_range;
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';

        % Baseline window per condition (for baselining in log space)
        dat2_base = ft_selectdata(struct('latency', t_base, 'trials', ind2), dataTFR);
        dat4_base = ft_selectdata(struct('latency', t_base, 'trials', ind4), dataTFR);
        dat6_base = ft_selectdata(struct('latency', t_base, 'trials', ind6), dataTFR);

        fooof2_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_base);
        fooof4_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_base);
        fooof6_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_base);

        % Full window (0-2) non-baselined
        dat2_full = ft_selectdata(struct('latency', t_full, 'trials', ind2), dataTFR);
        dat4_full = ft_selectdata(struct('latency', t_full, 'trials', ind4), dataTFR);
        dat6_full = ft_selectdata(struct('latency', t_full, 'trials', ind6), dataTFR);
        pow2_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_full);
        pow4_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_full);
        pow6_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_full);

        % Baselined full/early/late (absolute in log space)
        pow2_fooof_bl = pow2_fooof; pow2_fooof_bl.powspctrm = pow2_fooof.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl = pow4_fooof; pow4_fooof_bl.powspctrm = pow4_fooof.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl = pow6_fooof; pow6_fooof_bl.powspctrm = pow6_fooof.powspctrm - fooof6_base.powspctrm;

        dat2_early = ft_selectdata(struct('latency', t_early, 'trials', ind2), dataTFR);
        dat4_early = ft_selectdata(struct('latency', t_early, 'trials', ind4), dataTFR);
        dat6_early = ft_selectdata(struct('latency', t_early, 'trials', ind6), dataTFR);
        fooof2_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_early);
        fooof4_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_early);
        fooof6_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_early);
        pow2_fooof_bl_early = fooof2_early; pow2_fooof_bl_early.powspctrm = fooof2_early.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl_early = fooof4_early; pow4_fooof_bl_early.powspctrm = fooof4_early.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl_early = fooof6_early; pow6_fooof_bl_early.powspctrm = fooof6_early.powspctrm - fooof6_base.powspctrm;

        dat2_late = ft_selectdata(struct('latency', t_late, 'trials', ind2), dataTFR);
        dat4_late = ft_selectdata(struct('latency', t_late, 'trials', ind4), dataTFR);
        dat6_late = ft_selectdata(struct('latency', t_late, 'trials', ind6), dataTFR);
        fooof2_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_late);
        fooof4_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_late);
        fooof6_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_late);
        pow2_fooof_bl_late = fooof2_late; pow2_fooof_bl_late.powspctrm = fooof2_late.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl_late = fooof4_late; pow4_fooof_bl_late.powspctrm = fooof4_late.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl_late = fooof6_late; pow6_fooof_bl_late.powspctrm = fooof6_late.powspctrm - fooof6_base.powspctrm;

        % Hard guards against catastrophic numeric explosions in FOOOF output.
        cfg_guard = struct('r2_min', 0.90, 'hard_abs', 1e4, 'winsor_prc', [1 99]);
        fooof2_base = sanitize_fooof_struct(fooof2_base, cfg_guard);
        fooof4_base = sanitize_fooof_struct(fooof4_base, cfg_guard);
        fooof6_base = sanitize_fooof_struct(fooof6_base, cfg_guard);
        pow2_fooof = sanitize_fooof_struct(pow2_fooof, cfg_guard);
        pow4_fooof = sanitize_fooof_struct(pow4_fooof, cfg_guard);
        pow6_fooof = sanitize_fooof_struct(pow6_fooof, cfg_guard);
        pow2_fooof_bl = sanitize_fooof_struct(pow2_fooof_bl, cfg_guard);
        pow4_fooof_bl = sanitize_fooof_struct(pow4_fooof_bl, cfg_guard);
        pow6_fooof_bl = sanitize_fooof_struct(pow6_fooof_bl, cfg_guard);
        pow2_fooof_bl_early = sanitize_fooof_struct(pow2_fooof_bl_early, cfg_guard);
        pow4_fooof_bl_early = sanitize_fooof_struct(pow4_fooof_bl_early, cfg_guard);
        pow6_fooof_bl_early = sanitize_fooof_struct(pow6_fooof_bl_early, cfg_guard);
        pow2_fooof_bl_late = sanitize_fooof_struct(pow2_fooof_bl_late, cfg_guard);
        pow4_fooof_bl_late = sanitize_fooof_struct(pow4_fooof_bl_late, cfg_guard);
        pow6_fooof_bl_late = sanitize_fooof_struct(pow6_fooof_bl_late, cfg_guard);

        save('power_stern_fooof_windows.mat', ...
            'pow2_fooof','pow4_fooof','pow6_fooof', ...
            'pow2_fooof_bl','pow4_fooof_bl','pow6_fooof_bl', ...
            'pow2_fooof_bl_early','pow4_fooof_bl_early','pow6_fooof_bl_early', ...
            'pow2_fooof_bl_late','pow4_fooof_bl_late','pow6_fooof_bl_late')

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
for i = 1:length(pow2_full.label)
    label = pow2_full.label{i};
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
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {}, ...
    'AlphaPowerEarly', {}, 'AlphaPowerLate', {}, 'AlphaPowerFull', {}, ...
    'AlphaPowerEarlyBL', {}, 'AlphaPowerLateBL', {}, 'AlphaPowerFullBL', {}, ...
    'AlphaPower_FOOOF', {}, 'AlphaPower_FOOOF_bl', {}, 'AlphaPower_FOOOF_bl_early', {}, 'AlphaPower_FOOOF_bl_late', {});

for subj = 1:length(subjects)
    try
        clc
        disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath);
        load('power_stern_windows.mat');
        load('power_stern_fooof_windows.mat');

        % Channel selection
        channelIdx = find(ismember(pow2_full.label, channels));

        % Extract FULL-window power spectra for selected channels (used for IAF)
        powspctrm2 = mean(pow2_full.powspctrm(channelIdx, :), 1);
        powspctrm4 = mean(pow4_full.powspctrm(channelIdx, :), 1);
        powspctrm6 = mean(pow6_full.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(pow2_full.freq >= alphaRange(1) & pow2_full.freq <= alphaRange(2));

        % Calculate IAF for WM load 2
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = pow2_full.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(pow2_full.freq > (IAF2-4) & pow2_full.freq < (IAF2+2));
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
            IAF4 = pow4_full.freq(alphaIndices(locs(ind)));
            IAF_range4 = find(pow4_full.freq > (IAF4-4) & pow4_full.freq < (IAF4+2));
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
            IAF6 = pow6_full.freq(alphaIndices(locs(ind)));
            IAF_range6 = find(pow6_full.freq > (IAF6-4) & pow6_full.freq < (IAF6+2));
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
        powloads = {pow2_late_bl, pow4_late_bl, pow6_late_bl};
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

        AlphaPowerEarly    = [robust_roi_pow(pow2_early, channelIdx, IAF_band2);    robust_roi_pow(pow4_early, channelIdx, IAF_band4);    robust_roi_pow(pow6_early, channelIdx, IAF_band6)];
        AlphaPowerLate     = [robust_roi_pow(pow2_late,  channelIdx, IAF_band2);    robust_roi_pow(pow4_late,  channelIdx, IAF_band4);    robust_roi_pow(pow6_late,  channelIdx, IAF_band6)];
        AlphaPowerFull     = [robust_roi_pow(pow2_full,  channelIdx, IAF_band2);    robust_roi_pow(pow4_full,  channelIdx, IAF_band4);    robust_roi_pow(pow6_full,  channelIdx, IAF_band6)];
        AlphaPowerEarlyBL  = [robust_roi_pow(pow2_early_bl, channelIdx, IAF_band2); robust_roi_pow(pow4_early_bl, channelIdx, IAF_band4); robust_roi_pow(pow6_early_bl, channelIdx, IAF_band6)];
        AlphaPowerLateBL   = [robust_roi_pow(pow2_late_bl,  channelIdx, IAF_band2); robust_roi_pow(pow4_late_bl,  channelIdx, IAF_band4); robust_roi_pow(pow6_late_bl,  channelIdx, IAF_band6)];
        AlphaPowerFullBL   = [robust_roi_pow(pow2_full_bl,  channelIdx, IAF_band2); robust_roi_pow(pow4_full_bl,  channelIdx, IAF_band4); robust_roi_pow(pow6_full_bl,  channelIdx, IAF_band6)];
        AlphaPower_FOOOF         = [robust_roi_pow(pow2_fooof,          channelIdx, alphaRange); robust_roi_pow(pow4_fooof,          channelIdx, alphaRange); robust_roi_pow(pow6_fooof,          channelIdx, alphaRange)];
        AlphaPower_FOOOF_bl      = [robust_roi_pow(pow2_fooof_bl,       channelIdx, alphaRange); robust_roi_pow(pow4_fooof_bl,       channelIdx, alphaRange); robust_roi_pow(pow6_fooof_bl,       channelIdx, alphaRange)];
        AlphaPower_FOOOF_bl_early= [robust_roi_pow(pow2_fooof_bl_early, channelIdx, alphaRange); robust_roi_pow(pow4_fooof_bl_early, channelIdx, alphaRange); robust_roi_pow(pow6_fooof_bl_early, channelIdx, alphaRange)];
        AlphaPower_FOOOF_bl_late = [robust_roi_pow(pow2_fooof_bl_late,  channelIdx, alphaRange); robust_roi_pow(pow4_fooof_bl_late,  channelIdx, alphaRange); robust_roi_pow(pow6_fooof_bl_late,  channelIdx, alphaRange)];

        % Keep legacy AlphaPower as LATE raw (registered retention)
        AlphaPower = AlphaPowerLate;

        % Create a structure array for this subject (legacy fields preserved)
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'AlphaPower', num2cell(AlphaPower), 'IAF', num2cell([IAF2; IAF4; IAF6]), ...
            'Lateralization', num2cell([LatIdx2; LatIdx4; LatIdx6]));

        % Add windowed alpha power fields (raw + baselined, IAF band)
        tmp = num2cell(AlphaPowerEarly);   [subj_data_eeg.AlphaPowerEarly]   = tmp{:};
        tmp = num2cell(AlphaPowerLate);    [subj_data_eeg.AlphaPowerLate]    = tmp{:};
        tmp = num2cell(AlphaPowerFull);    [subj_data_eeg.AlphaPowerFull]    = tmp{:};
        tmp = num2cell(AlphaPowerEarlyBL); [subj_data_eeg.AlphaPowerEarlyBL] = tmp{:};
        tmp = num2cell(AlphaPowerLateBL);  [subj_data_eeg.AlphaPowerLateBL]  = tmp{:};
        tmp = num2cell(AlphaPowerFullBL);  [subj_data_eeg.AlphaPowerFullBL]  = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF);          [subj_data_eeg.AlphaPower_FOOOF]          = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl);       [subj_data_eeg.AlphaPower_FOOOF_bl]       = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_early); [subj_data_eeg.AlphaPower_FOOOF_bl_early] = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_late);  [subj_data_eeg.AlphaPower_FOOOF_bl_late]  = tmp{:};

        % Save
        if ispc == 1
            savepath = fullfile('W:\Students\Arne\AOC\data\features\', subjects{subj}, 'eeg');
        else
            savepath = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, 'eeg');
        end
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
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_sternberg eeg_data_sternberg
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_sternberg eeg_data_sternberg
end

function v = robust_roi_pow(S, channelIdx, band)
fmask = S.freq >= band(1) & S.freq <= band(2);
if ~any(fmask)
    v = NaN;
    return
end

function S = sanitize_fooof_struct(S, cfg)
if ~isfield(S, 'powspctrm') || isempty(S.powspctrm)
    return
end

X = S.powspctrm;
X(~isfinite(X)) = NaN;

if isfield(S, 'fooofparams') && ~isempty(S.fooofparams)
    fp = S.fooofparams;
    if iscell(fp), fp = fp{1}; end
    if isstruct(fp) && isfield(fp, 'r_squared')
        rsq = [fp.r_squared]';
        bad = ~isfinite(rsq) | rsq < cfg.r2_min;
        if isvector(bad) && numel(bad) == size(X, 1)
            X(bad, :) = NaN;
        end
    end
end

X(abs(X) > cfg.hard_abs) = NaN;

vals = X(isfinite(X));
if numel(vals) >= 20
    lo = prctile(vals, cfg.winsor_prc(1));
    hi = prctile(vals, cfg.winsor_prc(2));
    X = min(max(X, lo), hi);
end

S.powspctrm = X;
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