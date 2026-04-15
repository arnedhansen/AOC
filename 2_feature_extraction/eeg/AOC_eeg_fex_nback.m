%% AOC EEG Feature Extraction — N-Back
% Computes subject-level NON-FOOOF EEG features from preprocessed N-back EEG.
% This script is the non-FOOOF branch in the split pipeline and saves
% `eeg_data_nback` (MAT + CSV). FOOOF features are produced in
% AOC_eeg_fex_nback_TFR.m.
%
% Extracted features:
%   Power Spectrum windows (raw + baselined)
%   IAF (condition-wise from full window)
%   Alpha power in IAF band (early/late/full; raw + dB)
%   Lateralization index (late baselined)

%% POWSPCTRM (Baseline + Early/Late/Full) (subject-level)
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
scriptName = 'AOC_eeg_fex_nback';

for subj = 1:length(subjects)
    try
        clc
        disp(['Processing windowed POWSPCTRM (subject-level) for Subject AOC ', num2str(subjects{subj})])
        % Load data
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(dataTFR.trialinfo(:, 1) == 21);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind3 = find(dataTFR.trialinfo(:, 1) == 23);

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
        cfg.toi        = -1.25:0.05:2.25;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind1; tfr1 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind3; tfr3 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = t_base;
        cfgb.baselinetype = 'db';
        tfr1_bl = ft_freqbaseline(cfgb, tfr1);
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr3_bl = ft_freqbaseline(cfgb, tfr3);

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [3 30];
        pow1_early     = remove_time_dimension(select_data(t_early, freq_range, tfr1));
        pow2_early     = remove_time_dimension(select_data(t_early, freq_range, tfr2));
        pow3_early     = remove_time_dimension(select_data(t_early, freq_range, tfr3));
        pow1_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr1_bl));
        pow2_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr2_bl));
        pow3_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr3_bl));

        pow1_late      = remove_time_dimension(select_data(t_late, freq_range, tfr1));
        pow2_late      = remove_time_dimension(select_data(t_late, freq_range, tfr2));
        pow3_late      = remove_time_dimension(select_data(t_late, freq_range, tfr3));
        pow1_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr1_bl));
        pow2_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr2_bl));
        pow3_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr3_bl));

        pow1_full      = remove_time_dimension(select_data(t_full, freq_range, tfr1));
        pow2_full      = remove_time_dimension(select_data(t_full, freq_range, tfr2));
        pow3_full      = remove_time_dimension(select_data(t_full, freq_range, tfr3));
        pow1_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr1_bl));
        pow2_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr2_bl));
        pow3_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr3_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_nback_windows.mat', ...
            'pow1_early','pow2_early','pow3_early', ...
            'pow1_late','pow2_late','pow3_late', ...
            'pow1_full','pow2_full','pow3_full', ...
            'pow1_early_bl','pow2_early_bl','pow3_early_bl', ...
            'pow1_late_bl','pow2_late_bl','pow3_late_bl', ...
            'pow1_full_bl','pow2_full_bl','pow3_full_bl')

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
datapath = fullfile(path, subjects{1}, 'eeg');
cd(datapath);
load('power_nback_windows.mat');
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

% Load data and calculate power, IAF and lateralization index
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {}, ...
    'AlphaPowerEarly', {}, 'AlphaPowerLate', {}, 'AlphaPowerFull', {}, ...
    'AlphaPowerEarlyBL', {}, 'AlphaPowerLateBL', {}, 'AlphaPowerFullBL', {});

for subj = 1:length(subjects)
    try
        clc
        disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath);
        load('power_nback_windows.mat');

        % Channels selection based on CBPT
        channelIdx = find(ismember(pow1_full.label, channels));

        % FULL-window power spectra for IAF
        powspctrm1 = mean(pow1_full.powspctrm(channelIdx, :), 1);
        powspctrm2 = mean(pow2_full.powspctrm(channelIdx, :), 1);
        powspctrm3 = mean(pow3_full.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(pow1_full.freq >= alphaRange(1) & pow1_full.freq <= alphaRange(2));

        % Calculate IAF for 1-back
        alphaPower1 = powspctrm1(alphaIndices);
        [pks1,locs] = findpeaks(alphaPower1);
        if isempty(pks1)
            IAF1 = NaN;
            IAF_range1 = NaN;
            powerIAF1 = NaN;
        else
            [~, ind] = max(pks1);
            IAF1 = pow1_full.freq(alphaIndices(locs(ind)));
            IAF_range1 = find(pow1_full.freq > (IAF1-4) & pow1_full.freq < (IAF1+2));
            powerIAF1 = mean(powspctrm1(IAF_range1));
        end

        % Calculate IAF for 2-back
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

        % Calculate IAF for 3-back
        alphaPower3 = powspctrm3(alphaIndices);
        [pks3,locs] = findpeaks(alphaPower3);
        if isempty(pks3)
            IAF3 = NaN;
            IAF_range3 = NaN;
            powerIAF3 = NaN;
        else
            [~, ind] = max(pks3);
            IAF3 = pow3_full.freq(alphaIndices(locs(ind)));
            IAF_range3 = find(pow3_full.freq > (IAF3-4) & pow3_full.freq < (IAF3+2));
            powerIAF3 = mean(powspctrm3(IAF_range3));
        end

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF1 == alphaRange(1) || IAF1 == alphaRange(2)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF3 == alphaRange(1) || IAF3 == alphaRange(2)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF1 > max(pks1)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF3 > max(pks3)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Compute lateralization index on LATE BASELINED spectra (dB)
        powloads = {pow1_late_bl, pow2_late_bl, pow3_late_bl};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx1 = LatIdx(1);
        LatIdx2 = LatIdx(2);
        LatIdx3 = LatIdx(3);

        % Alpha power in IAF band (raw + baselined) for early/late/full
        IAF_band1 = [IAF1-4 IAF1+2]; if any(isnan(IAF_band1)); IAF_band1 = alphaRange; end
        IAF_band2 = [IAF2-4 IAF2+2]; if any(isnan(IAF_band2)); IAF_band2 = alphaRange; end
        IAF_band3 = [IAF3-4 IAF3+2]; if any(isnan(IAF_band3)); IAF_band3 = alphaRange; end

        AlphaPowerEarly    = [robust_roi_pow(pow1_early, channelIdx, IAF_band1);    robust_roi_pow(pow2_early, channelIdx, IAF_band2);    robust_roi_pow(pow3_early, channelIdx, IAF_band3)];
        AlphaPowerLate     = [robust_roi_pow(pow1_late,  channelIdx, IAF_band1);    robust_roi_pow(pow2_late,  channelIdx, IAF_band2);    robust_roi_pow(pow3_late,  channelIdx, IAF_band3)];
        AlphaPowerFull     = [robust_roi_pow(pow1_full,  channelIdx, IAF_band1);    robust_roi_pow(pow2_full,  channelIdx, IAF_band2);    robust_roi_pow(pow3_full,  channelIdx, IAF_band3)];
        AlphaPowerEarlyBL  = [robust_roi_pow(pow1_early_bl, channelIdx, IAF_band1); robust_roi_pow(pow2_early_bl, channelIdx, IAF_band2); robust_roi_pow(pow3_early_bl, channelIdx, IAF_band3)];
        AlphaPowerLateBL   = [robust_roi_pow(pow1_late_bl,  channelIdx, IAF_band1); robust_roi_pow(pow2_late_bl,  channelIdx, IAF_band2); robust_roi_pow(pow3_late_bl,  channelIdx, IAF_band3)];
        AlphaPowerFullBL   = [robust_roi_pow(pow1_full_bl,  channelIdx, IAF_band1); robust_roi_pow(pow2_full_bl,  channelIdx, IAF_band2); robust_roi_pow(pow3_full_bl,  channelIdx, IAF_band3)];
        % Keep legacy AlphaPower as LATE raw (registered retention)
        AlphaPower = AlphaPowerLate;

        % Create a structure array for this subject (legacy fields preserved)
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'AlphaPower', num2cell(AlphaPower), 'IAF', num2cell([IAF1; IAF2; IAF3]), 'Lateralization', num2cell([LatIdx1; LatIdx2; LatIdx3]));

        tmp = num2cell(AlphaPowerEarly);   [subj_data_eeg.AlphaPowerEarly]   = tmp{:};
        tmp = num2cell(AlphaPowerLate);    [subj_data_eeg.AlphaPowerLate]    = tmp{:};
        tmp = num2cell(AlphaPowerFull);    [subj_data_eeg.AlphaPowerFull]    = tmp{:};
        tmp = num2cell(AlphaPowerEarlyBL); [subj_data_eeg.AlphaPowerEarlyBL] = tmp{:};
        tmp = num2cell(AlphaPowerLateBL);  [subj_data_eeg.AlphaPowerLateBL]  = tmp{:};
        tmp = num2cell(AlphaPowerFullBL);  [subj_data_eeg.AlphaPowerFullBL]  = tmp{:};
        % Save
        if ispc == 1
            savepath = fullfile('W:\Students\Arne\AOC\data\features\', subjects{subj}, 'eeg');
        else
            savepath = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, 'eeg');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        save lateralization_nback LatIdx1 LatIdx2 LatIdx3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), ' ...
            '3-back: %f Hz (Power: %f) |Lateralization: %f %f %f \n'], subjects{subj}, IAF1, ...
            powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3, LatIdx1, LatIdx2, LatIdx3);
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_nback eeg_data_nback
    writetable(struct2table(eeg_data_nback), 'W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_nback.csv')
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_nback eeg_data_nback
    writetable(struct2table(eeg_data_nback), '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_nback.csv')
end

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
