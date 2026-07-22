%% AOC EEG Feature Extraction — N-Back
% Computes subject-level EEG features from preprocessed N-back EEG.
% Saves `eeg_data_nback` (MAT + CSV) and subject-level TFR products.
%
% Extracted features:
%   Power Spectrum windows (raw + baselined) via mtmconvol (2 Hz foi grid)
%   IAF (condition-wise): mtmfft+DPSS on retention window, trial-averaged (N-back [0 2] s), findpeaks [8 14] Hz
%   Alpha power in (IAF-4, IAF+2) Hz (early/late/full; raw + dB); NaN if no valid IAF
%   Lateralization index (late baselined)
%   ERSD_early / ERSD_late / ERSD_full ((IAF-4, IAF+2) Hz on baselined TFR, occipital ROI; fallback [8 14] if no valid IAF)

%% POWSPCTRM (Baseline + Early/Late/Full) (subject-level)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG FEX - NBACK] Windowed power spectrum extraction for Subject %d / %d \n', subj, length(subjects))
        % Load data
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(dataTFR.trialinfo(:, 1) == 21);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind3 = find(dataTFR.trialinfo(:, 1) == 23);

        % Retention windows [s] — suffix mapping for saved variables
        window.base  = [-1.5 -0.5];
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
        cfg.toi        = -1.5:0.05:2.25;  % full baseline [-1.5 -0.5] coverage
        cfg.pad        = 'nextpow2';
        cfg.keeptrials = 'no';

        cfg.trials = ind1; tfr1 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind3; tfr3 = ft_freqanalysis(cfg, dataTFR);

        cfgb              = [];
        cfgb.baseline     = window.base;
        cfgb.baselinetype = 'db';
        tfr1_bl = ft_freqbaseline(cfgb, tfr1);
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr3_bl = ft_freqbaseline(cfgb, tfr3);

        % Save trial-averaged TFR outputs
        save('tfr_nback.mat', 'tfr1', 'tfr2', 'tfr3', 'tfr1_bl', 'tfr2_bl', 'tfr3_bl')

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [2 40];
        pow1_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr1));
        pow2_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr2));
        pow3_raw_early  = remove_time_dimension(select_data(window.early, freq_range, tfr3));
        pow1_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr1_bl));
        pow2_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr2_bl));
        pow3_bl_early   = remove_time_dimension(select_data(window.early, freq_range, tfr3_bl));

        pow1_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr1));
        pow2_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr2));
        pow3_raw_late   = remove_time_dimension(select_data(window.late, freq_range, tfr3));
        pow1_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr1_bl));
        pow2_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr2_bl));
        pow3_bl_late    = remove_time_dimension(select_data(window.late, freq_range, tfr3_bl));

        pow1_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr1));
        pow2_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr2));
        pow3_raw_full   = remove_time_dimension(select_data(window.full, freq_range, tfr3));
        pow1_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr1_bl));
        pow2_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr2_bl));
        pow3_bl_full    = remove_time_dimension(select_data(window.full, freq_range, tfr3_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_nback_windows.mat', ...
            'pow1_raw_early','pow2_raw_early','pow3_raw_early', ...
            'pow1_raw_late','pow2_raw_late','pow3_raw_late', ...
            'pow1_raw_full','pow2_raw_full','pow3_raw_full', ...
            'pow1_bl_early','pow2_bl_early','pow3_bl_early', ...
            'pow1_bl_late','pow2_bl_late','pow3_bl_late', ...
            'pow1_bl_full','pow2_bl_full','pow3_bl_full')

    catch ME
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
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'IAF', {}, 'Lateralization', {}, ...
    'ERSD_early', {}, 'ERSD_late', {}, 'ERSD_full', {}, ...
    'AlphaPower_raw_early', {}, 'AlphaPower_raw_late', {}, 'AlphaPower_raw_full', {}, ...
    'AlphaPower_bl_early', {}, 'AlphaPower_bl_late', {}, 'AlphaPower_bl_full', {});

for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG FEX - NBACK] Alpha power, IAF and lateralization for Subject %d / %d \n', subj, length(subjects))
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath);
        load('power_nback_windows.mat');
        load dataEEG_TFR_nback

        ind1 = find(dataTFR.trialinfo(:, 1) == 21);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind3 = find(dataTFR.trialinfo(:, 1) == 23);

        % Channels selection based on CBPT
        channelIdx = find(ismember(pow1_raw_full.label, channels));
        chLabs = pow1_raw_full.label(channelIdx);

        % IAF: trial-averaged mtmfft on full retention [0 2] s, occipital ROI (~0.5 Hz resolution)
        winIAF = [0 2];
        [IAF1, powerIAF1] = iaf_from_retention_mtmfft(dataTFR, ind1, winIAF, chLabs, alphaRange);
        [IAF2, powerIAF2] = iaf_from_retention_mtmfft(dataTFR, ind2, winIAF, chLabs, alphaRange);
        [IAF3, powerIAF3] = iaf_from_retention_mtmfft(dataTFR, ind3, winIAF, chLabs, alphaRange);

        % ERSD from cached baselined TFR (avoid duplicate spectral transforms)
        % Band is subject IAF (IAF-4, IAF+2) per condition; fallback [8 14] when IAF undefined
        tfr_cache = load('tfr_nback.mat', 'tfr1_bl', 'tfr2_bl', 'tfr3_bl');
        [ERSD_early, ERSD_late, ERSD_full] = compute_ersd_scalars( ...
            {tfr_cache.tfr1_bl, tfr_cache.tfr2_bl, tfr_cache.tfr3_bl}, tfr_cache.tfr1_bl.label, ...
            [IAF1; IAF2; IAF3], alphaRange);

        % ERSD timecourse (occipital, IAF band with [8 14] fallback, dB) per condition for later figures
        ersd_timecourse = compute_ersd_timecourse( ...
            {tfr_cache.tfr1_bl, tfr_cache.tfr2_bl, tfr_cache.tfr3_bl}, tfr_cache.tfr1_bl.label, ...
            [1; 2; 3], [IAF1; IAF2; IAF3], alphaRange);

        % Compute lateralization index on LATE BASELINED spectra (dB)
        powloads = {pow1_bl_late, pow2_bl_late, pow3_bl_late};
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

        % Alpha power in subject IAF band (IAF-4, IAF+2); NaN when IAF undefined
        IAF_band1 = iaf_alpha_band(IAF1, alphaRange);
        IAF_band2 = iaf_alpha_band(IAF2, alphaRange);
        IAF_band3 = iaf_alpha_band(IAF3, alphaRange);

        AlphaPower_raw_early = [robust_roi_pow(pow1_raw_early, channelIdx, IAF_band1); robust_roi_pow(pow2_raw_early, channelIdx, IAF_band2); robust_roi_pow(pow3_raw_early, channelIdx, IAF_band3)];
        AlphaPower_raw_late  = [robust_roi_pow(pow1_raw_late,  channelIdx, IAF_band1); robust_roi_pow(pow2_raw_late,  channelIdx, IAF_band2); robust_roi_pow(pow3_raw_late,  channelIdx, IAF_band3)];
        AlphaPower_raw_full  = [robust_roi_pow(pow1_raw_full,  channelIdx, IAF_band1); robust_roi_pow(pow2_raw_full,  channelIdx, IAF_band2); robust_roi_pow(pow3_raw_full,  channelIdx, IAF_band3)];
        AlphaPower_bl_early  = [robust_roi_pow(pow1_bl_early, channelIdx, IAF_band1); robust_roi_pow(pow2_bl_early, channelIdx, IAF_band2); robust_roi_pow(pow3_bl_early, channelIdx, IAF_band3)];
        AlphaPower_bl_late   = [robust_roi_pow(pow1_bl_late,  channelIdx, IAF_band1); robust_roi_pow(pow2_bl_late,  channelIdx, IAF_band2); robust_roi_pow(pow3_bl_late,  channelIdx, IAF_band3)];
        AlphaPower_bl_full   = [robust_roi_pow(pow1_bl_full,  channelIdx, IAF_band1); robust_roi_pow(pow2_bl_full,  channelIdx, IAF_band2); robust_roi_pow(pow3_bl_full,  channelIdx, IAF_band3)];

        subID = str2double(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'IAF', num2cell([IAF1; IAF2; IAF3]), ...
            'Lateralization', num2cell([LatIdx1; LatIdx2; LatIdx3]), ...
            'ERSD_early', num2cell(ERSD_early), ...
            'ERSD_late', num2cell(ERSD_late), ...
            'ERSD_full', num2cell(ERSD_full), ...
            'AlphaPower_raw_early', num2cell(AlphaPower_raw_early), ...
            'AlphaPower_raw_late', num2cell(AlphaPower_raw_late), ...
            'AlphaPower_raw_full', num2cell(AlphaPower_raw_full), ...
            'AlphaPower_bl_early', num2cell(AlphaPower_bl_early), ...
            'AlphaPower_bl_late', num2cell(AlphaPower_bl_late), ...
            'AlphaPower_bl_full', num2cell(AlphaPower_bl_full));
        % Save
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        if ~isfolder(savepath)
            mkdir(savepath)
        end
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        save lateralization_nback LatIdx1 LatIdx2 LatIdx3
        save('ersd_nback_timecourse.mat', 'ersd_timecourse')
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), ' ...
            '3-back: %f Hz (Power: %f) | Lateralization: %f %f %f \n'], subjects{subj}, IAF1, ...
            powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3, LatIdx1, LatIdx2, LatIdx3);
    catch ME
        fprintf('Continuing to next subject...\n');
    end
end
save(fullfile(paths.features, 'AOC_eeg_matrix_nback.mat'), 'eeg_data_nback')
writetable(struct2table(eeg_data_nback), fullfile(paths.features, 'AOC_eeg_matrix_nback.csv'))

%%
function [IAF, powerIAF] = iaf_from_retention_mtmfft(dataTFR, trialinds, winSec, chLabs, alphaRange)
% Trial-averaged multitaper FFT on retention window; ROI-mean spectrum; legacy peak rules.
% Trials are averaged inside ft_freqanalysis (keeptrials = no), not concatenated.
IAF = NaN;
powerIAF = NaN;
if isempty(trialinds) || isempty(chLabs)
    return
end
try
    cfg_sel = [];
    cfg_sel.latency = winSec;
    cfg_sel.trials = trialinds(:)';
    cfg_sel.channel = chLabs(:);
    dw = ft_selectdata(cfg_sel, dataTFR);
    if isempty(dw.trial)
        return
    end
    cfgf = [];
    cfgf.method = 'mtmfft';
    cfgf.output = 'pow';
    cfgf.taper = 'dpss';
    cfgf.tapsmofrq = 2;
    cfgf.foilim = [6 18];
    cfgf.pad = 'nextpow2';
    cfgf.keeptrials = 'no';
    fr = ft_freqanalysis(cfgf, dw);
    ps = mean(fr.powspctrm, 1);
    [IAF, powerIAF] = iaf_peak_rules(fr.freq(:), ps(:), alphaRange);
catch
    IAF = NaN;
    powerIAF = NaN;
end
end

function [IAF, powerIAF] = iaf_peak_rules(freq, spec, alphaRange)
% Legacy IAF: tallest findpeaks local maximum in [8 14] Hz on ROI-mean spectrum.
% Border peak (first/last alpha bin) → keep IAF, powerIAF = NaN.
freq = freq(:);
spec = spec(:);
IAF = NaN;
powerIAF = NaN;
amask = freq >= alphaRange(1) & freq <= alphaRange(2);
alphaFreqs = freq(amask);
alphaSpec = spec(amask);
if numel(alphaSpec) < 3
    return
end
[pks, locs] = findpeaks(alphaSpec);
if isempty(pks)
    return
end
[~, ind] = max(pks);
IAF = alphaFreqs(locs(ind));
bandIdx = freq > (IAF - 4) & freq < (IAF + 2);
if any(bandIdx)
    powerIAF = mean(spec(bandIdx));
end
if locs(ind) == 1 || locs(ind) == numel(alphaFreqs)
    powerIAF = NaN;
end
end

function band = iaf_alpha_band(IAF, alphaRange)
if isnan(IAF)
    band = [NaN NaN];
else
    band = [IAF - 4, IAF + 2];
end
end

function band = ersd_alpha_band(IAF, alphaRange)
% ERSD band: (IAF-4, IAF+2) when IAF is valid; otherwise fixed alphaRange ([8 14]).
if ~isfinite(IAF)
    band = alphaRange;
else
    band = [IAF - 4, IAF + 2];
end
if any(~isfinite(band)) || band(1) >= band(2)
    band = alphaRange;
end
end

function v = robust_roi_pow(S, channelIdx, band)
fmask = S.freq > band(1) & S.freq < band(2);
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

function [ERSD_early, ERSD_late, ERSD_full] = compute_ersd_scalars(tfPack, labels, iafVals, alphaRange)
ERSD_early = nan(numel(tfPack), 1);
ERSD_late = nan(numel(tfPack), 1);
ERSD_full = nan(numel(tfPack), 1);
occ_ch = occ_channels_from_labels(labels);
latencyWins = {[0 1], [1 2], [0 2]};
for ic = 1:numel(tfPack)
    tf = tfPack{ic};
    chUse = occ_ch(ismember(occ_ch, tf.label));
    if isempty(chUse)
        continue
    end
    band = ersd_alpha_band(iafVals(ic), alphaRange);
    for iw = 1:3
        cfgE = [];
        cfgE.channel = chUse;
        cfgE.avgoverchan = 'yes';
        cfgE.frequency = band;
        cfgE.avgoverfreq = 'yes';
        cfgE.latency = latencyWins{iw};
        cfgE.avgovertime = 'yes';
        try
            outE = ft_selectdata(cfgE, tf);
            v = mean(outE.powspctrm(:), 'omitnan');
        catch
            v = NaN;
        end
        if iw == 1
            ERSD_early(ic) = v;
        elseif iw == 2
            ERSD_late(ic) = v;
        else
            ERSD_full(ic) = v;
        end
    end
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end + 1} = lab;
    end
end
if isempty(ch)
    ch = labels;
end
end

function ersd_tc = compute_ersd_timecourse(tfPack, labels, condVals, iafVals, alphaRange)
occ_ch = occ_channels_from_labels(labels);
nCond = numel(tfPack);
timeVec = tfPack{1}.time(:)';
nTime = numel(timeVec);
tc = nan(nCond, nTime);
for ic = 1:nCond
    tf = tfPack{ic};
    chUse = occ_ch(ismember(occ_ch, tf.label));
    if isempty(chUse)
        continue
    end
    band = ersd_alpha_band(iafVals(ic), alphaRange);
    cfgE = [];
    cfgE.channel = chUse;
    cfgE.avgoverchan = 'yes';
    cfgE.frequency = band;
    cfgE.avgoverfreq = 'yes';
    outE = ft_selectdata(cfgE, tf);
    tc(ic, :) = outE.powspctrm(:)';
end
ersd_tc = struct();
ersd_tc.time = timeVec;
ersd_tc.condition = condVals(:);
ersd_tc.ersd_occ_iaf_db = tc;
end
