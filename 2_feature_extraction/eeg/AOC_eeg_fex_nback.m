%% AOC EEG Feature Extraction — N-Back
% Computes subject-level EEG features from preprocessed N-back EEG.
% Saves `eeg_data_nback` (MAT + CSV). FOOOF TFR products are produced in
% AOC_eeg_fex_nback_TFR.m.
%
% Extracted features:
%   Power Spectrum windows (raw + baselined) via mtmconvol (2 Hz foi grid; unchanged)
%   IAF (condition-wise): concatenated [0 2]s retention data per condition, mtmfft + DPSS,
%       same peak rules as legacy; not read from the 2 Hz windowed spectra above
%   IAF_specParam: FOOOF alpha peak centre frequency (Hz) per condition (occ channels, median CF)
%   Alpha power in IAF band (early/late/full; raw + dB)
%   Lateralization index (late baselined)

%% POWSPCTRM (Baseline + Early/Late/Full) (subject-level)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

logDir = paths.logs;
scriptName = 'AOC_eeg_fex_nback';

% TEMP: set true to compute IAF_specParam via FOOOF (ft_freqanalysis_Arne_FOOOF).
RUN_FOOOF = false;

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
        cfg.toi        = -1.25:0.05:2.25;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind1; tfr1 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind3; tfr3 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = window.base;
        cfgb.baselinetype = 'db';
        tfr1_bl = ft_freqbaseline(cfgb, tfr1);
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr3_bl = ft_freqbaseline(cfgb, tfr3);

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

% Load data and calculate power, IAF and lateralization index
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'IAF', {}, 'IAF_specParam', {}, 'Lateralization', {}, ...
    'AlphaPower_raw_early', {}, 'AlphaPower_raw_late', {}, 'AlphaPower_raw_full', {}, ...
    'AlphaPower_bl_early', {}, 'AlphaPower_bl_late', {}, 'AlphaPower_bl_full', {});

winIAF = [0 2];  % full retention (IAF / IAF_specParam)

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

        % IAF: concatenated full retention [0 2]s + mtmfft + DPSS (high-res freq), same peak rules
        [IAF1, powerIAF1] = iaf_from_concat_dpss(dataTFR, ind1, winIAF, chLabs, alphaRange);
        [IAF2, powerIAF2] = iaf_from_concat_dpss(dataTFR, ind2, winIAF, chLabs, alphaRange);
        [IAF3, powerIAF3] = iaf_from_concat_dpss(dataTFR, ind3, winIAF, chLabs, alphaRange);

        % FOOOF alpha peak CF (median across occ channels), optional if ft_freqanalysis_Arne_FOOOF on path
        if RUN_FOOOF
            IAF_specParam1 = iaf_specparam_fooof(dataTFR, ind1, winIAF, chLabs, alphaRange);
            IAF_specParam2 = iaf_specparam_fooof(dataTFR, ind2, winIAF, chLabs, alphaRange);
            IAF_specParam3 = iaf_specparam_fooof(dataTFR, ind3, winIAF, chLabs, alphaRange);
        else
            IAF_specParam1 = NaN;
            IAF_specParam2 = NaN;
            IAF_specParam3 = NaN;
        end

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

        % Alpha power in IAF band (raw + baselined) for early/late/full
        IAF_band1 = [IAF1-4 IAF1+2]; if any(isnan(IAF_band1)); IAF_band1 = alphaRange; end
        IAF_band2 = [IAF2-4 IAF2+2]; if any(isnan(IAF_band2)); IAF_band2 = alphaRange; end
        IAF_band3 = [IAF3-4 IAF3+2]; if any(isnan(IAF_band3)); IAF_band3 = alphaRange; end

        AlphaPower_raw_early = [robust_roi_pow(pow1_raw_early, channelIdx, IAF_band1); robust_roi_pow(pow2_raw_early, channelIdx, IAF_band2); robust_roi_pow(pow3_raw_early, channelIdx, IAF_band3)];
        AlphaPower_raw_late  = [robust_roi_pow(pow1_raw_late,  channelIdx, IAF_band1); robust_roi_pow(pow2_raw_late,  channelIdx, IAF_band2); robust_roi_pow(pow3_raw_late,  channelIdx, IAF_band3)];
        AlphaPower_raw_full  = [robust_roi_pow(pow1_raw_full,  channelIdx, IAF_band1); robust_roi_pow(pow2_raw_full,  channelIdx, IAF_band2); robust_roi_pow(pow3_raw_full,  channelIdx, IAF_band3)];
        AlphaPower_bl_early  = [robust_roi_pow(pow1_bl_early, channelIdx, IAF_band1); robust_roi_pow(pow2_bl_early, channelIdx, IAF_band2); robust_roi_pow(pow3_bl_early, channelIdx, IAF_band3)];
        AlphaPower_bl_late   = [robust_roi_pow(pow1_bl_late,  channelIdx, IAF_band1); robust_roi_pow(pow2_bl_late,  channelIdx, IAF_band2); robust_roi_pow(pow3_bl_late,  channelIdx, IAF_band3)];
        AlphaPower_bl_full   = [robust_roi_pow(pow1_bl_full,  channelIdx, IAF_band1); robust_roi_pow(pow2_bl_full,  channelIdx, IAF_band2); robust_roi_pow(pow3_bl_full,  channelIdx, IAF_band3)];

        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'IAF', num2cell([IAF1; IAF2; IAF3]), ...
            'IAF_specParam', num2cell([IAF_specParam1; IAF_specParam2; IAF_specParam3]), ...
            'Lateralization', num2cell([LatIdx1; LatIdx2; LatIdx3]), ...
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
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        save lateralization_nback LatIdx1 LatIdx2 LatIdx3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), ' ...
            '3-back: %f Hz (Power: %f) | Lateralization: %f %f %f \n'], subjects{subj}, IAF1, ...
            powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3, LatIdx1, LatIdx2, LatIdx3);
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end
save(fullfile(paths.features, 'AOC_eeg_matrix_nback.mat'), 'eeg_data_nback')
writetable(struct2table(eeg_data_nback), fullfile(paths.features, 'AOC_eeg_matrix_nback.csv'))

function datc = local_ft_concat_trials(dw)
% One FieldTrip trial: all epochs of dw concatenated along time (per channel).
nTr = numel(dw.trial);
nCh = size(dw.trial{1}, 1);
Ttot = 0;
for k = 1:nTr
    Ttot = Ttot + size(dw.trial{k}, 2);
end
datc = struct();
datc.label = dw.label;
datc.fsample = dw.fsample;
datc.trial = {zeros(nCh, Ttot)};
t0 = 1;
for k = 1:nTr
    tk = size(dw.trial{k}, 2);
    datc.trial{1}(:, t0:t0 + tk - 1) = dw.trial{k};
    t0 = t0 + tk;
end
datc.time = {(0:Ttot - 1) ./ dw.fsample};
datc.sampleinfo = [1 Ttot];
if isfield(dw, 'hdr') && ~isempty(dw.hdr)
    datc.hdr = dw.hdr;
end
end

function [IAF, powerIAF] = iaf_from_concat_dpss(dataTFR, trialinds, winSec, chLabs, alphaRange)
% Concatenate retention window per trial, multitaper FFT, ROI-mean spectrum, legacy peak rules.
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
    datc = local_ft_concat_trials(dw);
    cfgf = [];
    cfgf.method = 'mtmfft';
    cfgf.output = 'pow';
    cfgf.taper = 'dpss';
    cfgf.tapsmofrq = 2;
    cfgf.foilim = [6 18];
    cfgf.pad = 5;
    fr = ft_freqanalysis(cfgf, datc);
    ps = mean(fr.powspctrm, 1);
    [IAF, powerIAF] = iaf_peak_rules(fr.freq(:), ps(:), alphaRange);
catch
    IAF = NaN;
    powerIAF = NaN;
end
end

function [IAF, powerIAF] = iaf_peak_rules(freq, spec, alphaRange)
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
df = median(diff(alphaFreqs));
if ~isfinite(df) || df <= 0
    df = min(diff(alphaFreqs));
    if ~isfinite(df) || df <= 0
        return
    end
end
if IAF <= alphaRange(1) + df || IAF >= alphaRange(2) - df
    IAF = NaN;
    return
end
bandIdx = freq > (IAF - 4) & freq < (IAF + 2);
if ~any(bandIdx)
    return
end
powerIAF = mean(spec(bandIdx));
if powerIAF > max(pks)
    IAF = NaN;
    powerIAF = NaN;
end
end

function cf = iaf_specparam_fooof(dataTFR, trialinds, winSec, chLabs, alphaRange)
% Median FOOOF peak centre frequency in alphaRange across selected channels (trial-averaged spectrum).
cf = NaN;
if isempty(trialinds) || isempty(chLabs) || exist('ft_freqanalysis_Arne_FOOOF', 'file') ~= 2
    return
end
try
    cfg_fooof = [];
    cfg_fooof.method = 'mtmfft';
    cfg_fooof.taper = 'hanning';
    cfg_fooof.foilim = [2 40];
    cfg_fooof.pad = 5;
    cfg_fooof.output = 'fooof';
    cfg_fooof.keeptrials = 'no';
    cfg_sel = [];
    cfg_sel.latency = winSec;
    cfg_sel.trials = trialinds(:)';
    cfg_sel.channel = chLabs(:);
    dat = ft_selectdata(cfg_sel, dataTFR);
    fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat);
    if iscell(fooof_out.fooofparams)
        rep = fooof_out.fooofparams{1};
    else
        rep = fooof_out.fooofparams;
    end
    cfs = [];
    for k = 1:numel(rep)
        if ~isfield(rep(k), 'peak_params') || isempty(rep(k).peak_params)
            continue
        end
        pk = rep(k).peak_params;
        inb = pk(:, 1) >= alphaRange(1) & pk(:, 1) <= alphaRange(2);
        if ~any(inb)
            continue
        end
        pk = pk(inb, :);
        [~, im] = max(pk(:, 2));
        cfs(end + 1) = pk(im, 1); %#ok<AGROW>
    end
    if ~isempty(cfs)
        cf = median(cfs, 'omitnan');
    end
catch
    cf = NaN;
end
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
