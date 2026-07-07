%% AOC EEG Feature Extraction — N-Back
% Computes subject-level EEG features from preprocessed N-back EEG.
% Saves `eeg_data_nback` (MAT + CSV) and subject-level TFR products.
%
% Extracted features:
%   Power Spectrum windows (raw + baselined) via mtmconvol (2 Hz foi grid)
%   IAF (condition-wise): mtmfft+DPSS on retention window, trial-averaged (N-back [0 2] s), findpeaks [8 14] Hz
%   IAF (condition-wise)
%   ERSD_early / ERSD_late / ERSD_full (fixed [8 14] Hz on baselined TFR, occipital ROI)

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

        % Save ERSD timecourse (occipital, 8 to 14 Hz, dB) per condition for later figures
        ersd_timecourse = compute_ersd_timecourse({tfr1_bl, tfr2_bl, tfr3_bl}, tfr1_bl.label, [1; 2; 3]);
        save('ersd_nback_timecourse.mat', 'ersd_timecourse')

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

%% IAF and ERSD
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

% Load data and calculate power, IAF and lateralization index
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'IAF', {}, ...
    'ERSD_early', {}, 'ERSD_late', {}, 'ERSD_full', {});

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
        tfr_cache = load('tfr_nback.mat', 'tfr1_bl', 'tfr2_bl', 'tfr3_bl');
        [ERSD_early, ERSD_late, ERSD_full] = compute_ersd_scalars( ...
            {tfr_cache.tfr1_bl, tfr_cache.tfr2_bl, tfr_cache.tfr3_bl}, tfr_cache.tfr1_bl.label);

        subID = str2double(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'IAF', num2cell([IAF1; IAF2; IAF3]), ...
            'ERSD_early', num2cell(ERSD_early), ...
            'ERSD_late', num2cell(ERSD_late), ...
            'ERSD_full', num2cell(ERSD_full));
        % Save
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        if ~isfolder(savepath)
            mkdir(savepath)
        end
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save IAF_nback IAF1 IAF2 IAF3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf('Subject %s IAF: 1-back %f Hz, 2-back %f Hz, 3-back %f Hz\n', subjects{subj}, IAF1, IAF2, IAF3);
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

function [ERSD_early, ERSD_late, ERSD_full] = compute_ersd_scalars(tfPack, labels)
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
    for iw = 1:3
        cfgE = [];
        cfgE.channel = chUse;
        cfgE.avgoverchan = 'yes';
        cfgE.frequency = [8 14];
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

function ersd_tc = compute_ersd_timecourse(tfPack, labels, condVals)
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
    cfgE = [];
    cfgE.channel = chUse;
    cfgE.avgoverchan = 'yes';
    cfgE.frequency = [8 14];
    cfgE.avgoverfreq = 'yes';
    outE = ft_selectdata(cfgE, tf);
    tc(ic, :) = outE.powspctrm(:)';
end
ersd_tc = struct();
ersd_tc.time = timeVec;
ersd_tc.condition = condVals(:);
ersd_tc.ersd_occ_8_14_db = tc;
end
