%% AOC EEG Feature Extraction — Sternberg
% Computes subject-level EEG features from preprocessed Sternberg EEG.
% Saves `eeg_data_sternberg` (MAT + CSV) and subject-level TFR products.
%
% Extracted features:
%   Power Spectrum (Early [0 1], Late [1 2], Full [0 2]) + Baseline [-1.5 -0.5]
%   Baselined spectra (dB) for each window (mtmconvol, 2 Hz foi grid)
%   IAF (condition-wise): mtmfft+DPSS on retention window, trial-averaged (Sternberg [1 2] s), findpeaks [8 14] Hz
%   IAF (condition-wise)
%   ERSD_early / ERSD_late / ERSD_full (fixed [8 14] Hz on baselined TFR, occipital ROI)

%% POWSPCTRM (Baseline + Early/Late/Full)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

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
        cfg.toi        = -1.5:0.05:3;
        cfg.pad        = 'nextpow2';
        cfg.keeptrials = 'no';

        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR); % chan_freq_time
        cfg.trials = ind4; tfr4 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind6; tfr6 = ft_freqanalysis(cfg, dataTFR);

        cfgb              = [];
        cfgb.baseline     = window.base;
        cfgb.baselinetype = 'db';
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr4_bl = ft_freqbaseline(cfgb, tfr4);
        tfr6_bl = ft_freqbaseline(cfgb, tfr6);

        % Save trial-averaged TFR outputs
        save('tfr_stern.mat', 'tfr2', 'tfr4', 'tfr6', 'tfr2_bl', 'tfr4_bl', 'tfr6_bl')

        % Save ERSD timecourse (occipital, 8 to 14 Hz, dB) per condition for later figures
        ersd_timecourse = compute_ersd_timecourse({tfr2_bl, tfr4_bl, tfr6_bl}, tfr2_bl.label, [2; 4; 6]);
        save('ersd_sternberg_timecourse.mat', 'ersd_timecourse')

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
        fprintf('Continuing to next subject...\n');
    end
end

%% IAF and ERSD
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
% Load data and calculate alpha power, IAF and lateralization index
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'IAF', {}, ...
    'ERSD_early', {}, 'ERSD_late', {}, 'ERSD_full', {});

for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG FEX - STERNBERG] Alpha power, IAF and lateralization for Subject %d / %d \n', subj, length(subjects))
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath);
        load('power_stern_windows.mat');
        load dataEEG_TFR_sternberg

        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind4 = find(dataTFR.trialinfo(:, 1) == 24);
        ind6 = find(dataTFR.trialinfo(:, 1) == 26);

        % Channel selection
        channelIdx = find(ismember(pow2_raw_full.label, channels));
        chLabs = pow2_raw_full.label(channelIdx);

        % IAF: trial-averaged mtmfft on late retention [1 2] s, occipital ROI (~0.5 Hz resolution)
        winIAF = [1 2];
        [IAF2, powerIAF2] = iaf_from_retention_mtmfft(dataTFR, ind2, winIAF, chLabs, alphaRange);
        [IAF4, powerIAF4] = iaf_from_retention_mtmfft(dataTFR, ind4, winIAF, chLabs, alphaRange);
        [IAF6, powerIAF6] = iaf_from_retention_mtmfft(dataTFR, ind6, winIAF, chLabs, alphaRange);

        % ERSD from cached baselined TFR (avoid duplicate spectral transforms)
        tfr_cache = load('tfr_stern.mat', 'tfr2_bl', 'tfr4_bl', 'tfr6_bl');
        [ERSD_early, ERSD_late, ERSD_full] = compute_ersd_scalars( ...
            {tfr_cache.tfr2_bl, tfr_cache.tfr4_bl, tfr_cache.tfr6_bl}, tfr_cache.tfr2_bl.label);

        subID = str2double(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'IAF', num2cell([IAF2; IAF4; IAF6]), ...
            'ERSD_early', num2cell(ERSD_early), ...
            'ERSD_late', num2cell(ERSD_late), ...
            'ERSD_full', num2cell(ERSD_full));
        % Save
        savepath = fullfile(paths.features, subjects{subj}, 'eeg');
        if ~isfolder(savepath)
            mkdir(savepath)
        end
        cd(savepath)
        save eeg_matrix_sternberg_subj subj_data_eeg
        save IAF_sternberg IAF2 IAF4 IAF6
        eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
        clc
        fprintf('Subject %s IAF: WM2 %f Hz, WM4 %f Hz, WM6 %f Hz\n', subjects{subj}, IAF2, IAF4, IAF6);
    catch ME
        fprintf('Continuing to next subject...\n');
    end
end
save(fullfile(paths.features, 'AOC_eeg_matrix_sternberg.mat'), 'eeg_data_sternberg')
writetable(struct2table(eeg_data_sternberg), fullfile(paths.features, 'AOC_eeg_matrix_sternberg.csv'))

function [IAF, powerIAF] = iaf_from_retention_mtmfft(dataTFR, trialinds, winSec, chLabs, alphaRange)
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
