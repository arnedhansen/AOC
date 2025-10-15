%% AOC EEG Feature Extraction Sternberg TRIAL-BY-TRIAL
%
% Extracted features:
%   Power Spectrum (Early, Late (= Registered Retention))  [trial-by-trial]
%   IAF (subject-level), Power at IAF (trial-wise), and Lateralization Index (trial-wise)
%   TFR (Raw, FOOOF and Baselined)

%% POWSPCTRM (Early, Late (= Registered Retention), Baseline Period) - TRIAL-BY-TRIAL
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing POWSPCTRM (TRIALS) for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg

        % Identify indices of trials belonging to conditions
        globalTrialID = dataTFR.trialinfo(:,2);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo(:, 1) == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo(:, 1) == 26); % WM load 6
        globalTrialID2 = globalTrialID(ind2);
        globalTrialID4 = globalTrialID(ind4);
        globalTrialID6 = globalTrialID(ind6);

        % ----------------------
        % Frequency analysis EARLY  (0-1 s) RAW
        % ----------------------
        % Select data
        cfg = [];                      % Empty configuration
        cfg.latency = [0 1];           % Segmentation for early retention interval
        dataEarly = ft_selectdata(cfg, dataTFR);

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'yes';        % Keep trial information
        cfg.pad = 5;                   % Zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind4;
        powload4_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind6;
        powload6_early = ft_freqanalysis(cfg, dataEarly);

       

        % ----------------------
        % Frequency analysis EARLY (0-1 s) BASELINED
        % ----------------------
        cfg = [];
        cfg.baseline = [-.5 -.25];
        cfg.baselinetype = 'db';
        powload2_early_bl = ft_freqbaseline(cfg, powload2_early);
        powload4_early_bl = ft_freqbaseline(cfg, powload4_early);
        powload6_early_bl = ft_freqbaseline(cfg, powload6_early);

        % Save early trial-wise power spectra
        cd(datapath)
        save power_stern_early_trials powload2_early powload4_early powload6_early ...
            powload2_early_bl powload4_early_bl powload6_early_bl

        % ----------------------
        % Frequency analysis LATE (1-2 s) RAW = Registered Retention
        % ----------------------
        % Select data
        cfg = [];
        cfg.latency = [1 2];               % Segmentation for registered retention interval
        dataLate = ft_selectdata(cfg, dataTFR); % Select data

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'yes';        % Keep trial information
        cfg.pad = 5;                   % Zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2_late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind4;
        powload4_late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind6;
        powload6_late = ft_freqanalysis(cfg, dataLate);

        % ----------------------
        % Frequency analysis LATE (1-2 s) BASELINED
        % ----------------------
        cfg = [];
        cfg.baseline = [-.5 -.25];
        cfg.baselinetype = 'db';
        powload2_late_bl  = ft_freqbaseline(cfg, powload2_late);
        powload4_late_bl  = ft_freqbaseline(cfg, powload4_late);
        powload6_late_bl  = ft_freqbaseline(cfg, powload6_late);

        % Save late trial-wise spectra
        cd(datapath)
        save power_stern_late_trials powload2_late powload4_late powload6_late ...
            powload2_late_bl powload4_late_bl powload6_late_bl

    catch ME
        ME.message
        error(['ERROR extracting trial-wise power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER (trial-wise), IAF (subject-level) and LATERALIZATION INDEX (trial-wise)
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define channels (from one subject-s labels)
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern_late_trials.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2_late.label)
    label = powload2_late.label{i};
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
        disp(['Midline or nonstandard channel: ', ch])
    end
end

% Load data and calculate alpha power (trial-wise), IAF (subject-level) and lateralization index (trial-wise)
alphaRange = [8 14];
epsP = 1e-12; % small epsilon for safe logs / denominators
ridgeFrac = 0.01; % ridge as fraction of typical (R+L)

eeg_data_sternberg_trials = struct('Trial', {}, 'ID', {}, 'Condition', {}, ...
    'AlphaPowerEarly', {}, 'AlphaPowerEarlyBL', {}, 'AlphaPowerLate', {}, ...
    'AlphaPowerLateBL', {}, 'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power (trial-wise), IAF (subject), Lateralization (trial-wise) for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);

        % Load trial-wise spectra
        load('power_stern_early_trials.mat');
        load('power_stern_late_trials.mat');
        load('power_stern_baseline_period_trials.mat');

        % Channel selection
        channelIdx = find(ismember(powload2_early.label, channels));

        % ----------------------
        % Subject-level IAF (from late retention, trial-averaged ROI)
        % ----------------------
        % Build subject-level ROI-averaged spectra (average across trials, then across ROI channels)
        % We average trials first to stabilise the IAF estimate.
        % Ensure dims: rpt x chan x freq
        S2 = squeeze(mean(powload2_late.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S4 = squeeze(mean(powload4_late.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S6 = squeeze(mean(powload6_late.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        subjSpec = nanmean([S2; S4; S6], 1);                                      % pooled across trials & loads
        freqs = powload2_late.freq(:)';

        % Find IAF with guards (smooth slightly, edges excluded)
        alphaMask  = (freqs >= alphaRange(1)) & (freqs <= alphaRange(2));
        alphaFreqs = freqs(alphaMask);
        alphaSpec  = subjSpec(alphaMask);
        if numel(alphaSpec) >= 5
            alphaSpec = movmean(alphaSpec, 3);
        end
        [pks, locs, w, p] = findpeaks(alphaSpec, alphaFreqs, 'MinPeakProminence', max(eps, 0.02*max(alphaSpec)));
        if isempty(pks)
            IAF_subj = NaN;
        else
            [~, idxMax] = max([p(:), pks(:)] * [1; 1e-3]);
            IAF_subj = locs(idxMax);
            df = median(diff(alphaFreqs));
            if (IAF_subj - alphaRange(1) <= df) || (alphaRange(2) - IAF_subj <= df)
                IAF_subj = NaN;
            end
        end

        % Define IAF integration band if valid (4/+2 Hz)
        if ~isnan(IAF_subj)
            IAF_band = [(IAF_subj-4) (IAF_subj+2)];
        else
            IAF_band = [NaN NaN];
        end        

        % ----------------------
        % Trial-wise Alpha Power (EARLY/LATE, RAW/BASELINED) at IAF band
        % ----------------------
        % EARLY RAW
        AlphaPowerEarly2   = bandpower_trials(powload2_early,  channelIdx, powload2_early.freq,  IAF_band);
        AlphaPowerEarly4   = bandpower_trials(powload4_early,  channelIdx, powload4_early.freq,  IAF_band);
        AlphaPowerEarly6   = bandpower_trials(powload6_early,  channelIdx, powload6_early.freq,  IAF_band);
        % EARLY BL (dB)
        AlphaPowerEarlyBL2 = bandpower_trials(powload2_early_bl, channelIdx, powload2_early_bl.freq, IAF_band);
        AlphaPowerEarlyBL4 = bandpower_trials(powload4_early_bl, channelIdx, powload4_early_bl.freq, IAF_band);
        AlphaPowerEarlyBL6 = bandpower_trials(powload6_early_bl, channelIdx, powload6_early_bl.freq, IAF_band);
        % LATE RAW
        AlphaPowerLate2    = bandpower_trials(powload2_late,   channelIdx,  powload2_late.freq,   IAF_band);
        AlphaPowerLate4    = bandpower_trials(powload4_late,   channelIdx,  powload4_late.freq,   IAF_band);
        AlphaPowerLate6    = bandpower_trials(powload6_late,   channelIdx,  powload6_late.freq,   IAF_band);
        % LATE BL (dB)
        AlphaPowerLateBL2  = bandpower_trials(powload2_late_bl, channelIdx,  powload2_late_bl.freq, IAF_band);
        AlphaPowerLateBL4  = bandpower_trials(powload4_late_bl, channelIdx,  powload4_late_bl.freq, IAF_band);
        AlphaPowerLateBL6  = bandpower_trials(powload6_late_bl, channelIdx,  powload6_late_bl.freq, IAF_band);

        % ----------------------
        % Trial-wise Lateralization Index (use LATE BL by default)
        % ----------------------
        [LI2_trials, ~] = lateralization_trials(powload2_late_bl, left_channels, right_channels, powload2_late_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI4_trials, ~] = lateralization_trials(powload4_late_bl, left_channels, right_channels, powload4_late_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI6_trials, ~] = lateralization_trials(powload6_late_bl, left_channels, right_channels, powload6_late_bl.freq, IAF_band, ridgeFrac, epsP);

        % ----------------------
        % Build subject trial-wise structure array (now with 4 alpha-power fields)
        % ----------------------
        subID = str2double(subjects{subj});
        n2 = size(powload2_late.powspctrm,1);
        n4 = size(powload4_late.powspctrm,1);
        n6 = size(powload6_late.powspctrm,1);

        IAFr2 = repmat(IAF_subj, n2, 1);
        IAFr4 = repmat(IAF_subj, n4, 1);
        IAFr6 = repmat(IAF_subj, n6, 1);

        subj_data_eeg_trials_2 = struct( ...
            'Trial',              num2cell(globalTrialID2), ...
            'ID',                 num2cell(repmat(subID, n2, 1)), ...
            'Condition',          num2cell(repmat(2, n2, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly2), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL2), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate2), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL2), ...
            'IAF',                num2cell(IAFr2), ...
            'Lateralization',     num2cell(LI2_trials) );

        subj_data_eeg_trials_4 = struct( ...
            'Trial',              num2cell(globalTrialID4), ...
            'ID',                 num2cell(repmat(subID, n4, 1)), ...
            'Condition',          num2cell(repmat(4, n4, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly4), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL4), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate4), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL4), ...
            'IAF',                num2cell(IAFr4), ...
            'Lateralization',     num2cell(LI4_trials) );

        subj_data_eeg_trials_6 = struct( ...
            'Trial',              num2cell(globalTrialID6), ...
            'ID',                 num2cell(repmat(subID, n6, 1)), ...
            'Condition',          num2cell(repmat(6, n6, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly6), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL6), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate6), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL6), ...
            'IAF',                num2cell(IAFr6), ...
            'Lateralization',     num2cell(LI6_trials) );

        subj_data_eeg_trials = [subj_data_eeg_trials_2; subj_data_eeg_trials_4; subj_data_eeg_trials_6];

        % Save (per subject + append to grand table)
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_sternberg_subj_trials subj_data_eeg_trials
        save IAF_sternberg_subject IAF_subj

        % Append to grand struct
        eeg_data_sternberg_trials = [eeg_data_sternberg_trials; subj_data_eeg_trials];

        % Console output
        clc
        fprintf(['Subject %s done...'], subjects{subj});

    catch ME
        ME.message
        error(['ERROR calculating trial-wise alpha power / LI for Subject ' num2str(subjects{subj}) '!'])
    end
end

% Save pooled table
if ispc == 1
    save W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg_trials eeg_data_sternberg_trials
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials eeg_data_sternberg_trials
end

%% FOOOF on FFT (trial-wise; EARLY and LATE) - optional add-on
% % Setup
% startup
% [subjects, path, ~ , ~] = setup('AOC');
%
% % FOOOF settings
% foo_settings = struct();
% foo_settings.peak_width_limits = [2 12];
% foo_settings.aperiodic_mode    = 'fixed';
% foo_settings.verbose           = false;
%
% foo_freq_range  = [3 30];     % fit range (Hz)
% foo_alpha_range = [8 14];     % peak selection range (can switch to [IAF_subj-1 IAF_subj+1] per subject)
%
% for subj = 1:length(subjects)
%     clc
%     disp(['FOOOF (trial-wise FFT) for Subject AOC ', num2str(subjects{subj})])
%
%     try
%         % Load per-condition FFT outputs
%         datapath = strcat(path, subjects{subj}, filesep, 'eeg');
%         cd(datapath)
%         load power_stern_early_trials
%         load power_stern_late_trials
%
%         % ROI (reuse 'channels' if it exists; else derive from labels)
%         if ~exist('channels','var') || isempty(channels)
%             occ_channels = {};
%             for i = 1:length(powload2_late.label)
%                 lbl = powload2_late.label{i};
%                 if contains(lbl, {'O'}) || contains(lbl, {'I'})
%                     occ_channels{end+1} = lbl;
%                 end
%             end
%             channels = occ_channels;
%         end
%
%         % Indices per object
%         ch_e2 = find(ismember(powload2_early.label, channels));
%         ch_e4 = find(ismember(powload4_early.label, channels));
%         ch_e6 = find(ismember(powload6_early.label, channels));
%         ch_l2 = find(ismember(powload2_late.label,  channels));
%         ch_l4 = find(ismember(powload4_late.label,  channels));
%         ch_l6 = find(ismember(powload6_late.label,  channels));
%
%         % Run FOOOF per trial on ROI-averaged FFT spectra
%         foo_early2 = fooof_trials_fft(powload2_early, ch_e2, foo_freq_range, foo_settings, foo_alpha_range);
%         foo_early4 = fooof_trials_fft(powload4_early, ch_e4, foo_freq_range, foo_settings, foo_alpha_range);
%         foo_early6 = fooof_trials_fft(powload6_early, ch_e6, foo_freq_range, foo_settings, foo_alpha_range);
%
%         foo_late2  = fooof_trials_fft(powload2_late,  ch_l2, foo_freq_range, foo_settings, foo_alpha_range);
%         foo_late4  = fooof_trials_fft(powload4_late,  ch_l4, foo_freq_range, foo_settings, foo_alpha_range);
%         foo_late6  = fooof_trials_fft(powload6_late,  ch_l6, foo_freq_range, foo_settings, foo_alpha_range);
%
%         % Save to disk (no change to your trial table)
%         save fooof_fft_trials ...
%             foo_early2 foo_early4 foo_early6 ...
%             foo_late2  foo_late4  foo_late6
%
%     catch ME
%         ME.message
%         error(['ERROR in trial-wise FFT FOOOF for Subject ' num2str(subjects{subj}) '!'])
%     end
% end

%% TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM
% % Setup
% startup
% [subjects, path, ~ , ~] = setup('AOC');
% runMode = askRunMode();
%
% % Read data, segment and convert to FieldTrip data structure
% for subj = 1:length(subjects)
%
%     % Only process new data
%     datapath = strcat(path, subjects{subj}, filesep, 'eeg');
%     newDataFolder = dir([datapath, filesep, 'power_stern_fooof_trials.mat']);
%
%     if strcmp(runMode,'all') || isempty(newDataFolder)
%         clc
%         disp(['Processing TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM for Subject AOC ', num2str(subjects{subj})])
%         try
%             cd(datapath)
%             close all
%             load dataEEG_TFR_sternberg
%
%             % Identify indices of trials belonging to conditions
%             globalTrialID = dataTFR.trialinfo(:,2);
%             ind2 = find(dataTFR.trialinfo(:, 1) == 22);
%             ind4 = find(dataTFR.trialinfo(:, 1) == 24);
%             ind6 = find(dataTFR.trialinfo(:, 1) == 26);
%
%             % ----------------------
%             % Time frequency analysis (averaged across trials for FOOOF stability)
%             % ----------------------
%             cfg              = [];
%             cfg.output       = 'pow';
%             cfg.method       = 'mtmconvol';
%             cfg.taper        = 'hanning';
%             cfg.foi          = 3:1:30;                         % 3 to 30 Hz in steps of 1 Hz
%             cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % 0.5 s windows
%             cfg.toi          = -1.5:0.1:3;                     % throttle to 100 ms steps to reduce FOOOF load
%             cfg.keeptrials   = 'no';                            % average for FOOOF robustness
%
%             cfg.trials = ind2;
%             tfr2 = ft_freqanalysis(cfg, dataTFR);
%             cfg.trials = ind4;
%             tfr4 = ft_freqanalysis(cfg, dataTFR);
%             cfg.trials = ind6;
%             tfr6 = ft_freqanalysis(cfg, dataTFR);
%
%             % ----------------------
%             % FOOOF (trial-wise)
%             % ----------------------
%             orig_freq = 3:1:30;
%             tfrs = {tfr2, tfr4, tfr6};
%             for tfr_conds = 1:3
%                 clc; disp('FOOOFing...')
%                 clear fspctrm
%                 tfr = tfrs{1, tfr_conds};
%
%                 % Pre-allocate
%                 nch = numel(tfr.label); nfr = numel(tfr.freq); nt = numel(tfr.time);
%                 fspctrm = nan(nch, nfr, nt, 'like', tfr.powspctrm);
%                 powspctrmff = nan(nch, nfr, 'like', tfr.powspctrm);
%
%                 % FOOOF settings
%                 settings = struct();
%                 settings.peak_width_limits = [2 12];
%                 settings.aperiodic_mode = 'fixed';
%                 settings.verbose = false;
%
%                 % Downsample time frames further for speed (every ~200 ms)
%                 timeIdx = 1:2:nt;
%
%                 for tt = timeIdx
%                     % Output progress
%                     clc
%                     disp(['subj      ' num2str(subj)])
%                     disp(['cond      ' num2str(tfr_conds)])
%                     disp(['timepnt   ' num2str(tt) ' / ' num2str(nt)])
%
%                     % Config
%                     cfgSel = [];
%                     cfgSel.latency = tfr.time(tt);
%                     % cfgSel.latency = [tfr.time(tt) tfr.time(tt)]
%                     tmp = ft_selectdata(cfgSel, tfr);             % dims: chan x freq  (since keeptrials='no')
%
%                     for chan = 1:length(tmp.label)
%                         % Prepare inputs
%                         freqs = orig_freq(:);                         % Equidistant freq vector (column)
%                         psd   = squeeze(tmp.powspctrm(chan, :))';     % row  flip to column below
%                         psd   = psd(:);
%
%                         % Keep only finite & positive bins; bail if too few
%                         good = isfinite(psd) & (psd > 0) & isfinite(freqs);
%                         if nnz(good) < 5
%                             powspctrmff(chan, :) = NaN;
%                             continue
%                         end
%                         freqs_use = freqs(good);
%                         psd_use   = psd(good);
%
%                         % Fit
%                         fooof_results = fooof(freqs_use, psd_use, [min(freqs_use), max(freqs_use)], settings, true);
%                         if isfield(fooof_results, 'fooofed_spectrum') && isfield(fooof_results, 'ap_fit')
%                             % Reconstruct on the original grid; simple nearest if needed
%                             % (Assumes equidistant freqs; here good==all typically)
%                             ff = nan(nfr,1);
%                             af = nan(nfr,1);
%                             ff(good) = fooof_results.fooofed_spectrum(:);
%                             af(good) = fooof_results.ap_fit(:);
%                             powspctrmff(chan,:) = (ff - af)'; % row
%                         else
%                             powspctrmff(chan,:) = NaN;
%                         end
%                     end
%                     fspctrm(:,:,tt) = powspctrmff;
%                 end
%
%                 % Hold last valid frame for skipped indices to keep size consistent
%                 if numel(timeIdx) < nt
%                     lastDone = timeIdx(end);
%                     for tt = setdiff(1:nt, timeIdx)
%                         fspctrm(:,:,tt) = fspctrm(:,:,lastDone);
%                     end
%                 end
%
%                 % Assign out
%                 if tfr_conds == 1
%                     tfr2_fooof = tfr;  tfr2_fooof.powspctrm = fspctrm;
%                 elseif tfr_conds == 2
%                     tfr4_fooof = tfr;  tfr4_fooof.powspctrm = fspctrm;
%                 elseif tfr_conds == 3
%                     tfr6_fooof = tfr;  tfr6_fooof.powspctrm = fspctrm;
%                 end
%             end
%             disp(upper('FOOOF done...'))
%
%             % ----------------------
%             % Baselined TFR
%             % ----------------------
%             % Raw powspctrm baselined (dB)
%             cfg                              = [];
%             cfg.baseline                     = [-.5 -.25];
%             cfg.baselinetype                 = 'db';
%             tfr2_bl                          = ft_freqbaseline(cfg, tfr2);
%             tfr4_bl                          = ft_freqbaseline(cfg, tfr4);
%             tfr6_bl                          = ft_freqbaseline(cfg, tfr6);
%
%             % FOOOFed powspctrm baselined (absolute)
%             cfg                              = [];
%             cfg.baseline                     = [-.5 -.25];
%             cfg.baselinetype                 = 'absolute';   % FOOOF already in log space; no dB here
%             tfr2_fooof_bl                    = ft_freqbaseline(cfg, tfr2_fooof);
%             tfr4_fooof_bl                    = ft_freqbaseline(cfg, tfr4_fooof);
%             tfr6_fooof_bl                    = ft_freqbaseline(cfg, tfr6_fooof);
%
%             % Save data
%             cd(datapath)
%             save tfr_stern ...
%                 tfr2 tfr4 tfr6 ...
%                 tfr2_fooof tfr4_fooof tfr6_fooof ...
%                 tfr2_bl tfr4_bl tfr6_bl ...
%                 tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl
%
%             % Convert TFR data to POWSPCTRM (channels x frequency) - using your helpers
%             analysisPeriodFull  = [0 2];
%             analysisPeriodEarly = [0 1];
%             analysisPeriodLate  = [1 2];
%             freq_range          = [3 30];
%
%             % Select data
%             pow2_fooof                  = select_data(analysisPeriodFull,  freq_range, tfr2_fooof);
%             pow2_fooof_bl               = select_data(analysisPeriodFull,  freq_range, tfr2_fooof_bl);
%             pow2_fooof_bl_early         = select_data(analysisPeriodEarly, freq_range, tfr2_fooof_bl);
%             pow2_fooof_bl_late          = select_data(analysisPeriodLate,  freq_range, tfr2_fooof_bl);
%
%             pow4_fooof                  = select_data(analysisPeriodFull,  freq_range, tfr4_fooof);
%             pow4_fooof_bl               = select_data(analysisPeriodFull,  freq_range, tfr4_fooof_bl);
%             pow4_fooof_bl_early         = select_data(analysisPeriodEarly, freq_range, tfr4_fooof_bl);
%             pow4_fooof_bl_late          = select_data(analysisPeriodLate,  freq_range, tfr4_fooof_bl);
%
%             pow6_fooof                  = select_data(analysisPeriodFull,  freq_range, tfr6_fooof);
%             pow6_fooof_bl               = select_data(analysisPeriodFull,  freq_range, tfr6_fooof_bl);
%             pow6_fooof_bl_early         = select_data(analysisPeriodEarly, freq_range, tfr6_fooof_bl);
%             pow6_fooof_bl_late          = select_data(analysisPeriodLate,  freq_range, tfr6_fooof_bl);
%
%             % Remove time dimension for POWSPCTRM (channels x frequency)
%             pow2_fooof                  = remove_time_dimension(pow2_fooof);
%             pow2_fooof_bl               = remove_time_dimension(pow2_fooof_bl);
%             pow2_fooof_bl_early         = remove_time_dimension(pow2_fooof_bl_early);
%             pow2_fooof_bl_late          = remove_time_dimension(pow2_fooof_bl_late);
%
%             pow4_fooof                  = remove_time_dimension(pow4_fooof);
%             pow4_fooof_bl               = remove_time_dimension(pow4_fooof_bl);
%             pow4_fooof_bl_early         = remove_time_dimension(pow4_fooof_bl_early);
%             pow4_fooof_bl_late          = remove_time_dimension(pow4_fooof_bl_late);
%
%             pow6_fooof                  = remove_time_dimension(pow6_fooof);
%             pow6_fooof_bl               = remove_time_dimension(pow6_fooof_bl);
%             pow6_fooof_bl_early         = remove_time_dimension(pow6_fooof_bl_early);
%             pow6_fooof_bl_late          = remove_time_dimension(pow6_fooof_bl_late);
%
%             save power_stern_fooof_trials ...
%                 pow2_fooof pow4_fooof pow6_fooof ...
%                 pow2_fooof_bl pow4_fooof_bl pow6_fooof_bl ...
%                 pow2_fooof_bl_early pow4_fooof_bl_early pow6_fooof_bl_early ...
%                 pow2_fooof_bl_late pow4_fooof_bl_late pow6_fooof_bl_late
%             clc
%         catch ME
%             ME.message
%             error(['ERROR extracting TFR for Subject ' num2str(subjects{subj}) '!'])
%         end
%     else
%         disp(['TFR and FOOOFed POWSPCTRM already exists for Subject AOC ', num2str(subjects{subj})])
%     end
% end
% disp('TFR and FOOOFed POWSPCTRM COMPUTED...');
