%% AOC EEG Feature Extraction N-back TRIAL-BY-TRIAL
%
% Extracted features:
%   Power Spectrum (Early, Late (= Registered Retention), Baseline Period)  [trial-by-trial]
%   IAF (subject-level), Power at IAF (trial-wise), and Lateralization Index (trial-wise, ridge stabilised)
%   TFR (Raw, FOOOF and Baselined)  [with safer, throttled FOOOF]

%% POWSPCTRM (Early, Late, Registered, Baseline Period) — TRIAL-BY-TRIAL
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
        load dataEEG_TFR_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(dataTFR.trialinfo == 21); % WM load 2
        ind2 = find(dataTFR.trialinfo == 22); % WM load 4
        ind3 = find(dataTFR.trialinfo == 23); % WM load 6

        % ----------------------
        % Frequency analysis EARLY  (0–1 s)
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
        cfg.trials = ind1;
        powload1_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind2;
        powload1_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind3;
        powload3_early = ft_freqanalysis(cfg, dataEarly);

        % Save early trial-wise power spectra
        cd(datapath)
        save power_nback_early_trials powload1_early powload1_early powload3_early

        % ----------------------
        % Frequency analysis LATE  (1–2 s)  = Registered Retention
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
        cfg.trials = ind1;
        powload1_late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind2;
        powload2_late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind3;
        powload3_late = ft_freqanalysis(cfg, dataLate);

        % Save late trial-wise spectra
        cd(datapath)
        save power_nback_late_trials powload1_late powload2_late powload3_late

        % ----------------------
        % Frequency analysis FULL  (0–2 s)
        % ----------------------
        % Select data
        cfg = [];                      % Empty configuration
        cfg.latency = [0 2];           % Segmentation for early retention interval
        dataFull = ft_selectdata(cfg, dataTFR);

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
        cfg.trials = ind1;
        powload1_full = ft_freqanalysis(cfg, dataFull);
        cfg.trials = ind2;
        powload2_full = ft_freqanalysis(cfg, dataFull);
        cfg.trials = ind3;
        powload3_full = ft_freqanalysis(cfg, dataFull);

        % Save early trial-wise power spectra
        cd(datapath)
        save power_nback_full_trials powload1_full powload2_full powload3_full

        % ----------------------
        % Frequency analysis BASELINE PERIOD  (-0.5–.25 s)
        % ----------------------
        % Select data
        cfg = [];                      % Empty configuration
        cfg.latency = [-.5 -.25];         % Baseline period
        dataBaselinePeriod = ft_selectdata(cfg, dataTFR);

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 4;             % Smoothening frequency around foi
        disp('BASELINED VALUES MAKE NO SENSE IN TRIALS ANALYSES BECAUSE FREQUENCY SMOOTHING IN BASELINE IS 4HZ');
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'yes';        % Keep trial information for trial-wise baseline
        cfg.pad = 5;                   % Zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind1;
        powload1_baseline_period = ft_freqanalysis(cfg, dataBaselinePeriod);
        cfg.trials = ind2;
        powload2_baseline_period = ft_freqanalysis(cfg, dataBaselinePeriod);
        cfg.trials = ind3;
        powload3_baseline_period = ft_freqanalysis(cfg, dataBaselinePeriod);

        % Save baseline trial-wise spectra
        cd(datapath)
        save power_nback_baseline_period_trials ...
            powload1_baseline_period powload2_baseline_period powload3_baseline_period

    catch ME
        ME.message
        error(['ERROR extracting trial-wise power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER (trial-wise), IAF (subject-level) and LATERALIZATION INDEX (trial-wise)
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define channels (from one subject’s labels)
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback_late_trials.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload1_late.label)
    label = powload1_late.label{i};
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

eeg_data_nback_trials = struct('Trial', {}, 'ID', {}, 'Condition', {}, ...
    'AlphaPowerEarly', {}, 'AlphaPowerEarlyBL', {}, 'AlphaPowerLate', {}, ...
    'AlphaPowerLateBL', {}, 'AlphaPowerFull', {}, 'AlphaPowerFullBL', {}, 'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power (trial-wise), IAF (subject), Lateralization (trial-wise) for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);

        % Load trial-wise spectra
        load('power_nback_early_trials.mat');
        load('power_nback_late_trials.mat');
        load('power_nback_full_trials.mat');
        load('power_nback_baseline_period_trials.mat');

        % Channel selection
        channelIdx = find(ismember(powload1_early.label, channels));

        % ----------------------
        % Subject-level IAF (from late retention, trial-averaged ROI)
        % ----------------------
        % Build subject-level ROI-averaged spectra (average across trials, then across ROI channels)
        % We average trials first to stabilise the IAF estimate.
        % Ensure dims: rpt x chan x freq
        S1 = squeeze(mean(powload1_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S2 = squeeze(mean(powload2_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S3 = squeeze(mean(powload3_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        subjSpec = nanmean([S1; S2; S3], 1);                                      % pooled across trials & loads
        freqs = powload1_full.freq(:)';

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

        % Define IAF integration band if valid (−4/+2 Hz)
        if ~isnan(IAF_subj)
            IAF_band = [(IAF_subj-4) (IAF_subj+2)];
        else
            IAF_band = [NaN NaN];
        end

        % ----------------------
        % Baseline-relative (per-trial) conversion for EARLY and LATE windows
        % ----------------------

        powload1_early_db = per_trial_db(powload1_early, powload1_baseline_period);
        powload2_early_db = per_trial_db(powload1_early, powload2_baseline_period);
        powload3_early_db = per_trial_db(powload3_early, powload3_baseline_period);

        powload1_late_db  = per_trial_db(powload1_late,  powload1_baseline_period);
        powload2_late_db  = per_trial_db(powload2_late,  powload2_baseline_period);
        powload3_late_db  = per_trial_db(powload3_late,  powload3_baseline_period);

        powload1_full_db  = per_trial_db(powload1_full,  powload1_baseline_period);
        powload2_full_db  = per_trial_db(powload2_full,  powload2_baseline_period);
        powload3_full_db  = per_trial_db(powload3_full,  powload3_baseline_period);

        % ----------------------
        % Trial-wise Alpha Power (EARLY/LATE, RAW/BASELINED) at IAF band
        % ----------------------
        % EARLY RAW
        AlphaPowerEarly2   = bandpower_trials(powload1_early,  channelIdx, powload1_early.freq,  IAF_band);
        AlphaPowerEarly4   = bandpower_trials(powload1_early,  channelIdx, powload1_early.freq,  IAF_band);
        AlphaPowerEarly6   = bandpower_trials(powload3_early,  channelIdx, powload3_early.freq,  IAF_band);
        % EARLY BL (dB)
        AlphaPowerEarlyBL2 = bandpower_trials(powload1_early_db, channelIdx, powload1_early_db.freq, IAF_band);
        AlphaPowerEarlyBL4 = bandpower_trials(powload2_early_db, channelIdx, powload2_early_db.freq, IAF_band);
        AlphaPowerEarlyBL6 = bandpower_trials(powload3_early_db, channelIdx, powload3_early_db.freq, IAF_band);
        % LATE RAW
        AlphaPowerLate2    = bandpower_trials(powload1_late,   channelIdx,  powload1_late.freq,   IAF_band);
        AlphaPowerLate4    = bandpower_trials(powload2_late,   channelIdx,  powload2_late.freq,   IAF_band);
        AlphaPowerLate6    = bandpower_trials(powload3_late,   channelIdx,  powload3_late.freq,   IAF_band);
        % LATE BL (dB)
        AlphaPowerLateBL2  = bandpower_trials(powload1_late_db, channelIdx,  powload1_late_db.freq, IAF_band);
        AlphaPowerLateBL4  = bandpower_trials(powload2_late_db, channelIdx,  powload2_late_db.freq, IAF_band);
        AlphaPowerLateBL6  = bandpower_trials(powload3_late_db, channelIdx,  powload3_late_db.freq, IAF_band);
        % FULL RAW
        AlphaPowerFull2    = bandpower_trials(powload1_full,   channelIdx,  powload1_full.freq,   IAF_band);
        AlphaPowerFull4    = bandpower_trials(powload2_full,   channelIdx,  powload2_full.freq,   IAF_band);
        AlphaPowerFull6    = bandpower_trials(powload3_full,   channelIdx,  powload3_full.freq,   IAF_band);
        % FULL BL (dB)
        AlphaPowerFullBL2  = bandpower_trials(powload1_full_db, channelIdx,  powload1_full_db.freq, IAF_band);
        AlphaPowerFullBL4  = bandpower_trials(powload2_full_db, channelIdx,  powload2_full_db.freq, IAF_band);
        AlphaPowerFullBL6  = bandpower_trials(powload3_full_db, channelIdx,  powload3_full_db.freq, IAF_band);

        % ----------------------
        % Trial-wise Lateralization Index (use LATE BL by default)
        % ----------------------
        [LI1_trials, ~] = lateralization_trials(powload1_full_db, left_channels, right_channels, powload1_late_db.freq, IAF_band, ridgeFrac, epsP);
        [LI2_trials, ~] = lateralization_trials(powload2_full_db, left_channels, right_channels, powload2_late_db.freq, IAF_band, ridgeFrac, epsP);
        [LI3_trials, ~] = lateralization_trials(powload3_full_db, left_channels, right_channels, powload3_late_db.freq, IAF_band, ridgeFrac, epsP);

        % ----------------------
        % Build subject trial-wise structure array (now with 4 alpha-power fields)
        % ----------------------
        subID = str2double(subjects{subj});
        n1 = size(powload1_full.powspctrm,1);
        n2 = size(powload2_full.powspctrm,1);
        n3 = size(powload3_full.powspctrm,1);

        trl1  = (1:n1)';  trl2  = (1:n2)';  trl3  = (1:n3)';
        IAFr1 = repmat(IAF_subj, n1, 1);
        IAFr2 = repmat(IAF_subj, n2, 1);
        IAFr3 = repmat(IAF_subj, n3, 1);

        subj_data_eeg_trials_1 = struct( ...
            'Trial',              num2cell(trl1), ...
            'ID',                 num2cell(repmat(subID, n1, 1)), ...
            'Condition',          num2cell(repmat(1, n1, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly2), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL2), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate2), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL2), ...
            'IAF',                num2cell(IAFr1), ...
            'Lateralization',     num2cell(LI1_trials) );

        subj_data_eeg_trials_2 = struct( ...
            'Trial',              num2cell(trl2), ...
            'ID',                 num2cell(repmat(subID, n2, 1)), ...
            'Condition',          num2cell(repmat(2, n2, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly4), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL4), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate4), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL4), ...
            'IAF',                num2cell(IAFr2), ...
            'Lateralization',     num2cell(LI2_trials) );

        subj_data_eeg_trials_3 = struct( ...
            'Trial',              num2cell(trl3), ...
            'ID',                 num2cell(repmat(subID, n3, 1)), ...
            'Condition',          num2cell(repmat(3, n3, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly6), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL6), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate6), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL6), ...
            'IAF',                num2cell(IAFr3), ...
            'Lateralization',     num2cell(LI3_trials) );

        subj_data_eeg_trials = [subj_data_eeg_trials_1; subj_data_eeg_trials_2; subj_data_eeg_trials_3];

        % Save (per subject + append to grand table)
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_nback_subj_trials subj_data_eeg_trials
        save IAF_nback_subject IAF_subj

        % Append to grand struct
        eeg_data_nback_trials = [eeg_data_nback_trials; subj_data_eeg_trials];

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
    save W:\Students\Arne\AOC\data\features\eeg_matrix_nback_trials eeg_data_nback_trials
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_nback_trials eeg_data_nback_trials
end

%% ADAPT STERNBERG FOOOFING
%
