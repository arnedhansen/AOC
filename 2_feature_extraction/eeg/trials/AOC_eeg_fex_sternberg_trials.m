%% AOC EEG Feature Extraction — Sternberg (Trial-Level)
% Trial-wise power spectra (early/late, retention, baseline), IAF, power at IAF, lateralization, and TFR (raw, FOOOF, baselined) from dataEEG_TFR_sternberg. Saves eeg_matrix_sternberg_trials and per-subject TFR/FOOOF.
%
% Extracted features:
%   Power Spectrum (Early, Late (= Registered Retention))  [trial-by-trial]
%   IAF (subject-level), Power at IAF (trial-wise), and Lateralization Index (trial-wise)
%   TFR (Raw, FOOOF and Baselined)

%% POWSPCTRM (Early, Late (= Registered Retention), Baseline Period) - TRIAL-BY-TRIAL
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

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
        % Common time-frequency transform
        % ----------------------
        cfg = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 3:1:30;                % 1-Hz bins
        cfg.t_ftimwin  = ones(size(cfg.foi))*1; % 1 s windows
        cfg.toi        = -1.5:0.05:3;           % as requested
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'yes';

        tfr_all = ft_freqanalysis(cfg, dataTFR); % dimord: rpt_chan_freq_time

        % Save raw trial-by-trial TFR data
        save tfr_stern_trials tfr_all

        % ----------------------
        % BASELINE (dB)
        % ----------------------
        cfgb = [];
        cfgb.baseline     = [-0.5 -0.25];
        cfgb.baselinetype = 'db';
        tfr_all_bl = ft_freqbaseline(cfgb, tfr_all);

        % ----------------------
        % EARLY (0-1 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [0 1];
        early_raw = ft_selectdata(cfgsel, tfr_all);
        early_db  = ft_selectdata(cfgsel, tfr_all_bl);

        % mean over time (4th dim), keep rpt-chan-freq
        early_raw.powspctrm = mean(early_raw.powspctrm, 4);
        early_db.powspctrm  = mean(early_db.powspctrm, 4);
        if isfield(early_raw, 'time'); early_raw = rmfield(early_raw,'time'); end
        if isfield(early_db,  'time'); early_db  = rmfield(early_db, 'time'); end
        early_raw.dimord = 'rpt_chan_freq';
        early_db.dimord  = 'rpt_chan_freq';

        % Split by condition (preserves trial order)
        powload2_early     = ft_selectdata(struct('trials', ind2), early_raw);
        powload4_early     = ft_selectdata(struct('trials', ind4), early_raw);
        powload6_early     = ft_selectdata(struct('trials', ind6), early_raw);
        powload2_early_bl  = ft_selectdata(struct('trials', ind2), early_db);
        powload4_early_bl  = ft_selectdata(struct('trials', ind4), early_db);
        powload6_early_bl  = ft_selectdata(struct('trials', ind6), early_db);

        % Save EARLY trial-wise power spectra
        cd(datapath)
        save power_stern_early_trials powload2_early powload4_early powload6_early ...
            powload2_early_bl powload4_early_bl powload6_early_bl

        % ----------------------
        % LATE (1-2 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [1 2];
        late_raw = ft_selectdata(cfgsel, tfr_all);
        late_db  = ft_selectdata(cfgsel, tfr_all_bl);

        late_raw.powspctrm = mean(late_raw.powspctrm, 4);
        late_db.powspctrm  = mean(late_db.powspctrm, 4);
        if isfield(late_raw, 'time'); late_raw = rmfield(late_raw,'time'); end
        if isfield(late_db,  'time'); late_db  = rmfield(late_db, 'time'); end
        late_raw.dimord = 'rpt_chan_freq';
        late_db.dimord  = 'rpt_chan_freq';

        % Split by condition
        powload2_late     = ft_selectdata(struct('trials', ind2), late_raw);
        powload4_late     = ft_selectdata(struct('trials', ind4), late_raw);
        powload6_late     = ft_selectdata(struct('trials', ind6), late_raw);
        powload2_late_bl  = ft_selectdata(struct('trials', ind2), late_db);
        powload4_late_bl  = ft_selectdata(struct('trials', ind4), late_db);
        powload6_late_bl  = ft_selectdata(struct('trials', ind6), late_db);

        % Save LATE trial-wise spectra
        cd(datapath)
        save power_stern_late_trials powload2_late powload4_late powload6_late ...
            powload2_late_bl powload4_late_bl powload6_late_bl

        % ----------------------
        % FULL (0-2 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [0 2];
        full_raw = ft_selectdata(cfgsel, tfr_all);
        full_db  = ft_selectdata(cfgsel, tfr_all_bl);

        full_raw.powspctrm = mean(full_raw.powspctrm, 4);
        full_db.powspctrm  = mean(full_db.powspctrm, 4);
        if isfield(full_raw, 'time'); full_raw = rmfield(full_raw,'time'); end
        if isfield(full_db,  'time'); full_db  = rmfield(full_db, 'time'); end
        full_raw.dimord = 'rpt_chan_freq';
        full_db.dimord  = 'rpt_chan_freq';

        % Split by condition
        powload2_full     = ft_selectdata(struct('trials', ind2), full_raw);
        powload4_full     = ft_selectdata(struct('trials', ind4), full_raw);
        powload6_full     = ft_selectdata(struct('trials', ind6), full_raw);
        powload2_full_bl  = ft_selectdata(struct('trials', ind2), full_db);
        powload4_full_bl  = ft_selectdata(struct('trials', ind4), full_db);
        powload6_full_bl  = ft_selectdata(struct('trials', ind6), full_db);

        % Save full trial-wise spectra
        cd(datapath)
        save power_stern_full_trials powload2_full powload4_full powload6_full ...
            powload2_full_bl powload4_full_bl powload6_full_bl

    catch ME
        ME.message
        error(['ERROR extracting trial-wise power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER (trial-wise), IAF (subject-level) and LATERALIZATION INDEX (trial-wise)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

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
    'AlphaPowerLateBL', {}, 'AlphaPowerFull', {}, 'AlphaPowerFullBL', {}, ...
    'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power (trial-wise), IAF (subject), Lateralization (trial-wise) for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);

        % Load trial-wise spectra
        load('power_stern_early_trials.mat');
        load('power_stern_late_trials.mat');
        load('power_stern_full_trials.mat');

        % Channel selection
        channelIdx = find(ismember(powload2_full.label, channels));

        % Rebuild global trial IDs for this subject from the saved trialinfo
        globalTrialID2 = powload2_full.trialinfo(:,2);
        globalTrialID4 = powload4_full.trialinfo(:,2);
        globalTrialID6 = powload6_full.trialinfo(:,2);

        % ----------------------
        % Subject-level IAF (from full retention, trial-averaged ROI)
        % ----------------------
        % Build subject-level ROI-averaged spectra (average across trials, then across ROI channels)
        % We average trials first to stabilise the IAF estimate.
        % Ensure dims: rpt x chan x freq
        S2 = squeeze(mean(powload2_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S4 = squeeze(mean(powload4_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        S6 = squeeze(mean(powload6_full.powspctrm(:, channelIdx, :), 2));   % [rpt x freq]
        subjSpec = nanmean([S2; S4; S6], 1);                                      % pooled across trials & loads
        freqs = powload2_full.freq(:)';

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
        if any(isnan(IAF_band))
            IAF_band = alphaRange;  % fallback 8-14 Hz if no clear IAF
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
        % FULL RAW
        AlphaPowerFull2    = bandpower_trials(powload2_full,   channelIdx,  powload2_full.freq,   IAF_band);
        AlphaPowerFull4    = bandpower_trials(powload4_full,   channelIdx,  powload4_full.freq,   IAF_band);
        AlphaPowerFull6    = bandpower_trials(powload6_full,   channelIdx,  powload6_full.freq,   IAF_band);
        % FULL BL (dB)full
        AlphaPowerFullBL2  = bandpower_trials(powload2_full_bl, channelIdx,  powload2_full_bl.freq, IAF_band);
        AlphaPowerFullBL4  = bandpower_trials(powload4_full_bl, channelIdx,  powload4_full_bl.freq, IAF_band);
        AlphaPowerFullBL6  = bandpower_trials(powload6_full_bl, channelIdx,  powload6_full_bl.freq, IAF_band);

        % ----------------------
        % Trial-wise Lateralization Index (with LATE BL)
        % ----------------------
        [LI2_trials, ~] = lateralization_trials(powload2_late_bl, left_channels, right_channels, powload2_late_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI4_trials, ~] = lateralization_trials(powload4_late_bl, left_channels, right_channels, powload4_late_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI6_trials, ~] = lateralization_trials(powload6_late_bl, left_channels, right_channels, powload6_late_bl.freq, IAF_band, ridgeFrac, epsP);

        % ----------------------
        % Build subject trial-wise structure array
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
            'AlphaPowerFull',     num2cell(AlphaPowerFull2), ...
            'AlphaPowerFullBL',   num2cell(AlphaPowerFullBL2), ...
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
            'AlphaPowerFull',     num2cell(AlphaPowerFull4), ...
            'AlphaPowerFullBL',   num2cell(AlphaPowerFullBL4), ...
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
            'AlphaPowerFull',     num2cell(AlphaPowerFull6), ...
            'AlphaPowerFullBL',   num2cell(AlphaPowerFullBL6), ...
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
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_sternberg_trials eeg_data_sternberg_trials
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_sternberg_trials eeg_data_sternberg_trials
end
