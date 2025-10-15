%% AOC EEG Feature Extraction N-back TRIAL-BY-TRIAL
%
% Extracted features:
%   Power Spectrum (Early, Late, Full)  [trial-by-trial]
%   IAF (subject-level), Power at IAF (trial-wise), and Lateralization Index (trial-wise, ridge stabilised)
%   TFR (Raw, FOOOF and Baselined)  [with safer, throttled FOOOF]

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
        load dataEEG_TFR_nback

        % Identify indices of trials belonging to conditions
        globalTrialID = dataTFR.trialinfo(:,2);
        ind1 = find(dataTFR.trialinfo(:, 1) == 21); % 1-back
        ind2 = find(dataTFR.trialinfo(:, 1) == 22); % 2-back
        ind3 = find(dataTFR.trialinfo(:, 1) == 23); % 3-back
        globalTrialID1 = globalTrialID(ind1);
        globalTrialID2 = globalTrialID(ind2);
        globalTrialID3 = globalTrialID(ind3);

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

        freq_all = ft_freqanalysis(cfg, dataTFR); % dimord: rpt_chan_freq_time

        % ----------------------
        % BASELINE (dB)
        % ----------------------
        cfgb = [];
        cfgb.baseline     = [-0.5 -0.25];
        cfgb.baselinetype = 'db';
        freq_all_bl = ft_freqbaseline(cfgb, freq_all);

        % ----------------------
        % EARLY (0-1 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [0 1];
        early_raw = ft_selectdata(cfgsel, freq_all);
        early_db  = ft_selectdata(cfgsel, freq_all_bl);

        % mean over time (4th dim), keep rpt-chan-freq
        early_raw.powspctrm = mean(early_raw.powspctrm, 4);
        early_db.powspctrm  = mean(early_db.powspctrm, 4);
        if isfield(early_raw, 'time'); early_raw = rmfield(early_raw,'time'); end
        if isfield(early_db,  'time'); early_db  = rmfield(early_db, 'time'); end
        early_raw.dimord = 'rpt_chan_freq';
        early_db.dimord  = 'rpt_chan_freq';

        % Split by condition (preserves trial order)
        powload1_early     = ft_selectdata(struct('trials', ind1), early_raw);
        powload2_early     = ft_selectdata(struct('trials', ind2), early_raw);
        powload3_early     = ft_selectdata(struct('trials', ind3), early_raw);
        powload1_early_bl  = ft_selectdata(struct('trials', ind1), early_db);
        powload2_early_bl  = ft_selectdata(struct('trials', ind2), early_db);
        powload3_early_bl  = ft_selectdata(struct('trials', ind3), early_db);

        % Save EARLY trial-wise power spectra
        cd(datapath)
        save power_nback_early_trials powload1_early powload2_early powload3_early ...
            powload1_early_bl powload2_early_bl powload3_early_bl

        % ----------------------
        % LATE (1-2 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [1 2];
        late_raw = ft_selectdata(cfgsel, freq_all);
        late_db  = ft_selectdata(cfgsel, freq_all_bl);

        late_raw.powspctrm = mean(late_raw.powspctrm, 4);
        late_db.powspctrm  = mean(late_db.powspctrm, 4);
        if isfield(late_raw, 'time'); late_raw = rmfield(late_raw,'time'); end
        if isfield(late_db,  'time'); late_db  = rmfield(late_db, 'time'); end
        late_raw.dimord = 'rpt_chan_freq';
        late_db.dimord  = 'rpt_chan_freq';

        % Split by condition
        powload1_late     = ft_selectdata(struct('trials', ind1), late_raw);
        powload2_late     = ft_selectdata(struct('trials', ind2), late_raw);
        powload3_late     = ft_selectdata(struct('trials', ind3), late_raw);
        powload1_late_bl  = ft_selectdata(struct('trials', ind1), late_db);
        powload2_late_bl  = ft_selectdata(struct('trials', ind2), late_db);
        powload3_late_bl  = ft_selectdata(struct('trials', ind3), late_db);

        % Save LATE trial-wise spectra
        cd(datapath)
        save power_nback_late_trials powload1_late powload2_late powload3_late ...
            powload1_late_bl powload2_late_bl powload3_late_bl

        % ----------------------
        % FULL (0-2 s)
        % ----------------------
        cfgsel = [];
        cfgsel.latency = [0 2];
        full_raw = ft_selectdata(cfgsel, freq_all);
        full_db  = ft_selectdata(cfgsel, freq_all_bl);

        full_raw.powspctrm = mean(full_raw.powspctrm, 4);
        full_db.powspctrm  = mean(full_db.powspctrm, 4);
        if isfield(full_raw, 'time'); full_raw = rmfield(full_raw,'time'); end
        if isfield(full_db,  'time'); full_db  = rmfield(full_db, 'time'); end
        full_raw.dimord = 'rpt_chan_freq';
        full_db.dimord  = 'rpt_chan_freq';

        % Split by condition
        powload1_full     = ft_selectdata(struct('trials', ind1), full_raw);
        powload2_full     = ft_selectdata(struct('trials', ind2), full_raw);
        powload3_full     = ft_selectdata(struct('trials', ind3), full_raw);
        powload1_full_bl  = ft_selectdata(struct('trials', ind1), full_db);
        powload2_full_bl  = ft_selectdata(struct('trials', ind2), full_db);
        powload3_full_bl  = ft_selectdata(struct('trials', ind3), full_db);

        % Save FULL trial-wise spectra
        cd(datapath)
        save power_nback_full_trials powload1_full powload2_full powload3_full ...
            powload1_full_bl powload2_full_bl powload3_full_bl

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

        % Channel selection
        channelIdx = find(ismember(powload1_early.label, channels));

        % Rebuild global trial IDs for this subject from the saved trialinfo
        globalTrialID1 = powload1_full.trialinfo(:,2);
        globalTrialID2 = powload2_full.trialinfo(:,2);
        globalTrialID3 = powload3_full.trialinfo(:,2);

        % ----------------------
        % Subject-level IAF (trial-averaged ROI)
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
        if any(isnan(IAF_band))
            IAF_band = alphaRange;  % fallback 8-14 Hz if no clear IAF
        end

        % ----------------------
        % Trial-wise Alpha Power (EARLY/LATE, RAW/BASELINED) at IAF band
        % ----------------------
        % EARLY RAW
        AlphaPowerEarly1   = bandpower_trials(powload1_early,  channelIdx, powload1_early.freq,  IAF_band);
        AlphaPowerEarly2   = bandpower_trials(powload2_early,  channelIdx, powload2_early.freq,  IAF_band);
        AlphaPowerEarly3   = bandpower_trials(powload3_early,  channelIdx, powload3_early.freq,  IAF_band);
        % EARLY BL (dB)
        AlphaPowerEarlyBL1 = bandpower_trials(powload1_early_bl, channelIdx, powload1_early_bl.freq, IAF_band);
        AlphaPowerEarlyBL2 = bandpower_trials(powload2_early_bl, channelIdx, powload2_early_bl.freq, IAF_band);
        AlphaPowerEarlyBL3 = bandpower_trials(powload3_early_bl, channelIdx, powload3_early_bl.freq, IAF_band);
        % LATE RAW
        AlphaPowerLate1    = bandpower_trials(powload1_late,   channelIdx,  powload1_late.freq,   IAF_band);
        AlphaPowerLate2    = bandpower_trials(powload2_late,   channelIdx,  powload2_late.freq,   IAF_band);
        AlphaPowerLate3    = bandpower_trials(powload3_late,   channelIdx,  powload3_late.freq,   IAF_band);
        % LATE BL (dB)
        AlphaPowerLateBL1  = bandpower_trials(powload1_late_bl, channelIdx,  powload1_late_bl.freq, IAF_band);
        AlphaPowerLateBL2  = bandpower_trials(powload2_late_bl, channelIdx,  powload2_late_bl.freq, IAF_band);
        AlphaPowerLateBL3  = bandpower_trials(powload3_late_bl, channelIdx,  powload3_late_bl.freq, IAF_band);
        % FULL RAW
        AlphaPowerFull1    = bandpower_trials(powload1_full,   channelIdx,  powload1_full.freq,   IAF_band);
        AlphaPowerFull2    = bandpower_trials(powload2_full,   channelIdx,  powload2_full.freq,   IAF_band);
        AlphaPowerFull3    = bandpower_trials(powload3_full,   channelIdx,  powload3_full.freq,   IAF_band);
        % FULL BL (dB)
        AlphaPowerFullBL1  = bandpower_trials(powload1_full_bl, channelIdx,  powload1_full_bl.freq, IAF_band);
        AlphaPowerFullBL2  = bandpower_trials(powload2_full_bl, channelIdx,  powload2_full_bl.freq, IAF_band);
        AlphaPowerFullBL3  = bandpower_trials(powload3_full_bl, channelIdx,  powload3_full_bl.freq, IAF_band);

        % ----------------------
        % Trial-wise Lateralization Index
        % ----------------------
        [LI1_trials, ~] = lateralization_trials(powload1_full_bl, left_channels, right_channels, powload1_full_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI2_trials, ~] = lateralization_trials(powload2_full_bl, left_channels, right_channels, powload2_full_bl.freq, IAF_band, ridgeFrac, epsP);
        [LI3_trials, ~] = lateralization_trials(powload3_full_bl, left_channels, right_channels, powload3_full_bl.freq, IAF_band, ridgeFrac, epsP);

        % ----------------------
        % Build subject trial-wise structure array (now with 4 alpha-power fields)
        % ----------------------
        subID = str2double(subjects{subj});
        n1 = size(powload1_full.powspctrm,1);
        n2 = size(powload2_full.powspctrm,1);
        n3 = size(powload3_full.powspctrm,1);

        IAFr1 = repmat(IAF_subj, n1, 1);
        IAFr2 = repmat(IAF_subj, n2, 1);
        IAFr3 = repmat(IAF_subj, n3, 1);

        subj_data_eeg_trials_1 = struct( ...
            'Trial',              num2cell(globalTrialID1), ...
            'ID',                 num2cell(repmat(subID, n1, 1)), ...
            'Condition',          num2cell(repmat(1, n1, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly1), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL1), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate1), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL1), ...
            'AlphaPowerFull',     num2cell(AlphaPowerFull1), ...
            'AlphaPowerFullBL',   num2cell(AlphaPowerFullBL1), ...
            'IAF',                num2cell(IAFr1), ...
            'Lateralization',     num2cell(LI1_trials) );

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

        subj_data_eeg_trials_3 = struct( ...
            'Trial',              num2cell(globalTrialID3), ...
            'ID',                 num2cell(repmat(subID, n3, 1)), ...
            'Condition',          num2cell(repmat(3, n3, 1)), ...
            'AlphaPowerEarly',    num2cell(AlphaPowerEarly3), ...
            'AlphaPowerEarlyBL',  num2cell(AlphaPowerEarlyBL3), ...
            'AlphaPowerLate',     num2cell(AlphaPowerLate3), ...
            'AlphaPowerLateBL',   num2cell(AlphaPowerLateBL3), ...
            'AlphaPowerFull',     num2cell(AlphaPowerFull3), ...
            'AlphaPowerFullBL',   num2cell(AlphaPowerFullBL3), ...
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
