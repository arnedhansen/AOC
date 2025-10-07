%% AOC EEG Feature Extraction Sternberg
%
% Extracted features:
%   Power Spectrum (Retention)
%   Power Spectrum (Baseline)
%   IAF, Power at IAF, and Lateralization Index
%   TFR (Raw, FOOOF and Baselined)

%% POWSPCTRM (Retention)
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing Retention POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_sternberg

        % Identify indices of trials belonging to conditions
        ind2 = find(dataEEG.trialinfo == 22); % WM load 2
        ind4 = find(dataEEG.trialinfo == 24); % WM load 4
        ind6 = find(dataEEG.trialinfo == 26); % WM load 6

        % Frequency analysis for retention interval = 1000 ms - 2000ms after stimulus presentation
        % Select data
        cfg = [];                          % Empty configuration
        cfg.latency = [1 2];               % Segmentation for retention interval 1000ms - 2000ms
        dat = ft_selectdata(cfg, dataEEG); % Select data

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'no';         % Discard trial information
        cfg.pad = 5;                   % Add zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind4;
        powload4 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind6;
        powload6 = ft_freqanalysis(cfg,dat);

        % Save data
        cd(datapath)
        save power_stern powload2 powload4 powload6

    catch ME
        ME.message
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% POWSPCTRM (Early and Late Retention SPLIT)
%  Early: 0-1000ms after stimulus presentation
%  Late: 1000-2000ms after stimulus presentation

% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing Early and Late Retention SPLIT POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo == 26); % WM load 6

        % Frequency analysis
        % Select EARLY data
        cfg = [];                      % Empty configuration
        cfg.latency = [0 1];
        dataEarly = ft_selectdata(cfg, dataTFR);

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'no';         % Discard trial information
        cfg.pad = 5;                   % Add zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind4;
        powload4_early = ft_freqanalysis(cfg, dataEarly);
        cfg.trials = ind6;
        powload6_early = ft_freqanalysis(cfg, dataEarly);

        % Save baselined power spectra
        cd(datapath)
        save power_stern_early powload2_early powload4_early powload6_early

        % Select LATE data
        cfg = [];                              % Empty configuration
        cfg.latency = [1 2];
        dataLate = ft_selectdata(cfg, dataTFR); % Select data

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'no';         % Discard trial information
        cfg.pad = 5;                   % Add zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind4;
        powload4late = ft_freqanalysis(cfg, dataLate);
        cfg.trials = ind6;
        powload6late = ft_freqanalysis(cfg, dataLate);

        % Save data
        cd(datapath)
        save power_stern_late powload2late powload4late powload6late

    catch ME
        ME.message
        error(['ERROR extracting baslined power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% POWSPCTRM (Baseline & Long Retention Interval)
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing Baseline POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo == 26); % WM load 6

        % Frequency analysis
        % Select data
        cfg = [];                      % Empty configuration
        cfg.latency = [-.5 0];     % Segmentation for retention interval
        datalong = ft_selectdata(cfg, dataTFR);

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'no';         % Discard trial information
        cfg.pad = 5;                   % Add zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2_baseline_period = ft_freqanalysis(cfg, datalong);
        cfg.trials = ind4;
        powload4_baseline_period = ft_freqanalysis(cfg, datalong);
        cfg.trials = ind6;
        powload6_baseline_period = ft_freqanalysis(cfg, datalong);

        % Save baselined power spectra
        cd(datapath)
        save power_stern_baseline_period powload2_baseline_period powload4_baseline_period powload6_baseline_period

        % Frequency analysis for retention interval = 200 ms - 2000ms after stimulus presentation
        disp(['Processing Long Retention Interval (200ms - 2000ms) POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

        % Select data
        cfg = [];                              % Empty configuration
        cfg.latency = [0.2 2];                 % Segmentation for retention interval 200ms - 2000ms
        datLong = ft_selectdata(cfg, dataTFR); % Select data

        % Analysis settings
        cfg = [];                      % Empty configuration
        cfg.output = 'pow';            % Estimate power only
        cfg.method = 'mtmfft';         % Multi-taper FFT method
        cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
        cfg.tapsmofrq = 1;             % Smoothening frequency around foi
        cfg.foilim = [3 30];           % Frequencies of interest
        cfg.keeptrials = 'no';         % Discard trial information
        cfg.pad = 5;                   % Add zero-padding

        % Conduct frequency analysis for each condition separately
        cfg.trials = ind2;
        powload2long = ft_freqanalysis(cfg, datLong);
        cfg.trials = ind4;
        powload4long = ft_freqanalysis(cfg, datLong);
        cfg.trials = ind6;
        powload6long = ft_freqanalysis(cfg, datLong);

        % Save data
        cd(datapath)
        save power_stern_late powload2long powload4long powload6long

    catch ME
        ME.message
        error(['ERROR extracting baslined power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
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
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_stern.mat');

        % Channel selection
        channelIdx = find(ismember(powload2.label, channels));

        % Extract power spectra for selected channels
        powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
        powspctrm4 = mean(powload4.powspctrm(channelIdx, :), 1);
        powspctrm6 = mean(powload6.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(powload2.freq >= alphaRange(1) & powload2.freq <= alphaRange(2));

        % Calculate IAF for WM load 2
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = powload2.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));
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
            IAF4 = powload4.freq(alphaIndices(locs(ind)));
            IAF_range4 = find(powload4.freq > (IAF4-4) & powload4.freq < (IAF4+2));
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
            IAF6 = powload6.freq(alphaIndices(locs(ind)));
            IAF_range6 = find(powload6.freq > (IAF6-4) & powload6.freq < (IAF6+2));
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

        % Compute lateralization index as done in Stroganova et al., 2007
        powloads = {powload2, powload4, powload6};
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

        % Create a structure array for this subject
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'AlphaPower', num2cell([powerIAF2; powerIAF4; powerIAF6]), 'IAF', num2cell([IAF2; IAF4; IAF6]), ...
            'Lateralization', num2cell([LatIdx2; LatIdx4; LatIdx6]));

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
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
        ME.message
        error(['ERROR calculating alpha power and IAF for Subject ' num2str(subjects{subj}) '!'])
    end
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg eeg_data_sternberg
else
    save /Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg eeg_data_sternberg
end

%% TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)

    % Check existing data
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    %if ~isfile([datapath, filesep, 'power_stern_fooof.mat'])
    clc
    disp(['Processing TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM for Subject AOC ', num2str(subjects{subj})])
    try
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo == 22);
        ind4 = find(dataTFR.trialinfo == 24);
        ind6 = find(dataTFR.trialinfo == 26);

        % Time frequency analysis
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:1:40;                       % 2..40 Hz (1-Hz step)
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; % 0.5 s window
        cfg.toi          = -2:0.05:3;
        cfg.keeptrials   = 'yes';                        % need to keep trials for FOOOF on single trials

        cfg.trials = ind2;
        tfr2 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind4;
        tfr4 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind6;
        tfr6 = ft_freqanalysis(cfg, dataTFR);

        % FOOOF
        orig_freq = 2:1:40;
        tfrs = {tfr2, tfr4, tfr6};
        tfrs_fooof_trials = cell(1,3);

        % FOOOF settings for 2–40 Hz
        settings = struct();
        settings.verbose = false;
        settings.peak_width_limits = [0.5 12];
        settings.max_n_peaks = 4;

        for cond_idx = 1:3
            clc
            disp('FOOOFing (keeptrials = yes)...')

            tfr = tfrs{cond_idx};
            freqs = tfr.freq(:);
            nChan = numel(tfr.label);
            nFreq = numel(freqs);
            nTime = numel(tfr.time);
            nTr   = size(tfr.powspctrm, 4);

            % Exclude time points too close to edges for the 0.5 s window
            t_win   = 0.5;
            t_valid = find( tfr.time >= (tfr.time(1)+t_win/2) & ...
                tfr.time <= (tfr.time(end)-t_win/2) );

            % Preallocate: chan x freq x time x trial (NaNs by default)
            fspctrm = nan(nChan, nFreq, nTime, nTr, 'single');

            f_range = [freqs(1) freqs(end)];
            for kk = 1:numel(t_valid)
                tt = t_valid(kk);

                % Select a single time point; FieldTrip will squeeze time to length 1
                cfgSel = [];
                cfgSel.latency = [tfr.time(tt) tfr.time(tt)];
                tmp = ft_selectdata(cfgSel, tfr);  % chan x freq x 1 x trials

                % reshape to chan x freq x trials for convenience
                tmpPSD = squeeze(tmp.powspctrm);   % chan x freq x trials

                % loop trials then channels
                for tr = 1:nTr
                    for ch = 1:nChan
                        psd = squeeze(tmpPSD(ch, :, tr)).';
                        % guard against non-positive / non-finite values
                        if any(~isfinite(psd)) || any(psd <= 0), continue, end

                        % Use tmp.freq to be maximally safe about spacing
                        fr = fooof(tmp.freq(:), psd, f_range, settings, true);

                        % periodic-only (log10 units): fooofed - aperiodic
                        fspctrm(ch, :, tt, tr) = fr.fooofed_spectrum - fr.ap_fit;
                    end
                end
            end

            % Build a FieldTrip object with TRIAL DIMENSION preserved
            tfr_fooof_trials = tfr;
            tfr_fooof_trials.powspctrm = fspctrm;

            % Store
            tfrs_fooof_trials{cond_idx} = tfr_fooof_trials;
        end
        disp(upper('FOOOF (keeptrials) done...'))

        % Baselining 
        % Raw power: baseline in dB, with trials kept
        cfgBL = [];
        cfgBL.baseline     = [-.75 -.25];
        cfgBL.baselinetype = 'db';

        tfr2_bl_trials = ft_freqbaseline(cfgBL, tfr2);
        tfr4_bl_trials = ft_freqbaseline(cfgBL, tfr4);
        tfr6_bl_trials = ft_freqbaseline(cfgBL, tfr6);

        % FOOOF-periodic (log10 units): baseline with subtraction (absolute), with trials kept
        cfgBLP = [];
        cfgBLP.baseline     = [-.75 -.25];
        cfgBLP.baselinetype = 'absolute';

        tfr2_fooof_bl_trials = ft_freqbaseline(cfgBLP, tfrs_fooof_trials{1});
        tfr4_fooof_bl_trials = ft_freqbaseline(cfgBLP, tfrs_fooof_trials{2});
        tfr6_fooof_bl_trials = ft_freqbaseline(cfgBLP, tfrs_fooof_trials{3});

        % Average across trials
        tfr2_bl = tfr2_bl_trials;
        tfr2_bl.powspctrm = mean(tfr2_bl_trials.powspctrm, 4, 'omitnan');

        tfr4_bl = tfr4_bl_trials;
        tfr4_bl.powspctrm = mean(tfr4_bl_trials.powspctrm, 4, 'omitnan');

        tfr6_bl = tfr6_bl_trials;
        tfr6_bl.powspctrm = mean(tfr6_bl_trials.powspctrm, 4, 'omitnan');

        tfr2_fooof_bl = tfr2_fooof_bl_trials;
        tfr2_fooof_bl.powspctrm = mean(tfr2_fooof_bl_trials.powspctrm, 4, 'omitnan');

        tfr4_fooof_bl = tfr4_fooof_bl_trials;
        tfr4_fooof_bl.powspctrm = mean(tfr4_fooof_bl_trials.powspctrm, 4, 'omitnan');

        tfr6_fooof_bl = tfr6_fooof_bl_trials;
        tfr6_fooof_bl.powspctrm = mean(tfr6_fooof_bl_trials.powspctrm, 4, 'omitnan');

        % Save
        cd(datapath)
        save tfr_stern ...
            tfr2 tfr4 tfr6 ...
            tfr2_bl tfr4_bl tfr6_bl ...
            tfr2_bl_trials tfr4_bl_trials tfr6_bl_trials

        save tfr_stern_fooof ...
            tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl ...
            tfr2_fooof_bl_trials tfr4_fooof_bl_trials tfr6_fooof_bl_trials

        clc
    catch ME
        ME.message
        error(['ERROR extracting TFR for Subject ' num2str(subjects{subj}) '!'])
    end
    %else
    %    disp(['TFR and FOOOFed POWSPCTRM already exists for Subject AOC ', num2str(subjects{subj})])
    %end
end
disp('TFR and FOOOFed POWSPCTRM COMPUTED...');
