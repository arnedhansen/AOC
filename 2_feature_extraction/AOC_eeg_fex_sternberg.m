%% AOC EEG Feature Extraction Sternberg
%
% Extracted features:
%   Power Spectrum
%   IAF  and Power at IAF
%   FOOOF Power
%   TFR

%% Setup
[subjects, path] = setup('AOC');

%% Extract POWER
for subj = 1:length(subjects)
    try
        datapath = strcat(path,subjects{subj}, '/eeg');
        cd(datapath)
        close all
        load dataEEG_sternberg
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        ind2=find(data.trialinfo==52);
        ind4=find(data.trialinfo==54);
        ind6=find(data.trialinfo==56);

        %% Frequency analysis
        cfg = [];
        cfg.latency = [1 2];% segment here only for retention interval
        dat = ft_selectdata(cfg,data);
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'no';
        cfg.pad = 10;
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind4;
        powload4 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind6;
        powload6 = ft_freqanalysis(cfg,dat);

        %% Save data
        cd(datapath)
        save power_stern powload2 powload4 powload6
    catch ME
        ME.message
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% Setup
[subjects, path] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data and calculate alpha power and IAF
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {});

for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, '/eeg');
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
        [pks,locs] = findpeaks(alphaPower2);
        [~, ind] = max(pks);
        IAF2 = powload2.freq(alphaIndices(locs(ind)));
        IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));

        % Calculate IAF for WM load 4
        alphaPower4 = powspctrm4(alphaIndices);
        [pks,locs] = findpeaks(alphaPower4);
        [~, ind] = max(pks);
        IAF4 = powload4.freq(alphaIndices(locs(ind)));
        IAF_range4 = find(powload4.freq > (IAF4-4) & powload4.freq < (IAF4+2));

        % Calculate IAF for WM load 6
        alphaPower6 = powspctrm6(alphaIndices);
        [pks,locs] = findpeaks(alphaPower6);
        [~, ind] = max(pks);
        IAF6 = powload6.freq(alphaIndices(locs(ind)));
        IAF_range6 = find(powload6.freq > (IAF6-4) & powload6.freq < (IAF6+2));

        % Store the power values at the calculated IAFs
        powerIAF2 = mean(powspctrm2(IAF_range2));
        powerIAF4 = mean(powspctrm4(IAF_range4));
        powerIAF6 = mean(powspctrm6(IAF_range6));

        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF2 == 8 || IAF2 == 14
            powerIAF2 = NaN;
        end
        if IAF4 == 8 || IAF4 == 14
            powerIAF4 = NaN;
        end
        if IAF6 == 8 || IAF6 == 14
            powerIAF6 = NaN;
        end

        %% Create a structure array for this subject
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'AlphaPower', num2cell([powerIAF2; powerIAF4; powerIAF6]), 'IAF', num2cell([IAF2; IAF4; IAF6]));

        %% Save
        savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_sternberg_subj subj_data_eeg
        save alpha_power_sternberg powerIAF2 powerIAF4 powerIAF6
        save IAF_sternberg IAF2 IAF4 IAF6
        eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
        clc
        fprintf('Subject %s IAF: load2: %f Hz (Power: %f), load4: %f Hz (Power: %f), load6: %f Hz (Power: %f) \n', subjects{subj}, IAF2, powerIAF2, IAF4, powerIAF4, IAF6, powerIAF6);
    catch ME
        ME.message
        error(['ERROR calculating alpha power and IAF for Subject ' num2str(subjects{subj}) '!'])
    end
end
save /Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg eeg_data_sternberg

%% Setup
[subjects, path] = setup('AOC');

%% Extract POWER WITH TRIAL INFO
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    try
        datapath = strcat(path,subjects{subj}, '/eeg');
        cd(datapath)
        close all
        load dataEEG_sternberg
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        ind2=find(data.trialinfo==52);
        ind4=find(data.trialinfo==54);
        ind6=find(data.trialinfo==56);

        %% Frequency analysis
        cfg = [];
        cfg.latency = [1 2];% segment here only for retention interval
        dat = ft_selectdata(cfg,data);
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'yes';
        cfg.pad = 10;
        cfg.trials = ind2;
        powload2_trials = ft_freqanalysis(cfg,dat);
        cfg.trials = ind4;
        powload4_trials = ft_freqanalysis(cfg,dat);
        cfg.trials = ind6;
        powload6_trials = ft_freqanalysis(cfg,dat);

        %% Save data
        cd(datapath)
        save power_stern_trials powload2_trials powload4_trials powload6_trials
    catch ME
        ME.message
        error(['ERROR extracting trial-by-trial power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% Setup
[subjects, path] = setup('AOC');

%% Extract FOOOF Power
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    try
        datapath = strcat(path,subjects{subj}, '/eeg');
        cd(datapath)
        close all
        load dataEEG_sternberg
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        ind2=find(data.trialinfo==52);
        ind4=find(data.trialinfo==54);
        ind6=find(data.trialinfo==56);

        %% Frequency analysis
        cfg = [];
        cfg.latency =[1 2]; % segment here only for retention interval
        dat = ft_selectdata(cfg,data);
        cfg = []; % empty config
        cfg.output = 'fooof_peaks'; % 1/f
        cfg.method = 'mtmfft'; % multi taper fft method
        cfg.taper = 'dpss'; % multiple tapers
        cfg.tapsmofrq = 1; % smoothening frequency around foi
        cfg.foilim = [3 30]; % frequencies of interest (foi)
        cfg.keeptrials = 'no';
        cfg.pad = 10;
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind4;
        powload4 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind6;
        powload6 = ft_freqanalysis(cfg,dat);

        %% Save data
        cd(datapath)
        save power_stern_fooof powload2 powload4 powload6
    catch ME
        ME.message
        error(['ERROR extracting FOOOF power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% Setup
[subjects, path] = setup('AOC');

%% Extract TFR
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    try
        datapath = strcat(path,subjects{subj}, '/eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg
        load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

        %% Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo == 52);
        ind4 = find(dataTFR.trialinfo == 54);
        ind6 = find(dataTFR.trialinfo == 56);

        %% Time frequency analysis
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:40;                         % analysis 2 to 40 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.toi          = -3:0.05:3;
        cfg.keeptrials = 'no';
        cfg.trials = ind2;
        tfr2 = ft_freqanalysis(cfg,dataTFR);
        cfg.trials = ind4;
        tfr4 = ft_freqanalysis(cfg,dataTFR);
        cfg.trials = ind6;
        tfr6 = ft_freqanalysis(cfg,dataTFR);

        %% Baseline
        cfg                              = [];
        cfg.baseline                     = [-1.5 -.5];
        cfg.baselinetype                 = 'db';
        tfr2_bl                          = ft_freqbaseline(cfg, tfr2);
        tfr4_bl                          = ft_freqbaseline(cfg, tfr4);
        tfr6_bl                          = ft_freqbaseline(cfg, tfr6);

        %% Save data
        cd(datapath)
        save tfr_stern tfr2 tfr4 tfr6 tfr2_bl tfr4_bl tfr6_bl
    catch ME
        ME.message
        error(['ERROR extracting TFR for Subject ' num2str(subjects{subj}) '!'])
    end
end