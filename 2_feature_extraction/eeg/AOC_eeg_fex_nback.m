%% AOC EEG Feature Extraction N-back
%
% Extracted features:
%   Power Spectrum
%   IAF  and Power at IAF
%   FOOOF Power
%   TFR

%% Setup
startup
[subjects, path, ~ , ant128lay] = setup('AOC');

%% POWER Spectrum
for subj = 1:length(subjects)
    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        %% Identify indices of trials belonging to conditions
        ind1=find(data.trialinfo==21);
        ind2=find(data.trialinfo==22);
        ind3=find(data.trialinfo==23);

        %% Frequency analysis
        cfg = [];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg,data);
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'no';% do not keep single trials in output
        cfg.pad = 10;

        cfg.trials = ind1;
        powload1 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind3;
        powload3 = ft_freqanalysis(cfg,dat);

        %% Save
        cd(datapath)
        save power_nback powload1 powload2 powload3
    catch ME
        ME.message
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER and IAF
% Setup
startup
[subjects, path, ~ , ant128lay] = setup('AOC');

% Define channels
datapath = strcat(path, subjects{1}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

% Load data and calculate power and IAF
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {});

for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_nback.mat');

        % Channels selection based on CBPT
        channelIdx = find(ismember(powload1.label, channels));

        % Extract power spectra for selected channels
        powspctrm1 = mean(powload1.powspctrm(channelIdx, :), 1);
        powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
        powspctrm3 = mean(powload3.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));

        % Calculate IAF for 1-back
        alphaPower1 = powspctrm1(alphaIndices);
        [pks1,locs] = findpeaks(alphaPower1);
        [~, ind] = max(pks1);
        IAF1 = powload1.freq(alphaIndices(locs(ind)));
        IAF_range1 = find(powload1.freq > (IAF1-4) & powload1.freq < (IAF1+2));

        % Calculate IAF for 2-back
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        [~, ind] = max(pks2);
        IAF2 = powload2.freq(alphaIndices(locs(ind)));
        IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));

        % Calculate IAF for 3-back
        alphaPower3 = powspctrm3(alphaIndices);
        [pks3,locs] = findpeaks(alphaPower3);
        [~, ind] = max(pks3);
        IAF3 = powload3.freq(alphaIndices(locs(ind)));
        IAF_range3 = find(powload3.freq > (IAF3-4) & powload3.freq < (IAF3+2));

        % Store the power values at the calculated IAFs
        powerIAF1 = mean(powspctrm1(IAF_range1));
        powerIAF2 = mean(powspctrm2(IAF_range2));
        powerIAF3 = mean(powspctrm3(IAF_range3));

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF1 == alphaRange(1) || IAF1 == alphaRange(2)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF3 == alphaRange(1) || IAF3 == alphaRange(2)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF1 > max(pks1)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF3 > max(pks3)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Create a structure array for this subject
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'AlphaPower', num2cell([powerIAF1; powerIAF2; powerIAF3]), 'IAF', num2cell([IAF1; IAF2; IAF3]));

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf('Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), 3-back: %f Hz (Power: %f)\n', subjects{subj}, IAF1, powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3);
    catch ME
        ME.message
        error(['ERROR calculating alpha power and IAF for Subject ' num2str(subjects{subj}) '!'])
    end
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\eeg_matrix_nback eeg_data_nback
else
    save /Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_nback eeg_data_nback
end

%% POWER WITH TRIAL INFO
% Setup
startup
[subjects, path, ~ , ant128lay] = setup('AOC');

for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        %% Identify indices of trials belonging to conditions
        ind1=find(data.trialinfo==21);
        ind2=find(data.trialinfo==22);
        ind3=find(data.trialinfo==23);

        %% Frequency analysis
        cfg=[];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg,data);
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'yes';% do not keep single trials in output
        cfg.pad = 10;
        cfg.trials = ind1;
        powload1_trials = ft_freqanalysis(cfg,dat);
        powload1_trials.trialinfo = ones(1,length(powload1_trials.powspctrm(:, 1, 1)));
        cfg.trials = ind2;
        powload2_trials = ft_freqanalysis(cfg,dat);
        powload2_trials.trialinfo = ones(1,length(powload2_trials.powspctrm(:, 1, 1)))*2;
        cfg.trials = ind3;
        powload3_trials = ft_freqanalysis(cfg,dat);
        powload3_trials.trialinfo = ones(1,length(powload1_trials.powspctrm(:, 1, 1)))*3;

        %% Save
        cd(datapath)
        save power_nback_trials powload1_trials powload2_trials powload3_trials
    catch ME
        ME.message
        error(['ERROR extracting trial-by-trial power for Subject ' num2str(subjects{subj}) '!'])
    end
end
%% FOOOF Power
% Setup
startup
[subjects, path, ~ , ant128lay] = setup('AOC');

% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        %% Identify indices of trials belonging to conditions
        ind1=find(data.trialinfo==21);
        ind2=find(data.trialinfo==22);
        ind3=find(data.trialinfo==23);

        %% Frequency analysis
        cfg = [];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg,data);
        cfg = []; % empty config
        cfg.output = 'fooof_peaks'; % 1/f
        cfg.method = 'mtmfft'; % multi taper fft method
        cfg.taper = 'dpss'; % multiple tapers
        cfg.tapsmofrq = 1; % smoothening frequency around foi
        cfg.foilim = [3 30]; % frequencies of interest (foi)
        cfg.keeptrials = 'no'; % do not keep single trials in output
        cfg.pad = 10;
        cfg.trials = ind1;
        powload1 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind3;
        powload3 = ft_freqanalysis(cfg,dat);

        %% Save data
        cd(datapath)
        save power_nback_fooof powload1 powload2 powload3
    catch ME
        ME.message
        error(['ERROR extracting FOOOF power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% TFR
% Setup
startup
[subjects, path, ~ , ant128lay] = setup('AOC');

% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_nback

        %% Identify indices of trials belonging to conditions
        ind1 = find(dataTFR.trialinfo==21);
        ind2 = find(dataTFR.trialinfo==22);
        ind3 = find(dataTFR.trialinfo==23);

        %% Time frequency analysis
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:40;                         % analysis 2 to 40 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.toi          = -1:0.05:2;
        cfg.keeptrials = 'no';
        cfg.trials = ind1;
        tfr1 = ft_freqanalysis(cfg,dataTFR);
        cfg.trials = ind2;
        tfr2 = ft_freqanalysis(cfg,dataTFR);
        cfg.trials = ind3;
        tfr3 = ft_freqanalysis(cfg,dataTFR);

        %% Baseline
        cfg                              = [];
        cfg.baseline                     = [-1.5 -0.5];
        cfg.baselinetype                 = 'db';
        tfr1_bl                          = ft_freqbaseline(cfg, tfr1);
        tfr2_bl                          = ft_freqbaseline(cfg, tfr2);
        tfr3_bl                          = ft_freqbaseline(cfg, tfr3);

        %% Save data
        cd(datapath)
        save tfr_nback tfr1 tfr2 tfr3 tfr1_bl tfr2_bl tfr3_bl
    catch ME
        ME.message
        error(['ERROR extracting TFR for Subject ' num2str(subjects{subj}) '!'])
    end
end