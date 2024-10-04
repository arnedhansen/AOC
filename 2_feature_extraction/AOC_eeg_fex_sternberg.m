%% AOC EEG Feature Extraction Sternberg
%
% Extracted features:
%   Power Spectrum
%   FOOOF Power
%   TFR

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract POWER
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
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
end

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract POWER WITH TRIAL INFO
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
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
end

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract FOOOF Power
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
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
end

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Extract TFR
% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG_TFR_sternberg
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind2 = find(dataTFR.trialinfo==52);
    ind4 = find(dataTFR.trialinfo==54);
    ind6 = find(dataTFR.trialinfo==56);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:1:31;                         % analysis 2 to 30 Hz in steps of 1 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    tfr2 = ft_freqanalysis(cfg,dataTFR);
    cfg.trials = ind4;
    tfr4 = ft_freqanalysis(cfg,dataTFR);
    cfg.trials = ind6;
    tfr6 = ft_freqanalysis(cfg,dataTFR);

    %% Save data
    cd(datapath)
    save tfr_stern tfr2 tfr4 tfr6
end