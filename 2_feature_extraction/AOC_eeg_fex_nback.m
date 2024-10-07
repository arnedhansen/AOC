%% AOC EEG Feature Extraction N-back
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
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
    
    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==1);
    ind2=find(data.trialinfo==2);
    ind3=find(data.trialinfo==3);

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
    cfg.keeptrials = 'no';% do not keep single trials in output
    cfg.pad = 10;
    cfg.trials = ind1;
    powload1= ft_freqanalysis(cfg,dat);
    cfg.trials = ind2;
    powload2= ft_freqanalysis(cfg,dat);
    cfg.trials = ind3;
    powload3= ft_freqanalysis(cfg,dat);

    %% Save
    cd(datapath)
    save power_nback powload1 powload2 powload3

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
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    close all
    load dataEEG_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
    
    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==1);
    ind2=find(data.trialinfo==2);
    ind3=find(data.trialinfo==3);

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
    load dataEEG_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==1);
    ind2=find(data.trialinfo==2);
    ind3=find(data.trialinfo==3);

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
    load dataEEG_TFR_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind1 = find(dataTFR.trialinfo==1);
    ind2 = find(dataTFR.trialinfo==2);
    ind3 = find(dataTFR.trialinfo==3);

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
    load1 = ft_freqanalysis(cfg,dataTFR);
    cfg.trials = ind2;
    load2 = ft_freqanalysis(cfg,dataTFR);
    cfg.trials = ind3;
    load3 = ft_freqanalysis(cfg,dataTFR);

    %% Save data
    cd(datapath)
    save tfr_nback load1 load2 load3
end