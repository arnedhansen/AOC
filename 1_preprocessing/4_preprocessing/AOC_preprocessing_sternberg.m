%% AOC Preprocessing for Sternberg task

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Read data, segment and convert to FieldTrip data structure
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)

    %% Read blocks
    for block = 1:6
        load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
        alleeg{block} = EEG;
        clear EEG
        fprintf('Subject %.3d/%.3d: Block %.1d loaded \n', subj, length(subjects), block)
    end

    %% Segment data by conditions
    for block = 1:6
        EEGload2 = pop_epoch(alleeg{block},{'52'},[-2 3.5]);
        EEGload4 = pop_epoch(alleeg{block},{'54'},[-2 3.5]);
        EEGload6 = pop_epoch(alleeg{block},{'56'},[-2 3.5]);
        data2{block} = eeglab2fieldtrip(EEGload2, 'raw');
        data4{block} = eeglab2fieldtrip(EEGload4, 'raw');
        data6{block} = eeglab2fieldtrip(EEGload6, 'raw');
    end

    %% Equalize labels
    for block = 1:6
        data2{block}.label = data2{1}.label;
        data4{block}.label = data4{1}.label;
        data6{block}.label = data6{1}.label;
    end

    %% Append data for conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data2 = ft_appenddata(cfg,data2{:});
    data4 = ft_appenddata(cfg,data4{:});
    data6 = ft_appenddata(cfg,data6{:});

    %% Add trialinfo data
    data2.trialinfo = ones(1,length(data2.trial))*52';
    data4.trialinfo = ones(1,length(data4.trial))*54';
    data6.trialinfo = ones(1,length(data6.trial))*56';

    %% Append
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data = ft_appenddata(cfg,data2,data4,data6);
    trialinfo = [data2.trialinfo,data4.trialinfo,data6.trialinfo];
    data.trialinfo = trialinfo';
    data.cfg =[ ];

    %% Get EyeTracking data
    cfg = [];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'};
    dataet = ft_selectdata(cfg,data);

    %% Get EEG data (excl. ET and EOG data)
    cfg = [];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA'};
    data = ft_selectdata(cfg,data);

    %% Resegment data to avoid filter ringing
    cfg = [];
    dataTFR = ft_selectdata(cfg,data); % TRF data long
    cfg = [];
    cfg.latency = [1 2]; % Time window for Sternberg task
    data = ft_selectdata(cfg,data); % EEG data
    dataet = ft_selectdata(cfg,dataet); % ET data

    %% Re-reference data to average or common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg,data);

    %% Save data
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save dataEEG_sternberg data 
    save dataEEG_TFR_sternberg dataTFR 
    savepathET = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepathET)
    cd(savepathET)
    save dataET_sternberg dataet
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])
end