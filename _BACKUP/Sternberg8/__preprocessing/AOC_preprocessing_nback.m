%% AOC Preprocessing for N-back task
% Add EEGLAB temporarly to segment data 
% later on it needs to be removed from matlab path to avoid colision with FT

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
        load(strcat(subjects{subj}, '_EEG_ET_Nback_block',num2str(block),'_merged.mat'))
        alleeg{block} = EEG;
        clear EEG
        fprintf('Subject %.3d/%.3d: Block %.1d loaded \n', subj, length(subjects), block)
    end

    %% Segment data by conditions
    for block = 1:6
        % 1-back
        try
            EEG1back = pop_epoch(alleeg{block},{'21'},[-1.5 2.5]);
            % Exclude matching trials to avoid including motor responses
            matching_trials = find(strcmp({EEG1back.event.type}, '4'));
            exclude_epochs = unique([EEG1back.event(matching_trials).epoch]);
            EEG1back = pop_select(EEG1back, 'notrial', exclude_epochs);
            data1{block} = eeglab2fieldtrip(EEG1back, 'raw');
            fprintf('1-back data processed for block %d.\n', block);
        catch ME
            fprintf('Error processing 1-back for block %d: %s\n', block, ME.message);
            data1{block} = struct;
        end

        % 2-back
        try
            EEG2back = pop_epoch(alleeg{block},{'22'},[-1.5 2.5]);
            % Exclude matching trials to avoid including motor responses
            matching_trials = find(strcmp({EEG2back.event.type}, '4'));
            exclude_epochs = unique([EEG2back.event(matching_trials).epoch]);
            EEG2back = pop_select(EEG2back, 'notrial', exclude_epochs);
            data2{block} = eeglab2fieldtrip(EEG2back, 'raw');
            fprintf('2-back data processed for block %d.\n', block);
        catch ME
            fprintf('Error processing 2-back for block %d: %s\n', block, ME.message);
            data2{block} = struct;
        end

        % 3-back
        try
            EEG3back = pop_epoch(alleeg{block},{'23'},[-1.5 2.5]);
            % Exclude matching trials to avoid including motor responses
            matching_trials = find(strcmp({EEG3back.event.type}, '4'));
            exclude_epochs = unique([EEG3back.event(matching_trials).epoch]);
            EEG3back = pop_select(EEG3back, 'notrial', exclude_epochs);
            data3{block} = eeglab2fieldtrip(EEG3back, 'raw');
            fprintf('3-back data processed for block %d.\n', block);
        catch ME
            fprintf('Error processing 3-back for block %d: %s\n', block, ME.message);
            data3{block} = struct;
        end
    end

    %% Equalize labels
    % for block = 1:6
    %     data1{block}.label = data1{1}.label;
    %     data2{block}.label = data2{1}.label;
    %     data3{block}.label = data3{1}.label;
    % end

    %% Delete empty blocks
    block = 1;
    while block <= length(data1)
        if isfield(data1{block}, 'fsample')
            block = block + 1;
        else
            data1(block) = [];
        end
    end
    block = 1;
    while block <= length(data2)
        if isfield(data2{block}, 'fsample')
            block = block + 1;
        else
            data2(block) = [];
        end
    end
        block = 1;
    while block <= length(data3)
        if isfield(data3{block}, 'fsample')
            block = block + 1;
        else
            data3(block) = [];
        end
    end

    %% Append data for conditions
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data1 = ft_appenddata(cfg,data1{:});
    data2 = ft_appenddata(cfg,data2{:});
    data3 = ft_appenddata(cfg,data3{:});

    %% Convert to FT data structure
    data1.trialinfo = ones(1,length(data1.trial));
    data2.trialinfo = ones(1,length(data2.trial))*2;
    data3.trialinfo = ones(1,length(data2.trial))*3;
    data1.trialinfo = data1.trialinfo';
    data2.trialinfo = data2.trialinfo';
    data3.trialinfo = data3.trialinfo';

    %% Append
    cfg = [];
    cfg.keepsampleinfo = 'no';
    data = ft_appenddata(cfg,data1,data2,data3);
    trialinfo = [data1.trialinfo,data2.trialinfo,data3.trialinfo];
    data.trialinfo = trialinfo';
    data.cfg =[ ];

    %% Get EyeTracking data
    cfg = [];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'}; % NEEDS TO BE ADJUSTED FOR DISTANCE
    dataet = ft_selectdata(cfg,data);

    %% Get EEG data (excl. ET and EOG data)
    cfg = [];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA'}; % NEEDS TO BE ADJUSTED FOR DISTANCE
    data = ft_selectdata(cfg,data);

    %% Resegment data to avoid filter ringing
    cfg = [];
    dataTFR = ft_selectdata(cfg,data); % TRF data long
    dataTFR.trialinfo = trialinfo';
    cfg = [];
    cfg.latency = [0 2]; % Time window for N-back task
    data = ft_selectdata(cfg,data);
    dataet = ft_selectdata(cfg,dataet);

    %% Re-reference data to average or common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg,data);
    data.trialinfo = trialinfo';

    %% Save to disk
    savepathEEG = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepathEEG)
    cd(savepathEEG)
    save dataEEG_nback data
    save dataEEG_TFR_nback dataTFR 
    savepathET = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/gaze/');
    mkdir(savepathET)
    cd(savepathET)
    save dataET_nback dataet
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])
end