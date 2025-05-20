%% AOC Alpha Power Sternberg (trial-by-trial)
% Outputs: eeg_data_sternberg_trials (struct array)

%% Startup
startup;           % initialise paths, EEGLAB, etc.
clear; clc; close all;
addEEGLab;         % make sure pop_epoch etc. are on the path

%% Define paths & subjects
if ispc
    basepath = 'W:\Students\Arne\AOC\data\merged\';
else
    basepath = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
end

% find subject folders
dirs    = dir(basepath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name},{'.','..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');

% prepare output
eeg_data_sternberg_trials = struct('ID',{},'Trial',{},'Condition',{},'AlphaPower',{},'IAF',{});

%% Loop over subjects
for isub = 1:numel(subjects)
    clc
    subj = subjects{isub};
    fprintf('Subject %s (%d/%d)\n', subj, isub, numel(subjects));

    % Load all blocks
    alleeg = cell(1,6);
    for block = 1:6
        fname = fullfile(basepath, subj, sprintf('%s_EEG_ET_Sternberg_block%d_merged.mat', subj, block));
        if exist(fname,'file')
            tmp = load(fname);
            alleeg{block} = tmp.EEG;
            fprintf('Subject %s: Block %.1d loaded \n', subj, block)
        else
            warning('  Block %d missing, skipping', block);
        end
    end
    alleeg = alleeg(~cellfun(@isempty,alleeg));
    if isempty(alleeg)
        warning('  No data for subject %s, skipping\n', subj);
        continue;
    end

    % Epoch all Sternberg stimuli in each block
    epoch_window = [-2 3.5];  % seconds
    dataBlocks = {};
    for b = 1:numel(alleeg)
        EEGb = alleeg{b};
        try
            EEGepo = pop_epoch(EEGb, {'22','24','26'}, epoch_window, 'epochinfo','yes' );
        catch ME
            warning('  Error epoching block %d: %s', b, ME.message);
            continue;
        end

        % convert to FieldTrip raw data structure
        ftb = eeglab2fieldtrip(EEGepo, 'raw', 'none' );

        % extract trialinfo = trigger code for each epoch, in order
        nTr = EEGepo.trials;
        conds = nan(nTr,1);
        for t = 1:nTr
            for nTypes = 1:length(EEGepo.epoch(t).eventtype)
                et = EEGepo.epoch(t).eventtype{nTypes};
                disp(et)
                if strcmp(et, '22')|| strcmp(et, '24')|| strcmp(et, '26')
                    conds(t)= str2double(et);
                end
            end
        end
        ftb.trialinfo = conds;
        % ftb.trialOrder = b*length(conds)-(length(conds)-1):b*length(conds);

        dataBlocks{end+1} = ftb;  %#ok<AGROW>
    end

    % if no epoched blocks, skip
    if isempty(dataBlocks)
        warning('  No epochs for subject %s, skipping\n', subj);
        continue;
    end

    % Append blocks so trials remain in order
    cfg        = [];
    cfg.keepsampleinfo = 'yes';
    dataAll    = ft_appenddata(cfg, dataBlocks{:});

    % EEG data for Sternberg retention interval
    cfg = [];
    cfg.latency = [1 2]; % Time window for Sternberg task
    dataAll = ft_selectdata(cfg, dataAll); % EEG data

    % Re-reference data to average or common reference
    cfg = [];
    cfg.reref   = 'yes';
    cfg.refchannel = 'all';
    dataAll = ft_preprocessing(cfg, dataAll);
    dataAll.trialinfo
    % dataAll.trialOrder = ftb.trialOrder;

    % Single-trial frequency analysis (3�30Hz)
    cfg         = [];
    cfg.output  = 'pow';
    cfg.method  = 'mtmfft';
    cfg.taper   = 'dpss';
    cfg.tapsmofrq = 1;
    cfg.foilim  = [3 30];
    cfg.pad     = 5;
    cfg.keeptrials = 'yes';
    freqAll     = ft_freqanalysis(cfg, dataAll);

    % carry over trialinfo
    freqAll.trialinfo = dataAll.trialinfo;

    % Identify occipital channels
    occIdx = find(contains(freqAll.label, 'O')& ~contains(freqAll.label, {'EOG'})| contains(freqAll.label, {'I'} ));

    % Set up alpha band indices
    freqs      = freqAll.freq;
    alphaRange = [8 14];
    alphaIdx   = find(freqs>=alphaRange(1)& freqs<=alphaRange(2));

    % Loop through trials and extract IAF & alpha power
    nTrials = size(freqAll.powspctrm,1);
    sid      = str2double(subj);
    for t = 1:nTrials
        % average over occ channels, then take only alpha band
        spec = squeeze(mean(freqAll.powspctrm(t,occIdx,alphaIdx), 2 ));

        % find peaks within 8�14Hz
        [pks, locs] = findpeaks(spec);
        if isempty(pks)
            IAF      = NaN;
            alphaPow = NaN;
        else
            [~, imax] = max(pks);
            IAF      = freqs(alphaIdx(locs(imax)));
            alphaPow = pks(imax);
        end

        % build struct entry
        eeg_data_sternberg_trials(end+1)= struct(...
            'ID',         sid, ...
            'Trial',      t, ...
            'Condition',  freqAll.trialinfo(t)-20, ...  % 222, 244, 266
            'AlphaPower', alphaPow, ...
            'IAF',        IAF );
    end
end

%% Save results
if ispc
    save('W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
else
    save('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
end
fprintf('AOC STERNBERG TRIAL-BY-TRIAL ALPHA POWER COMPUTED');
