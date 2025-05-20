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
subj_eeg_data_sternberg_trials = struct('ID',{},'Trial',{},'Condition',{},'AlphaPower',{},'IAF',{});
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
                % disp(et)
                if strcmp(et, '22')|| strcmp(et, '24')|| strcmp(et, '26')
                    conds(t)= str2double(et);
                end
            end
        end
        ftb.trialinfo = conds;
        % ftb.trialOrder = b*length(conds)-(length(conds)-1):b*length(conds);

        dataBlocks{end+1} = ftb;
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
    % dataAll.trialOrder = ftb.trialOrder;

    % Single-trial frequency analysis (330Hz)
    cfg            = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmfft';
    cfg.taper      = 'dpss';
    cfg.tapsmofrq  = 1;
    cfg.foilim     = [3 30];
    cfg.pad        = 5;
    cfg.keeptrials = 'yes';
    PowAll         = ft_freqanalysis(cfg, dataAll);

    % carry over trialinfo
    PowAll.trialinfo = dataAll.trialinfo;

    % Identify occipital channels
    occIdx = find(contains(PowAll.label, 'O')& ~contains(PowAll.label, {'EOG'}) | contains(PowAll.label, {'I'} ));

    % Set up alpha band indices
    freqs      = PowAll.freq;
    alphaRange = [8 14];
    alphaIdx   = find(freqs>=alphaRange(1)& freqs<=alphaRange(2));

    % Loop through trials and extract IAF & alpha power
    nTrials = size(PowAll.powspctrm,1);
    sid      = str2double(subj);
    for t = 1:nTrials
        % average over occ channels, then take only alpha band
        spec = squeeze(mean(PowAll.powspctrm(t,occIdx,alphaIdx), 2));

        % find peaks within 814Hz
        [pks, locs] = findpeaks(spec);
        if isempty(pks)
            IAF      = NaN;
            powerIAF = NaN;
        else
            [~, imax] = max(pks);
            IAF       = freqs(alphaIdx(locs(imax)));
            IAF_range = find(freqs > (IAF-4) & freqs < (IAF+2));
            specIAF   = squeeze(mean(PowAll.powspctrm(t,occIdx,IAF_range), 2));
            powerIAF  = mean(specIAF);
            disp(['Alpha Power ' num2str(powerIAF) ' at IAF range ' num2str(freqs(IAF_range(1))) ' to ' num2str(freqs(IAF_range(end))) ' Hz'])
        end

        % build struct entry
        subj_eeg_data_sternberg_trials(end+1) = struct(...
            'ID',         sid, ...
            'Trial',      t, ...
            'Condition',  PowAll.trialinfo(t)-20, ...
            'AlphaPower', powerIAF, ...
            'IAF',        IAF );
        eeg_data_sternberg_trials(end+1) = struct(...
            'ID',         sid, ...
            'Trial',      t, ...
            'Condition',  PowAll.trialinfo(t)-20, ...
            'AlphaPower', powerIAF, ...
            'IAF',        IAF );
    end
    % Save subject results
    if ispc == 1
        save(['W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\eeg_matrix_sternberg_subj_trials.mat'], 'subj_eeg_data_sternberg_trials');
        save(['W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\power_nback_trials.mat'], 'PowAll');
    else
        save(['/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/eeg_matrix_sternberg_subj_trials.mat'], 'subj_eeg_data_sternberg_trials');
        save(['/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/power_nback_trials.mat'], 'PowAll');
    end
end

%% Save results
if ispc
    save('W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
else
    save('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
end
fprintf('AOC STERNBERG TRIAL-BY-TRIAL ALPHA POWER COMPUTED');
