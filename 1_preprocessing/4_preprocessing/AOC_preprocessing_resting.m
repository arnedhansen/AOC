%% AOC Preprocessing — Resting State EEG to FieldTrip
% Reads merged EEG–ET resting state file (RestingEO), converts to
% FieldTrip, separates EEG and ET channels, re-references EEG, and
% extracts segment boundaries (fixation cross vs blank periods).
%
% Input:
%   {subjectID}_EEG_ET_RestingEO_merged.mat (EEGLAB format, from merge step)
%
% Key outputs per subject:
%   dataEEG_resting.mat — FieldTrip EEG (continuous, avg-referenced) + segInfo
%   dataET_resting.mat  — FieldTrip ET  (continuous) + segInfo
%
% segInfo fields:
%   .fs        — sampling rate (Hz)
%   .fixcross  — [start, end] sample indices for fixation-cross period
%   .blank     — [start, end] sample indices for blank-screen period

%% Setup
startup
clear
[~, ~, ~, ~] = setup('AOC');
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\merged\';
else
    path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/merged/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');
runMode = askRunMode();

%% Process subjects
tic;
for subj = 1:length(subjects)
    try
        datapath = fullfile(path, subjects{subj});
        cd(datapath)

        % Check if already processed
        if ispc == 1
            checkFile = dir(['W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\dataEEG_resting.mat']);
        else
            checkFile = dir(['/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/dataEEG_resting.mat']);
        end

        if strcmp(runMode, 'all') || isempty(checkFile)
            clc
            disp(['Preprocessing Resting State — Subject AOC ', num2str(subjects{subj})])

            %% Load merged EEGLAB file
            restingFile = fullfile(datapath, [subjects{subj}, '_EEG_ET_RestingEO_merged.mat']);
            if ~exist(restingFile, 'file')
                fprintf('No resting state file for subject %s. SKIPPING...\n', subjects{subj});
                continue;
            end
            load(restingFile);  % → EEG (EEGLAB struct)
            fprintf('  Loaded: %s (%d channels, %d samples, %.1f s)\n', ...
                restingFile, EEG.nbchan, EEG.pnts, EEG.pnts / EEG.srate);

            %% Extract trigger events for segment boundaries
            segInfo       = struct();
            segInfo.fs    = EEG.srate;

            % Find trigger events (handle both string and numeric types)
            evTypes     = {EEG.event.type};
            evLatencies = [EEG.event.latency];

            % Helper: find trigger by value (handles string or numeric)
            findTrig = @(val) find(cellfun(@(x) ...
                (ischar(x) && strcmp(x, num2str(val))) || ...
                (isnumeric(x) && x == val), evTypes));

            startIdx = findTrig(10);   % START trigger
            flipIdx  = findTrig(50);   % FLIPTONOFIXCROSS trigger
            endIdx   = findTrig(90);   % END trigger

            if ~isempty(startIdx) && ~isempty(flipIdx)
                segInfo.fixcross = [round(evLatencies(startIdx(1))), round(evLatencies(flipIdx(1)))];
            else
                segInfo.fixcross = [1, round(EEG.pnts / 2)];
                fprintf('  WARNING: trigger 10 or 50 not found, using midpoint split.\n');
            end

            if ~isempty(flipIdx) && ~isempty(endIdx)
                segInfo.blank = [round(evLatencies(flipIdx(1))), round(evLatencies(endIdx(1)))];
            else
                segInfo.blank = [round(EEG.pnts / 2) + 1, EEG.pnts];
                fprintf('  WARNING: trigger 50 or 90 not found, using midpoint split.\n');
            end

            fprintf('  Fixation cross: samples %d – %d  (%.1f s)\n', ...
                segInfo.fixcross(1), segInfo.fixcross(2), diff(segInfo.fixcross) / segInfo.fs);
            fprintf('  Blank period:   samples %d – %d  (%.1f s)\n', ...
                segInfo.blank(1), segInfo.blank(2), diff(segInfo.blank) / segInfo.fs);

            % Sanity check: segment durations should be ~150 s (2.5 min)
            fixDur   = diff(segInfo.fixcross) / segInfo.fs;
            blankDur = diff(segInfo.blank)    / segInfo.fs;
            if fixDur < 100 || fixDur > 200
                fprintf('  WARNING: fixcross duration %.1f s outside expected range [100, 200] s.\n', fixDur);
            end
            if blankDur < 100 || blankDur > 200
                fprintf('  WARNING: blank duration %.1f s outside expected range [100, 200] s.\n', blankDur);
            end

            %% Convert to FieldTrip (continuous, single "trial")
            data = eeglab2fieldtrip(EEG, 'raw');
            clear EEG

            %% Separate Eye-Tracking channels
            cfg         = [];
            cfg.channel = {'L-GAZE-X', 'L-GAZE-Y', 'L-AREA', ...
                           'R-GAZE-X', 'R-GAZE-Y', 'R-AREA'};
            dataET = ft_selectdata(cfg, data);

            %% Separate EEG channels (exclude ET, EOG, unused)
            cfg         = [];
            cfg.channel = {'all', '-B*', '-HEOGR', '-HEOGL', '-VEOGU', '-VEOGL', ...
                '-L-GAZE-X', '-L-GAZE-Y', '-L-AREA', ...
                '-R-GAZE-X', '-R-GAZE-Y', '-R-AREA'};
            dataEEG = ft_selectdata(cfg, data);
            clear data

            %% Re-reference EEG to common average
            cfg            = [];
            cfg.reref      = 'yes';
            cfg.refchannel = 'all';
            dataEEG = ft_preprocessing(cfg, dataEEG);

            fprintf('  EEG: %d channels, %d samples\n', numel(dataEEG.label), numel(dataEEG.time{1}));
            fprintf('  ET:  %d channels, %d samples\n', numel(dataET.label),  numel(dataET.time{1}));

            %% Save EEG data
            clc
            if ispc == 1
                savepathEEG = ['W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\'];
            else
                savepathEEG = ['/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/'];
            end
            mkdir(savepathEEG)
            cd(savepathEEG)
            disp('SAVING dataEEG_resting...')
            save('dataEEG_resting', 'dataEEG', 'segInfo', '-v7.3')

            %% Save ET data
            if ispc == 1
                savepathET = ['W:\Students\Arne\AOC\data\features\', subjects{subj}, '\gaze\'];
            else
                savepathET = ['/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/'];
            end
            mkdir(savepathET)
            cd(savepathET)
            disp('SAVING dataET_resting...')
            save('dataET_resting', 'dataET', 'segInfo', '-v7.3')

            clc
            if subj == length(subjects)
                disp(['Subject AOC ', num2str(subjects{subj}), ' (', num2str(subj), '/', num2str(length(subjects)), ') done. RESTING PREPROCESSING FINALIZED.'])
            else
                disp(['Subject AOC ', num2str(subjects{subj}), ' (', num2str(subj), '/', num2str(length(subjects)), ') done. Loading next subject...'])
            end
        else
            disp(['Subject ', num2str(subjects{subj}), ' already done. SKIPPING...'])
        end

    catch ME
        fprintf('[ERROR] Subject %s (iteration %d/%d): %s\n', ...
            subjects{subj}, subj, length(subjects), ME.message);
    end
end
disp('RESTING STATE PREPROCESSING COMPLETE.')
toc;
