%% AOC Merge — EEG and Eyelink After Automagic
% Loads Automagic-cleaned EEG and Eyelink .mat, merges them with EYE-EEG (pop_importeyetracker), and saves one merged file per block. Runs: merge over subjects/blocks.
%
% Key outputs:
%   *_EEG_ET_*_blockN_merged.mat per subject and block (Sternberg, N-back)
%   Merged dataset with co-registered EEG and eye-tracking

%% Setup
startup
clear
[~, paths, ~, ~] = setup('AOC');
addpath(paths.eeglab_share);
eeglab
clc
close all

path = paths.automagic_nohp;
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjectIDs = {folders.name};

%% Merge data
tic;
for subjects = 1 : length(subjectIDs)
    subjectID = subjectIDs{subjects};
    fprintf('Processing Subject %s\n', subjectID)

    mergedPattern = fullfile(paths.merged, subjectID, [subjectID, '*_merged.mat']);
    if isempty(dir(mergedPattern))
        filePathET = fullfile(paths.raw_occ, subjectID);
        filePathEEG = fullfile(paths.automagic_nohp, subjectID);
        resultFolder = fullfile(paths.merged, subjectID);
        mkdir(resultFolder)
        dEEG = dir(fullfile(filePathEEG, '*ip*EEG.mat'));
        dET = dir(fullfile(filePathET, '*ET.mat'));
        for files = [size(dEEG, 1), 1 : size(dEEG, 1)-1]
            close all
            try
                ETnameShort = dET(files).name(1:end-7);
                ETname = dET(files).name;

                idxEEG = contains({dEEG.name}, ETnameShort);

                EEGname = dEEG(idxEEG).name;

                load(fullfile(dEEG(idxEEG).folder, EEGname));
                ETfile = fullfile(dET(1).folder, ETname);

                fileTaskName = strsplit(EEGname, '_');
                task = sprintf('%s', char(fileTaskName(3)), '_', char(fileTaskName(4)));
                if strcmp(task, 'Resting_EEG.mat')
                    task = 'Resting';
                end
                if ~strcmp(task, 'Resting')
                    block = sprintf('%s', char(fileTaskName(5)));
                    block = block(6:end);
                end

                %% Define start and end triggers
                % Resting
                if strcmp(task, 'Resting')
                    startTrigger = 10;
                    endTrigger = 90;
                % Sternberg & Nback
                else
                    startTriggers = [31:38, 61:66];
                    endTriggers = [41:48, 71:76];
                    startTriggersCell = arrayfun(@num2str, [31:38, 61:66], 'UniformOutput', 0);

                    startTrigger = startTriggers (ismember(startTriggersCell, {EEG.event.type}));
                    endTrigger = endTriggers(ismember(startTriggersCell, {EEG.event.type}));
                end

                %% Merge files
                EEG = pop_importeyetracker(EEG, ETfile,[startTrigger endTrigger],[2 3 4 5 6 7],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA', 'R_GAZE_X', 'R_GAZE_Y', 'R_AREA'},1,1,1,1);

                %% Save to disk
                if strcmp(task, 'Resting') == 1
                    fileName = [char(subjectID) '_EEG_ET_RestingEO_merged'];
                elseif strcmp(task, 'AOC_Sternberg') == 1
                    fileName = [char(subjectID) '_EEG_ET_Sternberg_block' num2str(block) '_merged'];
                elseif strcmp(task, 'AOC_Nback') == 1
                    fileName = [char(subjectID) '_EEG_ET_Nback_block' num2str(block) '_merged'];
                end
                save(fullfile(resultFolder, fileName), 'EEG', '-v7.3')
                if strcmp(task, 'Resting') == 1
                    disp(['AOC' char(subjectID) ': Resting done' ])
                else
                    step = sprintf('%s', char(fileTaskName(4)), '_', char(fileTaskName(5)));
                    disp(['AOC' char(subjectID) ': ' step ' done' ])
                end
            catch ME
                ME.message;
                disp(ME.message)
                disp(['ERROR merging file ' num2str(files) '!'])
            end
        end
    end
end
disp('SYNCHRONIZATION COMPLETE')
toc
