%% Merge ET and EEG data
%
% Use data after AUTOMAGIC

%% Setup
startup
clear
if ispc == 1
    addpath W:\4marius_bdf\eeglab
else
    addpath /Volumes/methlab/4marius_bdf/eeglab % for pop_importeyetracker (EYE-EEG)
end
eeglab
clc
close all
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\automagic';
else
    path = '/Volumes/methlab/Students/Arne/AOC/data/automagic/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjectIDs = {folders.name};

%% Merge data
tic;
for subjects = 19 : length(subjectIDs) %%% 1 : length(subjectIDs)
    subjectID = subjectIDs{subjects};
    fprintf('Processing Subject %s\n', subjectID)

    % Check if subject files have already been merged
    %if isempty(dir(['/Volumes/methlab/Students/Arne/AOC/data/merged/', char(subjectID), filesep, char(subjectID), '*_merged.mat']))
    % Set up data paths
    if ispc == 1
        filePathET = ['V:\OCC\AOC\data\', char(subjectID)];
        filePathEEG = ['W:\Students\Arne\AOC\data\automagic\',  char(subjectID)];
        resultFolder = ['W:\Students\Arne\AOC\data\merged\', char(subjectID)];
    else
        filePathET = ['/Volumes/methlab_data/OCC/AOC/data/', char(subjectID)];
        filePathEEG = ['/Volumes/methlab/Students/Arne/AOC/data/automagic/',  char(subjectID)];
        resultFolder = ['/Volumes/methlab/Students/Arne/AOC/data/merged/', char(subjectID)];
    end
    mkdir(resultFolder)
    dEEG = dir([filePathEEG, filesep, '*ip*EEG.mat']);
    dET = dir([filePathET, filesep, '*ET.mat']);
    for files = 1 : size(dEEG, 1)
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
            EEG = pop_importeyetracker(EEG, ETfile,[startTrigger endTrigger],[2 3 4],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA', 'R_GAZE_X', 'R_GAZE_Y', 'R_AREA'},1,1,1,1);

            %% Save merge info as image
            set(gcf, "Position", [0 0 1200 800], "Color", "W")
            if ispc == 1
                savepath = strcat('W:\Students\Arne\AOC\data\controls\', sujectID);
            else
                savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/controls/', subjectID);
            end
            mkdir(savepath)
            if strcmp(task, 'Resting') == 1
                taskName = 'Resting';
            elseif strcmp(task, 'AOC_Sternberg') == 1
                taskName = ['Sternberg_block' num2str(block)];
            elseif strcmp(task, 'AOC_Nback') == 1
                taskName = ['Nback_block' num2str(block)];
            end
            saveName = [savepath, filesep, num2str(subjectID) '_mergeInfo_', taskName, '.png'];
            saveas(gcf, saveName);

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
    %end
end
disp('SYNCHRONIZATION COMPLETE')
toc
