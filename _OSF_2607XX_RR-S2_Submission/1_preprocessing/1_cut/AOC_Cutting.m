%% AOC Cutting — Segment Raw EEG by Task
% Finds unconverted .cnt files, cuts each into task-specific segments via AOC_DataCuttingFunction, then loads and synchronizes EEG with Eyelink. Runs: cutting, .mat conversion, EEG–ET sync.
%
% Key outputs:
%   Task-specific .mat files (Resting, Sternberg, N-back, etc.) per subject
%   Synchronized EEG and Eyelink in subject folders

%% EEGlab
startup
clear
clc
[~, paths, ~, ~] = setup('AOC', 0);
p = fileparts(mfilename('fullpath'));
cd(paths.eeglab_share)
eeglab
close()
cd(p)
clc

%% Define path, cut data and convert asc to mat files.
tic
d = dir(fullfile(paths.raw_occ, '*', '*cnt'));
disp(upper([num2str(size(d, 1)), ' subjects to compute']))
ids = {};
for f = 1 : size(d, 1)
    filePath = fullfile( d(f).folder, d(f).name);
    AOC_DataCuttingFunction(filePath)
    id = strsplit(d(f).name, '_');
    ids{f} = id{1};
    clc
    fprintf('Cutting of data for Subject AOC %.3s done \n', ids{f})
end

%% Load and synchronize EEG & Eyelink
clc
disp(upper('Loading and Synchronizing EEG & Eyelink'));
for id = 1 : length(ids)
    ID = ids{id};

    filePath = fullfile(paths.raw_occ, ID);
    d = dir(fullfile(filePath, '*.asc'));

    if not(isempty(d))
        fprintf('Synchronizing data for Subject AOC %.3s \n', ID)
        for f = 1 : size(d, 1)

            % Convert ET asc to mat
            inputFile = fullfile(d(f).folder, d(f).name);

            % Rename the files (EyeLink can't deal with long filenames, edf filenames has to be short)
            x = strsplit(d(f).name, '_');
            name = x{2};
            id = x{end-1};
            block = str2double(name(1));

            if isnan(block)
                if name(1) == 'R'
                    task = 'Resting';
                    eegfile = [id, '_Resting_EEG.mat'];
                    etfile = [id, '_Resting_ET.mat'];
                else
                    task = 'Training';
                    eegfile = '';
                    if name(1) == 'N'
                        nameTr = 'AOC_Nback';
                        etfile = [id, '_', nameTr, '_Training.mat'];
                    end
                end
            else
                if name(2) == 'S'
                    task = 'AOC_Sternberg';
                elseif name(2) == 'N'
                    task = 'AOC_Nback';
                end

                etfile = [id, '_' task, '_block', num2str(block), '_task_ET.mat'];
            end

            outputFile = fullfile(d(f).folder, etfile);
            ET = parseeyelink(inputFile, outputFile);
        end
    end

    % Move asc files to archive
    d = dir(fullfile(filePath, '*.asc'));
    for f = 1 : size(d, 1)
        source = fullfile(d(f).folder, d(f).name);
        destination = fullfile(filePath, 'archive');
        movefile(source,destination)
    end

    % Move edf files to archive
    d = dir(fullfile(filePath, '*.edf'));
    for f = 1 : size(d, 1)
        source = fullfile(d(f).folder, d(f).name);
        destination = fullfile(filePath, 'archive');
        movefile(source,destination)
    end
end
disp('SYNCHRONIZATION COMPLETE')
toc;
%finishedScriptMail;
