function AOC_DataCuttingFunction(filePath)
% ANT EEG data comes in 1 file containing all tasks (e.g. Resting, CDA, etc.).
% This function cuts the data into distinct tasks, saves the tasks to .mat files
% and moves the original file to filePath/archive
%
% EEGlab with 'pop_loadeep_v4' required

%% Load the data
[EEGorig, command] = pop_loadeep_v4(filePath);

%% Extract folder, base file name and extension using fileparts
[folder, baseFile, ext] = fileparts(filePath);
fileName = [baseFile, ext];
subjectParts = strsplit(baseFile, '_');
subjectID = subjectParts{1};

%% Remove photodiode data and save to a file
diode = pop_select(EEGorig, 'channel', 129);
savename_pd = [subjectID, '_Photodiode.mat'];
save(fullfile(folder, savename_pd), 'diode', '-v7.3')

% Exclude photodiode data if present
if EEGorig.nbchan > 128
    EEGorig = pop_select(EEGorig, 'nochannel', 129:EEGorig.nbchan);
end

%% Add reference channel and data, and load channel location file
EEGorig.data(129, :) = 0;
EEGorig.nbchan = 129;
EEGorig.chanlocs(129).labels = 'CPz';
locspath = 'standard_1005.elc';
EEGorig = pop_chanedit(EEGorig, 'lookup', locspath);

%% Find start and end triggers of the data recording:
% Resting
i10 = find(ismember({EEGorig.event.type}, '10')); % Start
i90 = find(ismember({EEGorig.event.type}, '90')); % End

% Sternberg block 1
i31 = find(ismember({EEGorig.event.type}, '31'));
i41 = find(ismember({EEGorig.event.type}, '41'));

% Sternberg block 2
i32 = find(ismember({EEGorig.event.type}, '32'));
i42 = find(ismember({EEGorig.event.type}, '42'));

% Sternberg block 3
i33 = find(ismember({EEGorig.event.type}, '33'));
i43 = find(ismember({EEGorig.event.type}, '43'));

% Sternberg block 4
i34 = find(ismember({EEGorig.event.type}, '34'));
i44 = find(ismember({EEGorig.event.type}, '44'));

% Sternberg block 5
i35 = find(ismember({EEGorig.event.type}, '35'));
i45 = find(ismember({EEGorig.event.type}, '45'));

% Sternberg block 6
i36 = find(ismember({EEGorig.event.type}, '36'));
i46 = find(ismember({EEGorig.event.type}, '46'));

% NBack block 1
i61 = find(ismember({EEGorig.event.type}, '61'));
i71 = find(ismember({EEGorig.event.type}, '71'));

% NBack block 2
i62 = find(ismember({EEGorig.event.type}, '62'));
i72 = find(ismember({EEGorig.event.type}, '72'));

% NBack block 3
i63 = find(ismember({EEGorig.event.type}, '63'));
i73 = find(ismember({EEGorig.event.type}, '73'));

% NBack block 4
i64 = find(ismember({EEGorig.event.type}, '64'));
i74 = find(ismember({EEGorig.event.type}, '74'));

% NBack block 5
i65 = find(ismember({EEGorig.event.type}, '65'));
i75 = find(ismember({EEGorig.event.type}, '75'));

% NBack block 6
i66 = find(ismember({EEGorig.event.type}, '66'));
i76 = find(ismember({EEGorig.event.type}, '76'));

%% Cut the data into tasks

% Resting
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10(1)).latency, EEGorig.event(i90(1)).latency]);
    task = [subjectID, '_Resting_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Resting is missing...')
end

% Sternberg Block 1
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i31(end)).latency, EEGorig.event(i41(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block1_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 1 is missing...')
end

% Sternberg Block 2
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i32(end)).latency, EEGorig.event(i42(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block2_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 2 is missing...')
end

% Sternberg Block 3
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i33(end)).latency, EEGorig.event(i43(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block3_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 3 is missing...')
end

% Sternberg Block 4
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i34(end)).latency, EEGorig.event(i44(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block4_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 4 is missing...')
end

% Sternberg Block 5
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i35(end)).latency, EEGorig.event(i45(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block5_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 5 is missing...')
end

% Sternberg Block 6
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i36(end)).latency, EEGorig.event(i46(end)).latency]);
    task = [subjectID, '_AOC_Sternberg_block6_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('Block 6 is missing...')
end

% NBack Block 1
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i61(end)).latency, EEGorig.event(i71(end)).latency]);
    task = [subjectID, '_AOC_Nback_block1_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 1 is missing...')
end

% NBack Block 2
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i62(end)).latency, EEGorig.event(i72(end)).latency]);
    task = [subjectID, '_AOC_Nback_block2_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 2 is missing...')
end

% NBack Block 3
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i63(end)).latency, EEGorig.event(i73(end)).latency]);
    task = [subjectID, '_AOC_Nback_block3_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 3 is missing...')
end

% NBack Block 4
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i64(end)).latency, EEGorig.event(i74(end)).latency]);
    task = [subjectID, '_AOC_Nback_block4_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 4 is missing...')
end

% NBack Block 5
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i65(end)).latency, EEGorig.event(i75(end)).latency]);
    task = [subjectID, '_AOC_Nback_block5_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 5 is missing...')
end

% NBack Block 6
try
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i66(end)).latency, EEGorig.event(i76(end)).latency]);
    task = [subjectID, '_AOC_Nback_block6_task_EEG.mat'];
    save(fullfile(folder, task), 'EEG', '-v7.3')
catch ME
    disp(ME.message)
    warning('NBack Block 6 is missing...')
end

%% Create archive folder and move the original files there
archiveFolder = fullfile(folder, 'archive');
if ~exist(archiveFolder, 'dir')
    mkdir(archiveFolder)
end

% Move the .cnt file
source = fullfile(folder, fileName);
destination = fullfile(archiveFolder);
movefile(source, destination)

% Move the .evt file
source = fullfile(folder, [baseFile, '.evt']);
movefile(source, destination)

% Move the .seg file
source = fullfile(folder, [baseFile, '.seg']);
movefile(source, destination)

end
