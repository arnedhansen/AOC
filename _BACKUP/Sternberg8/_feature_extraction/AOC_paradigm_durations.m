%% AOC Durations of Paradigms

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
durationNback = [];
durationSternberg = [];

%% Read N-back data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    
    %% Read blocks
    for block = 1:6
        fileNback = dir(strcat(subjects{subj}, '_AOC_Nback_block', num2str(block), '_*_task.mat'));
        load(fileNback.name)
        durationNback = [durationNback; saves.data.timing.duration'];
    end
end

%% Read Sternberg data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    
    %% Read blocks
    for block = 1:8
        load(strcat(subjects{subj}, '_AOC_Sternberg_block', num2str(block), '_task.mat'))
        durationSternberg = [durationSternberg; saves.timing.duration'];
    end
end

%% Display
clc
disp('Nback')
disp(durationNback)
disp(['Total: ', num2str(sum(durationNback)) 's = ' num2str((sum(durationNback))/60), 'min'])

disp('Sternberg')
disp(durationSternberg)
disp(['Total: ', num2str(sum(durationSternberg)) 's = ' num2str((sum(durationSternberg))/60), 'min'])
