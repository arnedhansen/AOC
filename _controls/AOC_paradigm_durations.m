%% AOC Paradigm Durations
% Reads *_AOC_Nback_block*_task.mat and *_AOC_Sternberg_block*_task.mat, extracts trial durations. Computes and displays duration stats (e.g. mean, SD) per task.
%
% Key outputs:
%   durationNback, durationSternberg arrays; console or figure summary

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
durationNback = [];
durationSternberg = [];

%% Read N-back data
for subj = 1:length(subjects)
    try
    datapath = strcat(path, subjects{subj});
    cd(datapath)

    %% Read blocks
    for block = 1:6
        fileNback = dir(strcat(subjects{subj}, '_AOC_Nback_block', num2str(block), '_*_task.mat'));
        load(fileNback.name)
        try
            durationNback = [durationNback; saves.data.timing.duration'];
        end
        try
            durationNback = [durationNback; saves.timing.duration'];
        end
    end
    end
end

%% Read Sternberg data
for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj});
        cd(datapath)

        %% Read blocks
        for block = 1:6
            load(strcat(subjects{subj}, '_AOC_Sternberg_block', num2str(block), '_task.mat'))
            durationSternberg = [durationSternberg; saves.timing.duration'];
        end
    end
end

%% Display
clc
disp('Nback')
%disp(durationNback)
disp(['Total: ', num2str(sum(durationNback)/length(subjects)) 's = ' num2str((sum(durationNback))/60/length(subjects)), 'min'])

disp(' ')
disp('Sternberg')
%disp(durationSternberg)
disp(['Total: ', num2str(sum(durationSternberg)/length(subjects)) 's = ' num2str((sum(durationSternberg))/60/length(subjects)), 'min'])
