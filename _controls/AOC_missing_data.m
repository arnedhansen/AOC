%% AOC check for missing data

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Check missing data
missingFiles = {};
for subj = 1 : length(subjects)
    for type = {'Sternberg', 'Nback'}
        for block = 1 : 6
            filePath = ['/Volumes/methlab/Students/Arne/AOC/data/merged/', char(subjects(subj)), filesep, char(subjects(subj)), '_EEG_ET_', char(type), '_block', num2str(block), '_merged.mat'];
            fileName = [char(subjects(subj)), '_EEG_ET_', char(type), '_block', num2str(block), '_merged.mat'];
            if isfile(filePath)
                disp([fileName '...'])
            elseif ~isfile(filePath)
                disp([fileName '... MISSING.'])
                missingFiles{end+1} = fileName; %#ok<SAGROW><
            end
        end
    end
end

%% Display all missing files
if ~isempty(missingFiles)
    disp(' ');
    disp(' ');
    disp('MISSING FILES:');
    disp(' ');
    for i = 1:length(missingFiles)
        disp(missingFiles{i});
    end
else
    disp('No files are missing.');
end