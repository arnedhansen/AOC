%% AOC check for missing data

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/AOC/data';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Check missing data
missingFiles = {};
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
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
    disp(upper('No files are missing.'));
end

missingFilesPercentage = size(missingFiles, 2) / (length(subjects)*12)*100;
fprintf('%.2f%% of files are missing\n', missingFilesPercentage);

%% bei wie vielen fehlt eine condition
