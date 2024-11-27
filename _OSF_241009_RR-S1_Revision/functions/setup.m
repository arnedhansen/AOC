%% Function to set up data processing for any project

function setup(projectName)
    % Clear environment
    clearvars -except projectName;
    eeglab;
    clc;
    close all;

    % Set the base path according to the provided project name
    baseDir = '/Volumes/methlab/Students/Arne/';
    path = fullfile(baseDir, projectName, 'data/features/');

    % Check if the path exists
    if ~isfolder(path)
        error('The specified project path does not exist: %s', path);
    end

    % List directories in the selected path
    dirs = dir(path);
    folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {folders.name};

    % Display the loaded subjects
    disp('Loaded subjects:');
    disp(subjects);
end
