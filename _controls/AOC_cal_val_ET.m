%% Computation of EyeLink CALIBRATION and VALIDATION data

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load gaze data
% Define the number of files for each pattern (N and S)
numFiles = 6;

% Loop through each subject
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/archive');
    cd(datapath);

    if str2double(subjects{subj}) == 319 || str2double(subjects{subj}) == 320 || str2double(subjects{subj}) == 378 || str2double(subjects{subj}) == 412
        continue;
    end
    % Pre-allocate CAL and VAL for the current subject
    CAL{subj} = cell(1, numFiles * 2); % Pre-allocating for both N and S files
    VAL{subj} = cell(1, numFiles * 2); % Pre-allocating for both N and S files

    % Initialize file counter
    fileCounter = 1;

    % Loop through all 'N' and 'S' files
    for type = {'Ntk', 'Stk'}
        for fileNum = 1:numFiles
            % Create filename
            filePath = sprintf('%s_%d%s.asc', subjects{subj}, fileNum, type{1});

            % Load data file
            [lastValidationOffsets] = parseASCFile(filePath);

            % Save positions
            VAL{subj}{fileCounter} = lastValidationOffsets;

            % Increment file counter
            fileCounter = fileCounter + 1;
        end
    end
end

%% VISUALIZE for each subject
numSubjects = numel(VAL);
for subjIdx = 1:numSubjects
    try
        close all
        figure;
        set(gcf, 'Position', [300, 200, 1000, 600], 'Color', 'w');

        % Get data
        validationData = VAL{subjIdx};

        % Plot
        hold on;
        for fileIdx = 1:numel(validationData)
            % Extract validation values for the current file
            values = validationData{fileIdx};

            % Plot the values
            plot(fileIdx * ones(size(values)), values, 'bo');
        end

        % Add a red line at y=1
        yline(1, 'r', 'LineWidth', 2);
        yline(1.25, 'k', 'LineWidth', 0.5, 'LineStyle', '--');

        % Set the x-axis and y-axis labels
        xlabel('');
        ylabel('Validation Values [dva]');

        % Set the x-axis limits
        xlim([0, 13]);

        % Set the x-axis ticks and labels
        xticks(1:12);
        xticklabels({'1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6'});

        % Set the y-axis ticks
        yticks(0:0.1:10);

        % Add a title
        title(['Validation Data for Subject ' num2str(subjects{subjIdx})], 'FontSize', 20);
        hold off;

        % Add subtitles
        annotation('textbox', [0.3, 0.03, 0.1, 0.05], 'String', 'N-back', 'EdgeColor', 'none', 'FontSize', 15);
        annotation('textbox', [0.665, 0.03, 0.1, 0.05], 'String', 'Sternberg', 'EdgeColor', 'none', 'FontSize', 15);

        % Save
        savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/controls/ET_validations/');
        mkdir(savepath)
        cd(savepath)
        saveName = [savepath, filesep, num2str(subjects{subjIdx}) '_validations.png'];
        saveas(gcf, saveName)
    catch
        disp(['Could not create VALIDATION overview for subject ' num2str(subjects{subjIdx})])
    end
end

%% Correlate VL with VAL data


%% FUNCTION
function [lastValidationOffsets] = parseASCFile(filePath)
% Initialize cell arrays to store calibration and validation data
calibrationData = {};
validationData = {};

% Open the file for reading
fid = fopen(filePath, 'rt');

if fid == -1
    error('Cannot open file: %s', filePath);
end

% Read the file line by line
line = fgetl(fid);

currentCalibration = [];
currentValidation = [];

while ~contains(line, '!MODE RECORD CR')

    % Check for validation offset lines
    if contains(line, 'VALIDATE') && contains(line, 'POINT')
        tokens = regexp(line, 'OFFSET (\d+\.\d+) deg.', 'tokens');
        if ~isempty(tokens)
            offset = str2double(tokens{1}{1});
            currentValidation = [currentValidation, offset];
        end
    end

    % Read the next line
    line = fgetl(fid);
end

fclose(fid);
try
    lastValidationOffsets = currentValidation(end-8:end);
catch
    lastValidationOffsets = currentValidation;
end
end