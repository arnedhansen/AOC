%% AOC check for missing data
%
% AOC 319 - Sternberg only block 1
% AOC 320 - N-back only block 1-4; NO ET DATA
% AOC 322 - N-back only block 1-2
% AOC 333 - NO DATA
% AOC 347 - no Sternberg data (triggers bug)
% AOC 351 - No N-back data
% AOC 358 - no data (triggers bug)
% AOC 360 - N-back only block 1, no Sternberg data
% AOC 378 - N-back block 1-5, Sternberg only block 1
% AOC 380 - N-back only block 1
% AOC 381 - no Sternberg data 
% AOC 387 - Only Sternberg block 1
% AOC 405 - no data (triggers bug)
% AOC 407 - N-back only block 1
% AOC 408 - N-back only block 1, no Sternberg data
% AOC 412 - N-back only block 1-3, no Sternberg data
% AOC 414 - 

%% Setup
startup
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

%% Heatmap Visualization
close all

% Define paths
data_path = '/Volumes/methlab_data/OCC/AOC/data';
check_path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';

% Get subject list
dirs = dir(data_path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

% Define conditions and blocks
conditions = {'Sternberg', 'Nback'};
num_blocks = 6;

% Preallocate matrix for missing files (1 = present, 0 = missing)
missingMatrix = ones(length(subjects), length(conditions) * num_blocks);

% Check for missing data
missingFiles = {};
for subj = 1:length(subjects)
    for c = 1:length(conditions)
        for block = 1:num_blocks
            % Construct expected file path
            fileName = sprintf('%s_EEG_ET_%s_block%d_merged.mat', subjects{subj}, conditions{c}, block);
            filePath = fullfile(check_path, subjects{subj}, fileName);
            
            % Check file existence
            if ~isfile(filePath)
                missingFiles{end+1} = fileName; %#ok<SAGROW>
                missingMatrix(subj, (c-1)*num_blocks + block) = 0; % Mark as missing
            end
        end
    end
end

% Compute missing files percentage from the precomputed matrix
missingPercentage = sum(missingMatrix(:) == 0) / numel(missingMatrix) * 100;

% Heatmap figure
figure('Position', [100, 100, 2000, 1200], 'Color', 'W');
imagesc(missingMatrix');
colormap([1 0 0; 0 1 0]); % Red for missing, Green for present
caxis([0 1]);

% Custom X-tick labels (Subjects)
xticks(1:length(subjects));
xticklabels(subjects);

% Load exclusion list
fid = fopen('/Volumes/methlab/Students/Arne/AOC/data/controls/AOC_exclusion_participants.rtf','rt');
raw = fread(fid, '*char')';
fclose(fid);
exclIDs = regexp(raw, '\d+', 'match');  % cell array of strings

% Determine which subjects are excluded
isExcluded = ismember(subjects, exclIDs);

% Replace X‑tick labels with coloured text
hAx = gca;
hAx.XTickLabel = [];   % remove default labels
hAx.XTick = 1:length(subjects);

for i = 1:length(subjects)
    txt = text(i, hAx.YLim(2) + 0.1, subjects{i}, ...
               'Rotation', 90, ...
               'HorizontalAlignment', 'right', ...
               'FontSize', 12);
    if isExcluded(i)
        txt.Color = 'red';
    else
        txt.Color = 'black';
    end
end

% Custom Y-tick labels (Conditions & Blocks)
yticks(1:length(conditions)*num_blocks);
yticklabels(arrayfun(@(i) sprintf('%s Block %d', conditions{ceil(i/num_blocks)}, mod(i-1, num_blocks)+1), 1:length(conditions)*num_blocks, 'UniformOutput', false));

% Labels and title
ylabel('Task & Block');
title(sprintf('AOC Missing Data Heatmap (%.d / %.d (%.2f%%) files missing)', sum(missingMatrix(:) == 0), numel(missingMatrix), missingPercentage), 'FontSize', 50, 'FontWeight', 'bold');

% Add grid lines (grid lines removed and replaced with explicit lines)
hold on;
[numRows, numCols] = size(missingMatrix');
numRows = numRows+1;
numCols = numCols+1;
for i = 1:numCols
    % Vertical lines between columns
    plot([i+0.5, i+0.5], [0, numRows], 'k-', 'LineWidth', 1.5);
end
for j = 1:numRows
    % Horizontal lines between rows
    plot([0, numCols], [j+0.5, j+0.5], 'k-', 'LineWidth', 1.5);
end
hold off;

% Set font sizes for labels and title
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
titleHandle = get(gca, 'Title');
set(titleHandle, 'FontSize', 30);
xlabelHandle = get(gca, 'XLabel');
set(xlabelHandle, 'FontSize', 25);
ylabelHandle = get(gca, 'YLabel');
set(ylabelHandle, 'FontSize', 25);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/data/controls/AOC_missing_data_heatmap.png')
