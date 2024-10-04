%% Pupil Size per BLOCK for AOC STERNBERG

%% Setup
clear
clc
close all
% path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
% dirs = dir(path);
% folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
% subjects = {folders.name};
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
base_path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Preallocate pupil size
pupil_size = cell(length(subjects), 1);

%% Load and process data
for subj = 1:length(subjects)
    pupil_size_subj = cell(1, 8); % Preallocate for this subject's 8 blocks
    for block = 1:8
        datapath = strcat(base_path, subjects{subj});
        load([datapath, filesep, num2str(subjects{subj}), '_EEG_ET_Sternberg_block', num2str(block), '_merged.mat']);

        %% Select only ET data
        data = EEG.data(130:end, :);

        %% Remove data points that contain zeros (blinks)
        window_size = 50;
        data = remove_blink_window(data, window_size);

        %% Only keep pupil size data
        pupil_data = data(3, :);

        %% Clean data
        meanData = mean(pupil_data);
        stdData = std(pupil_data);
        threshold = 3;
        outliers = (pupil_data < (meanData - threshold * stdData)) | (pupil_data > (meanData + threshold * stdData));
        cleanData = pupil_data(~outliers);

        %% Save clean pupil size data for this block
        pupil_size_subj{block} = cleanData;
    end
    pupil_size{subj} = pupil_size_subj; % Store this subject's cleaned data
    fprintf('pupil size calculated for subject %d/%d \n', subj, length(subjects))
end

% Save pupil size data for all subjects at once
save([base_path, 'pupil_size.mat'], 'pupil_size');
fprintf('Pupil size saved for all subjects. \n');

%% Load data
clc
clear
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
datapath = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

load([datapath, filesep, 'pupil_size']); % Load the 'pupil_size' variable from the file
for subj = 1:length(subjects)
    pupil_size_allsubs{subj} = pupil_size{subj}; % Assign to the corresponding cell in pupil_size_allsubs
end
%% Plot of individual subjects
clc
numSubjects = length(pupil_size_allsubs);
numBlocks = 8;
figure;
set(gcf, 'Color', 'w', 'Position', [0, 0, 2000, 1200]);

% Determining the global y-axis limits
allData = [];
for subj = 1:numSubjects
    for block = 1:numBlocks
        allData = [allData pupil_size_allsubs{1, subj}{1, block}];
    end
end
globalYLim = [min(allData), max(allData)];

% Plot each subject's data
for subj = 1:numSubjects
    subplot(ceil(numSubjects/2), 2, subj); % Adjust the subplot grid as needed
    hold on;
    
    startTime = 1;
    for block = 1:numBlocks
        blockData = pupil_size_allsubs{1, subj}{1, block};
        time = (startTime:(startTime + length(blockData) - 1));
        plot(time, blockData);
        plot([time(end) time(end)], globalYLim, 'k--'); % Block separator        
        startTime = time(end) + 1; % Update start time for the next block
    end
    
    ylim(globalYLim);
    xlim([0 1650000])
    title(sprintf('Subject %d', subj));
    xlabel('Samples');
    ylabel('Pupil Size');
    hold off;
end

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/pupil/pupil_size_blocks/pupil_size_individual.png');

%% Plot average
clc
close all
% Find the maximum length of any block across all subjects
maxLength = 0;
for subj = 1:numSubjects
    for block = 1:numBlocks
        blockLength = length(pupil_size_allsubs{1, subj}{1, block});
        if blockLength > maxLength
            maxLength = blockLength;
        end
    end
end

% Initialize an array to hold padded data
paddedData = NaN(maxLength, numBlocks, numSubjects);

% Pad data with NaNs where necessary
for subj = 1:numSubjects
    for block = 1:numBlocks
        currentLength = length(pupil_size_allsubs{1, subj}{1, block});
        paddedData(1:currentLength, block, subj) = pupil_size_allsubs{1, subj}{1, block};
    end
end

% Calculate the mean across subjects, ignoring NaNs
averageData = nanmean(paddedData, 3); % The third dimension is subjects

% Plot the averaged data
figure; set(gcf, 'Color', 'w', 'Position', [100, 100, 2000, 650]);
hold on;

% Plot blocks
for block = 1:numBlocks
    blockData = averageData(:, block);
    blockData = blockData/1000;
    time = ((1:length(blockData)) + (block - 1) * maxLength - 1) * 2;

    % Plot and get the handle of the line
    if block == 8
        hLine = plot(time, blockData, 'LineWidth', 1, 'Color', [0.91, 0.76, 0.65]); % Light brown for block 8
    else
        hLine = plot(time, blockData, 'LineWidth', 1); % Default color for other blocks
    end

    % Retrieve the colour of the plotted line
    lineColor = get(hLine, 'Color');

    % Use the colour for text. Here, placing text at the start of each block
    text(time(1)+120000, 4.500, ['Block ' num2str(block)], 'Color', lineColor, 'FontSize', 15);

    % Draw a separator
    plot([time(end), time(end)], [0 5], 'k--');
end

ylim([2.300 4.650]);
title('SternbergSIM average pupil size across all subjects');
xlabel('Time [ms]');
ylabel('Pupil Size [a.u.]');
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/pupil/pupil_size_blocks/pupil_size_avg.png');

%% Plot smoothed average
clc
close all

meanData = nanmean(paddedData, 3); 
stdData = nanstd(paddedData, 0, 3); 

% Choose a window size and order for the Savitzky-Golay filter
windowSize = 5001; % Auf 10s
polyOrder = 3; % Polynomial order

% Smoothed data and standard deviation
smoothedData = zeros(size(meanData));
smoothedStd = zeros(size(stdData));
for block = 1:numBlocks
    smoothedData(:, block) = sgolayfilt(meanData(:, block), polyOrder, windowSize);
    smoothedStd(:, block) = sgolayfilt(stdData(:, block), polyOrder, windowSize);
end

% Plotting
figure; set(gcf, 'Color', 'w', 'Position', [100, 100, 2000, 650]);
hold on;

for block = 1:numBlocks
    blockData = smoothedData(:, block);
    blockData = blockData/1000;
    blockStd = smoothedStd(:, block);
    blockStd = blockStd/1000;

    time = ((1:length(blockData)) + (block - 1) * maxLength - 1) * 2;

    % Plot and get the handle of the line
    if block == 8
        hLine = plot(time, blockData, 'LineWidth', 1.5, 'Color', [0.91, 0.76, 0.65]); % Light brown for block 8
    else
        hLine = plot(time, blockData, 'LineWidth', 1.5); % Default color for other blocks
    end

    % Retrieve the colour of the plotted line
    lineColor = get(hLine, 'Color');

    % Use the colour for text. Here, placing text at the start of each block
    text(time(1)+120000, 4.500, ['Block ' num2str(block)], 'Color', lineColor, 'FontSize', 20);

    % Draw a separator
    plot([time(end), time(end)], [0 5], 'k--');
end

ylim([3.1 4.6]);
title('SternbergSIM average pupil size across all subjects (smoothed with 10s window)');
xlabel('Time [ms]');
ylabel('Pupil Size [a.u.]');
hold off;

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/pupil/pupil_size_blocks/pupil_size_avg_smoothed.png');

%% Plot all data vs. skipping first 3 minutes
close all;
sample_rate = 500; % 500 Hz
skip_samples = 2 * 60 * sample_rate; % Number of samples to skip (2 minutes)

% Find the maximum length of any block across all subjects (adjusting for skipping 2 minutes)
maxLength = 0;
for subj = 1:numSubjects
    for block = 1:numBlocks
        blockLength = length(pupil_size_allsubs{1, subj}{1, block});
        if blockLength > maxLength
            maxLength = blockLength;
        end
    end
end

% Initialize an array to hold padded data
paddedData = NaN(maxLength, numBlocks, numSubjects);
paddedData_noAdapt = NaN(maxLength, numBlocks, numSubjects); % Data excluding the first 2 minutes

% Pad data with NaNs where necessary
for subj = 1:numSubjects
    for block = 1:numBlocks
        currentLength = length(pupil_size_allsubs{1, subj}{1, block});
        paddedData(1:currentLength, block, subj) = pupil_size_allsubs{1, subj}{1, block};
        paddedData_noAdapt(skip_samples+1:currentLength, block, subj) = pupil_size_allsubs{1, subj}{1, block}(skip_samples+1:end);
    end
end

% Calculate the mean across subjects, ignoring NaNs
averageData = nanmean(paddedData, 3); % Full data
averageData_noAdapt = nanmean(paddedData_noAdapt, 3); % Data excluding the first 2 minutes

% Smooth both datasets
smoothedData = zeros(size(averageData));
smoothedData_noAdapt = zeros(size(averageData_noAdapt));
windowSize = 5001; % Adjust window size for Savitzky-Golay filter
polyOrder = 3; % Polynomial order for smoothing

for block = 1:numBlocks
    smoothedData(:, block) = sgolayfilt(averageData(:, block), polyOrder, windowSize);
    smoothedData_noAdapt(:, block) = sgolayfilt(averageData_noAdapt(:, block), polyOrder, windowSize);
end

% Plotting the averaged and smoothed data
figure; set(gcf, 'Color', 'w', 'Position', [100, 100, 2000, 650]);
hold on;

for block = 1:numBlocks
    % Full data (blue plot)
    blockData = smoothedData(:, block) / 1000; % Scale for plot readability
    time = ((1:length(blockData)) + (block - 1) * maxLength - 1) * 2; % Time in milliseconds

    % Data excluding first 3 minutes (red plot)
    blockData_noAdapt = smoothedData_noAdapt(:, block) / 1000;

    % Plot the full data
    plot(time, blockData, 'b', 'LineWidth', 1.5);

    % Plot the data excluding the first three minutes 
    plot(time, blockData_noAdapt, 'r', 'LineWidth', 1.5);

    % Add block separators
    plot([time(end), time(end)], [0 5], 'k--');
    
    % Add block label
    text(time(1)+120000, 4.500, ['Block ' num2str(block)], 'Color', 'k', 'FontSize', 15);
end

ylim([3.2 4.6]);
set(gca, 'FontSize', 20)
title('Sternberg Average Pupil Size (smoothed) with and without first 2 minutes excluded', 'FontSize', 25);
xlabel('Time [ms]');
ylabel('Pupil Size [a.u.]');
legend('Full Data', 'Excluding first 2 Minutes', 'Location', 'eastoutside');
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/pupil/pupil_size_blocks/pupil_size_avg_with_vs_without_3min.png');

%% Define function for blink removal
function cleaned_data = remove_blink_window(data, window_size)
blink_indices = find(all(data(1:2, :) == 0, 1));
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end
data(:, removal_indices) = NaN;
cleaned_data = data;
end