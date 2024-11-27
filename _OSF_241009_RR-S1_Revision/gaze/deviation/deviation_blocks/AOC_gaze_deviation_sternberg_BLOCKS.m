%% Gaze Deviation per BLOCK for AOC STERNBERG

%% Setup
clear
clc
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
base_path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Preallocate gaze deviation
gaze_deviation = cell(length(subjects), 1);

%% Load and process data
for subj = 1:length(subjects)
    gaze_deviation_subj = cell(1, 8); % Preallocate for this subject's 8 block
    for block = 1:8
        datapath = strcat(base_path, subjects{subj});
        load([datapath, filesep, num2str(subjects{subj}), '_EEG_ET_Sternberg_block', num2str(block), '_merged.mat']);

        %% Select only ET data
        data = EEG.data(130:end, :);

        %% Remove data points that contain zeros (blinks)
        window_size = 50;
        data = remove_blink_window(data, window_size);

        %% Extract gaze data
        gaze_x = data(1, :);
        gaze_y = data(2, :);
        
        %% Compute gaze deviation with Euclidean distances
        num_points = length(gaze_x);
        gaze_euclidean_dev = zeros(1, num_points - 1);
        % Calculate the Euclidean distances between successive points
        for i = 1:num_points - 1
            dx = gaze_x(i + 1) - gaze_x(i);
            dy = gaze_y(i + 1) - gaze_y(i);
            gaze_euclidean_dev(i) = sqrt(dx^2 + dy^2);
        end
        % Calculate the mean Euclidean distance
        euclidean_dist{subj, block} = gaze_euclidean_dev;
    end
    fprintf('Gaze deviation calculated for subject %d/%d \n', subj, length(subjects))
end

%% Plot average
clc
close all

% Find the maximum length of any block across all subjects
maxLength = 0;
numSubjects = length(subjects);
numblock = 8;
for subj = 1:numSubjects
    for block = 1:numblock
        blockLength{block} = length(euclidean_dist{subj, block});
        if blockLength{block} > maxLength
            maxLength = blockLength{block};
        end
    end
end

% Initialize an array to hold padded data
paddedData = NaN(maxLength, numblock, numSubjects);

% Pad data with NaNs where necessary
for subj = 1:numSubjects
    for block = 1:numblock
        currentLength = length(euclidean_dist{subj, block});
        paddedData(1:currentLength, block, subj) = euclidean_dist{subj, block};
    end
end

% Calculate the mean across subjects, ignoring NaNs
averageData = nanmean(paddedData, 3); % The third dimension is subjects

% Get rid of outliers
filteredData = zeros(size(averageData));
% Loop through each block (column) to identify and remove outliers
for col = 1:size(averageData, 2)
    % Calculate the first (Q1) and third (Q3) quartiles
    Q1 = quantile(averageData(:, col), 0.25);
    Q3 = quantile(averageData(:, col), 0.75);
    
    % Calculate the Interquartile Range (IQR)
    IQR = Q3 - Q1;
    
    % Determine the lower and upper bounds for outlier detection
    lowerBound = 0.075; % Q1 - 1.5 * IQR
    upperBound = Q3 + 1.5 * IQR;
    
    % Remove outliers by replacing them with NaN 
    filteredData(:, col) = averageData(:, col);
    filteredData(filteredData(:, col) < lowerBound | filteredData(:, col) > upperBound, col) = NaN; % Or replace with mean/median
end
averageData = filteredData;

% Smooth data
windowSize = 51; % Choose appropriate window size
polyOrder = 2; % Polynomial order
for block = 1:numblock
    smoothedData(:, block) = sgolayfilt(averageData(:, block), polyOrder, windowSize);
end
averageData = smoothedData;

% Plotting
figure; set(gcf, 'Color', 'w', 'Position', [100, 100, 2000, 650]);
hold on;
sample_rate = 500; % Samples per second
time_per_sample = 1000 / sample_rate; % Time in ms per sample
skip_samples = 2 * 60 * sample_rate; % Number of samples to skip (2 minutes)

for block = 1:numblock
    time = ((1:size(averageData(:, block))) + (block - 1) * maxLength - 1) * time_per_sample; % Time in milliseconds

    % Plot the data of the first two minutes
    plot(time(1:skip_samples), averageData(1:skip_samples, block), 'Color', 'b');

    % Plot the data excluding the first two minutes
    if skip_samples < length(averageData(:, block))
        plot(time(skip_samples+1:end), averageData(skip_samples+1:end, block), 'Color', 'r'); 
    end

    % Add block separators
    plot([time(end), time(end)], [0 max(averageData(:)) * 1.1], 'k--');

    % Add block label
    text(time(end) - (time(end) / 2), 7.5, ['Block ' num2str(block)], 'Color', 'k', 'FontSize', 15);
end

set(gca, 'FontSize', 20);
ylim([0, max(averageData(:)) * 1.05]);
title('Sternberg Average Gaze Deviation (smoothed) with and without the first 2 minutes excluded');
xlabel('Time [ms]');
ylabel('Gaze Deviation [px]');
legend('Full Data', 'Excluding first 2 Minutes', 'Location', 'southeast');
hold off;

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/gaze/deviation/deviation_blocks/gaze_deviation_avg_smoothed.png');

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
