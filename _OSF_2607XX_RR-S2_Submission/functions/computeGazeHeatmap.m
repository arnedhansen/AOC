function dataOut = computeGazeHeatmap(data, num_bins, smoothing_factor)
%% Filter data for out-of-screen data points and zeros from blinks

% Filter out data points outside the screen boundaries
valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
valid_data = data(1:3, valid_data_indices);

% Remove blinks with a window of 100ms (= 50 samples/timepoints)
win_size = 50;
data = remove_blinks(valid_data, win_size);

x_positions = data(1, :);
y_positions = data(2, :);

%% Create scatterplot for data check
% figure;
% scatterhist(x_positions, y_positions, 'Location', 'SouthEast', 'Color', 'k', 'Marker', '.');
%
% % Calculate mean values
% mean_x = mean(x_positions);
% mean_y = mean(y_positions);
%
% % Add mean markers and labels
% hold on;
% plot(mean_x, mean_y, 'ro', 'MarkerSize', 10);
%
% % Set axis labels
% xlabel('X Position');
% ylabel('Y Position');
% title('Scatterhist of Eye Tracker Data');
% % Invert y-axis to match the typical eye-tracking coordinate system
% set(gca, 'YDir','reverse')
% xlim([0 800]);
% ylim([0 600]);

%% Bin and smooth data
% Create custom grid for heatmap in pixels
x_grid_pixels = linspace(0, 800, num_bins);
y_grid_pixels = linspace(0, 600, num_bins);

% Bin data
binned_data_pixels = histcounts2(x_positions, y_positions, x_grid_pixels, y_grid_pixels); % unit: coutns/samples
binned_time = binned_data_pixels / 500; % Divide by sampling rate (500Hz) -> unit: seconds per bin

% Normalise by window duration
window_length = size(data, 2) / 500;
binned_rate = binned_time / window_length; % unit: proportion of time spent in that spatial bin

% Apply gaussian smoothing
smoothed_data_pixels = imgaussfilt(binned_rate, smoothing_factor);

% Treat ET data as TFR for stats
freq = [];
freq.powspctrm(1,:,:) = squeeze(smoothed_data_pixels');
freq.time       = x_grid_pixels(2:end); % x-axis
freq.freq       = y_grid_pixels(2:end); % y-axis
freq.label      = {'et'};
freq.dimord     = 'chan_freq_time';

% Save for output
dataOut = freq;
