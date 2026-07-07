%% Function to remove blinks with sliding window
function cleaned_data = remove_blinks(data, window_size)
blink_indices = find(data(2, :) < 100 | data(2, :) > 500);
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end
data(:, removal_indices) = NaN;
cleaned_data = data;
end