function [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength)
% Detect microsaccades and compute microsaccade rate
% Based on Engbert & Kliegl, 2003

% Define convolution kernel (Engbert et al., 2003)
kernel = [1 1 0 -1 -1] * (fsample / 6); 
velthres = 6; % Threshold: 6x median SD of velocity
mindur = 6;   % Minimum duration of 6 samples = 12ms at 500 Hz

% Padding and convolution
pad = ceil(length(kernel) / 2);
dat = ft_preproc_padding(velData, 'localmean', pad);
vel = convn(dat, kernel, 'same');
vel = ft_preproc_padding(vel, 'remove', pad);

% Compute velocity thresholds (Engbert et al., 2003, Eqn. 2)
medianstd = sqrt(median(vel.^2, 2, 'omitnan') - (median(vel, 2, 'omitnan')).^2);
radius = velthres * medianstd;

% Threshold crossing for microsaccade detection
test = sum((vel ./ radius(:, ones(1, size(vel, 2)))).^2, 1);
sacsmp = find(test > 1);

% Initialise microsaccades array
microsaccades = [];

% Identify consecutive microsaccade samples
if ~isempty(sacsmp)
    sacdiff = diff(sacsmp);
    breaks = [0, find(sacdiff > 1), length(sacsmp)];
    for b = 1:(length(breaks) - 1)
        duration = sacsmp((breaks(b) + 1):breaks(b + 1));
        if length(duration) >= mindur
            onset = duration(1); 
            offset = duration(end);
            peak = duration(round(length(duration) / 2));
            microsaccades = [microsaccades; onset, offset, peak];
        end
    end
end

% Compute microsaccade rate and details
if isempty(microsaccades)
    microsaccade_rate = 0;
    microsaccade_details = struct('Onset', [], 'Offset', [], 'Peak', [], 'Rate', microsaccade_rate);
else
    microsaccade_rate = size(microsaccades, 1) / (trlLength / fsample);
    microsaccade_details = struct('Onset', microsaccades(:, 1), ...
                                   'Offset', microsaccades(:, 2), ...
                                   'Peak', microsaccades(:, 3), ...
                                   'Rate', microsaccade_rate);
end

end
