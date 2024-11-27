%% Function to detect microsaccades and calculate microsaccade rate
% Based on Engbert & Kliegl, 2003

function microsaccade_rate = detect_microsaccades(fsample, velData, trlLength)
% Convolution kernel as per Engbert et al., 2003
kernel = [1 1 0 -1 -1].*(fsample/6);
velthres = 6; % Threshold as per Engbert et al. (2003): 6x median SD of velocity
mindur = 6; % Minimum duration of microsaccades of 6 samples = 12ms at 500 Hz

% Padding and convolution
n = size(kernel, 2);
pad = ceil(n/2);
dat = ft_preproc_padding(velData, 'localmean', pad);
vel = convn(dat, kernel, 'same');
vel = ft_preproc_padding(vel, 'remove', pad);

% Compute velocity thresholds (Engbert et al. 2003, Eqn. 2)
medianstd = sqrt(median(vel.^2, 2, 'omitnan') - (median(vel, 2, 'omitnan')).^2);
radius = velthres * medianstd;

% Microsaccade detection based on threshold crossing
test = sum((vel ./ radius(:, ones(1, size(vel, 2)))).^2, 1);
sacsmp = find(test > 1); % Find samples where velocity exceeds threshold

% Initialize microsaccades array
microsaccades = [];

% Find consecutive samples that are microsaccades
j = find(diff(sacsmp) == 1);
j1 = [j; j + 1];
com = intersect(j, j + 1);
cut = ~ismember(j1, com);
sacidx = reshape(j1(cut), 2, []);

% Loop through detected saccades
for k = 1:size(sacidx, 2)
    duration = sacidx(1, k):sacidx(2, k);
    if length(duration) >= mindur
        % Store the onset, offset, and peak of each microsaccade
        onset = sacsmp(duration(1)); % Onset of the microsaccade
        offset = sacsmp(duration(end)); % Offset of the microsaccade
        peak = sacsmp(duration(round(length(duration) / 2))); % Peak of the microsaccade
        microsaccades = [microsaccades; onset, offset, peak]; % Append to the microsaccades array
    end
end

% Compute microsaccade rate (microsaccades per second)
microsaccade_rate = length(microsaccades) / ((trlLength-1) / fsample);
end
