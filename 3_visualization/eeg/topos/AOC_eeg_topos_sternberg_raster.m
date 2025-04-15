%% AOC Alpha Power Sternberg

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_stern
    tfr2_all{subj} = tfr2_bl;
    tfr4_all{subj} = tfr4_bl;
    tfr6_all{subj} = tfr6_bl;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

% Compute grand average
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr4 = ft_freqgrandaverage([],tfr4_all{:});
gatfr6 = ft_freqgrandaverage([],tfr6_all{:});

%% Plot alpha power TOPOS
% Set up frequency bands and timepoints
freq_bands = {
    'Theta', [4 8];
    'Alpha', [8 14];
    'Beta', [14 30];
    'Gamma', [30 90];
};

% Define timepoints every 200 ms starting from 0
timepoints = 0:0.2:0.4 %%%2;

% Loop over frequency bands and timepoints
figure('Color', 'w', 'Position', [100, 100, 300 * length(timepoints), 300 * size(freq_bands, 1)]);
nRows = size(freq_bands, 1);
nCols = length(timepoints);
fig = figure('Color', 'w', 'Position', [100, 100, 300 * nCols, 300 * nRows]);
t = tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for f = 1:nRows
    band_name = freq_bands{f, 1};
    band_range = freq_bands{f, 2};
    
    for tp = 1:nCols
        nexttile;
        cfg = [];
        cfg.layout = ant128lay;
        cfg.highlight = 'on';
        cfg.highlightchannel = channels;
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 10;
        cfg.marker = 'off';
        cfg.comment = 'no';
        cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
        cfg.colormap = cmap;
        cfg.gridscale = 300;
        cfg.ylim = band_range;
        cfg.zlim = 'maxabs';
        startTmpnt = find(gatfr6.time == timepoints(tp-1))
        endTmpnt = find(gatfr6.time == timepoints(tp))
        cfg.xlim = gatfr6.time(startTmpnt:endTmpnt) 

        % Make FieldTrip plot into current tile
        ft_topoplotTFR(cfg, gatfr6);
        title(sprintf('%.1f s', timepoints(tp)), 'FontSize', 10);
        
        if tp == 1
            ylabel(band_name, 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end

% Save figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_sternberg_rasterplot_WM6.png');
