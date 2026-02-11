%% AOC TFR Raster — Sternberg (Baselined)
% Loads tfr_stern (tfr*_bl), grand-averages, plots raster TFR topoplots. Saves figures.
%
% Key outputs:
%   TFR raster figures (baselined, per condition)

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

% Compute grand average per condition
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr4 = ft_freqgrandaverage([],tfr4_all{:});
gatfr6 = ft_freqgrandaverage([],tfr6_all{:});

% Collapse across all WM loads
gatfr = gatfr2;
gatfr.powspctrm = (gatfr2.powspctrm + gatfr4.powspctrm + gatfr6.powspctrm) / 3;

%% Plot alpha power TOPOS
close all

% Set up frequency bands and timepoints
freq_bands = {
    'Theta', [4 8];
    'Alpha', [8 14];
    'Beta', [14 30];
    };
freqs = size(freq_bands, 1);
timepoints = 0:0.25:2;
timepnts = length(timepoints)-1;
nCols = timepnts + 1;

% Colormap (shared across all subplots)
cmap = flipud(cbrewer('div', 'RdBu', 64));

% Plot
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W');
sgtitle('Sternberg Theta, Alpha, and Beta Topoplots over Time', 'FontSize', 20, 'FontWeight', 'bold')

for f = 1:freqs
    freqband_name = freq_bands{f, 1};
    freqband_range = freq_bands{f, 2};

    % Find max power value for this frequency band (EEG channels only, collapsed data, all timepoints 0-2s)
    freq_idx = find(gatfr.freq >= freqband_range(1) & gatfr.freq <= freqband_range(2));
    timepnts_idx = find(gatfr.time >= 0 & gatfr.time <= 2);
    dat = gatfr.powspctrm(1:125, freq_idx, timepnts_idx);
    max_spctrm = prctile(abs(dat(:)), 95);

    for tp = 1:timepnts
        ax = subplot(freqs, nCols, (f-1)*nCols + tp);
        cfg = [];
        cfg.layout = headmodel.layANThead;
        allchannels = cfg.layout.label;
        cfg.figure = ax;
        cfg.channel = allchannels(1:end-2);
        cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
        cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
        cfg.highlight = 'on';
        cfg.highlightchannel = channels;
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 5;
        cfg.marker = 'off';
        cfg.comment = 'no';
        cfg.colormap = cmap;
        cfg.gridscale = 300;
        cfg.ylim = freqband_range;
        cfg.zlim = [-max_spctrm max_spctrm];
        startTmpnt = find(gatfr.time == timepoints(tp));
        endTmpnt = find(gatfr.time == timepoints(tp+1));
        cfg.xlim = [gatfr.time(startTmpnt) gatfr.time(endTmpnt)];

        ft_topoplotTFR(cfg, gatfr);
        title(sprintf('%.0f - %.0f ms', 1000*timepoints(tp), 1000*timepoints(tp+1)), 'FontSize', 12);
    end

    % Colorbar for this frequency band
    ax_cb = subplot(freqs, nCols, f*nCols);
    set(ax_cb, 'Visible', 'off');
    cb = colorbar(ax_cb, 'Location', 'eastoutside');
    colormap(ax_cb, cmap);
    clim(ax_cb, [-max_spctrm max_spctrm]);
    set(cb, 'FontSize', 15);
    ylabel(cb, 'Power [dB]', 'FontSize', 15);
end

% Frequency band row labels
ax_label = axes('Position', [0 0 1 1], 'Visible', 'off');
text(ax_label, 0.02, 0.80, 'Theta (4-8Hz)', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', 'Units', 'normalized');
text(ax_label, 0.02, 0.48, 'Alpha (8-14Hz)', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', 'Units', 'normalized');
text(ax_label, 0.02, 0.17, 'Beta (14-30Hz)', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', 'Units', 'normalized');

% Save figure
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/raster/AOC_eeg_topos_sternberg_raster.png');
