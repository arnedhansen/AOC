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

% Plot 
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'W'); 
title('Sternberg WM load 6 Theta, Alpha, and Beta Topolots over Time', 'FontSize', 25)
for f = 1:freqs
    freqband_name = freq_bands{f, 1};
    freqband_range = freq_bands{f, 2};

    % Find max power value for frequency band
    freq_idx = find(gatfr6.freq >= freqband_range(1) & gatfr6.freq <= freqband_range(2));
    timepnts_idx = find(gatfr6.time == 0) : find(gatfr6.time == 2);
    dat = gatfr6.powspctrm(1:125, freq_idx, timepnts_idx);
    dat = abs(dat(:));
    max_spctrm = prctile(dat, 95);

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
        cmap = flipud(cbrewer('div', 'RdBu', 64));
        cfg.colormap = cmap;
        cfg.gridscale = 300;
        cfg.ylim = freqband_range;
        cfg.zlim = [-max_spctrm max_spctrm];
        startTmpnt = find(gatfr6.time == timepoints(tp));
        endTmpnt = find(gatfr6.time == timepoints(tp+1));
        cfg.xlim = [gatfr6.time(startTmpnt) gatfr6.time(endTmpnt)];

        ft_topoplotTFR(cfg, gatfr6);
        title(sprintf('%.0f - %.0f ms', 1000*timepoints(tp), 1000*timepoints(tp+1)), 'FontSize', 12);
    end

    % Colorbar for frequency band
    ax_cb = subplot(freqs, nCols, f*nCols);
    set(ax_cb, 'Visible', 'off');
    cb = colorbar(ax_cb, 'Location', 'eastoutside');
    colormap(cmap);
    clim(cfg.zlim);
    set(cb, 'FontSize', 15);
    ylabel(cb, 'Power [dB]', 'FontSize', 15);
end

% Freq bands
text(-20, 2.85, "Theta (4-8Hz)", 'FontSize', 20, 'FontWeight', 'bold', 'Rotation', 90);
text(-20, 1.35, "Alpha (8-14Hz)", 'FontSize', 20, 'FontWeight', 'bold', 'Rotation', 90);
text(-20, -.05, "Beta (14-30Hz)", 'FontSize', 20, 'FontWeight', 'bold', 'Rotation', 90);

% Save figure
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/raster/AOC_alpha_power_sternberg_rasterplot_WM6.png');
