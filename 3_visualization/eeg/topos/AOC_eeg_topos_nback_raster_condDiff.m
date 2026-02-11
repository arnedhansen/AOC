%% AOC TFR Raster Condition Difference — N-Back (Baselined)
% Loads tfr_nback (tfr*_bl), grand-averages, computes difference
% (3-back - 1-back), and plots raster TFR topoplots. Saves figures.
%
% Key outputs:
%   TFR raster figures showing condition difference over time

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');
% Exclude subject 361 because of faulty data in electrode Cz
subjects = setdiff(subjects, {'361'});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload1.label)
    label = powload1.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_nback
    tfr1_all{subj} = tfr1_bl;
    tfr3_all{subj} = tfr3_bl;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

% Compute grand average per condition
gatfr1 = ft_freqgrandaverage([],tfr1_all{:});
gatfr3 = ft_freqgrandaverage([],tfr3_all{:});

% Compute condition difference (highest WM load - lowest WM load)
gatfr_diff = gatfr3;
gatfr_diff.powspctrm = gatfr3.powspctrm - gatfr1.powspctrm;

%% Plot condition difference TOPOS
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
sgtitle('N-back Condition Difference (3-back - 1-back) over Time', 'FontSize', 20, 'FontWeight', 'bold')

for f = 1:freqs
    freqband_name = freq_bands{f, 1};
    freqband_range = freq_bands{f, 2};

    % Find max power value for this frequency band (EEG channels only, difference data, all timepoints 0-2s)
    freq_idx = find(gatfr_diff.freq >= freqband_range(1) & gatfr_diff.freq <= freqband_range(2));
    timepnts_idx = find(gatfr_diff.time >= 0 & gatfr_diff.time <= 2);
    dat = gatfr_diff.powspctrm(1:125, freq_idx, timepnts_idx);
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
        startTmpnt = find(gatfr_diff.time == timepoints(tp));
        endTmpnt = find(gatfr_diff.time == timepoints(tp+1));
        cfg.xlim = [gatfr_diff.time(startTmpnt) gatfr_diff.time(endTmpnt)];

        ft_topoplotTFR(cfg, gatfr_diff);
        title(sprintf('%.0f - %.0f ms', 1000*timepoints(tp), 1000*timepoints(tp+1)), 'FontSize', 12);
    end

    % Colorbar for this frequency band
    ax_cb = subplot(freqs, nCols, f*nCols);
    set(ax_cb, 'Visible', 'off');
    cb = colorbar(ax_cb, 'Location', 'eastoutside');
    colormap(ax_cb, cmap);
    clim(ax_cb, [-max_spctrm max_spctrm]);
    set(cb, 'FontSize', 15);
    ylabel(cb, '\Delta Power [dB]', 'FontSize', 15);
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
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/raster/AOC_eeg_topos_nback_raster_condDiff.png');
