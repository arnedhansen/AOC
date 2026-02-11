%% AOC Topos — Sternberg (Power)
% Loads power_stern_raw, grand-averages powspctrm, plots topographical maps per condition. Saves figures.
%
% Key outputs:
%   Topographic maps (occipital alpha power, per condition)

startup
[subjects, path, colors, headmodel] = setup('AOC');

%% Load raw powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    clc
    disp('LOADING DATA...')
    disp(subj)
    load power_stern_raw
    powl2{subj} = powload2;
    powl4{subj} = powload4;
    powl6{subj} = powload6;
end

% Compute grand avg of raw powspctrm data
gapow2 = ft_freqgrandaverage([],powl2{:});
gapow4 = ft_freqgrandaverage([],powl4{:});
gapow6 = ft_freqgrandaverage([],powl6{:});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern_raw.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot alpha power TOPOS
close all;
clc;
fontSize = 50;

cfg = [];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.figure = 'gcf';
cfg.marker = 'off';
%cmap = cbrewer('seq', 'Reds', 64); % 'RdBu' for blue to red diverging color map
cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];

% Global alpha zlim across all channels & conditions
[~, channel_idx] = ismember(channels, gapow2.label);
freq_idx = find(gapow2.freq >= 8 & gapow2.freq <= 14);
A2 = mean(gapow2.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
A4 = mean(gapow4.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
A6 = mean(gapow6.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
all_alpha = [A2(:); A4(:); A6(:)];
global_max = prctile(all_alpha,99);
cfg.zlim = [0 global_max];

% Plot WM load 2
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow2);
cb = colorbar;
set(gca, 'FontSize', fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('WM load 2');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_sternberg_load2.png');

% Plot WM load 4
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow4);
cb = colorbar;
set(gca, 'FontSize', fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('WM load 4');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_sternberg_load4.png');

% Plot WM load 6
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow6);
cb = colorbar;
set(gca, 'FontSize', fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('WM load 6');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_sternberg_load6.png');

%% Plot alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (WM load 6 - WM load 2)
ga_diff = gapow6;
ga_diff.powspctrm = gapow6.powspctrm - gapow2.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 1000, 800]);
cfg = [];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.figure = 'gcf';
cfg.marker = 'off';
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);
cfg.colormap = cmap;
cb = colorbar;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
cfg.zlim = 'maxabs';
set(cb, 'FontSize', 20);
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25);
title('Sternberg Task Alpha Power Difference (WM load 6 - WM load 2)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_sternberg_diff.png');
