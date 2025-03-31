%% AOC Alpha Power Sternberg

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

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
% Load power at IAF
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath);
    load('alpha_power_sternberg.mat');
    powIAF2(subj) = powerIAF2;
    powIAF4(subj) = powerIAF4;
    powIAF6(subj) = powerIAF6;
end

% Load powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern
    powl2{subj} = powload2;
    powl4{subj} = powload4;
    powl6{subj} = powload6;
end

% Compute grand avg of raw powspctrm data
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow4_raw = ft_freqgrandaverage([],powl4{:});
gapow6_raw = ft_freqgrandaverage([],powl6{:});

% Load baselined powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern_bl
    powl2_bl{subj} = powload2_bl;
    powl4_bl{subj} = powload4_bl;
    powl6_bl{subj} = powload6_bl;
end

% Compute grand avg of baselined powspctrm data
gapow2_bl = ft_freqgrandaverage([], powl2_bl{:});
gapow4_bl = ft_freqgrandaverage([], powl4_bl{:});
gapow6_bl = ft_freqgrandaverage([], powl6_bl{:});

%% Plot alpha power TOPOS

gapow2 = gapow2_raw;
gapow4 = gapow4_raw;
gapow6 = gapow6_raw;

close all;
clc;
cfg = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat');
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
addpath('/Users/Arne/Documents/matlabtools/customcolormap/')
cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
[~, channel_idx] = ismember(channels, gapow2.label);
freq_idx = find(gapow2.freq >= 8 & gapow2.freq <= 14);
max_spctrm = max([mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow4.powspctrm(channel_idx, freq_idx), 2); mean(gapow6.powspctrm(channel_idx, freq_idx), 2)]);
cfg.zlim = [0 max_spctrm];

% Plot WM load 2
figure('Color', 'w');
set(gcf, 'Position', [0, 300, 800, 600]);
ft_topoplotER(cfg, gapow2);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 2', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_sternberg_topo2.png');

% Plot WM load 4
figure('Color', 'w');
set(gcf, 'Position', [700, 300, 800, 600]);
ft_topoplotER(cfg, gapow4);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 4', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_sternberg_topo4.png');

% Plot WM load 6
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 800, 600]);
ft_topoplotER(cfg, gapow6);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 6', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_sternberg_topo6.png');

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
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat');
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
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_sternberg_topo_diff.png');

%% INDIVIDUAL alpha power TOPOS
for subj = 1:length(subjects)
    close all;
    subjectID = subjects{subj};
    outFolder = fullfile('/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos', subjectID);
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end
    indivPow2 = powl2{subj};
    indivPow4 = powl4{subj};
    indivPow6 = powl6{subj};
    
    % Set up the configuration for topographic plotting
    cfg = [];
    cfg.layout = layANThead;
    allchannels = cfg.layout.label;
    cfg.channel = allchannels(1:end-2);
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
    cfg.marker = 'off';
    cfg.highlight = 'on';
    cfg.highlightchannel = channels;  
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 10;
    cfg.figure = 'gcf';
    
    % Define custom colormap
    addpath('/Users/Arne/Documents/matlabtools/customcolormap/')
    cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
    cfg.colormap = cmap;
    cfg.gridscale = 300;
    cfg.comment = 'no';
    cfg.xlim = [8 14];
    
    % Determine the maximum power over the occipital channels and frequency range
    channel_idx = ismember(indivPow2.label, channels);
    freq_idx = find(indivPow2.freq >= 8 & indivPow2.freq <= 14);
    max2 = max(mean(indivPow2.powspctrm(channel_idx, freq_idx), 2));
    max4 = max(mean(indivPow4.powspctrm(channel_idx, freq_idx), 2));
    max6 = max(mean(indivPow6.powspctrm(channel_idx, freq_idx), 2));
    max_spctrm = max([max2, max4, max6]);
    cfg.zlim = [0 max_spctrm];
    
    % 1-back
    figure('Color', 'w');
    set(gcf, 'Position', [0, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow2);
    title(['Sub ', num2str(subjectID), ' WM load 2'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_sternberg_topo_WM2_%s.png', subjectID)));
    
    % 2-back
    figure('Color', 'w');
    set(gcf, 'Position', [300, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow4);
    title(['Sub ', num2str(subjectID), ' WM load 4'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_sternberg_topo_WM4_%s.png', subjectID)));
    
    % 3-back
    figure('Color', 'w');
    set(gcf, 'Position', [600, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow6);
    title(['Sub ', num2str(subjectID), ' WM load 6'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_sternberg_topo_WM6_%s.png', subjectID)));

    % Difference topo
    cfg = [];
    cfg.layout = layANThead;
    allchannels = cfg.layout.label;
    cfg.channel = allchannels(1:end-2);
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
    cfg.marker = 'off';
    cfg.highlight = 'on';
    cfg.highlightchannel = channels;  
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 10;
    cfg.gridscale = 300;
    cfg.comment = 'no';
    cfg.xlim = [8 14];
    cfg.figure = 'gcf';
    cmap = cbrewer('div', 'RdBu', 100);
    cmap = max(min(cmap, 1), 0);
    cmap = flipud(cmap);
    cfg.colormap = cmap;
    figure('Color', 'w');
    set(gcf, 'Position', [600, 300, 800, 600]);
    indivDiff = indivPow6;
    indivDiff.powspctrm = indivPow6.powspctrm - indivPow2.powspctrm;
    ft_topoplotER(cfg, indivDiff);
    cfg.zlim = 'maxabs';
    title(['Sub ', num2str(subjectID), ' Difference'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_sternberg_topo_diff_%s.png', subjectID))); 
end

