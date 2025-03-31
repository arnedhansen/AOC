%% AOC Alpha Power N-back

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
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
% Load power at IAF and IAF data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath);
    load('alpha_power_nback.mat');
    powIAF1(subj) = powerIAF1;
    powIAF2(subj) = powerIAF2;
    powIAF3(subj) = powerIAF3;
end

% Load powspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_nback
    powl1{subj} = powload1;
    powl2{subj} = powload2;
    powl3{subj} = powload3;
end

% Compute grand avg of raw powspctrm data
gapow1_raw = ft_freqgrandaverage([],powl1{:});
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow3_raw = ft_freqgrandaverage([],powl3{:});

% Load baselined powspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_nback_bl
    powl1_bl{subj} = powload1_bl;
    powl2_bl{subj} = powload2_bl;
    powl3_bl{subj} = powload3_bl;
end

% Compute grand avg of baselined powspctrm data
gapow1_bl = ft_freqgrandaverage([], powl1_bl{:});
gapow2_bl = ft_freqgrandaverage([], powl2_bl{:});
gapow3_bl = ft_freqgrandaverage([], powl3_bl{:});

%% Plot grand average alpha power TOPOS
gapow1 = gapow1_raw;
gapow2 = gapow2_raw;
gapow3 = gapow3_raw;

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
[~, channel_idx] = ismember(channels, gapow1.label);
freq_idx = find(gapow1.freq >= 8 & gapow1.freq <= 14);
max_spctrm = max([mean(gapow1.powspctrm(channel_idx, freq_idx), 2); mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow3.powspctrm(channel_idx, freq_idx), 2)]);
cfg.zlim = [0 max_spctrm];

% Plot 1-back
figure('Color', 'w');
set(gcf, 'Position', [0, 300, 800, 600]);
ft_topoplotER(cfg, gapow1);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('1-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_nback_topo1.png');

% Plot 2-back
figure('Color', 'w');
set(gcf, 'Position', [300, 300, 800, 600]);
ft_topoplotER(cfg, gapow2);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('2-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_nback_topo2.png');

% Plot 3-back
figure('Color', 'w');
set(gcf, 'Position', [600, 300, 800, 600]);
ft_topoplotER(cfg, gapow3);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('3-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_nback_topo3.png');

%% Plot grand average alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (3-back - 1-back)
ga_diff = gapow3;
ga_diff.powspctrm = gapow3.powspctrm - gapow1.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]);
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
title('N-back Alpha Power Difference (3-back - 1-back)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/AOC_alpha_power_nback_topo_diff.png');

%% INDIVIDUAL alpha power TOPOS
for subj = 1:length(subjects)
    close all;
    subjectID = subjects{subj};
    outFolder = fullfile('/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos', subjectID);
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end
    indivPow1 = powl1{subj};
    indivPow2 = powl2{subj};
    indivPow3 = powl3{subj};
    
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
    channel_idx = ismember(indivPow1.label, channels);
    freq_idx = find(indivPow1.freq >= 8 & indivPow1.freq <= 14);
    max1 = max(mean(indivPow1.powspctrm(channel_idx, freq_idx), 2));
    max2 = max(mean(indivPow2.powspctrm(channel_idx, freq_idx), 2));
    max3 = max(mean(indivPow3.powspctrm(channel_idx, freq_idx), 2));
    max_spctrm = max([max1, max2, max3]);
    cfg.zlim = [0 max_spctrm];
    
    % 1-back
    figure('Color', 'w');
    set(gcf, 'Position', [0, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow1);
    title(['Sub ', num2str(subjectID), ' 1-back'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_nback_topo1_%s.png', subjectID)));
    
    % 2-back
    figure('Color', 'w');
    set(gcf, 'Position', [300, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow2);
    title(['Sub ', num2str(subjectID), ' 2-back'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_nback_topo2_%s.png', subjectID)));
    
    % 3-back
    figure('Color', 'w');
    set(gcf, 'Position', [600, 300, 800, 600]);
    ft_topoplotER(cfg, indivPow3);
    title(['Sub ', num2str(subjectID), ' 3-back'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_nback_topo3_%s.png', subjectID)));

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
    indivDiff = indivPow3;
    indivDiff.powspctrm = indivPow3.powspctrm - indivPow1.powspctrm;
    ft_topoplotER(cfg, indivDiff);
    cfg.zlim = 'maxabs';
    title(['Sub ', num2str(subjectID), ' Difference'], 'FontSize', 40);
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
    saveas(gcf, fullfile(outFolder, sprintf('AOC_alpha_power_nback_topo_diff_%s.png', subjectID))); 
end
