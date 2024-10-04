%% AOC Alpha Power Sternberg

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data and calculate alpha power and IAF
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {});

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load('power_stern.mat');
    
    % Channel selection
    channelIdx = find(ismember(powload2.label, channels));
    
    % Extract power spectra for selected channels
    powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
    powspctrm4 = mean(powload4.powspctrm(channelIdx, :), 1);
    powspctrm6 = mean(powload6.powspctrm(channelIdx, :), 1);

    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload2.freq >= alphaRange(1) & powload2.freq <= alphaRange(2));
    
    % Calculate IAF for WM load 2
    alphaPower2 = powspctrm2(alphaIndices);
    [pks,locs] = findpeaks(alphaPower2);
    [~, ind] = max(pks);
    IAF2 = powload2.freq(alphaIndices(locs(ind)));
    IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));

    % Calculate IAF for WM load 4
    alphaPower4 = powspctrm4(alphaIndices);
    [pks,locs] = findpeaks(alphaPower4);
    [~, ind] = max(pks);
    IAF4 = powload4.freq(alphaIndices(locs(ind)));
    IAF_range4 = find(powload4.freq > (IAF4-4) & powload4.freq < (IAF4+2));

    % Calculate IAF for WM load 6
    alphaPower6 = powspctrm6(alphaIndices);
    [pks,locs] = findpeaks(alphaPower6);
    [~, ind] = max(pks);
    IAF6 = powload6.freq(alphaIndices(locs(ind)));
    IAF_range6 = find(powload6.freq > (IAF6-4) & powload6.freq < (IAF6+2));

    % Store the power values at the calculated IAFs
    powerIAF2 = mean(powspctrm2(IAF_range2));
    powerIAF4 = mean(powspctrm4(IAF_range4));
    powerIAF6 = mean(powspctrm6(IAF_range6));

    % Check if any IAF is 8 or 14 and set the corresponding power to NaN
    if IAF2 == 8 || IAF2 == 14
        powerIAF2 = NaN;
    end
    if IAF4 == 8 || IAF4 == 14
        powerIAF4 = NaN;
    end
    if IAF6 == 8 || IAF6 == 14
        powerIAF6 = NaN;
    end

    %% Create a structure array for this subject
    subID = str2num(subjects{subj});
    subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
        'AlphaPower', num2cell([powerIAF2; powerIAF4; powerIAF6]), 'IAF', num2cell([IAF2; IAF4; IAF6]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_sternberg_subj subj_data_eeg
    save alpha_power_sternberg powerIAF2 powerIAF4 powerIAF6
    save IAF_sternberg IAF2 IAF4 IAF6
    eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
    clc
    fprintf('Subject %s IAF: load2: %f Hz (Power: %f), load4: %f Hz (Power: %f), load6: %f Hz (Power: %f) \n', subjects{subj}, IAF2, powerIAF2, IAF4, powerIAF4, IAF6, powerIAF6);
end
save /Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg eeg_data_sternberg

%% Plot alpha power BOXPLOT
close all
figure;
set(gcf, 'Position', [0, 0, 1200, 900], 'Color', 'w');
colors = {'b', 'k', 'r'};
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
numSubjects = length(subjects);

% Collect data into a matrix for plotting
dataAlphaPower = [powerIAF2; powerIAF4; powerIAF6]';

% Plot for alpha power
hold on;
% Boxplots
boxplot(dataAlphaPower, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for subj = 1:numSubjects
    plot(1:length(conditions), dataAlphaPower(subj, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end
% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, numSubjects), dataAlphaPower(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end
xlabel('Condition', 'FontName', 'Arial', 'FontSize', 25);
ylabel('Alpha Power at IAF [\muV^2/Hz]', 'FontName', 'Arial', 'FontSize', 25);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 20, 'Location', 'best');
title('Sternberg Alpha Power at IAF', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'XLim', [0.5 3.5]);
hold off;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/boxplot/AOC_alpha_power_sternberg_boxplot.png');

%% Plot alpha power POWERSPECTRUM
close all
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
colors = {'b', 'r', 'k'};
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};
numSubjects = length(subjects);

% Load powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load power_stern
    powl2{subj} = powload2;
    powl4{subj} = powload4;
    powl6{subj} = powload6;
end
% Compute grand avg
load('/Volumes/methlab/Students/Arne/AOC/toolboxes/headmodel/layANThead.mat');
gapow2 = ft_freqgrandaverage([],powl2{:});
gapow4 = ft_freqgrandaverage([],powl4{:});
gapow6 = ft_freqgrandaverage([],powl6{:});

% Plot 
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linecolor = 'brk';
cfg.linewidth = 1;
ft_singleplotER(cfg,gapow2,gapow4,gapow6);
hold on;

% Add shadedErrorBar
addpath('/Volumes/methlab/Students/Arne/AOC/toolboxes')
channels_seb = ismember(gapow2.label, cfg.channel);
l2ebar = shadedErrorBar(gapow2.freq, mean(gapow2.powspctrm(channels_seb, :), 1), std(gapow2.powspctrm(channels_seb, :))/sqrt(size(gapow2.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
l4ebar = shadedErrorBar(gapow4.freq, mean(gapow4.powspctrm(channels_seb, :), 1), std(gapow4.powspctrm(channels_seb, :))/sqrt(size(gapow4.powspctrm(channels_seb, :), 1)), {'g', 'markerfacecolor', 'g'});
l6ebar = shadedErrorBar(gapow6.freq, mean(gapow6.powspctrm(channels_seb, :), 1), std(gapow6.powspctrm(channels_seb, :))/sqrt(size(gapow6.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.25;
set(l2ebar.patch, 'FaceAlpha', transparency);
set(l4ebar.patch, 'FaceAlpha', transparency);
set(l6ebar.patch, 'FaceAlpha', transparency);

% Adjust plotting
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow2.label);
freq_idx = find(gapow2.freq >= 8 & gapow2.freq <= 14);
max_spctrm = max([mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow4.powspctrm(channel_idx, freq_idx), 2); mean(gapow6.powspctrm(channel_idx, freq_idx), 2)]);
ylim([0 max_spctrm*1.4])
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
legend([l2ebar.mainLine, l4ebar.mainLine, l6ebar.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontName', 'Arial', 'FontSize', 20);
title('Sternberg Power Spectrum', 'FontSize', 30);
hold off;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/AOC_alpha_power_sternberg_powspctrm.png');

%% Plot alpha power TOPOS
close all;
clc;
cfg = [];
load('/Volumes/methlab/Students/Arne/AOC/toolboxes/headmodel/layANThead.mat');
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
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo2.png');

% Plot WM load 4 
figure('Color', 'w'); 
set(gcf, 'Position', [700, 300, 800, 600]); 
ft_topoplotER(cfg, gapow4); 
title('');
cb = colorbar;
set(cb, 'FontSize', 20); 
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); 
title('WM load 4', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo4.png'); 

% Plot WM load 6 
figure('Color', 'w'); 
set(gcf, 'Position', [0, 0, 800, 600]);
ft_topoplotER(cfg, gapow6); 
title('');
cb = colorbar;
set(cb, 'FontSize', 20); 
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 6', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo6.png');

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
load('/Volumes/methlab/Students/Arne/AOC/toolboxes/headmodel/layANThead.mat');
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
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo_diff.png');
