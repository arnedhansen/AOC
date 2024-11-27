%% AOC Alpha Power N-back

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
load('power_nback.mat');
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
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {});

for numSubjects = 1:length(subjects)
    datapath = strcat(path, subjects{numSubjects}, '/eeg');
    cd(datapath);
    load('power_nback.mat');
    
    % Channels selection based on CBPT    
    channelIdx = find(ismember(powload1.label, channels));
    
    % Extract power spectra for selected channels
    powspctrm1 = mean(powload1.powspctrm(channelIdx, :), 1);
    powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
    powspctrm3 = mean(powload3.powspctrm(channelIdx, :), 1);

    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));
    
    % Calculate IAF for 1-back
    alphaPower1 = powspctrm1(alphaIndices);
    [pks,locs] = findpeaks(alphaPower1);
    [~, ind] = max(pks);
    IAF1 = powload1.freq(alphaIndices(locs(ind)));
    IAF_range1 = find(powload1.freq > (IAF1-4) & powload1.freq < (IAF1+2));

    % Calculate IAF for 2-back
    alphaPower2 = powspctrm2(alphaIndices);
    [pks,locs] = findpeaks(alphaPower2);
    [~, ind] = max(pks);
    IAF2 = powload2.freq(alphaIndices(locs(ind)));
    IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));

    % Calculate IAF for 3-back
    alphaPower3 = powspctrm3(alphaIndices);
    [pks,locs] = findpeaks(alphaPower3);
    [~, ind] = max(pks);
    IAF3 = powload3.freq(alphaIndices(locs(ind)));
    IAF_range3 = find(powload3.freq > (IAF3-4) & powload3.freq < (IAF3+2));

    % Store the power values at the calculated IAFs
    powerIAF1 = mean(powspctrm1(IAF_range1));
    powerIAF2 = mean(powspctrm2(IAF_range2));
    powerIAF3 = mean(powspctrm3(IAF_range3));

    % Check if any IAF is 8 or 14 and set the corresponding power to NaN
    if IAF1 == 8 || IAF1 == 14
        powerIAF1 = NaN;
    end
    if IAF2 == 8 || IAF2 == 14
        powerIAF2 = NaN;
    end
    if IAF3 == 8 || IAF3 == 14
        powerIAF3 = NaN;
    end

    %% Create a structure array for this subject
    subID = str2num(subjects{subj});
    subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
        'AlphaPower', num2cell([powerIAF1; powerIAF2; powerIAF3]), 'IAF', num2cell([IAF1; IAF2; IAF3]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_nback_subj subj_data_eeg
    save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
    save IAF_nback IAF1 IAF2 IAF3
    eeg_data_nback = [eeg_data_nback; subj_data_eeg];
    clc
    fprintf('Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), 3-back: %f Hz (Power: %f)\n', subjects{numSubjects}, IAF1, powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3);
end
save /Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_nback eeg_data_nback

%% Plot alpha power BOXPLOT
close all

% Collect data into a matrix for plotting
dataAlphaPower = [powerIAF1; powerIAF2; powerIAF3]';

% Plot
figure;
set(gcf, 'Position', [0, 0, 1200, 900], 'Color', 'w');
colors = {'b', 'g', 'r'};
conditions = {'1-back', '2-back', '3-back'};
numSubjects = length(subjects);

% Boxplot
hold on;
boxplot(dataAlphaPower, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
for numSubjects = 1:numSubjects
    plot(1:length(conditions), dataAlphaPower(numSubjects, :), '-o', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor', 'w');
end

% Scatter plot for individual points
scatterHandles = gobjects(1, length(conditions));
for condIdx = 1:length(conditions)
    scatterHandles(condIdx) = scatter(repelem(condIdx, numSubjects), dataAlphaPower(:, condIdx), 100, colors{condIdx}, 'filled', 'MarkerEdgeColor', 'k');
end

% Adjust plot
xlabel('Conditions', 'FontName', 'Arial', 'FontSize', 25);
ylabel('Alpha Power at IAF [\muV^2/Hz]', 'FontName', 'Arial', 'FontSize', 25);
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
legend(scatterHandles, conditions, 'FontName', 'Arial', 'FontSize', 20, 'Location', 'northeast');
set(gca, 'XLim', [0.5 3.5]);
set(gca, 'YLim', [0 max(dataAlphaPower(:))*1.15]);
title('N-back Alpha Power at IAF', 'FontName', 'Arial', 'FontSize', 30);
hold off;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/boxplot/AOC_alpha_power_nback_boxplot.png');

%% Plot alpha power POWERSPECTRUM
close all
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
colors = {'b', 'g', 'r'};
conditions = {'1-back', '2-back', '3-back'};
numSubjects = length(subjects);

% Load powerspctrm data
for subj= 1:length(subjects)
    datapath = strcat(path, subjects{numSubjects}, '/eeg');
    cd(datapath)
    load power_nback
    powl1{subj} = powload1;
    powl2{subj} = powload2;
    powl3{subj} = powload3;
end
% Compute grand avg
load('/Volumes/methlab/Students/Arne/AOC/toolboxes/headmodel/layANThead.mat');
gapow1 = ft_freqgrandaverage([],powl1{:});
gapow2 = ft_freqgrandaverage([],powl2{:});
gapow3 = ft_freqgrandaverage([],powl3{:});

% Plot 
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linecolor = 'bgr';
cfg.linewidth = 1;
ft_singleplotER(cfg,gapow1,gapow2,gapow3);
hold on;

% Add shadedErrorBar
addpath('/Volumes/methlab/Students/Arne/AOC/toolboxes')
channels_seb = ismember(gapow1.label, cfg.channel);
l1ebar = shadedErrorBar(gapow1.freq, mean(gapow1.powspctrm(channels_seb, :), 1), std(gapow1.powspctrm(channels_seb, :))/sqrt(size(gapow1.powspctrm(channels_seb, :), 1)), {'b', 'markerfacecolor', 'b'});
l2ebar = shadedErrorBar(gapow2.freq, mean(gapow2.powspctrm(channels_seb, :), 1), std(gapow2.powspctrm(channels_seb, :))/sqrt(size(gapow2.powspctrm(channels_seb, :), 1)), {'g', 'markerfacecolor', 'g'});
l3ebar = shadedErrorBar(gapow3.freq, mean(gapow3.powspctrm(channels_seb, :), 1), std(gapow3.powspctrm(channels_seb, :))/sqrt(size(gapow3.powspctrm(channels_seb, :), 1)), {'r', 'markerfacecolor', 'r'});
transparency = 0.5;
set(l1ebar.patch, 'FaceAlpha', transparency);
set(l2ebar.patch, 'FaceAlpha', transparency);
set(l3ebar.patch, 'FaceAlpha', transparency);

% Adjust plotting
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow1.label);
freq_idx = find(gapow1.freq >= 8 & gapow1.freq <= 14);
max_spctrm = max([mean(gapow1.powspctrm(channel_idx, freq_idx), 2); mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow3.powspctrm(channel_idx, freq_idx), 2)]);
ylim([0 max_spctrm*1.25])
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
legend([l1ebar.mainLine, l2ebar.mainLine, l3ebar.mainLine], {'1 back', '2 back', '3 back'}, 'FontName', 'Arial', 'FontSize', 20);
title('N-back Power Spectrum', 'FontSize', 30)
hold off;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/AOC_alpha_power_nback_powspctrm.png');

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
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo1.png'); 

% Plot 2-back 
figure('Color', 'w'); 
set(gcf, 'Position', [300, 300, 800, 600]); 
ft_topoplotER(cfg, gapow2); 
title('');
cb = colorbar;
set(cb, 'FontSize', 20); 
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); 
title('2-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo2.png'); 

% Plot 3-back 
figure('Color', 'w'); 
set(gcf, 'Position', [600, 300, 800, 600]);
ft_topoplotER(cfg, gapow3); 
title('');
cb = colorbar;
set(cb, 'FontSize', 20); 
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); 
title('3-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo3.png'); 

%% Plot alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (3-back - 1-back)
ga_diff = gapow3;
ga_diff.powspctrm = gapow3.powspctrm - gapow1.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]);
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
title('N-back Alpha Power Difference (3-back - 1-back)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo_diff.png');
