%% Alpha Power Time Frequency Analysis for AOC Sternberg data

%% Setup
setup('AOC');

%% Compute grand average time and frequency data GATFR
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_stern
    tfr2_all{subj} = tfr2;
    tfr4_all{subj} = tfr4;
    tfr6_all{subj} = tfr6;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR loaded.'])
end

% Compute grand average
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr4 = ft_freqgrandaverage([],tfr4_all{:});
gatfr6 = ft_freqgrandaverage([],tfr6_all{:});

%% Apply percentage change for baseline correction
baseline_period = [-Inf -0.5];  
gatfr2 = baseline_corr_powspctrm(gatfr2, baseline_period);
gatfr4 = baseline_corr_powspctrm(gatfr4, baseline_period);
gatfr6 = baseline_corr_powspctrm(gatfr6, baseline_period);

%% Define channels
load('tfr_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(tfr2.label)
    label = tfr2.label{i};
    if contains(label, {'O'}) 
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for each individual condition
close all

% Common configuration
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
cfg.showlabels = 'yes'; % show channel labels
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = ([-0.5 2]);
cfg.ylim = [5 20];
% clim = ([-1.25 1.25]);
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Wm load 2
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
% set(gca, 'CLim', clim);
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Sternberg WM load 2 TFR');
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_sternberg_2.png');

% WM load 4
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr4);
colormap(color_map);
% set(gca, 'CLim', clim); 
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Sternberg WM load 4 TFR');
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_sternberg_4.png');

% WM load 6
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr6);
colormap(color_map);
% set(gca, 'CLim', clim); 
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Sternberg WM load 6 TFR');
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_sternberg_6.png');

%% Plot the grand averages for the difference between condition 3 and condition 1
close all
% Plot the difference
diff = gatfr6;
diff.powspctrm = gatfr6.powspctrm - gatfr2.powspctrm;

% Define configuration for multiplot
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead; % your specific layout
cfg.channel = channels; % specify the channels to include
cfg.showlabels = 'yes'; % show channel labels
cfg.colorbar = 'yes'; % include color bar
cfg.xlim = [1 2]; 
cfg.ylim = [5 20];

% Plot: Difference Time-Frequency Response
figure;
ft_singleplotTFR(cfg, diff);
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
xlabel('Time [s]');
ylabel('Frequency [Hz]');
colorbar;
colormap(flipud(cbrewer('div', 'RdBu', 64))); 
% set(gca, 'CLim', [-50, 50]);
set(gca, 'FontSize', 25);
title('Sternberg Time-Frequency Response Difference (WM load 6 minus WM load 2)', 'FontName', 'Arial', 'FontSize', 30);

% Save 
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_sternberg_diff.png');