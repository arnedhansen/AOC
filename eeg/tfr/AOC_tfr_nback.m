%% Alpha Power Time Frequency Analysis for AOC Nback data
clear
clc
close all
run startup

path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Compute grand average time and frequency data GATFR
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_nback
    l1{subj} = load1;
    l2{subj} = load2;
    l3{subj} = load3;
    disp(['Subject ' num2str(subj) '/10 TRF loaded.'])
end

% Compute grand average
gatfr1 = ft_freqgrandaverage([],l1{:});
gatfr2 = ft_freqgrandaverage([],l2{:});
gatfr3 = ft_freqgrandaverage([],l3{:});

%% Apply percentage change for baseline correction
baseline_period = [-Inf -0.5];  
gatfr1 = baseline_corr_powspctrm(gatfr1, baseline_period);
gatfr2 = baseline_corr_powspctrm(gatfr2, baseline_period);
gatfr3 = baseline_corr_powspctrm(gatfr3, baseline_period);

%% Define channels
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, 'O')
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for each individual condition
close all

% Define the common configuration
cfg = [];
% cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in secon
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map
% Find maximum deviation across conditions
% [~, channel_idx] = ismember(channels, gatfr1.label);
% max_spctrm = max([
%     max(abs(gatfr1.powspctrm(channel_idx, :, :)), [], 'all'), ...
%     max(abs(gatfr2.powspctrm(channel_idx, :, :)), [], 'all'), ...
%     max(abs(gatfr3.powspctrm(channel_idx, :, :)), [], 'all')
% ]);
% clim = [-max_spctrm, max_spctrm]
clim = [-5.5, 5.5];

% 1-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr1);
colormap(color_map);
set(gca, 'CLim', clim); 
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', 25);
title('1-back - Time-Frequency Response', 'FontSize', 30);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_1back.png');

% 2-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
set(gca, 'CLim', clim); 
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('2-back - Time-Frequency Response', 'FontSize', 30);
set(gca, 'FontSize', 25);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_2back.png');

% 3-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr3);
colormap(color_map);
set(gca, 'CLim', clim); 
colorbar;
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('3-back - Time-Frequency Response', 'FontSize', 30);
set(gca, 'FontSize', 25);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_3back.png');


%% Compute the difference between condition 3 and condition 1
close all
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');

% Plot the difference
diff = gatfr3;
diff.powspctrm = gatfr3.powspctrm - gatfr1.powspctrm;

% Define configuration for multiplot
cfg = [];
cfg.layout = layANThead; % your specific layout
cfg.channel = channels; % specify the channels to include
% cfg.baseline = [-Inf -0.5]; % baseline correction window in seconds
% cfg.baselinetype = 'absolute'; % type of baseline correction
cfg.colorbar = 'yes'; % include color bar
cfg.xlim = [0 2]; % Time axis limits in secon
cfg.ylim = [5 20];

% Plot: Difference (Condition 3 minus Condition 1) - Time-Frequency Response
figure;
ft_singleplotTFR(cfg, diff);
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
xlabel('Time [s]', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Frequency [Hz]', 'FontName', 'Arial', 'FontSize', 20);
colorbar;
colormap(flipud(cbrewer('div', 'RdBu', 64))); % 'RdBu' for blue to red diverging color map
caxis([-max(abs(diff.powspctrm(:))), max(abs(diff.powspctrm(:)))]); % Symmetric limits
set(gca, 'XLim', [0 2]); % Set time limits directly on axes
set(gca, 'CLim', [-3, 3]) % Set color limits directly on axes
set(gca, 'FontSize', 25);
title('Time-Frequency Response Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_nback_diff.png');