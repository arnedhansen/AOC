%% Alpha Power Time Frequency Analysis for AOC Nback data

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Compute grand average time and frequency data GATFR
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_nback
    l1{subj} = tfr1_bl;
    l2{subj} = tfr2_bl;
    l3{subj} = tfr3_bl;
    disp(['Subject ' num2str(subjects{subj}) ' TFR data loaded.'])
end

% Compute grand average
gatfr1 = ft_freqgrandaverage([],l1{:});
gatfr2 = ft_freqgrandaverage([],l2{:});
gatfr3 = ft_freqgrandaverage([],l3{:});

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
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in secon
cfg.ylim = [4 20];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr1.label);
max_spctrm = max([
    max(abs(gatfr1.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(gatfr2.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(gatfr3.powspctrm(channel_idx, :, :)), [], 'all')]);
clim = double([-max_spctrm * 0.9, max_spctrm * 0.9]);

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
title('1-back TFR', 'FontSize', 30);
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
title('2-back TFR', 'FontSize', 30);
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
title('3-back TFR', 'FontSize', 30);
set(gca, 'FontSize', 25);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_ga_3back.png');


%% Compute the difference between condition 3 and condition 1
close all
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');

% Plot the difference
diff = gatfr3;
diff.powspctrm = gatfr3.powspctrm - gatfr1.powspctrm;

% Define configuration 
cfg = [];
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in secon
cfg.ylim = [4 20];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation
[~, channel_idx] = ismember(channels, gatfr1.label);
max_spctrm = max(abs(diff.powspctrm(channel_idx, :, :)), [], 'all');
max_spctrm = 0.75
clim = double([-max_spctrm, max_spctrm]);

% Plot: Difference (Condition 3 minus Condition 1) - Time-Frequency Response
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim); 
colorbar;
line([0 0], [2 2], 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
%line(14, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('N-back TFR Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', 25);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_nback_diff.png');