%% Alpha Power Time Frequency Analysis for AOC Nback data

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Compute grand average time and frequency data GATFR
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_nback
    l1{subj} = tfr1%%%%%_bl;
    l2{subj} = tfr2%%%%%_bl;
    l3{subj} = tfr3%%%%%_bl;
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
    if contains(label, 'O') || contains(label, 'I')
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for each individual condition
close all
fontSize = 40;

% Define the common configuration
cfg = [];
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
%cfg.xlim = [-.5 2]; % Time axis limits in seconds
cfg.xlim = [0 2]
cfg.ylim = [7 20];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat');
cfg.layout = layANThead; % your specific layout
%color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map
color_map = cbrewer('seq', 'Reds', 64); % 'RdBu' for blue to red diverging color map

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr1.label);
max_spctrm = max([
    max(abs(gatfr1.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(gatfr2.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(gatfr3.powspctrm(channel_idx, :, :)), [], 'all')]);
max_spctrm = 10
%clim = [-max_spctrm, max_spctrm];
clim = [0 max_spctrm];
 
% 1-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr1);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
%ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('1-back TFR', 'FontSize', 30);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_1back.png');

% 2-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
%ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('2-back TFR', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_2back.png');

% 3-back
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, gatfr3);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
%ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('3-back TFR', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_3back.png');

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
cfg.ylim = [7 20];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation
[~, channel_idx] = ismember(channels, gatfr1.label);
max_spctrm = max(abs(diff.powspctrm(channel_idx, :, :)), [], 'all');
max_spctrm = 3
clim = double([-max_spctrm, max_spctrm]);

% Plot: Difference (Condition 3 minus Condition 1) - Time-Frequency Response
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w'); 
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
%ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
line([0 0], [2 2], 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
%line(14, 'LineWidth', 5, 'LineStyle', '-', 'Color', 'r');
xlabel('Time [s]');
ylabel('Frequency [Hz]');
ylim([6 20]);
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('N-back TFR Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_nback_diff.png');