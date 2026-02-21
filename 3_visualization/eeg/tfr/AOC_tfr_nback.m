%% AOC TFR — N-Back (Raw)
% Loads tfr_nback per subject, grand-averages TFR across conditions (0/1/2-back), plots occipital TFRs. Saves figures.
%
% Key outputs:
%   TFR figures (occipital, per condition)

%% Setup
startup
[subjects, path, ~, headmodel] = setup('AOC');
if ispc
    figpath = 'W:\Students\Arne\AOC\figures\eeg\tfr\';
else
    figpath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/';
end

%% Compute grand average time and frequency data GATFR
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_nback
    l1{subj} = tfr1;
    l2{subj} = tfr2;
    l3{subj} = tfr3;
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
fontSize = 50;

% Define the common configuration
cfg = [];
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in seconds
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;
color_map = interp1(linspace(0,1,5), [1 0.96 0.94; 0.99 0.73 0.63; 0.98 0.42 0.29; 0.80 0.09 0.11; 0.40 0 0.05], linspace(0,1,64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr1.label);
freq_idx = gatfr1.freq >= 5 & gatfr1.freq <= 30;
time_idx = gatfr1.time >= -0.5 & gatfr1.time <= 2;
max_spctrm = max([ ...
    max(mean(gatfr1.powspctrm(channel_idx, freq_idx, time_idx), 1), [], 'all'), ...
    max(mean(gatfr2.powspctrm(channel_idx, freq_idx, time_idx), 1), [], 'all'), ...
    max(mean(gatfr3.powspctrm(channel_idx, freq_idx, time_idx), 1), [], 'all')]);
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
title('1-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_1back.png']);

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
title('2-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_2back.png']);

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
title('3-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_3back.png']);

%% Compute the difference between condition 3 and condition 1
close all

% Plot the difference
diff = gatfr3;
diff.powspctrm = gatfr3.powspctrm - gatfr1.powspctrm;

% Define configuration 
cfg = [];
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in secon
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;
color_map = interp1(linspace(0,1,5), [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], linspace(0,1,64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation
[~, channel_idx] = ismember(channels, diff.label);
freq_idx = diff.freq >= 5 & diff.freq <= 30;
time_idx = diff.time >= -0.5 & diff.time <= 2;
max_spctrm = max(abs(mean(diff.powspctrm(channel_idx, freq_idx, time_idx), 1)), [], 'all');
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
ylim([5 30]);
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('N-back TFR Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

% Save the figure
saveas(gcf, [figpath 'AOC_tfr_nback_diff.png']);