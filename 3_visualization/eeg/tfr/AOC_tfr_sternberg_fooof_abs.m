%% AOC TFR — Sternberg (FOOOF, Absolute)
% Loads FOOOFed tfr_stern (tfr*_fooof) per subject, grand-averages, plots occipital TFRs. Saves figures.
%
% Key outputs:
%   TFR figures (FOOOF, absolute power, per condition)

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Compute grand average time and frequency data GATFR
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_stern
    tfr2_all{subj} = tfr2_fooof;
    tfr4_all{subj} = tfr4_fooof;
    tfr6_all{subj} = tfr6_fooof;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

% Compute grand average
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr4 = ft_freqgrandaverage([],tfr4_all{:});
gatfr6 = ft_freqgrandaverage([],tfr6_all{:});

%% Define channels
occ_channels = {};
for i = 1:length(tfr2_fooof.label)
    label = tfr2_fooof.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
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
cfg.xlim = [-.5 2];
cfg.ylim = [4 30];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64)); % Red diverging color map

% Baseline
cfg.baseline      = [-.5 -0.25];   % example baseline window
cfg.baselinetype  = 'absolute';   % options: 'absolute', 'relative', 'relchange', 'db'

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr2.label);
freq_idx = find(gatfr2.freq >= 8 & gatfr2.freq <= 14);
time_idx = find(gatfr2.time >= -1 & gatfr2.time <= 2);
max_spctrm = max([mean(gatfr2.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(gatfr4.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(gatfr6.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan')]);
max_spctrm = max(max(abs(max_spctrm)));
max_spctrm = 0.25
clim = [-max_spctrm, max_spctrm];

% WM load 2
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
%rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 2 TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_sternberg_2_fooof_bl_absX.png');

% WM load 4
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr4);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
%rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 4 TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_sternberg_4_fooof_bl_absX.png');

% WM load 6
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr6);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
%rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 6 TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_sternberg_6_fooof_bl_absX.png');

%% Plot the grand averages for the difference between condition 3 and condition 1
close all

% Plot the difference
diff = gatfr6;
diff.powspctrm = gatfr6.powspctrm - gatfr2.powspctrm;

% Define configuration 
cfg = [];
cfg.channel = channels; % specify the channels to include
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxabs'; % color limits
cfg.xlim = [-.5 2]; % Time axis limits in secon
cfg.ylim = [4 20];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation
[~, channel_idx] = ismember(channels, gatfr2.label);
time_idx = find(diff.time >= -0.5 & diff.time <= 2);
freq_idx = find(gatfr2.freq >= 8 & gatfr2.freq <= 14);
max_spctrm = max(abs(diff.powspctrm(channel_idx, freq_idx, time_idx)), [], 'all');
max_spctrm = .035
clim = double([-max_spctrm max_spctrm]);

% Plot: Difference Time-Frequency Response
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('Sternberg TFR Difference (WM load 6 minus WM load 2)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_sternberg_diff_fooof_bl_absX.png');