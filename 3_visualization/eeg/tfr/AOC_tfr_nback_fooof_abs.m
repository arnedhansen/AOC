%% AOC TFR — N-Back (FOOOF, Absolute)
% Loads FOOOFed tfr_nback (tfr*_fooof) per subject, grand-averages, plots occipital TFRs. Saves figures.
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
    load tfr_nback
    tfr1_all{subj} = tfr1_fooof;
    tfr2_all{subj} = tfr2_fooof;
    tfr3_all{subj} = tfr3_fooof;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

% Compute grand average
gatfr1 = ft_freqgrandaverage([],tfr1_all{:});
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr3 = ft_freqgrandaverage([],tfr3_all{:});

%% Define channels
occ_channels = {};
for i = 1:length(tfr1.label)
    label = tfr1_fooof.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
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
cfg.xlim = [-.5 2];
cfg.ylim = [4 30];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64)); % Red diverging color map

% Baseline
cfg.baseline      = [-.5 -0.25];   % example baseline window
cfg.baselinetype  = 'absolute';   % options: 'absolute', 'relative', 'relchange', 'db'

% Find maximum deviation across conditions
[~, channel_idx] = ismember(channels, gatfr1.label);
freq_idx = find(gatfr1.freq >= 8 & gatfr1.freq <= 14);
time_idx = find(gatfr1.time >= 0 & gatfr1.time <= 2);
max_spctrm = max([mean(gatfr1.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(gatfr2.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(gatfr3.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan')]);
max_spctrm = max(max(abs(max_spctrm)));
max_spctrm = 0.25
clim = [-max_spctrm, max_spctrm];

% 1-back
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr1);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
%rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('1-back TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_1back_fooof_bl_abs.png');

% 2-back
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
title('2-back TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_2back_fooof_bl_abs.png');

% 3-back
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr3);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
%rectangle('Position', [0, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('3-back TFR');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_3back_fooof_bl_abs.png');

%% Plot the grand averages for the difference between condition 3 and condition 1
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
cfg.ylim = [7 20];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat'); % Load layout
cfg.layout = layANThead; % your specific layout
color_map = flipud(cbrewer('div', 'RdBu', 64)); % 'RdBu' for blue to red diverging color map

% Find maximum deviation
[~, channel_idx] = ismember(channels, gatfr1.label);
time_idx = find(diff.time >= -0.5 & diff.time <= 2);
freq_idx = find(gatfr1.freq >= 8 & gatfr1.freq <= 14);
max_spctrm = max(abs(diff.powspctrm(channel_idx, freq_idx, time_idx)), [], 'all');
max_spctrm = .5
clim = double([-max_spctrm max_spctrm]);

% Plot: Difference Time-Frequency Response
figure;
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim); 
cb = colorbar;
ylabel(cb, 'Absolute Power Change [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
%rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'r', 'LineWidth', 5);
title('nback TFR Difference (WM load 6 minus WM load 2)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/AOC_tfr_nback_diff_fooof_bl_abs.png');