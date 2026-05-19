%% AOC TFR — Sternberg (Baselined)
% Loads raw tfr_stern per subject, applies ft_freqbaseline (db), 
% grand-averages, plots occipital TFRs.
%
% Data source: tfr_stern.mat from AOC_eeg_fex_sternberg_TFR.m
%
% Key outputs:
%   TFR figures (dB baseline, per condition)

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('AOC');
path = paths.features;
figpath = fullfile(paths.figures, 'eeg', 'tfr');

%% Compute grand average time and frequency data GATFR
baseline_window = [-.5 -0.25];
cfgb = [];
cfgb.baseline = baseline_window;
cfgb.baselinetype = 'db';

for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath)
    load tfr_stern
    tfr2_all{subj} = ft_freqbaseline(cfgb, tfr2);
    tfr4_all{subj} = ft_freqbaseline(cfgb, tfr4);
    tfr6_all{subj} = ft_freqbaseline(cfgb, tfr6);
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

gatfr2 = ft_freqgrandaverage([], tfr2_all{:});
gatfr4 = ft_freqgrandaverage([], tfr4_all{:});
gatfr6 = ft_freqgrandaverage([], tfr6_all{:});

%% Define channels
occ_channels = {};
for i = 1:length(gatfr2.label)
    label = gatfr2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Plot TFR for each individual condition
close all
fontSize = 40;
cfg = [];
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;
color_map = interp1(linspace(0,1,5), [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], linspace(0,1,64));

% Get max power for clim
[~, channel_idx] = ismember(channels, gatfr2.label);
freq_idx = gatfr2.freq >= 8 & gatfr2.freq <= 14;
time_idx = gatfr2.time >= 0 & gatfr2.time <= 2;
avg2 = squeeze(mean(gatfr2.powspctrm(channel_idx, :, :), 1));
avg4 = squeeze(mean(gatfr4.powspctrm(channel_idx, :, :), 1));
avg6 = squeeze(mean(gatfr6.powspctrm(channel_idx, :, :), 1));
max_spctrm = max([ ...
    max(abs(avg2(freq_idx, time_idx)), [], 'all'), ...
    max(abs(avg4(freq_idx, time_idx)), [], 'all'), ...
    max(abs(avg6(freq_idx, time_idx)), [], 'all')]);
clim = [-max_spctrm, max_spctrm];

% WM load 2
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 2');
drawnow; pause(0.05);
exportgraphics(gcf, fullfile(figpath, 'AOC_tfr_sternberg_2_bl.png'), 'Resolution', 600);

% WM load 4
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr4);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 4');
drawnow; pause(0.05);
exportgraphics(gcf, fullfile(figpath, 'AOC_tfr_sternberg_4_bl.png'), 'Resolution', 600);

% WM load 6
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr6);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 6');
drawnow; pause(0.05);
exportgraphics(gcf, fullfile(figpath, 'AOC_tfr_sternberg_6_bl.png'), 'Resolution', 600);

%% Plot difference (WM load 6 minus WM load 2)
close all

diff = gatfr6;
diff.powspctrm = gatfr6.powspctrm - gatfr2.powspctrm;

cfg = [];
cfg.channel = channels;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [5 30];
cfg.layout = headmodel.layANThead;
color_map = interp1(linspace(0,1,5), [0.02 0.19 0.58; 0.40 0.67 0.87; 0.97 0.97 0.97; 0.94 0.50 0.36; 0.40 0 0.05], linspace(0,1,64));

[~, channel_idx] = ismember(channels, diff.label);
freq_idx = diff.freq >= 5 & diff.freq <= 30;
time_idx = diff.time >= -0.5 & diff.time <= 2;
max_spctrm = max(abs(mean(diff.powspctrm(channel_idx, freq_idx, time_idx), 1)), [], 'all');
clim = double([-max_spctrm max_spctrm]);

figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('Difference (WM load 6 - WM load 2)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

drawnow; pause(0.05);
exportgraphics(gcf, fullfile(figpath, 'AOC_tfr_sternberg_diff_bl.png'), 'Resolution', 600);
