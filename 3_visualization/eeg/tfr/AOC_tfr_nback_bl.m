%% AOC TFR — N-Back (Baselined, Raw)
% Loads raw tfr_nback (tfr1/2/3, not FOOOF) per subject, applies FieldTrip
% ft_freqbaseline (db), grand-averages, plots occipital TFRs. Saves figures.
%
% Data source: tfr_nback.mat from AOC_eeg_fex_nback_TFR.m
%
% Key outputs:
%   TFR figures (dB baseline, per condition)

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('AOC');
path = paths.features;
figpath = fullfile(paths.figures, 'eeg', 'tfr');
baseline_window = [-.5 -0.25];
if ~isfolder(figpath)
    mkdir(figpath);
end

cfgb = [];
cfgb.baseline = baseline_window;
cfgb.baselinetype = 'db';

%% Compute grand average time and frequency data GATFR
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath)
    load tfr_nback
    tfr1_all{subj} = ft_freqbaseline(cfgb, tfr1);
    tfr2_all{subj} = ft_freqbaseline(cfgb, tfr2);
    tfr3_all{subj} = ft_freqbaseline(cfgb, tfr3);
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

gatfr1 = ft_freqgrandaverage([], tfr1_all{:});
gatfr2 = ft_freqgrandaverage([], tfr2_all{:});
gatfr3 = ft_freqgrandaverage([], tfr3_all{:});

%% Define channels
occ_channels = {};
for i = 1:length(gatfr1.label)
    label = gatfr1.label{i};
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

[~, channel_idx] = ismember(channels, gatfr1.label);
freq_idx = gatfr1.freq >= 5 & gatfr1.freq <= 30;
time_idx = gatfr1.time >= -0.5 & gatfr1.time <= 2;
avg1 = squeeze(mean(gatfr1.powspctrm(channel_idx, :, :), 1));
avg2 = squeeze(mean(gatfr2.powspctrm(channel_idx, :, :), 1));
avg3 = squeeze(mean(gatfr3.powspctrm(channel_idx, :, :), 1));
max_spctrm = max([ ...
    max(abs(avg1(freq_idx, time_idx)), [], 'all'), ...
    max(abs(avg2(freq_idx, time_idx)), [], 'all'), ...
    max(abs(avg3(freq_idx, time_idx)), [], 'all')]);
clim = [-max_spctrm, max_spctrm];

% 1-back
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr1);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('1-back');
drawnow;
print(gcf, [figpath 'AOC_tfr_1back_bl.png'], '-dpng', '-r300');

% 2-back
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('2-back');
drawnow;
print(gcf, [figpath 'AOC_tfr_2back_bl.png'], '-dpng', '-r300');

% 3-back
figure('Position', [0, 0, 1512, 982], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
ft_singleplotTFR(cfg, gatfr3);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('3-back');
drawnow;
print(gcf, [figpath 'AOC_tfr_3back_bl.png'], '-dpng', '-r300');

%% Plot difference (3-back minus 1-back)
close all

cfgm = [];
cfgm.operation = 'subtract';
diff = ft_freqmath(cfgm, gatfr3, gatfr1);

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
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('Difference (3-back - 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);

drawnow;
print(gcf, [figpath 'AOC_tfr_nback_diff_bl.png'], '-dpng', '-r300');
