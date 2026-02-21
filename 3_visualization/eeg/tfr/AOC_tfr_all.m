%% AOC TFR — All Tasks (Raw + FOOOF Baselined)
% Single script that generates all TFR plots:
%   1. Sternberg Raw (WM load 2/4/6 + diff)
%   2. Sternberg FOOOF + absolute baseline (WM load 2/4/6 + diff)
%   3. N-back Raw (1/2/3-back + diff)
%   4. N-back FOOOF + absolute baseline (1/2/3-back + diff)

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat');
figpath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/tfr/';

%% ========================================================================
%  STERNBERG — Load data
%  ========================================================================
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_stern
    stern_raw2{subj} = tfr2;
    stern_raw4{subj} = tfr4;
    stern_raw6{subj} = tfr6;
    stern_fooof2{subj} = tfr2_fooof;
    stern_fooof4{subj} = tfr4_fooof;
    stern_fooof6{subj} = tfr6_fooof;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' Sternberg TFR loaded.'])
end

% Grand averages — Raw
ga_s_raw2 = ft_freqgrandaverage([], stern_raw2{:});
ga_s_raw4 = ft_freqgrandaverage([], stern_raw4{:});
ga_s_raw6 = ft_freqgrandaverage([], stern_raw6{:});

% Grand averages — FOOOF
ga_s_f2 = ft_freqgrandaverage([], stern_fooof2{:});
ga_s_f4 = ft_freqgrandaverage([], stern_fooof4{:});
ga_s_f6 = ft_freqgrandaverage([], stern_fooof6{:});

% Occipital channels (Sternberg)
occ_channels = {};
for i = 1:length(tfr2.label)
    label = tfr2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
stern_ch = occ_channels;

%% ========================================================================
%  STERNBERG RAW — Condition TFRs
%  ========================================================================
close all
fontSize = 50;

cfg = [];
cfg.channel   = stern_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = cbrewer('seq', 'Reds', 64);

[~, channel_idx] = ismember(stern_ch, ga_s_raw2.label);
freq_idx = find(ga_s_raw2.freq >= 8 & ga_s_raw2.freq <= 14);
max_spctrm = max([mean(ga_s_raw2.powspctrm(channel_idx, freq_idx), 'omitnan'); ...
                  mean(ga_s_raw4.powspctrm(channel_idx, freq_idx), 'omitnan'); ...
                  mean(ga_s_raw6.powspctrm(channel_idx, freq_idx), 'omitnan')]);
max_spctrm = max(abs(max_spctrm));
max_spctrm = 6.75;
clim = [0 max_spctrm];

% WM load 2
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_raw2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 2');
saveas(gcf, [figpath 'AOC_tfr_sternberg_2.png']);

% WM load 4
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_raw4);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 4');
saveas(gcf, [figpath 'AOC_tfr_sternberg_4.png']);

% WM load 6
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_raw6);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'k', 'LineWidth', 5);
set(gca, 'FontSize', fontSize);
title('WM load 6');
saveas(gcf, [figpath 'AOC_tfr_sternberg_6.png']);

%% Sternberg Raw — Difference (WM load 6 minus WM load 2)
close all

diff = ga_s_raw6;
diff.powspctrm = ga_s_raw6.powspctrm - ga_s_raw2.powspctrm;

cfg = [];
cfg.channel   = stern_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(stern_ch, ga_s_raw2.label);
time_idx = find(diff.time >= -0.5 & diff.time <= 2);
freq_idx = find(ga_s_raw2.freq >= 8 & ga_s_raw2.freq <= 14);
max_spctrm = max(abs(diff.powspctrm(channel_idx, freq_idx, time_idx)), [], 'all');
max_spctrm = .5;
clim = double([-max_spctrm max_spctrm]);

figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [1, 8, 1, 6], 'EdgeColor', 'r', 'LineWidth', 5);
title('Sternberg TFR Difference (WM load 6 minus WM load 2)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_sternberg_diff.png']);

%% ========================================================================
%  STERNBERG FOOOF — Condition TFRs (absolute baseline)
%  ========================================================================
close all
fontSize = 40;

cfg = [];
cfg.channel      = stern_ch;
cfg.colorbar     = 'yes';
cfg.zlim         = 'maxabs';
cfg.xlim         = [-.5 2];
cfg.ylim         = [0 30];
cfg.layout       = layANThead;
cfg.baseline     = [-.5 -0.25];
cfg.baselinetype = 'absolute';
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(stern_ch, ga_s_f2.label);
freq_idx = find(ga_s_f2.freq >= 8 & ga_s_f2.freq <= 14);
time_idx = find(ga_s_f2.time >= -1 & ga_s_f2.time <= 2);
max_spctrm = max([mean(ga_s_f2.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(ga_s_f4.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(ga_s_f6.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan')]);
max_spctrm = max(max(abs(max_spctrm)));
max_spctrm = 0.25;
clim = [-max_spctrm, max_spctrm];

% WM load 2
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_f2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 2 TFR');
saveas(gcf, [figpath 'AOC_tfr_sternberg_2_fooof_bl_abs.png']);

% WM load 4
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_f4);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 4 TFR');
saveas(gcf, [figpath 'AOC_tfr_sternberg_4_fooof_bl_abs.png']);

% WM load 6
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_s_f6);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('Sternberg WM load 6 TFR');
saveas(gcf, [figpath 'AOC_tfr_sternberg_6_fooof_bl_abs.png']);

%% Sternberg FOOOF — Difference (WM load 6 minus WM load 2)
close all

diff = ga_s_f6;
diff.powspctrm = ga_s_f6.powspctrm - ga_s_f2.powspctrm;

cfg = [];
cfg.channel   = stern_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(stern_ch, ga_s_f2.label);
time_idx = find(diff.time >= -0.5 & diff.time <= 2);
freq_idx = find(ga_s_f2.freq >= 8 & ga_s_f2.freq <= 14);
max_spctrm = max(abs(diff.powspctrm(channel_idx, freq_idx, time_idx)), [], 'all');
max_spctrm = .035;
clim = double([-max_spctrm max_spctrm]);

figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
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
saveas(gcf, [figpath 'AOC_tfr_sternberg_diff_fooof_bl_abs.png']);

%% ========================================================================
%  N-BACK — Load data
%  ========================================================================
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_nback
    nb_raw1{subj} = tfr1;
    nb_raw2{subj} = tfr2;
    nb_raw3{subj} = tfr3;
    nb_fooof1{subj} = tfr1_fooof;
    nb_fooof2{subj} = tfr2_fooof;
    nb_fooof3{subj} = tfr3_fooof;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' N-back TFR loaded.'])
end

% Grand averages — Raw
ga_n_raw1 = ft_freqgrandaverage([], nb_raw1{:});
ga_n_raw2 = ft_freqgrandaverage([], nb_raw2{:});
ga_n_raw3 = ft_freqgrandaverage([], nb_raw3{:});

% Grand averages — FOOOF
ga_n_f1 = ft_freqgrandaverage([], nb_fooof1{:});
ga_n_f2 = ft_freqgrandaverage([], nb_fooof2{:});
ga_n_f3 = ft_freqgrandaverage([], nb_fooof3{:});

% Occipital channels (N-back)
occ_channels = {};
for i = 1:length(tfr1.label)
    label = tfr1.label{i};
    if contains(label, 'O') || contains(label, 'I')
        occ_channels{end+1} = label;
    end
end
nb_ch = occ_channels;

%% ========================================================================
%  N-BACK RAW — Condition TFRs
%  ========================================================================
close all
fontSize = 50;

cfg = [];
cfg.channel   = nb_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = cbrewer('seq', 'Reds', 64);

[~, channel_idx] = ismember(nb_ch, ga_n_raw1.label);
max_spctrm = max([
    max(abs(ga_n_raw1.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(ga_n_raw2.powspctrm(channel_idx, :, :)), [], 'all'), ...
    max(abs(ga_n_raw3.powspctrm(channel_idx, :, :)), [], 'all')]);
max_spctrm = 9.5;
clim = [0 max_spctrm];

% 1-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_raw1);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('1-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_1back.png']);

% 2-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_raw2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('2-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_2back.png']);

% 3-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_raw3);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('3-back', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_3back.png']);

%% N-back Raw — Difference (3-back minus 1-back)
close all

diff = ga_n_raw3;
diff.powspctrm = ga_n_raw3.powspctrm - ga_n_raw1.powspctrm;

cfg = [];
cfg.channel   = nb_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(nb_ch, ga_n_raw1.label);
max_spctrm = max(abs(diff.powspctrm(channel_idx, :, :)), [], 'all');
max_spctrm = 2;
clim = double([-max_spctrm, max_spctrm]);

figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 5);
title('N-back TFR Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_nback_diff.png']);

%% ========================================================================
%  N-BACK FOOOF — Condition TFRs (absolute baseline)
%  ========================================================================
close all
fontSize = 40;

cfg = [];
cfg.channel      = nb_ch;
cfg.colorbar     = 'yes';
cfg.zlim         = 'maxabs';
cfg.xlim         = [-.5 2];
cfg.ylim         = [0 30];
cfg.layout       = layANThead;
cfg.baseline     = [-.5 -0.25];
cfg.baselinetype = 'absolute';
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(nb_ch, ga_n_f1.label);
freq_idx = find(ga_n_f1.freq >= 8 & ga_n_f1.freq <= 14);
time_idx = find(ga_n_f1.time >= 0 & ga_n_f1.time <= 2);
max_spctrm = max([mean(ga_n_f1.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(ga_n_f2.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan'); ...
                  mean(ga_n_f3.powspctrm(channel_idx, freq_idx, time_idx), 'omitnan')]);
max_spctrm = max(max(abs(max_spctrm)));
max_spctrm = 0.25;
clim = [-max_spctrm, max_spctrm];

% 1-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_f1);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('1-back TFR');
saveas(gcf, [figpath 'AOC_tfr_1back_fooof_bl_abs.png']);

% 2-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_f2);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('2-back TFR');
saveas(gcf, [figpath 'AOC_tfr_2back_fooof_bl_abs.png']);

% 3-back
figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, ga_n_f3);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'FontSize', fontSize);
title('3-back TFR');
saveas(gcf, [figpath 'AOC_tfr_3back_fooof_bl_abs.png']);

%% N-back FOOOF — Difference (3-back minus 1-back)
close all

diff = ga_n_f3;
diff.powspctrm = ga_n_f3.powspctrm - ga_n_f1.powspctrm;

cfg = [];
cfg.channel   = nb_ch;
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.xlim      = [-.5 2];
cfg.ylim      = [0 30];
cfg.layout    = layANThead;
color_map = flipud(cbrewer('div', 'RdBu', 64));

[~, channel_idx] = ismember(nb_ch, ga_n_f1.label);
time_idx = find(diff.time >= -0.5 & diff.time <= 2);
freq_idx = find(ga_n_f1.freq >= 8 & ga_n_f1.freq <= 14);
max_spctrm = max(abs(diff.powspctrm(channel_idx, freq_idx, time_idx)), [], 'all');
max_spctrm = .5;
clim = double([-max_spctrm max_spctrm]);

figure('Position', [0 0 1512 982]);
set(gcf, 'Color', 'w');
ft_singleplotTFR(cfg, diff);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar;
ylabel(cb, 'Power [dB]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('N-back TFR Difference (3-back minus 1-back)', 'FontName', 'Arial', 'FontSize', 30);
set(gca, 'FontSize', fontSize);
saveas(gcf, [figpath 'AOC_tfr_nback_diff_fooof_bl_abs.png']);
