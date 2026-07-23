%% AOC Power Spectrum — N-back (Hanning/mtmconvol, baselined dB, no FOOOF)
% Loads power_nback_windows.mat fields pow*_bl_full, grand-averages across 1/2/3-back,
% plots baseline-relative power (dB), and saves figures.
%
% Data source: power_nback_windows.mat from AOC_eeg_fex_nback.m (ft_freqbaseline, db).
% For raw µV^2/Hz full-window spectra, use AOC_eeg_powspctrm_nback.m (pow*_raw_full).

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end

%% Define channels (reference: first subject)
datapath = fullfile(path, subjects{1}, 'eeg');
D0 = load(fullfile(datapath, 'power_nback_windows.mat'), 'pow1_bl_full');
pow1_bl_full = D0.pow1_bl_full;

occ_channels = {};
for i = 1:length(pow1_bl_full.label)
    label = pow1_bl_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
clc
disp('LOADING BASELINED DATA...')
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    clc
    disp('LOADING BASELINED DATA...')
    disp(subj)
    D = load(fullfile(datapath, 'power_nback_windows.mat'), ...
        'pow1_bl_full', 'pow2_bl_full', 'pow3_bl_full');
    powl1{subj} = D.pow1_bl_full;
    powl2{subj} = D.pow2_bl_full;
    powl3{subj} = D.pow3_bl_full;
end

gapow1 = ft_freqgrandaverage([], powl1{:});
gapow2 = ft_freqgrandaverage([], powl2{:});
gapow3 = ft_freqgrandaverage([], powl3{:});

%% Plot grand-average power spectrum
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 5;
hold on;
yline(0, '--')

channels_seb = ismember(gapow1.label, cfg.channel);
freqs = gapow1.freq;
freq_mask = freqs >= 5 & freqs <= 20;

pow1_e = gapow1.powspctrm(channels_seb, :);
pow2_e = gapow2.powspctrm(channels_seb, :);
pow3_e = gapow3.powspctrm(channels_seb, :);

m1 = mean(pow1_e, 1, 'omitnan');
n1 = sum(isfinite(pow1_e), 1);
se1 = std(pow1_e, 0, 1, 'omitnan') ./ sqrt(n1);
m2 = mean(pow2_e, 1, 'omitnan');
n2 = sum(isfinite(pow2_e), 1);
se2 = std(pow2_e, 0, 1, 'omitnan') ./ sqrt(n2);
m3 = mean(pow3_e, 1, 'omitnan');
n3 = sum(isfinite(pow3_e), 1);
se3 = std(pow3_e, 0, 1, 'omitnan') ./ sqrt(n3);

se1(n1 < 2) = NaN;
se2(n2 < 2) = NaN;
se3(n3 < 2) = NaN;

eb1 = shadedErrorBar(freqs, m1, se1, 'lineProps', {'-','Color',colors(1,:)});
eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(2,:)});
eb3 = shadedErrorBar(freqs, m3, se3, 'lineProps', {'-','Color',colors(3,:)});
ebs = [eb1, eb2, eb3];

for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(k,:));
    set(ebs(k).patch, 'FaceColor', colors(k,:), 'FaceAlpha', 0.20);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end

plot_vals = [m1(freq_mask), m1(freq_mask) + se1(freq_mask), m1(freq_mask) - se1(freq_mask), ...
    m2(freq_mask), m2(freq_mask) + se2(freq_mask), m2(freq_mask) - se2(freq_mask), ...
    m3(freq_mask), m3(freq_mask) + se3(freq_mask), m3(freq_mask) - se3(freq_mask)];
ymax = max(abs(plot_vals), [], 'omitnan');

set(gca, 'FontSize', 20);
box off
xlim([5 20]);
ylim([-ymax * 1.15, ymax * 1.15]);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p3 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1, leg_p2, leg_p3], {' 1-back', ' 2-back', ' 3-back'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');

drawnow;
exportgraphics(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_nback_bl_full.png'), 'Resolution', 600);
