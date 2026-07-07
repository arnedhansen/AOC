%% AOC Power Spectrum — Sternberg (Hanning/mtmconvol, baselined dB, late [1 2]s)
% Loads power_stern_windows.mat fields pow*_bl_late, grand-averages across WM load 2/4/6,
% plots baseline-relative power (dB), and saves figures.
%
% Data source: power_stern_windows.mat from AOC_eeg_fex_sternberg.m (ft_freqbaseline, db).
% For raw µV^2/Hz late-window spectra, use AOC_eeg_powspctrm_sternberg.m (pow*_raw_late).
% For FOOOF + baselined late-window spectra, use AOC_eeg_powspctrm_sternberg_fooof_bl.m
% (power_stern_fooof_TFR.mat).

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end

%% Define channels (reference: first subject)
datapath = fullfile(path, subjects{1}, 'eeg');
D0 = load(fullfile(datapath, 'power_stern_windows.mat'), 'pow2_bl_late');
pow2_bl_late = D0.pow2_bl_late;

occ_channels = {};
for i = 1:length(pow2_bl_late.label)
    label = pow2_bl_late.label{i};
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
    D = load(fullfile(datapath, 'power_stern_windows.mat'), ...
        'pow2_bl_late', 'pow4_bl_late', 'pow6_bl_late');
    powl2{subj} = D.pow2_bl_late;
    powl4{subj} = D.pow4_bl_late;
    powl6{subj} = D.pow6_bl_late;
end

gapow2 = ft_freqgrandaverage([], powl2{:});
gapow4 = ft_freqgrandaverage([], powl4{:});
gapow6 = ft_freqgrandaverage([], powl6{:});

%% Plot grand-average power spectrum
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel   = channels;
cfg.figure    = 'gcf';
cfg.linewidth = 3;
hold on;
yline(0, '--')

elecs = ismember(gapow2.label, cfg.channel);
freqs = gapow2.freq;
freq_mask = freqs >= 5 & freqs <= 20;

pow2_e = gapow2.powspctrm(elecs, :);
pow4_e = gapow4.powspctrm(elecs, :);
pow6_e = gapow6.powspctrm(elecs, :);

m2 = mean(pow2_e, 1, 'omitnan');
n2 = sum(isfinite(pow2_e), 1);
se2 = std(pow2_e, 0, 1, 'omitnan') ./ sqrt(n2);
m4 = mean(pow4_e, 1, 'omitnan');
n4 = sum(isfinite(pow4_e), 1);
se4 = std(pow4_e, 0, 1, 'omitnan') ./ sqrt(n4);
m6 = mean(pow6_e, 1, 'omitnan');
n6 = sum(isfinite(pow6_e), 1);
se6 = std(pow6_e, 0, 1, 'omitnan') ./ sqrt(n6);

se2(n2 < 2) = NaN;
se4(n4 < 2) = NaN;
se6(n6 < 2) = NaN;

eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(1,:)});
eb4 = shadedErrorBar(freqs, m4, se4, 'lineProps', {'-','Color',colors(2,:)});
eb6 = shadedErrorBar(freqs, m6, se6, 'lineProps', {'-','Color',colors(3,:)});
ebs = [eb2, eb4, eb6];

for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(k,:));
    set(ebs(k).patch, 'FaceColor', colors(k,:), 'FaceAlpha', 0.20);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end

plot_vals = [m2(freq_mask), m2(freq_mask) + se2(freq_mask), m2(freq_mask) - se2(freq_mask), ...
    m4(freq_mask), m4(freq_mask) + se4(freq_mask), m4(freq_mask) - se4(freq_mask), ...
    m6(freq_mask), m6(freq_mask) + se6(freq_mask), m6(freq_mask) - se6(freq_mask)];
ymax = max(abs(plot_vals), [], 'omitnan');

set(gca, 'FontSize', 20);
box off
xlim([5 20]);
ylim([-ymax * 1.15, ymax * 1.15]);
ylabel('Power [dB]');
xlabel('Frequency [Hz]');
leg_p2 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p6 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], {' WM load 2', ' WM load 4', ' WM load 6'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');

drawnow;
exportgraphics(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_sternberg_bl_late.png'), 'Resolution', 600);
