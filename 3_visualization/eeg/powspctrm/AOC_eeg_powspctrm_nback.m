%% AOC Power Spectrum N Back Raw
% Plot grand average raw power spectra across N back conditions in the full retention window.
% Output: `AOC_powspctrm_nback_raw_full.png`.

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end
fig_dir_ind = fullfile(paths.figures, 'eeg', 'alpha_power', 'powspctrm');

%% Define channels (reference: first subject)
datapath = fullfile(path, subjects{1}, 'eeg');
D0 = load(fullfile(datapath, 'power_nback_windows.mat'), 'pow1_raw_full');
pow1_raw_full = D0.pow1_raw_full;
% Occipital channels
occ_channels = {};
for i = 1:length(pow1_raw_full.label)
    label = pow1_raw_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
% Load power at IAF and IAF data
% for subj = 1:length(subjects)
%     datapath = fullfile(path, subjects{subj}, 'eeg');
%     cd(datapath);
%     load('alpha_power_nback.mat');
%     powIAF1(subj) = powerIAF1;
%     powIAF2(subj) = powerIAF2;
%     powIAF3(subj) = powerIAF3;
% end

% Load powspctrm data (raw full-window spectra; same MAT as baselined script, different fields)
clc
for subj = 1:length(subjects)
    clc; fprintf('[VIZ POWSPCTRM NBACK RAW] Loading spectra for Subject %d/%d (%s)\n', subj, length(subjects), subjects{subj})
    datapath = fullfile(path, subjects{subj}, 'eeg');
    D = load(fullfile(datapath, 'power_nback_windows.mat'), ...
        'pow1_raw_full', 'pow2_raw_full', 'pow3_raw_full');
    powl1{subj} = D.pow1_raw_full;
    powl2{subj} = D.pow2_raw_full;
    powl3{subj} = D.pow3_raw_full;
end

% Compute grand avg of powspctrm data
gapow1_raw = ft_freqgrandaverage([],powl1{:});
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow3_raw = ft_freqgrandaverage([],powl3{:});

%% Plot alpha power grand average POWERSPECTRUM
gapow1 = gapow1_raw;
gapow2 = gapow2_raw;
gapow3 = gapow3_raw;
yLabel = 'Power [\muV^2/Hz]';

% Create figure
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
conditions = {' 1-back', ' 2-back', ' 3-back'};

% Plot
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 5;
hold on;

% Add shadedErrorBar
channels_seb = ismember(gapow1.label, cfg.channel);
m1  = mean(gapow1.powspctrm(channels_seb, :), 1);
se1 = std(gapow1.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow1.powspctrm(channels_seb, :), 1));
m2  = mean(gapow2.powspctrm(channels_seb, :), 1);
se2 = std(gapow2.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow2.powspctrm(channels_seb, :), 1));
m3  = mean(gapow3.powspctrm(channels_seb, :), 1);
se3 = std(gapow3.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow3.powspctrm(channels_seb, :), 1));

% Light display-only smoothing (does not alter source data)
freqs = gapow1.freq;
freqs_plot = linspace(min(freqs), max(freqs), 250);
gauss_win = 20;
m1s  = smoothdata(interp1(freqs, m1,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se1s = smoothdata(interp1(freqs, se1, freqs_plot, 'pchip'), 'gaussian', gauss_win);
m2s  = smoothdata(interp1(freqs, m2,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se2s = smoothdata(interp1(freqs, se2, freqs_plot, 'pchip'), 'gaussian', gauss_win);
m3s  = smoothdata(interp1(freqs, m3,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se3s = smoothdata(interp1(freqs, se3, freqs_plot, 'pchip'), 'gaussian', gauss_win);

eb1 = shadedErrorBar(freqs_plot, m1s, se1s, 'lineProps', {'-','Color',colors(1,:)});
eb2 = shadedErrorBar(freqs_plot, m2s, se2s, 'lineProps', {'-','Color',colors(1,:)});
eb3 = shadedErrorBar(freqs_plot, m3s, se3s, 'lineProps', {'-','Color',colors(1,:)});
eb1.mainLine.Color = colors(1, :);
eb2.mainLine.Color = colors(2, :);
eb3.mainLine.Color = colors(3, :);
eb1.patch.FaceColor = colors(1, :);
eb2.patch.FaceColor = colors(2, :);
eb3.patch.FaceColor = colors(3, :);
set(eb1.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
set(eb3.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
set(eb1.edge(1), 'Color', 'none');
set(eb1.edge(2), 'Color', 'none');
set(eb2.edge(1), 'Color', 'none');
set(eb2.edge(2), 'Color', 'none');
set(eb3.edge(1), 'Color', 'none');
set(eb3.edge(2), 'Color', 'none');
set(eb1.patch, 'FaceAlpha', 0.20);
set(eb2.patch, 'FaceAlpha', 0.20);
set(eb3.patch, 'FaceAlpha', 0.20);

% Adjust plotting
set(gcf,'color','w');
set(gca,'Fontsize',20);
box off
xlim([5 20]);
ylim([0 6.1])
xlabel('Frequency [Hz]');
ylabel(yLabel);
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p3 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1, leg_p2, leg_p3], conditions, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');
hold off;

% Save
drawnow;
exportgraphics(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_nback_raw_full.png'), 'Resolution', 600);