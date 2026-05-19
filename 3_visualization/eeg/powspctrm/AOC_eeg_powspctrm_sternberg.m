%% AOC Power Spectrum — Sternberg
% Loads power_stern (powspctrm), grand-averages across load 2/4/6. Plots power spectra. Saves figures.
%
% Key outputs:
%   Power spectrum figures (per condition)

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
fig_dir_ind = fullfile(paths.figures, 'eeg', 'alpha_power', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end
if ~isfolder(fig_dir_ind), mkdir(fig_dir_ind); end

%% Define channels
subj = 1;
datapath = fullfile(path, subjects{subj}, 'eeg');
cd(datapath);
load('power_stern_windows.mat');
% Occipital channels
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
disp('LOADING DATA...')

% Load raw LATE-window powerspctrm data
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath)
    clc
    disp('LOADING DATA...')
    disp(subj)
    load('power_stern_windows.mat')
    powl2{subj} = pow2_bl_late;
    powl4{subj} = pow4_bl_late;
    powl6{subj} = pow6_bl_late;
end

% Compute grand avg of raw LATE-window powspctrm data
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow4_raw = ft_freqgrandaverage([],powl4{:});
gapow6_raw = ft_freqgrandaverage([],powl6{:});

%% Plot alpha power grand average POWERSPECTRUM
close all
gapow2 = gapow2_raw;
gapow4 = gapow4_raw;
gapow6 = gapow6_raw;

% Configs
cfg = [];
set(gcf, 'Position', [0, 0, 1512*0.4, 982], 'Color', 'w');
cfg.channel   = channels;
cfg.figure    = 'gcf';
cfg.linewidth = 3;

% Plot
hold on;

% Shaded Error Bars
elecs = ismember(gapow2.label, cfg.channel);
freqs  = gapow2.freq;
m2     = mean(gapow2.powspctrm(elecs,:),1);
se2    = std (gapow2.powspctrm(elecs,:),0,1) / sqrt(sum(elecs));
m4     = mean(gapow4.powspctrm(elecs,:),1);
se4    = std (gapow4.powspctrm(elecs,:),0,1) / sqrt(sum(elecs));
m6     = mean(gapow6.powspctrm(elecs,:),1);
se6    = std (gapow6.powspctrm(elecs,:),0,1) / sqrt(sum(elecs));

% Light display-only smoothing (does not alter source data)
freqs_plot = linspace(min(freqs), max(freqs), 250);
gauss_win = 20;
m2s  = smoothdata(interp1(freqs, m2,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se2s = smoothdata(interp1(freqs, se2, freqs_plot, 'pchip'), 'gaussian', gauss_win);
m4s  = smoothdata(interp1(freqs, m4,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se4s = smoothdata(interp1(freqs, se4, freqs_plot, 'pchip'), 'gaussian', gauss_win);
m6s  = smoothdata(interp1(freqs, m6,  freqs_plot, 'pchip'), 'gaussian', gauss_win);
se6s = smoothdata(interp1(freqs, se6, freqs_plot, 'pchip'), 'gaussian', gauss_win);

eb2 = shadedErrorBar(freqs_plot, m2s, se2s, 'lineProps', {'-','Color',colors(1,:)});
eb4 = shadedErrorBar(freqs_plot, m4s, se4s, 'lineProps', {'-','Color',colors(2,:)});
eb6 = shadedErrorBar(freqs_plot, m6s, se6s, 'lineProps', {'-','Color',colors(3,:)});

ebs = [eb2, eb4, eb6];

% loop and set both line and patch properties
for k = 1:numel(ebs)
    set( ebs(k).mainLine, ...
        'LineWidth', cfg.linewidth, ...
        'Color',     colors(k,:) );
    set( ebs(k).patch,    ...
        'FaceColor', colors(k,:), ...
        'FaceAlpha', 0.25 );
    set( ebs(k).edge(1), 'Color', 'none' );
    set( ebs(k).edge(2), 'Color', 'none' );
end

% Aesthetics
set(gca,'FontSize',20);
box off
xlim([5 20]);
ylim([0 6.1])
ylabel('Power [\muV^2/Hz]');
xlabel('Frequency [Hz]');
leg_p2 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p6 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], ...
    {'WM load 2','WM load 4','WM load 6'}, ...
    'FontName','Arial','FontSize', 20, 'Box', 'off');
title('');

% Save
drawnow;
exportgraphics(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_sternberg_bl_late.png'), 'Resolution', 600);