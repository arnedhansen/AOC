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
for i = 1:length(pow2_early.label)
    label = pow2_early.label{i};
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
    powl2{subj} = pow2_late;
    powl4{subj} = pow4_late;
    powl6{subj} = pow6_late;
end

% Compute grand avg of raw LATE-window powspctrm data
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow4_raw = ft_freqgrandaverage([],powl4{:});
gapow6_raw = ft_freqgrandaverage([],powl6{:});

%% Plot alpha power grand average POWERSPECTRUM
close all

% Prepare your data-sets
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
        'FaceAlpha', 0.20 );
    set( ebs(k).edge(1), 'Color', 'none' );
    set( ebs(k).edge(2), 'Color', 'none' );
end

% Aesthetics
set(gca,'FontSize',20);
box off
xlim([5 20]);
%ylim([0 1])
ylabel('Power [\muV^2/Hz]');
xlabel('Frequency [Hz]');
leg_p2 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p6 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], ...
    {'WM load 2','WM load 4','WM load 6'}, ...
    'FontName','Arial','FontSize',20, 'Box', 'off');
title('Sternberg Power Spectrum', 'FontSize', 25);

% Save
saveas(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_sternberg_raw.png'));

%% Plot INDIVIDUAL power spectra
% close all
% output_dir = [fig_dir_ind filesep];
% load(fullfile(paths.features, 'AOC_eeg_matrix_sternberg.mat'))
% 
% for subj = 1:length(subjects)
%     close all;
%     figure;
%     set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
% 
%     % Extract participant data
%     pow2 = powl2{subj};
%     pow4 = powl4{subj};
%     pow6 = powl6{subj};
% 
%     % Figure common config
%     cfg = [];
%     cfg.channel = channels;
%     cfg.figure = 'gcf';
%     cfg.linewidth = 1;
% 
%     % Plot power spectrum for low and high contrast
%     ft_singleplotER(cfg, pow2, pow4, pow6);
%     hold on;
% 
%     % Add shaded error bars
%     channels_seb = ismember(pow2.label, cfg.channel);
%     eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
%         std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     eb4 = shadedErrorBar(pow4.freq, mean(pow4.powspctrm(channels_seb, :), 1), ...
%         std(pow4.powspctrm(channels_seb, :)) / sqrt(size(pow4.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     eb6 = shadedErrorBar(pow6.freq, mean(pow6.powspctrm(channels_seb, :), 1), ...
%         std(pow6.powspctrm(channels_seb, :)) / sqrt(size(pow6.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     eb2.mainLine.Color = colors(1, :);
%     eb4.mainLine.Color = colors(2, :);
%     eb6.mainLine.Color = colors(3, :);
%     eb2.patch.FaceColor = colors(1, :);
%     eb4.patch.FaceColor = colors(2, :);
%     eb6.patch.FaceColor = colors(3, :);
%     set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
%     set(eb4.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
%     set(eb6.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
%     set(eb2.edge(1), 'Color', colors(1, :));
%     set(eb2.edge(2), 'Color', colors(1, :));
%     set(eb4.edge(1), 'Color', colors(2, :));
%     set(eb4.edge(2), 'Color', colors(2, :));
%     set(eb6.edge(1), 'Color', colors(3, :));
%     set(eb6.edge(2), 'Color', colors(3, :));
%     set(eb2.patch, 'FaceAlpha', 0.5);
%     set(eb4.patch, 'FaceAlpha', 0.5);
%     set(eb6.patch, 'FaceAlpha', 0.5);
% 
%     % Add lines for IAF and AlphaPower for each condition
%     currentSubj = str2double(subjects{subj});
%     C1 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 2);
%     if ~isempty(C1)
%         xIAF = C1.IAF;
%         yAlpha = C1.AlphaPower;
%         line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
%         line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
%     end
%     C2 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 4);
%     if ~isempty(C2)
%         xIAF = C2.IAF;
%         yAlpha = C2.AlphaPower;
%         line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
%         line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
%     end
%     C3 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 6);
%     if ~isempty(C3)
%         xIAF = C3.IAF;
%         yAlpha = C3.AlphaPower;
%         line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
%         line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
%     end
% 
%     % Adjust plot aesthetics
%     set(gca, 'FontSize', 20);
%     max_spctrm = max([max(max(pow2.powspctrm(channels_seb, :))), max(max(pow4.powspctrm(channels_seb, :))), max(max(pow6.powspctrm(channels_seb, :)))]);
%     ylim([0 max_spctrm*0.75]);
%     xlim([5 30]);
%     xlabel('Frequency [Hz]');
%     ylabel('Power [a.u.]');
%     legend([eb2.mainLine, eb4.mainLine eb6.mainLine], {'WM Load 2', 'WM Load 4', 'WM Load 6'}, 'FontName', 'Arial', 'FontSize', 20);
%     title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
%     hold off;
% 
%     % Save individual plot
%     save_path = fullfile(output_dir, sprintf('AOC_powspctrm_sternberg_subj%s.png', subjects{subj}));
%     saveas(gcf, save_path);
%end

%% Plot SUBPLOT of INDIVIDUAL powerspectra
% close all;
% output_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/';
% num_subj = length(subjects);
%
% % Determine subplot grid size
% default_cols = 5;
% nrows = ceil(num_subj / default_cols);
% ncols = min(num_subj, default_cols);
%
% % Create figure
% figure;
% set(gcf, 'Color', 'w', 'Position', [0, 0, 300 * ncols, 300 * nrows]);
%
% for subj = 1:num_subj
%     % Extract participant data
%     pow2 = powl2{subj};
%     pow4 = powl4{subj};
%     pow6 = powl6{subj};
%
%     % Select subplot position
%     subplot(nrows, ncols, subj);
%     hold on;
%
%     % Figure common config
%     cfg = [];
%     cfg.channel = channels;
%     cfg.figure = 'gcf';
%     cfg.linewidth = 1;
%
%     % Plot power spectrum for low and high contrast
%     ft_singleplotER(cfg, pow2, pow4, pow6);
%
%     % Add shaded error bars
%     channels_seb = ismember(pow2.label, cfg.channel);
%     eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
%         std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     eb4 = shadedErrorBar(pow4.freq, mean(pow4.powspctrm(channels_seb, :), 1), ...
%         std(pow4.powspctrm(channels_seb, :)) / sqrt(size(pow4.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%     eb6 = shadedErrorBar(pow6.freq, mean(pow6.powspctrm(channels_seb, :), 1), ...
%         std(pow6.powspctrm(channels_seb, :)) / sqrt(size(pow6.powspctrm(channels_seb, :), 1)), {'-'}, 0);
%
%     eb2.mainLine.Color = colors(1, :);
%     eb4.mainLine.Color = colors(2, :);
%     eb6.mainLine.Color = colors(3, :);
%     eb2.patch.FaceColor = colors(1, :);
%     eb4.patch.FaceColor = colors(2, :);
%     eb6.patch.FaceColor = colors(3, :);
%     set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
%     set(eb4.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
%     set(eb6.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
%     set(eb2.patch, 'FaceAlpha', 0.5);
%     set(eb4.patch, 'FaceAlpha', 0.5);
%     set(eb6.patch, 'FaceAlpha', 0.5);
%
%     % Adjust plot aesthetics
%     set(gca, 'FontSize', 12);
%     max_spctrm = max([max(max(pow2.powspctrm(channels_seb, :))), max(max(pow4.powspctrm(channels_seb, :))), max(max(pow6.powspctrm(channels_seb, :))) ]);
%     ylim([0 max_spctrm * 0.75]);
%     xlim([5 30]);
%     xlabel('Frequency [Hz]');
%     ylabel('Power [a.u.]');
%     if subj == 1
%         legend([eb2.mainLine, eb4.mainLine eb6.mainLine], {'WM Load 2', 'WM Load 4', 'WM Load 6'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
%     end
%     title(sprintf('Subject %s', subjects{subj}), 'FontSize', 14);
%     hold off;
% end
%
% % Save combined figure
% save_path = fullfile(output_dir, 'AOC_powspctrm_sternberg_all_subs.png');
% saveas(gcf, save_path);

