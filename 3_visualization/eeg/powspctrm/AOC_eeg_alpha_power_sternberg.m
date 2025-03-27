%% AOC Alpha Power Sternberg

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
% Load power at IAF
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath);
    load('alpha_power_sternberg.mat');
    powIAF2(subj) = powerIAF2;
    powIAF4(subj) = powerIAF4;
    powIAF6(subj) = powerIAF6;
end

% Load powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern
    powl2{subj} = powload2;
    powl4{subj} = powload4;
    powl6{subj} = powload6;
end

% Compute grand avg of raw powspctrm data
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow4_raw = ft_freqgrandaverage([],powl4{:});
gapow6_raw = ft_freqgrandaverage([],powl6{:});

% Load baselined powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern_bl
    powl2_bl{subj} = powload2_bl;
    powl4_bl{subj} = powload4_bl;
    powl6_bl{subj} = powload6_bl;
end

% Compute grand avg of baselined powspctrm data
gapow2_bl = ft_freqgrandaverage([], powl2_bl{:});
gapow4_bl = ft_freqgrandaverage([], powl4_bl{:});
gapow6_bl = ft_freqgrandaverage([], powl6_bl{:});

%% Plot alpha power grand average POWERSPECTRUM
close all
figure;
set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
conditions = {'WM load 2', 'WM load 4', 'WM load 6'};

% Plot
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 1;
ft_singleplotER(cfg, gapow2, gapow4, gapow6);
hold on;

% Add shadedErrorBar
channels_seb = ismember(gapow2.label, cfg.channel);
eb2 = shadedErrorBar(gapow2.freq, mean(gapow2.powspctrm(channels_seb, :), 1), ...
    std(gapow2.powspctrm(channels_seb, :)) / sqrt(size(gapow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
eb4 = shadedErrorBar(gapow4.freq, mean(gapow4.powspctrm(channels_seb, :), 1), ...
    std(gapow4.powspctrm(channels_seb, :)) / sqrt(size(gapow4.powspctrm(channels_seb, :), 1)), {'-'}, 0);
eb6 = shadedErrorBar(gapow6.freq, mean(gapow6.powspctrm(channels_seb, :), 1), ...
    std(gapow6.powspctrm(channels_seb, :)) / sqrt(size(gapow6.powspctrm(channels_seb, :), 1)), {'-'}, 0);
eb2.mainLine.Color = colors(1, :);
eb4.mainLine.Color = colors(2, :);
eb6.mainLine.Color = colors(3, :);
eb2.patch.FaceColor = colors(1, :);
eb4.patch.FaceColor = colors(2, :);
eb6.patch.FaceColor = colors(3, :);
set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
set(eb4.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
set(eb6.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
set(eb2.edge(1), 'Color', colors(1, :));
set(eb2.edge(2), 'Color', colors(1, :));
set(eb4.edge(1), 'Color', colors(2, :));
set(eb4.edge(2), 'Color', colors(2, :));
set(eb6.edge(1), 'Color', colors(3, :));
set(eb6.edge(2), 'Color', colors(3, :));
set(eb2.patch, 'FaceAlpha', 0.25);
set(eb4.patch, 'FaceAlpha', 0.25);
set(eb6.patch, 'FaceAlpha', 0.25);

% Adjust plotting
set(gcf,'color','w');
set(gca,'Fontsize',20);
[~, channel_idx] = ismember(channels, gapow2.label);
freq_idx = find(gapow2.freq >= 8 & gapow2.freq <= 14);
max_spctrm = max([mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow4.powspctrm(channel_idx, freq_idx), 2); mean(gapow6.powspctrm(channel_idx, freq_idx), 2)]);
ylim([0 max_spctrm*1.4])
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
legend([eb2.mainLine, eb4.mainLine, eb6.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontName', 'Arial', 'FontSize', 20);
title('Sternberg Power Spectrum', 'FontSize', 30);
hold off;

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/AOC_alpha_power_sternberg_powspctrm.png');

%% Plot INDIVIDUAL power spectra
close all
output_dir = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/';
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg.mat')

for subj = 1:length(subjects)
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

    % Extract participant data
    pow2 = powl2{subj};
    pow4 = powl4{subj};
    pow6 = powl6{subj};

    % Figure common config
    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow2, pow4, pow6);
    hold on;

    % Add shaded error bars
    channels_seb = ismember(pow2.label, cfg.channel);
    eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
        std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb4 = shadedErrorBar(pow4.freq, mean(pow4.powspctrm(channels_seb, :), 1), ...
        std(pow4.powspctrm(channels_seb, :)) / sqrt(size(pow4.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb6 = shadedErrorBar(pow6.freq, mean(pow6.powspctrm(channels_seb, :), 1), ...
        std(pow6.powspctrm(channels_seb, :)) / sqrt(size(pow6.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb2.mainLine.Color = colors(1, :);
    eb4.mainLine.Color = colors(2, :);
    eb6.mainLine.Color = colors(3, :);
    eb2.patch.FaceColor = colors(1, :);
    eb4.patch.FaceColor = colors(2, :);
    eb6.patch.FaceColor = colors(3, :);
    set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(eb4.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(eb6.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
    set(eb2.edge(1), 'Color', colors(1, :));
    set(eb2.edge(2), 'Color', colors(1, :));
    set(eb4.edge(1), 'Color', colors(2, :));
    set(eb4.edge(2), 'Color', colors(2, :));
    set(eb6.edge(1), 'Color', colors(3, :));
    set(eb6.edge(2), 'Color', colors(3, :));
    set(eb2.patch, 'FaceAlpha', 0.5);
    set(eb4.patch, 'FaceAlpha', 0.5);
    set(eb6.patch, 'FaceAlpha', 0.5);

    % Add lines for IAF and AlphaPower for each condition
    currentSubj = str2double(subjects{subj});
    C1 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 2);
    if ~isempty(C1)
        xIAF = C1.IAF;
        yAlpha = C1.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
    end
    C2 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 4);
    if ~isempty(C2)
        xIAF = C2.IAF;
        yAlpha = C2.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
    end
    C3 = eeg_data_sternberg([eeg_data_sternberg.ID] == currentSubj & [eeg_data_sternberg.Condition] == 6);
    if ~isempty(C3)
        xIAF = C3.IAF;
        yAlpha = C3.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
    end

    % Adjust plot aesthetics
    set(gca, 'FontSize', 20);
    max_spctrm = max([max(max(pow2.powspctrm(channels_seb, :))), max(max(pow4.powspctrm(channels_seb, :))), max(max(pow6.powspctrm(channels_seb, :)))]);
    ylim([0 max_spctrm*0.75]);
    xlim([5 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    legend([eb2.mainLine, eb4.mainLine eb6.mainLine], {'WM Load 2', 'WM Load 4', 'WM Load 6'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    save_path = fullfile(output_dir, sprintf('AOC_powspctrm_sternberg_subj%s.png', subjects{subj}));
    saveas(gcf, save_path);
end

%% Plot SUBPLOT of INDIVIDUAL powerspectra
close all;
output_dir = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/';
num_subj = length(subjects);

% Determine subplot grid size
default_cols = 5;
nrows = ceil(num_subj / default_cols);
ncols = min(num_subj, default_cols);

% Create figure
figure;
set(gcf, 'Color', 'w', 'Position', [0, 0, 300 * ncols, 300 * nrows]);

for subj = 1:num_subj
    % Extract participant data
    pow2 = powl2{subj};
    pow4 = powl4{subj};
    pow6 = powl6{subj};

    % Select subplot position
    subplot(nrows, ncols, subj);
    hold on;

    % Figure common config
    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow2, pow4, pow6);

    % Add shaded error bars
    channels_seb = ismember(pow2.label, cfg.channel);
    eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
        std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb4 = shadedErrorBar(pow4.freq, mean(pow4.powspctrm(channels_seb, :), 1), ...
        std(pow4.powspctrm(channels_seb, :)) / sqrt(size(pow4.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb6 = shadedErrorBar(pow6.freq, mean(pow6.powspctrm(channels_seb, :), 1), ...
        std(pow6.powspctrm(channels_seb, :)) / sqrt(size(pow6.powspctrm(channels_seb, :), 1)), {'-'}, 0);

    eb2.mainLine.Color = colors(1, :);
    eb4.mainLine.Color = colors(2, :);
    eb6.mainLine.Color = colors(3, :);
    eb2.patch.FaceColor = colors(1, :);
    eb4.patch.FaceColor = colors(2, :);
    eb6.patch.FaceColor = colors(3, :);
    set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(eb4.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(eb6.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
    set(eb2.patch, 'FaceAlpha', 0.5);
    set(eb4.patch, 'FaceAlpha', 0.5);
    set(eb6.patch, 'FaceAlpha', 0.5);

    % Adjust plot aesthetics
    set(gca, 'FontSize', 12);
    max_spctrm = max([max(max(pow2.powspctrm(channels_seb, :))), max(max(pow4.powspctrm(channels_seb, :))), max(max(pow6.powspctrm(channels_seb, :))) ]);
    ylim([0 max_spctrm * 0.75]);
    xlim([5 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    if subj == 1
        legend([eb2.mainLine, eb4.mainLine eb6.mainLine], {'WM Load 2', 'WM Load 4', 'WM Load 6'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
    end
    title(sprintf('Subject %s', subjects{subj}), 'FontSize', 14);
    hold off;
end

% Save combined figure
save_path = fullfile(output_dir, 'AOC_powspctrm_sternberg_all_subs.png');
saveas(gcf, save_path);

%% Plot alpha power TOPOS
close all;
clc;
cfg = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.figure = 'gcf';
cfg.marker = 'off';
addpath('/Users/Arne/Documents/matlabtools/customcolormap/')
cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
[~, channel_idx] = ismember(channels, gapow2.label);
freq_idx = find(gapow2.freq >= 8 & gapow2.freq <= 14);
max_spctrm = max([mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow4.powspctrm(channel_idx, freq_idx), 2); mean(gapow6.powspctrm(channel_idx, freq_idx), 2)]);
cfg.zlim = [0 max_spctrm];

% Plot WM load 2
figure('Color', 'w');
set(gcf, 'Position', [0, 300, 800, 600]);
ft_topoplotER(cfg, gapow2);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 2', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo2.png');

% Plot WM load 4
figure('Color', 'w');
set(gcf, 'Position', [700, 300, 800, 600]);
ft_topoplotER(cfg, gapow4);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 4', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo4.png');

% Plot WM load 6
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 800, 600]);
ft_topoplotER(cfg, gapow6);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('WM load 6', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo6.png');

%% Plot alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (WM load 6 - WM load 2)
ga_diff = gapow6;
ga_diff.powspctrm = gapow6.powspctrm - gapow2.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 1000, 800]);
cfg = [];
load('/Volumes/methlab/Students/Arne/toolboxes/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = channels;
cfg.highlightsymbol = '.';
cfg.highlightsize = 10;
cfg.figure = 'gcf';
cfg.marker = 'off';
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);
cfg.colormap = cmap;
cb = colorbar;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
cfg.zlim = 'maxabs';
set(cb, 'FontSize', 20);
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25);
title('Sternberg Task Alpha Power Difference (WM load 6 - WM load 2)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_sternberg_topo_diff.png');
