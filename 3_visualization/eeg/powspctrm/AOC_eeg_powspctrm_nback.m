%% AOC Alpha Power N-back

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
% Load power at IAF and IAF data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath);
    load('alpha_power_nback.mat');
    powIAF1(subj) = powerIAF1;
    powIAF2(subj) = powerIAF2;
    powIAF3(subj) = powerIAF3;
end

% Load powspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_nback
    powl1{subj} = powload1;
    powl2{subj} = powload2;
    powl3{subj} = powload3;
end

% Compute grand avg of raw powspctrm data
gapow1_raw = ft_freqgrandaverage([],powl1{:});
gapow2_raw = ft_freqgrandaverage([],powl2{:});
gapow3_raw = ft_freqgrandaverage([],powl3{:});

% Load baselined powspctrm data
% for subj = 1:length(subjects)
%     datapath = strcat(path, subjects{subj}, filesep, 'eeg');
%     cd(datapath)
%     load power_nback_bl
%     powl1_bl{subj} = powload1_bl;
%     powl2_bl{subj} = powload2_bl;
%     powl3_bl{subj} = powload3_bl;
% end
% 
% % Compute grand avg of baselined powspctrm data
% gapow1_bl = ft_freqgrandaverage([], powl1_bl{:});
% gapow2_bl = ft_freqgrandaverage([], powl2_bl{:});
% gapow3_bl = ft_freqgrandaverage([], powl3_bl{:});

%% Plot alpha power grand average POWERSPECTRUM
% Loop over analysis conditions
analysis_conditions = {'raw'}%, 'bl'};

for i = 1:length(analysis_conditions)
    for electrodes = {'occ_cluster'}%, 'POz'}
        switch analysis_conditions{i}
            case 'raw'
                gapow1 = gapow1_raw;
                gapow2 = gapow2_raw;
                gapow3 = gapow3_raw;
                yLabel = 'Power [\muV^2/Hz]';
                figTitle = 'N-back Power Spectrum';
                if strcmp(electrodes, 'occ_cluster')
                    savePath = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_powspctrm_nback_raw.png';
                elseif strcmp(electrodes, 'POz')
                    savePath = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_powspctrm_nback_raw_elecPOz.png';
                end
            % case 'bl'
            %     gapow1 = gapow1_bl;
            %     gapow2 = gapow2_bl;
            %     gapow3 = gapow3_bl;
            %     yLabel = 'Power [dB]';
            %     figTitle = 'Baselined N-back Power Spectrum';
            %     if strcmp(electrodes, 'occ_cluster')
            %         savePath = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_alpha_power_nback_powspctrm_bl.png';
            %     elseif strcmp(electrodes, 'POz')
            %         savePath = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_alpha_power_nback_powspctrm_bl_elecPOz.png';
            %     end
        end
        % Create figure
        close all
        figure;
        set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
        conditions = {'1-back', '2-back', '3-back'};

        % Plot
        cfg = [];
        if strcmp(electrodes, 'occ_cluster')
            cfg.channel = channels;
        elseif strcmp(electrodes, 'POz')
            cfg.channel = 'POz';
        end
        cfg.figure = 'gcf';
        cfg.linewidth = 3;
        ft_singleplotER(cfg, gapow1, gapow2, gapow3);
        hold on;

        % Add shadedErrorBar
        channels_seb = ismember(gapow1.label, cfg.channel);
        eb1 = shadedErrorBar(gapow1.freq, mean(gapow1.powspctrm(channels_seb, :), 1), ...
            std(gapow1.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow1.powspctrm(channels_seb, :), 1)), 'lineProps', {'-','Color',colors(1,:)});
        eb2 = shadedErrorBar(gapow2.freq, mean(gapow2.powspctrm(channels_seb, :), 1), ...
            std(gapow2.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow2.powspctrm(channels_seb, :), 1)), 'lineProps', {'-','Color',colors(1,:)});
        eb3 = shadedErrorBar(gapow3.freq, mean(gapow3.powspctrm(channels_seb, :), 1), ...
            std(gapow3.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow3.powspctrm(channels_seb, :), 1)), 'lineProps', {'-','Color',colors(1,:)});
        eb1.mainLine.Color = colors(1, :);
        eb2.mainLine.Color = colors(2, :);
        eb3.mainLine.Color = colors(3, :);
        eb1.patch.FaceColor = colors(1, :);
        eb2.patch.FaceColor = colors(2, :);
        eb3.patch.FaceColor = colors(3, :);
        set(eb1.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
        set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
        set(eb3.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
        set(eb1.edge(1), 'Color', colors(1, :));
        set(eb1.edge(2), 'Color', colors(1, :));
        set(eb2.edge(1), 'Color', colors(2, :));
        set(eb2.edge(2), 'Color', colors(2, :));
        set(eb3.edge(1), 'Color', colors(3, :));
        set(eb3.edge(2), 'Color', colors(3, :));
        set(eb1.patch, 'FaceAlpha', 0.5);
        set(eb2.patch, 'FaceAlpha', 0.5);
        set(eb3.patch, 'FaceAlpha', 0.5);

        % Adjust plotting
        set(gcf,'color','w');
        set(gca,'Fontsize',20);
        [~, channel_idx] = ismember(cfg.channel, gapow1.label);
        freq_idx = find(gapow1.freq >= 8 & gapow1.freq <= 14);
        max_spctrm = max([mean(gapow1.powspctrm(channel_idx, freq_idx), 2); mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow3.powspctrm(channel_idx, freq_idx), 2)]);
        if strcmp(electrodes, 'occ_cluster')
            ylim([0 max_spctrm*1.25])
        elseif strcmp(electrodes, 'POz')
            ylim([0 0.45])
        end
        box on
        xlim([5 20]);
        xlabel('Frequency [Hz]');
        ylabel(yLabel);
        legend([eb1.mainLine, eb2.mainLine, eb3.mainLine], {'1 back', '2 back', '3 back'}, 'FontName', 'Arial', 'FontSize', 20);
        title(figTitle, 'FontSize', 30)
        if strcmp(electrodes, 'POz')
            text(25, max_spctrm*1.5, ['Electrodes: ', electrodes]);
        end
        hold off;

        % Save
        saveas(gcf, savePath);
    end
end

%% Plot INDIVIDUAL power spectra
close all
output_dir = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/';
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_nback.mat')

for subj = 1:length(subjects)
    clear pow1 pow2 pow3
    close all;
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');

    % Extract participant data
    pow1 = powl1{subj};
    pow2 = powl2{subj};
    pow3 = powl3{subj};

    % Figure common config
    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow1, pow2, pow3);
    hold on;

    % Add shaded error bars
    channels_seb = ismember(pow1.label, cfg.channel);
    eb1 = shadedErrorBar(pow1.freq, mean(pow1.powspctrm(channels_seb, :), 1), ...
        std(pow1.powspctrm(channels_seb, :)) / sqrt(size(pow1.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
        std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb3 = shadedErrorBar(pow3.freq, mean(pow3.powspctrm(channels_seb, :), 1), ...
        std(pow3.powspctrm(channels_seb, :)) / sqrt(size(pow3.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb1.mainLine.Color = colors(1, :);
    eb2.mainLine.Color = colors(2, :);
    eb3.mainLine.Color = colors(3, :);
    eb1.patch.FaceColor = colors(1, :);
    eb2.patch.FaceColor = colors(2, :);
    eb3.patch.FaceColor = colors(3, :);
    set(eb1.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(eb3.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
    set(eb1.edge(1), 'Color', colors(1, :));
    set(eb1.edge(2), 'Color', colors(1, :));
    set(eb2.edge(1), 'Color', colors(2, :));
    set(eb2.edge(2), 'Color', colors(2, :));
    set(eb3.edge(1), 'Color', colors(3, :));
    set(eb3.edge(2), 'Color', colors(3, :));
    set(eb1.patch, 'FaceAlpha', 0.5);
    set(eb2.patch, 'FaceAlpha', 0.5);
    set(eb3.patch, 'FaceAlpha', 0.5);

    % Add lines for IAF and AlphaPower for each condition
    currentSubj = str2double(subjects{subj});
    C1 = eeg_data_nback([eeg_data_nback.ID] == currentSubj & [eeg_data_nback.Condition] == 1);
    if ~isempty(C1)
        xIAF = C1.IAF;
        yAlpha = C1.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(1, :), 'LineWidth', 2);
    end
    C2 = eeg_data_nback([eeg_data_nback.ID] == currentSubj & [eeg_data_nback.Condition] == 2);
    if ~isempty(C2)
        xIAF = C2.IAF;
        yAlpha = C2.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(2, :), 'LineWidth', 2);
    end
    C3 = eeg_data_nback([eeg_data_nback.ID] == currentSubj & [eeg_data_nback.Condition] == 3);
    if ~isempty(C3)
        xIAF = C3.IAF;
        yAlpha = C3.AlphaPower;
        line([xIAF, xIAF], [0, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
        line([5, xIAF], [yAlpha, yAlpha], 'LineStyle', '--', 'Color', colors(3, :), 'LineWidth', 2);
    end

    % Adjust plot aesthetics
    set(gca, 'FontSize', 20);
    max_spctrm = max([max(max(pow1.powspctrm(channels_seb, :))), max(max(pow2.powspctrm(channels_seb, :))), max(max(pow3.powspctrm(channels_seb, :)))]);
    ylim([0 max_spctrm*0.75]);
    xlim([5 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    legend([eb1.mainLine, eb2.mainLine eb3.mainLine], {'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 20);
    title(sprintf('Subject %s: Power Spectrum', subjects{subj}), 'FontSize', 30);
    hold off;

    % Save individual plot
    save_path = fullfile(output_dir, sprintf('AOC_powspctrm_nback_subj%s.png', subjects{subj}));
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
    pow1 = powl1{subj};
    pow2 = powl2{subj};
    pow3 = powl3{subj};

    % Select subplot position
    subplot(nrows, ncols, subj);
    hold on;

    % Figure common config
    cfg = [];
    cfg.channel = channels;
    cfg.figure = 'gcf';
    cfg.linewidth = 1;

    % Plot power spectrum for low and high contrast
    ft_singleplotER(cfg, pow1, pow2, pow3);

    % Add shaded error bars
    channels_seb = ismember(pow1.label, cfg.channel);
    eb1 = shadedErrorBar(pow1.freq, mean(pow1.powspctrm(channels_seb, :), 1), ...
        std(pow1.powspctrm(channels_seb, :)) / sqrt(size(pow1.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb2 = shadedErrorBar(pow2.freq, mean(pow2.powspctrm(channels_seb, :), 1), ...
        std(pow2.powspctrm(channels_seb, :)) / sqrt(size(pow2.powspctrm(channels_seb, :), 1)), {'-'}, 0);
    eb3 = shadedErrorBar(pow3.freq, mean(pow3.powspctrm(channels_seb, :), 1), ...
        std(pow3.powspctrm(channels_seb, :)) / sqrt(size(pow3.powspctrm(channels_seb, :), 1)), {'-'}, 0);

    eb1.mainLine.Color = colors(1, :);
    eb2.mainLine.Color = colors(2, :);
    eb3.mainLine.Color = colors(3, :);
    eb1.patch.FaceColor = colors(1, :);
    eb2.patch.FaceColor = colors(2, :);
    eb3.patch.FaceColor = colors(3, :);
    set(eb1.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
    set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
    set(eb3.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
    set(eb1.patch, 'FaceAlpha', 0.5);
    set(eb2.patch, 'FaceAlpha', 0.5);
    set(eb3.patch, 'FaceAlpha', 0.5);

    % Adjust plot aesthetics
    set(gca, 'FontSize', 12);
    max_spctrm = max([max(max(pow1.powspctrm(channels_seb, :))), max(max(pow2.powspctrm(channels_seb, :))), max(max(pow3.powspctrm(channels_seb, :))) ]);
    ylim([0 max_spctrm * 0.75]);
    xlim([5 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    if subj == 1
        legend([eb1.mainLine, eb2.mainLine eb3.mainLine], {'1-back', '2-back', '3-back'}, 'FontName', 'Arial', 'FontSize', 15, 'Location', 'best');
    end
    title(sprintf('Subject %s', subjects{subj}), 'FontSize', 14);
    hold off;
end

% Save combined figure
save_path = fullfile(output_dir, 'AOC_powspctrm_nback_all_subs.png');
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
[~, channel_idx] = ismember(channels, gapow1.label);
freq_idx = find(gapow1.freq >= 8 & gapow1.freq <= 14);
max_spctrm = max([mean(gapow1.powspctrm(channel_idx, freq_idx), 2); mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow3.powspctrm(channel_idx, freq_idx), 2)]);
cfg.zlim = [0 max_spctrm];

% Plot 1-back
figure('Color', 'w');
set(gcf, 'Position', [0, 300, 800, 600]);
ft_topoplotER(cfg, gapow1);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('1-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo1.png');

% Plot 2-back
figure('Color', 'w');
set(gcf, 'Position', [300, 300, 800, 600]);
ft_topoplotER(cfg, gapow2);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('2-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo2.png');

% Plot 3-back
figure('Color', 'w');
set(gcf, 'Position', [600, 300, 800, 600]);
ft_topoplotER(cfg, gapow3);
title('');
cb = colorbar;
set(cb, 'FontSize', 20);
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25);
title('3-back', 'FontSize', 40);
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo3.png');

%% Plot alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (3-back - 1-back)
ga_diff = gapow3;
ga_diff.powspctrm = gapow3.powspctrm - gapow1.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]);
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
title('N-back Alpha Power Difference (3-back - 1-back)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/topos/AOC_alpha_power_nback_topo_diff.png');
