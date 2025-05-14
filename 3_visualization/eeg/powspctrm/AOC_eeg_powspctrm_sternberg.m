%% AOC POWSPCTRM Sternberg

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
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;
%channels(25:31) = [{'Pz'}, {'P1'}, {'P2'}, {'P3'}, {'P4'}, {'P5'}, {'P6'}]

% Left and right channels
left_channels = {};
right_channels = {};
for i = 1:length(channels)
    try
        ch = channels{i};
        % Find the first numeric part in the channel name
        numStr = regexp(ch, '\d+', 'match');
        % Convert the first numerical token to a number
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch;
        else
            right_channels{end+1} = ch;
        end
    catch ME
       ME.message
       disp(['Midline channel: ', ch])
    end
end

%% Load data
disp('LOADING DATA...')

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

% Load baseline period powerspctrm data
% for subj = 1:length(subjects)
%     datapath = strcat(path,subjects{subj}, filesep, 'eeg');
%     cd(datapath)
%     load power_stern_baseline_period
%     powl2_blperiod{subj} = powload2_baseline_period;
%     powl4_blperiod{subj} = powload4_baseline_period;
%     powl6_blperiod{subj} = powload6_baseline_period;
% end

%% Compute baseline-corrected POWERSPECTRUM
% Relative change: (retention - baseline) / baseline
for subj = 1:length(subjects)
    pow2_bl{subj} = powl2{subj};
    pow2_bl{subj}.powspctrm = 100 * (powl2{subj}.powspctrm - powl2_blperiod{subj}.powspctrm) ./ powl2_blperiod{subj}.powspctrm;
    pow4_bl{subj} = powl4{subj};
    pow4_bl{subj}.powspctrm = 100 * (powl4{subj}.powspctrm - powl4_blperiod{subj}.powspctrm) ./ powl4_blperiod{subj}.powspctrm;
    pow6_bl{subj} = powl6{subj};
    pow6_bl{subj}.powspctrm = 100 * (powl6{subj}.powspctrm - powl6_blperiod{subj}.powspctrm) ./ powl6_blperiod{subj}.powspctrm;
end

% Compute grand average of baseline-corrected power spectra
gapow2_bl = ft_freqgrandaverage([], pow2_bl{:});
gapow4_bl = ft_freqgrandaverage([], pow4_bl{:});
gapow6_bl = ft_freqgrandaverage([], pow6_bl{:});

%% Plot alpha power grand average POWERSPECTRUM
gapow2 = gapow2_raw;
gapow4 = gapow4_raw;
gapow6 = gapow6_raw;
% gapow2 = gapow2_bl;
% gapow4 = gapow4_bl;
% gapow6 = gapow6_bl;

for electrodes = {'occ_cluster'} %, 'POz'} %, 'right_hemisphere'}
    close all
    cfg = [];
    figure;
    set(gcf, 'Position', [0, 0, 800, 1600], 'Color', 'w');
    conditions = {'WM load 2', 'WM load 4', 'WM load 6'};

    % Plot
    cfg = [];
    if strcmp(electrodes, 'occ_cluster')
        cfg.channel = channels;
    elseif strcmp(electrodes, 'POz')
        cfg.channel = 'POz';
    elseif strcmp(electrodes, 'right_hemisphere')
        cfg.channel = right_channels;
    end
    cfg.figure = 'gcf';
    cfg.linewidth = 1;
    ft_singleplotER(cfg, gapow2, gapow4, gapow6);
    hold on;

    % Add shadedErrorBar
    channels_seb = ismember(gapow2.label, cfg.channel);
    lp2     = {'-','Color', colors(1,:)};
    lp4     = {'-','Color', colors(2,:)};
    lp6     = {'-','Color', colors(3,:)};
    eb2 = shadedErrorBar(gapow2.freq, mean(gapow2.powspctrm(channels_seb, :), 1), ...
        std(gapow2.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow2.powspctrm(channels_seb, :), 1)), 'lineProps', lp2);
    eb4 = shadedErrorBar(gapow4.freq, mean(gapow4.powspctrm(channels_seb, :), 1), ...
        std(gapow4.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow4.powspctrm(channels_seb, :), 1)), 'lineProps', lp4);
    eb6 = shadedErrorBar(gapow6.freq, mean(gapow6.powspctrm(channels_seb, :), 1), ...
        std(gapow6.powspctrm(channels_seb, :), 0, 1) / sqrt(size(gapow6.powspctrm(channels_seb, :), 1)), 'lineProps', lp6);
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
    [~, channel_idx] = ismember(cfg.channel, gapow2.label);
    freq_idx = find(gapow2.freq >= 4 & gapow2.freq <= 30);
    max_spctrm = max([mean(gapow2.powspctrm(channel_idx, freq_idx), 2); mean(gapow4.powspctrm(channel_idx, freq_idx), 2); mean(gapow6.powspctrm(channel_idx, freq_idx), 2)]);
    %ylim([0 max_spctrm*1.25])
    %ylim([-200 200])
    ylim([0 2.6])
    xlim([4 30])
    if strcmp(electrodes, 'POz')
        ylim([0 5])
    elseif strcmp(electrodes, 'right_hemisphere')
        ylim([0 1])
    end
    box on
    yline(0, '--')
    xlabel('Frequency [Hz]');
    %ylabel('Relative Change [%]');
    ylabel('Power [\muV^2/Hz]');
    legend([eb2.mainLine, eb4.mainLine, eb6.mainLine], {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontName', 'Arial', 'FontSize', 20);
    title('Sternberg Power Spectrum', 'FontSize', 30);
    if ~strcmp(electrodes, 'occ_cluster')
        text(25, max_spctrm*1.5, ['Electrodes: ', electrodes]);
    end
    hold off;

    % Save
    if strcmp(electrodes, 'occ_cluster')
        saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_alpha_power_sternberg_powspctrm.png');
    elseif strcmp(electrodes, 'POz')
        saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_alpha_power_sternberg_powspctrm_elecPOz.png');
    elseif strcmp(electrodes, 'right_hemisphere')
        saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm/AOC_alpha_power_sternberg_powspctrm_elecRH.png');
    end
end

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
% close all;
% output_dir = '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_power/powspctrm/';
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

