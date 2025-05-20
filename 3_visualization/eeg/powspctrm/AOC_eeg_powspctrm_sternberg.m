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
        disp(['Midline channel: ', ch])
    end
end

%% Load data
clc
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

% Load raw powerspctrm data
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

% Load long time window (200ms - 2000ms) powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern_long
    powl2long{subj} = powload2long;
    powl4long{subj} = powload4long;
    powl6long{subj} = powload6long;
end

% Compute grand avg of long time window (200ms - 2000ms) powerspctrm data
gapow2_long = ft_freqgrandaverage([],powl2long{:});
gapow4_long = ft_freqgrandaverage([],powl4long{:});
gapow6_long = ft_freqgrandaverage([],powl6long{:});

% Load baseline period powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern_baseline_period
    powl2_blperiod{subj} = powload2_baseline_period;
    powl4_blperiod{subj} = powload4_baseline_period;
    powl6_blperiod{subj} = powload6_baseline_period;
end

% Compute baseline-corrected POWERSPECTRUM
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

% Load FOOOF powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_stern_fooof

    % FOOOF
    powl2_fooof{subj} = pow2_fooof;
    powl4_fooof{subj} = pow4_fooof;
    powl6_fooof{subj} = pow6_fooof;

    % FOOOF Baselined
    powl2_fooof_bl{subj} = pow2_fooof_bl;
    powl4_fooof_bl{subj} = pow4_fooof_bl;
    powl6_fooof_bl{subj} = pow6_fooof_bl;

    % FOOOF Baselined Long
    powl2_fooof_bl_long{subj} = pow2_fooof_bl_long;
    powl4_fooof_bl_long{subj} = pow4_fooof_bl_long;
    powl6_fooof_bl_long{subj} = pow6_fooof_bl_long;
end

% Compute grand avg of long time window (200ms - 2000ms) powerspctrm data
gapow2_fooof = ft_freqgrandaverage([], powl2_fooof{:});
gapow4_fooof = ft_freqgrandaverage([], powl4_fooof{:});
gapow6_fooof = ft_freqgrandaverage([], powl6_fooof{:});
gapow2_fooof_bl = ft_freqgrandaverage([], powl2_fooof_bl{:});
gapow4_fooof_bl = ft_freqgrandaverage([], powl4_fooof_bl{:});
gapow6_fooof_bl = ft_freqgrandaverage([], powl6_fooof_bl{:});
gapow2_fooof_bl_long = ft_freqgrandaverage([], powl2_fooof_bl_long{:});
gapow4_fooof_bl_long = ft_freqgrandaverage([], powl4_fooof_bl_long{:});
gapow6_fooof_bl_long = ft_freqgrandaverage([], powl6_fooof_bl_long{:});

%% Plot alpha power grand average POWERSPECTRUM
% Prepare your data-sets
datasets = { ...
    struct('name','raw', 'pow2', gapow2_raw,  'pow4', gapow4_raw,  'pow6', gapow6_raw), ...
    struct('name','baselined', 'pow2', gapow2_bl,   'pow4', gapow4_bl,   'pow6', gapow6_bl), ...
    struct('name','long', 'pow2', gapow2_long, 'pow4', gapow4_long, 'pow6', gapow6_long),  ...
    struct('name','fooof', 'pow2', gapow2_fooof_bl_long, 'pow4', gapow4_fooof_bl_long, 'pow6', gapow6_fooof_bl_long)  ...
    };

% Prepare your electrode clusters
electrodeSets = struct( ...
    'occ_cluster', {channels}); %, ...
    %'POz',         {{'POz'}}), ...
    % 'right hemisphere',{right_channels} ...
    %);

% Loop over data types and electrode sets
for d = 1%%%%:numel(datasets)
    D = datasets{d};
    for fn = fieldnames(electrodeSets)'
        close all
        elecName = fn{1};
        elecChans = electrodeSets.(elecName);

        % extract the three conditions
        gapow2 = D.pow2;
        gapow4 = D.pow4;
        gapow6 = D.pow6;

        % now do exactly your plotting block, but using elecChans …
        cfg = [];
        figure('Position',[0,0,800,1600],'Color','w');
        cfg.channel   = elecChans;
        cfg.figure    = 'gcf';
        cfg.linewidth = 3;

        ft_singleplotER(cfg, gapow2, gapow4, gapow6);
        hold on;

        %–– shaded error bars ––
        isElec = ismember(gapow2.label, cfg.channel);
        freqs  = gapow2.freq;
        m2     = mean(gapow2.powspctrm(isElec,:),1);
        se2    = std (gapow2.powspctrm(isElec,:),0,1) / sqrt(sum(isElec));
        m4     = mean(gapow4.powspctrm(isElec,:),1);
        se4    = std (gapow4.powspctrm(isElec,:),0,1) / sqrt(sum(isElec));
        m6     = mean(gapow6.powspctrm(isElec,:),1);
        se6    = std (gapow6.powspctrm(isElec,:),0,1) / sqrt(sum(isElec));

        eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(1,:)});
        eb4 = shadedErrorBar(freqs, m4, se4, 'lineProps', {'-','Color',colors(2,:)});
        eb6 = shadedErrorBar(freqs, m6, se6, 'lineProps', {'-','Color',colors(3,:)});

        % collect them into a 1×3 struct array
        ebs = [eb2, eb4, eb6];

        % loop and set both line and patch properties
        for k = 1:numel(ebs)
            set( ebs(k).mainLine, ...
                'LineWidth', cfg.linewidth, ...
                'Color',     colors(k,:) );
            set( ebs(k).patch,    ...
                'FaceColor', colors(k,:), ...
                'FaceAlpha', 0.25 );
        end

        set(gca,'FontSize',20);
        box on
        %ylim([-0.175 0.175])
        ylim([0 2.75])
        if d == 4
            ylim([-0.25 0.25])
        end
        xlim([4 30]);
        ylabel('Power [\muV^2/Hz]');
        xlabel('Frequency [Hz]');
        legend([eb2.mainLine, eb4.mainLine, eb6.mainLine], ...
            {'WM load 2','WM load 4','WM load 6'}, ...
            'FontName','Arial','FontSize',20);
        if strcmp(D.name, 'raw') & strcmp(elecName, 'occ_cluster')
            title('Sternberg Power Spectrum', 'FontSize', 30);
        else
        title(sprintf('Sternberg Power Spectrum — %s — %s', D.name, elecName), 'FontSize',30);
        end

        %–– save out with a descriptive filename ––
        outfn = sprintf('AOC_powspctrm_sternberg_%s_%s.png', D.name, elecName);
        saveas(gcf, fullfile('/Volumes/methlab/Students/Arne/AOC/figures/eeg/powspctrm', outfn));
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

