%% AOC TFR — Sternberg (Baselined, Baseline Period)
% Loads tfr_stern (tfr*_bl), grand-averages, plots TFR with focus on baseline period. Saves figures.
%
% Key outputs:
%   TFR figures (baselined, baseline window)

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

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

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, '/eeg');
    cd(datapath)
    load tfr_stern
    tfr2_all{subj} = tfr2_bl;
    tfr4_all{subj} = tfr4_bl;
    tfr6_all{subj} = tfr6_bl;
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' TFR data loaded.'])
end

% Compute grand average
gatfr2 = ft_freqgrandaverage([],tfr2_all{:});
gatfr4 = ft_freqgrandaverage([],tfr4_all{:});
gatfr6 = ft_freqgrandaverage([],tfr6_all{:});

%% Plot BASELINED alpha power TOPOS
close all
alpha_range = [8 14];
conds = 1;
for gatfr = {gatfr2, gatfr4, gatfr6}
    gatfr = gatfr{1};

    % Plot
    fig = figure('Color', 'w', 'Position', [100, 100, 1200, 1200]);
    title(['Sternberg WM load ', num2str(conds*2)], 'FontSize', 25)

    % Set up configuration
    cfg = [];
    cfg.layout = headmodel.layANThead;
    allchannels = cfg.layout.label;
    cfg.figure = fig;
    cfg.channel = allchannels(1:end-2);
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
    cfg.highlight = 'on';
    cfg.highlightchannel = channels;
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 5;
    cfg.marker = 'off';
    cfg.comment = 'no';
    cmap = flipud(cbrewer('div', 'RdBu', 64));
    cfg.colormap = cmap;
    cfg.gridscale = 300;
    cfg.ylim = alpha_range;
    cb = colorbar;
    set(cb, 'FontSize', 20);
    ylabel(cb, 'Power [dB]', 'FontSize', 25);

    % Set time window
    startTmpnt = find(gatfr.time == 1);
    endTmpnt = find(gatfr.time == 2);
    cfg.xlim = [gatfr.time(startTmpnt) gatfr.time(endTmpnt)];

    % Find max power value for frequency band
    freq_idx = find(gatfr.freq >= alpha_range(1) & gatfr.freq <= alpha_range(2));
    timepnts_idx = find(gatfr.time == 0) : find(gatfr.time == 2); % Timepoints over entire interval
    mat.powspctrm = gatfr.powspctrm(1:end-3, freq_idx, timepnts_idx); % [chan x freq x time]
    max_spctrm = max(mat.powspctrm(:), [], 'omitnan');
    max_spctrm = abs(max_spctrm);
    max_spctrm = 0.575
    cfg.zlim = [-max_spctrm max_spctrm];

    % Create Topoplot
    ft_topoplotTFR(cfg, gatfr);

    % Save figure
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/eeg/topos/baseline/AOC_alpha_power_sternberg_topo', num2str(conds*2) '_bl.png']);
    conds = conds+1;
end
