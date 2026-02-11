%% AOC Topos — N-Back (Power)
% Loads power_nback, grand-averages powspctrm, plots topographical maps per condition. Saves figures.
%
% Key outputs:
%   Topographic maps (occipital alpha power, per condition)

startup
[subjects, path, colors, headmodel] = setup('AOC');
% exclude subject 364 because of faulty data in electrode Cz
subjects = setdiff(subjects, {'361'});

%% Load raw powerspctrm data
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj}, filesep, 'eeg');
    cd(datapath)
    clc
    disp('LOADING DATA...')
    disp(subj)
    load power_nback
    powl1{subj} = powload1;
    powl2{subj} = powload2;
    powl3{subj} = powload3;
end

% Compute grand avg of raw powspctrm data
gapow1 = ft_freqgrandaverage([],powl1{:});
gapow2 = ft_freqgrandaverage([],powl2{:});
gapow3 = ft_freqgrandaverage([],powl3{:});

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload1.label)
    label = powload1.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Data check: Histogram / bar plot of CZ alpha power per subject
% close all
% clc
% 
% % Which condition to inspect (powl3)
% data_cond = powl3; 
% 
% % Electrode to check
% target_ch = 'Cz';
% 
% % Frequency window
% f_low = 8; 
% f_high = 14;
% 
% alpha_vals = nan(length(subjects),1);
% 
% for subj = 1:length(subjects)
%     % extract subject freq data
%     freq = data_cond{subj};
% 
%     % find channel index
%     ch_idx = find(strcmp(freq.label, target_ch));
%     if isempty(ch_idx)
%         warning(['Channel ', target_ch, ' not found for subject ', num2str(subjects{subj})]);
%         continue
%     end
% 
%     % find freq indices
%     f_idx = freq.freq >= f_low & freq.freq <= f_high;
% 
%     % average alpha power (averaged across freq and time collapsed already)
%     % If powspctrm is chans x freqs, use:
%     alpha_vals(subj) = mean(freq.powspctrm(ch_idx, f_idx), 'omitnan');
% 
%     if alpha_vals(subj) > 1
%         disp(subj)
%         disp(subjects{subj})
%         disp(alpha_vals(subj))
%     end
% end
% 
% % Plot
% figure('Color','w');
% set(gcf, 'Position', [200 200 1400 500])
% bar(alpha_vals, 'FaceColor', [0.4 0.4 0.9])
% xlabel('Subject')
% ylabel(['Alpha Power at ', target_ch, ' (', num2str(f_low), '-', num2str(f_high), ' Hz)'])
% title(['Cz Alpha Power Across Subjects (Condition = data\_cond)'])
% set(gca, 'FontSize', 18)
% set(gca, "YScale", "log")

%% Plot alpha power TOPOS
close all;
clc;
fontSize = 75;

cfg = [];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/layANThead.mat');
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
%cmap = cbrewer('seq', 'Reds', 64); % 'RdBu' for blue to red diverging color map
cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 14];
%cfg.fontsize = fontSize*100;

% Global alpha zlim across all channels & conditions (robust)
[~, channel_idx] = ismember(channels, gapow1.label);
freq_idx = find(gapow1.freq >= 8 & gapow1.freq <= 14);
A1 = mean(gapow1.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
A2 = mean(gapow2.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
A3 = mean(gapow3.powspctrm(channel_idx, freq_idx), 2, 'omitnan');
all_alpha = [A1(:); A2(:); A3(:)];
global_max = prctile(all_alpha,99);
cfg.zlim = [0 global_max];

% Plot 1-back
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow1);
cb = colorbar;
set(findall(gcf,'Type','axes'),'FontSize',fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('1-back');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_nback_load1.png');

% Plot 2-back
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow2);
cb = colorbar;
set(findall(gcf,'Type','axes'),'FontSize',fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('2-back');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_nback_load2.png');

% Plot 3-back
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 2000, 2000]);
ft_topoplotER(cfg, gapow3);
cb = colorbar;
set(findall(gcf,'Type','axes'),'FontSize',fontSize)
ylabel(cb, 'Power [\muV^2/Hz]');
title('3-back');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_nback_load3.png');

%% Plot alpha power TOPOS DIFFERENCE
close all
clc

% Calculate the difference between the conditions (3-back - 1-back)
ga_diff = gapow3;
ga_diff.powspctrm = gapow3.powspctrm - gapow1.powspctrm;

% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 1000, 800]);
cfg = [];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/layANThead.mat');
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
title('N-back Task Alpha Power Difference (3-back - 1-back)', 'FontSize', 25);
ft_topoplotER(cfg, ga_diff);

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/topos/AOC_eeg_topos_nback_diff.png');
