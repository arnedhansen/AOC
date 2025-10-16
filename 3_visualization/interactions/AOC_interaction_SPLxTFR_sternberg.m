%% AOC Sternberg — Median-split by within-subject Scan Path Length (SPL) to compare TFRs

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat') % FieldTrip layout
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

colors = color_def('AOC');
fontSize = 36;

% Plot colormap
color_map = cbrewer('seq', 'Reds', 64);



subjects = subjects(1:10)




%% Preallocate holders for subject-level TFRs
low_tfr_subs = cell(1, length(subjects)); % per-subject LOW-SPL TFR (avg over trials)
high_tfr_subs = cell(1, length(subjects)); % per-subject HIGH-SPL TFR (avg over trials)

%% Per-subject split (by SPL) and aggregation of EEG TFRs
for s = 1:length(subjects)
clc
subjID = str2double(subjects{s});
fprintf('Subject %s (%d/%d)\n', subjects{s}, s, length(subjects));

% Get per-trial SPL and trial numbers for this subject
rows = merged_data_sternberg_trials(merged_data_sternberg_trials.ID == subjID, :);
spl  = rows.ScanPathLengthFull;  % total scan path length per trial (Full window)
trlN = rows.Trial;

good = isfinite(spl) & isfinite(trlN);
if ~any(good)
    warning('No finite SPL for subject %s. Skipping subject.', subjects{s})
    continue
end

% Random tie-breaking median split (within-subject, by SPL)
spl_sub  = spl(good);
trl_sub  = trlN(good);
rp       = randperm(numel(spl_sub));  % random permutation to break ties stably
spl_rand = spl_sub(rp);
trl_rand = trl_sub(rp);

[spl_sorted, idx_sorted] = sort(spl_rand, 'ascend'); %#ok<ASGLU> % keep order for index
trl_sorted = trl_rand(idx_sorted);

nHalf = floor(numel(trl_sorted)/2);
lowTrials  = trl_sorted(1:nHalf);
highTrials = trl_sorted(nHalf+1:end);

% Load EEG TFR and average LOW/HIGH-SPL trials over repetitions (occipital channels)
tfr_all = [];
try
    datapath_eeg = fullfile(path, subjects{s}, 'eeg');
    cd(datapath_eeg)
    load tfr_stern_trials   % -> tfr_all (dimord 'rpt_chan_freq_time'), trialinfo(:,2) = Trial
catch
    warning('Missing EEG TFR for subject %s, skipping.', subjects{s})
    continue
end

if isempty(tfr_all)
    warning('Empty tfr_all for subject %s, skipping.', subjects{s})
    continue
end

% Occipital channel heuristic (labels containing 'O' or 'I'); fallback to all
occ_channels = {};
for i = 1:length(tfr_all.label)
    lab = tfr_all.label{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        occ_channels{end+1} = lab; %#ok<AGROW>
    end
end
if isempty(occ_channels)
    occ_channels = tfr_all.label;
end

if size(tfr_all.trialinfo,2) < 2
    warning('tfr_all.trialinfo missing Trial column for subject %s. Skipping.', subjects{s})
    continue
end

eegTrials = tfr_all.trialinfo(:,2);

idxLow  = ismember(eegTrials,  lowTrials);
idxHigh = ismember(eegTrials, highTrials);

low_tfr  = [];
high_tfr = [];

if any(idxLow)
    cfgS = [];
    cfgS.trials     = find(idxLow);
    cfgS.channel    = occ_channels;
    cfgS.avgoverrpt = 'yes';
    low_tfr  = ft_selectdata(cfgS, tfr_all);
end

if any(idxHigh)
    cfgS = [];
    cfgS.trials     = find(idxHigh);
    cfgS.channel    = occ_channels;
    cfgS.avgoverrpt = 'yes';
    high_tfr = ft_selectdata(cfgS, tfr_all);
end

% Store per-subject results (skip if either side is empty)
if ~isempty(low_tfr)
    low_tfr_subs{s} = low_tfr;
end
if ~isempty(high_tfr)
    high_tfr_subs{s} = high_tfr;
end


end

%% Grand-average TFRs (LOW vs HIGH SPL)
low_tfr_subs = low_tfr_subs(~cellfun(@isempty, low_tfr_subs));
high_tfr_subs = high_tfr_subs(~cellfun(@isempty, high_tfr_subs));

if isempty(low_tfr_subs) || isempty(high_tfr_subs)
error('No subject TFRs available for grand average (LOW/HIGH SPL).')
end

gatfr_low = ft_freqgrandaverage([], low_tfr_subs{:});
gatfr_high = ft_freqgrandaverage([], high_tfr_subs{:});

%% Plot settings (FieldTrip singleplotTFR)
cfg = [];
cfg.channel = gatfr_low.label; % already restricted per subject to occipital
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [0 2];
cfg.ylim = [4 30];
cfg.layout = layANThead;

% Harmonise colour scaling across the two plots, using 0–2 s × 8–14 Hz window
alpha_idx = gatfr_low.freq >= 8 & gatfr_low.freq <= 14;
time_idx = gatfr_low.time >= 0 & gatfr_low.time <= 2;

[~, ch_low_idx] = ismember(cfg.channel, gatfr_low.label);
[~, ch_high_idx] = ismember(cfg.channel, gatfr_high.label);

low_alpha_power = mean(gatfr_low.powspctrm(ch_low_idx, alpha_idx, time_idx), 1:3, 'omitnan');
high_alpha_power = mean(gatfr_high.powspctrm(ch_high_idx, alpha_idx, time_idx), 1:3, 'omitnan');
max_spctrm = max([low_alpha_power(:); high_alpha_power(:)]);
clim = [0 max_spctrm];

%% Plot LOW-SPL TFR
close all
figure
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr_low);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4); % 0–2 s × 8–14 Hz
set(gca, 'FontSize', fontSize);
title('Sternberg TFR — LOW SPL trials');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_LOW_SPL_trials.png');

%% Plot HIGH-SPL TFR
figure
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfg, gatfr_high);
colormap(color_map);
set(gca, 'CLim', clim);
cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4);
set(gca, 'FontSize', fontSize);
title('Sternberg TFR — HIGH SPL trials');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_HIGH_SPL_trials.png');

%% Done
disp('Completed LOW/HIGH SPL median-split TFRs.')