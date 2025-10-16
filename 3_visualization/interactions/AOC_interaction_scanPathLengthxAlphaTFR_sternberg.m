%% AOC Sternberg — Median-split by within-subject Scan Path Length to compare with Alpha TFRs

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat') % FieldTrip layout
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

colors = color_def('AOC');
fontSize = 36;

% Common reference grid for gaze step-length series
fs = 500; % Hz
t_series = -0.5:1/fs:2; % reference grid
T = numel(t_series) - 1; % step series align to t_series(2:end)

% Preallocate
low_tfr_subs = cell(1, length(subjects)); % per-subject LOW-SPL TFR (avg over trials)
high_tfr_subs = cell(1, length(subjects)); % per-subject HIGH-SPL TFR (avg over trials)
scan_low = nan(length(subjects), T); % per-subject LOW-SPL scan-path series (avg over trials)
scan_high = nan(length(subjects), T); % per-subject HIGH-SPL scan-path series (avg over trials)

%% Per-subject split and aggregation
for s = 1:10%%%%%length(subjects)
subjID = str2double(subjects{s});
fprintf('Subject %s (%d/%d)\n', subjects{s}, s, length(subjects));

% Total Scan Path Length per trial = Early + Late (row-wise, ignore NaNs)
rows = merged_data_sternberg_trials(merged_data_sternberg_trials.ID == subjID, :);
if isempty(rows)
    warning('No merged rows for subject %s. Skipping.', subjects{s})
    continue
end

% Row-wise sum with NaN omission (requires numeric matrix)
spl_mat = [rows.ScanPathLengthEarly, rows.ScanPathLengthLate];
total_spl = sum(spl_mat, 2, 'omitnan');     % px
trlN      = rows.Trial;

good = isfinite(total_spl) & isfinite(trlN);
if ~any(good)
    warning('No finite total SPL for subject %s. Skipping subject.', subjects{s})
    continue
end

% Random tie-breaking median split (within-subject)
spl_sub  = total_spl(good);
trl_sub  = trlN(good);

rp        = randperm(numel(spl_sub));      % random permutation to break ties stably
spl_rand  = spl_sub(rp);
trl_rand  = trl_sub(rp);

[spl_sorted, idx_sorted] = sort(spl_rand, 'ascend'); %#ok<NASGU>
trl_sorted = trl_rand(idx_sorted);

nHalf      = floor(numel(trl_sorted)/2);
lowTrials  = trl_sorted(1:nHalf);
highTrials = trl_sorted(nHalf+1:end);

% EEG TFR — average LOW/HIGH SPL trials over repetitions (occipital channels)
tfr_all = [];
try
    datapath_eeg = fullfile(path, subjects{s}, 'eeg');
    cd(datapath_eeg)
    load tfr_stern_trials   % -> tfr_all
catch
    warning('Missing EEG TFR for subject %s, skipping EEG part.', subjects{s})
end

if ~isempty(tfr_all)
    % Occipital set by label name (contains 'O' or 'I'); fallback to all
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
        warning('tfr_all.trialinfo missing Trial column for subject %s. Skipping EEG part.', subjects{s})
    else
        eegTrials = tfr_all.trialinfo(:,2);

        idxLow  = ismember(eegTrials,  lowTrials);
        idxHigh = ismember(eegTrials, highTrials);

        if any(idxLow)
            cfgS = [];
            cfgS.trials     = find(idxLow);
            cfgS.channel    = occ_channels;
            cfgS.avgoverrpt = 'yes';
            low_tfr  = ft_selectdata(cfgS, tfr_all);
        else
            low_tfr = [];
        end

        if any(idxHigh)
            cfgS = [];
            cfgS.trials     = find(idxHigh);
            cfgS.channel    = occ_channels;
            cfgS.avgoverrpt = 'yes';
            high_tfr = ft_selectdata(cfgS, tfr_all);
        else
            high_tfr = [];
        end

        low_tfr_subs{s}  = low_tfr;
        high_tfr_subs{s} = high_tfr;
    end
end

% Gaze scan-path series — average LOW/HIGH SPL trials
ScanPathSeries  = {};
ScanPathSeriesT = {};
trialinfo = [];
try
    datapath_gaze = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features', subjects{s}, 'gaze', 'gaze_series_sternberg_trials.mat');
    load(datapath_gaze, 'ScanPathSeries', 'ScanPathSeriesT', 'trialinfo')
catch
    warning('Missing gaze series for subject %s, skipping gaze part.', subjects{s})
end

if ~isempty(ScanPathSeries)
    if size(trialinfo,2) < 2
        warning('Unexpected gaze trialinfo shape for subject %s. Skipping gaze part.', subjects{s})
    else
        % Interpolate each trial to the common grid
        subj_trials = nan(numel(ScanPathSeries), T);
        for trl = 1:numel(ScanPathSeries)
            srl = ScanPathSeries{trl};
            tt  = ScanPathSeriesT{trl};
            if isempty(srl) || numel(tt) ~= numel(srl)
                continue
            end
            try
                subj_trials(trl,:) = interp1(tt, srl, t_series(2:end), 'linear', NaN);
            catch
                % leave as NaN
            end
        end

        gazeTrials = trialinfo(:,2);  % column 2 = Trial
        lowMask  = ismember(gazeTrials,  lowTrials);
        highMask = ismember(gazeTrials, highTrials);

        if any(lowMask)
            scan_low(s,:)  = nanmean(subj_trials(lowMask,:), 1);
        end
        if any(highMask)
            scan_high(s,:) = nanmean(subj_trials(highMask,:), 1);
        end
    end
end


end

%% Grand-average TFRs (LOW vs HIGH SPL)
low_tfr_subs = low_tfr_subs(~cellfun(@isempty, low_tfr_subs));
high_tfr_subs = high_tfr_subs(~cellfun(@isempty, high_tfr_subs));

if isempty(low_tfr_subs) || isempty(high_tfr_subs)
error('No subject TFRs available for grand average.')
end

gatfr_low = ft_freqgrandaverage([], low_tfr_subs{:});
gatfr_high = ft_freqgrandaverage([], high_tfr_subs{:});

% Plot config (FieldTrip singleplotTFR)
cfg = [];
cfg.channel = gatfr_low.label; % already restricted to occipital in per-subject step
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [0 2];
cfg.ylim = [4 30];
cfg.layout = layANThead;

color_map = cbrewer('seq', 'Reds', 64);

% Common colour scaling from alpha band (8–14 Hz) and 0–2 s
alpha_idx = gatfr_low.freq >= 8 & gatfr_low.freq <= 14;
time_idx = gatfr_low.time >= 0 & gatfr_low.time <= 2;

[~, ch_low_idx] = ismember(cfg.channel, gatfr_low.label);
[~, ch_high_idx] = ismember(cfg.channel, gatfr_high.label);

low_alpha_power = mean(gatfr_low.powspctrm(ch_low_idx, alpha_idx, time_idx), 1:3, 'omitnan');
high_alpha_power = mean(gatfr_high.powspctrm(ch_high_idx, alpha_idx, time_idx), 1:3, 'omitnan');
max_spctrm = max([low_alpha_power(:); high_alpha_power(:)]);
clim = [0 max_spctrm];

%% Plot TFR — LOW SPL trials
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
title('Sternberg TFR — LOW Scan Path Length trials');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_TFR_LOWscanpath_trials.png');

%% Plot TFR — HIGH SPL trials
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
title('Sternberg TFR — HIGH Scan Path Length trials');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_TFR_HIGHscanpath_trials.png');

%% Grand-average Scan Path Length series (LOW vs HIGH SPL) with SEM
grand_low = nanmean(scan_low, 1);
grand_high = nanmean(scan_high, 1);
sem_low = nanstd(scan_low, [], 1) ./ sqrt(sum(isfinite(scan_low), 1));
sem_high = nanstd(scan_high, [], 1) ./ sqrt(sum(isfinite(scan_high), 1));
t_plot = t_series(2:end);

figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200]); hold on
shadedErrorBar(t_plot, grand_low, sem_low, 'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
shadedErrorBar(t_plot, grand_high, sem_high, 'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);
xlabel('Time (s)')
ylabel('Scan Path Length (px)')
title('Sternberg — Scan Path Length over time (LOW vs HIGH total SPL trials)')
xlim([t_series(1) t_series(end)])
box on
set(gca, 'FontSize', 25)
legend({'LOW total SPL ± SEM','HIGH total SPL ± SEM'}, 'Location','northwest')
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_scanPathSeries_LOWvsHIGHtotalSPL.png')

%% Done
disp('Completed LOW/HIGH Scan Path Length median-split TFRs and scan path plots.')