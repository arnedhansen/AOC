%% AOC Sternberg — Median-split by within-subject Scan Path Length (SPL) to compare TFRs & SPL time-course

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
%subjects = subjects(1:20)

load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat') % FieldTrip layout
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

colors = color_def('AOC');
fontSize = 36;

% Plot colormap
color_map = cbrewer('seq', 'Reds', 64);

% Common reference grid for gaze step-length series
t_series = linspace(-0.5, 2, 51); % 51 points -> 50 steps
T = numel(t_series) - 1; % step series aligns to t_series(2:end)

%% Preallocate holders for subject-level TFRs and SPL time-courses
low_tfr_subs = cell(1, length(subjects)); % per-subject LOW-SPL TFR (avg over trials)
high_tfr_subs = cell(1, length(subjects)); % per-subject HIGH-SPL TFR (avg over trials)
scan_low = nan(length(subjects), T); % per-subject LOW-SPL scan-path series (avg over trials)
scan_high = nan(length(subjects), T); % per-subject HIGH-SPL scan-path series (avg over trials)

%% Per-subject split (by SPL) and aggregation of EEG TFRs + SPL series
for s = 1:length(subjects)
    clc
    subjID = str2double(subjects{s});
    fprintf('Subject %s (%d/%d)\n', subjects{s}, s, length(subjects));

    % Get per-trial SPL and trial numbers for this subject
    rows = merged_data_sternberg_trials(merged_data_sternberg_trials.ID == subjID, :);
    spl  = rows.ScanPathLengthFull;  % total scan path length per trial (Full window)
    trlN = rows.Trial;

    good = isfinite(spl) & isfinite(trlN) & spl < 1200;
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

    [~, idx_sorted] = sort(spl_rand, 'ascend');
    trl_sorted = trl_rand(idx_sorted);

    nHalf      = floor(numel(trl_sorted)/2);
    lowTrials  = trl_sorted(1:nHalf);
    highTrials = trl_sorted(nHalf+1:end);

    %% Check median split
    % close all
    % figure; set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200]); 
    % hold on
    % [spl_sorted, ord] = sort(spl_sub, 'ascend');
    % x = 1:numel(spl_sorted);
    % plot(x, spl_sorted, '-', 'Color', 0.85*[1 1 1], 'LineWidth', 3)
    % scatter(x(1:nHalf), spl_sorted(1:nHalf), 60, colors(1,:), 'filled')
    % scatter(x(nHalf+1:end), spl_sorted(nHalf+1:end), 60, colors(3,:), 'filled')
    % yl = ylim; plot([nHalf+0.5 nHalf+0.5], yl, '--k', 'LineWidth', 2)
    % m = median(spl_sub, 'omitnan'); 
    % plot([1 numel(spl_sorted)], [m m], ':k', 'LineWidth', 2)
    % xlabel('Trials (ranked by SPL)'); ylabel('Scan Path Length (Full) [px]')
    % title(sprintf('Subject %s — within-subject median split (LOW=%d, HIGH=%d)', subjects{s}, nHalf, numel(spl_sorted)-nHalf))
    % set(gca, 'FontSize', fontSize-6), box on
    % saveas(gcf, fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/subjects/', ['AOC_sternberg_SPL_binned_medianSplit_' subjects{s} '.png']))

    %% EEG: TFR LOW/HIGH-SPL
    % ----------------------
    tfr_all = [];
    try
        datapath_eeg = fullfile(path, subjects{s}, 'eeg');
        cd(datapath_eeg)
        load tfr_stern_trials   % -> tfr_all (dimord 'rpt_chan_freq_time'), trialinfo(:,2) = Trial
    catch
        warning('Missing EEG TFR for subject %s, skipping EEG part.', subjects{s})
    end

    if ~isempty(tfr_all)
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

    % -----------------------------------
    % Gaze: SPL time-series LOW/HIGH-SPL
    % -----------------------------------
    ScanPathSeriesBins = {};
    ScanPathSeriesT    = {};
    trialinfo          = [];
    try
        datapath_gaze = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features', subjects{s}, 'gaze', 'gaze_series_sternberg_trials.mat');
        load(datapath_gaze, 'ScanPathSeriesBins', 'ScanPathSeriesT', 'trialinfo')
    catch
        warning('Missing gaze series for subject %s, skipping gaze part.', subjects{s})
    end

    if ~isempty(ScanPathSeriesBins)
        if size(trialinfo,2) < 2
            warning('Unexpected gaze trialinfo shape for subject %s. Skipping gaze part.', subjects{s})
        else
            % Interpolate each trial to the common grid of step times (t_series(2:end))
            subj_trials = nan(numel(ScanPathSeriesBins), T);
            ScanPathSeriesT = linspace(-0.5,2,50);
            for trl = 1:numel(ScanPathSeriesBins)
                srl = ScanPathSeriesBins{trl};
                tt  = ScanPathSeriesT;
                if isempty(srl) || isempty(tt) || numel(tt) ~= numel(srl)
                    continue
                end
                try
                    subj_trials(trl,:) = interp1(tt, srl, t_series(2:end), 'linear', NaN);
                catch
                    % leave as NaN
                end
            end

            gazeTrials = trialinfo(:,2);
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
    error('No subject TFRs available for grand average (LOW/HIGH SPL).')
end

gatfr_low = ft_freqgrandaverage([], low_tfr_subs{:});
gatfr_high = ft_freqgrandaverage([], high_tfr_subs{:});

%% Plot Alpha Power TFRs
cfg = [];
cfg.channel = gatfr_low.label; % already restricted per subject to occipital
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
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
max_spctrm = 4.15
clim = [0 max_spctrm];

% Plot LOW-SPL TFR
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

% Plot HIGH-SPL TFR
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

%% Grand-average Scan Path Length time-course (LOW vs HIGH SPL) with SEM
grand_low = nanmean(scan_low, 1);
grand_high = nanmean(scan_high, 1);
sem_low = nanstd(scan_low, [], 1) ./ sqrt(sum(isfinite(scan_low), 1));
sem_high = nanstd(scan_high, [], 1) ./ sqrt(sum(isfinite(scan_high), 1));
t_plot = t_series(2:end);

figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200]); hold on
shadedErrorBar(t_plot, grand_low, sem_low, 'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
shadedErrorBar(t_plot, grand_high, sem_high, 'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);
xlabel('Time [s]')
ylabel('Scan Path Length [px]')
title('Sternberg — Scan Path Length over time (LOW vs HIGH SPL trials)')
xlim([t_series(1) t_series(end)])
box on
set(gca, 'FontSize', 25)
legend({'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northwest')
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_scanPathLength_LOWvsHIGH_SPL.png')

%% Done
disp('Completed LOW/HIGH SPL median-split TFRs and SPL time-course plots.')