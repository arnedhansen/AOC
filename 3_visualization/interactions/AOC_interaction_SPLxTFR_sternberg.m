%% AOC Sternberg — Median-split by within-subject Scan Path Length (SPL) to compare TFRs & SPL time-course

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
subjects = subjects(1:10)

% aesthetics
color_map = cbrewer('seq', 'Reds', 64);
colors = color_def('AOC');
fontSize = 36;

% Common reference grid for gaze step-length series
time_series = linspace(-0.5, 2, 51); % 51 points -> 50 steps
T = numel(time_series) - 1; % step series aligns to time_series(2:end)

%% Load data
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

%% Preallocate holders for subject-level TFRs and SPL time-courses
low_tfr_subs = cell(1, length(subjects)); % per-subject LOW-SPL TFR (avg over trials)
high_tfr_subs = cell(1, length(subjects)); % per-subject HIGH-SPL TFR (avg over trials)
scan_low = nan(length(subjects), T); % per-subject LOW-SPL scan-path series (avg over trials)
scan_high = nan(length(subjects), T); % per-subject HIGH-SPL scan-path series (avg over trials)

%% Per-subject split (by SPL) and aggregation of EEG TFRs + SPL series
for subj = 1:length(subjects)
    clc
    subjID = str2double(subjects{subj});
    fprintf('Subject %s (%d/%d)\n', subjects{subj}, subj, length(subjects));

    % Get per-trial SPL and trial numbers for this subject
    disp('Extracting SPL')
    subjTable = merged_data_sternberg_trials(merged_data_sternberg_trials.ID == subjID, :);
    spl  = subjTable.ScanPathLengthFull;  % total scan path length per trial (Full window)
    trlN = subjTable.Trial;

    good_idx = isfinite(spl) & isfinite(trlN) & spl < 1200;
    if ~any(good_idx)
        warning('No finite SPL for subject %subj. Skipping subject.', subjects{subj})
        continue
    end

    % Random tie-breaking median split (within-subject, by SPL)
    spl_sub  = spl(good_idx);
    trl_sub  = trlN(good_idx);
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
    % title(sprintf('Subject %subj — within-subject median split (LOW=%d, HIGH=%d)', subjects{subj}, nHalf, numel(spl_sorted)-nHalf))
    % set(gca, 'FontSize', fontSize-6), box on
    % saveas(gcf, fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/subjects/', ['AOC_sternberg_SPL_binned_medianSplit_' subjects{subj} '.png']))

    %% EEG: TFR LOW/HIGH-SPL
    
    tfr_all = [];
    disp('Extracting EEG: TFR LOW/HIGH-SPL')
    try
        datapath_eeg = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath_eeg)
        load tfr_stern_trials   % tfr_all (dimord 'rpt_chan_freq_time'), trialinfo(:,2) = Trial
    catch
        warning('Missing EEG TFR for subject %subj, skipping EEG part.', subjects{subj})
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
            warning('tfr_all.trialinfo missing Trial column for subject %subj. Skipping EEG part.', subjects{subj})
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

            low_tfr_subs{subj}  = low_tfr;
            high_tfr_subs{subj} = high_tfr;
        end
    end

    % Gaze: SPL time-series LOW/HIGH-SPL
    ScanPathSeriesBins = {};
    ScanPathSeriesT    = {};
    trialinfo          = [];
    disp('Extracting Gaze: SPL time-series LOW/HIGH-SPL')
    try
        datapath_gaze = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features', subjects{subj}, 'gaze', 'gaze_series_sternberg_trials.mat');
        load(datapath_gaze, 'ScanPathSeriesBins', 'ScanPathSeriesT', 'trialinfo')
    catch
        warning('Missing gaze series for subject %subj, skipping gaze part.', subjects{subj})
    end

    if ~isempty(ScanPathSeriesBins)
        if size(trialinfo,2) < 2
            warning('Unexpected gaze trialinfo shape for subject %subj. Skipping gaze part.', subjects{subj})
        else
            % Interpolate each trial to the common grid of step times (time_series(2:end))
            subj_trials = nan(numel(ScanPathSeriesBins), T);
            ScanPathSeriesT = linspace(-0.5,2,50);
            for trl = 1:numel(ScanPathSeriesBins)
                srl = ScanPathSeriesBins{trl};
                tt  = ScanPathSeriesT;
                if isempty(srl) || isempty(tt) || numel(tt) ~= numel(srl)
                    continue
                end
                try
                    subj_trials(trl,:) = interp1(tt, srl, time_series(2:end), 'linear', NaN);
                catch
                    % leave as NaN
                end
            end

            gazeTrials = trialinfo(:,2);
            lowMask  = ismember(gazeTrials,  lowTrials);
            highMask = ismember(gazeTrials, highTrials);

            if any(lowMask)
                scan_low(subj,:)  = nanmean(subj_trials(lowMask,:), 1);
            end
            if any(highMask)
                scan_high(subj,:) = nanmean(subj_trials(highMask,:), 1);
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
cfg.channel = gatfr_low.label;
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
cfg.xlim = [-.5 2];
cfg.ylim = [4 30];
load('/Volumes/g_psyplafor_methlab$/Students/Arne/MA/headmodel/layANThead.mat')
cfg.layout = layANThead;

% Maxspctrm for plotting
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
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4); % 0–2 subj × 8–14 Hz
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

% Plot DIFF TFR
gatfr_diff = gatfr_high;
gatfr_diff.powspctrm = gatfr_high.powspctrm - gatfr_low.powspctrm;
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
title('Sternberg TFR — DIFF (HIGH SPL - LOW SPL trials)');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_DIFF_SPL_trials.png');

%% TFR DIFF STATS (CBPT)
N = numel(low_tfr_subs);
assert(N == numel(high_tfr_subs), 'LOW/HIGH cell arrays must be same length.');

% Neighbours (EEG): derive from electrode positions
cfg_n              = [];
cfg_n.method       = 'distance';          % safe default if you have .elec
cfg_n.neighbourdist= 0.12;                % ~12 cm; tweak to your cap geometry
cfg_n.elec         = low_tfr_subs{1}.elec;
neighbours         = ft_prepare_neighbours(cfg_n);

% Design (within-subject)
design = zeros(2, 2*N);
design(1, :) = [1:N,           1:N          ]; % subject id
design(2, :) = [ones(1, N),    2*ones(1, N) ]; % condition code
cfgs               = [];
cfgs.parameter     = 'powspctrm';
cfgs.method        = 'montecarlo';
cfgs.statistic     = 'ft_statfun_depsamplesT';
cfgs.correctm      = 'cluster';
cfgs.clusteralpha  = 0.05;
cfgs.clusterstatistic = 'maxsum';
cfgs.minnbchan     = 2;                    % require ≥2 neighbouring chans in a cluster
cfgs.neighbours    = neighbours;
cfgs.tail          = 0;                    % two-sided
cfgs.clustertail   = 0;
cfgs.alpha         = 0.025;                % report threshold (two-sided)
cfgs.numrandomization = 2000;              % increase if you want more precision
cfgs.design        = design;
cfgs.uvar          = 1;                    % unit of observation = subject
cfgs.ivar          = 2;                    % independent variable = condition
cfgs.channel       = 'all';                % or e.g., occipital ROI
cfgs.latency       = [ 0  2 ];           % optionally restrict time window
cfgs.frequency     = [ 4 30 ];            % optionally restrict frequency range

% Run stats across the full TFR (space × freq × time)
stat_tfr = ft_freqstatistics(cfgs, low_tfr_subs{:}, high_tfr_subs{:});

% Grand averages for plotting (and DIFF)
cfgGA          = [];
cfgGA.parameter= 'powspctrm';
ga_low         = ft_freqgrandaverage(cfgGA,  low_tfr_subs{:});
ga_high        = ft_freqgrandaverage(cfgGA, high_tfr_subs{:});
ga_diff        = ga_high; % copy metadata
ga_diff.powspctrm = ga_high.powspctrm - ga_low.powspctrm;

% Attach mask from stats (significant clusters)
ga_diff.mask   = stat_tfr.mask; % logical mask time×freq×chan re-ordered internally by FT

% Plot masked DIFF TFR (average over chosen channels)
cfgp                  = [];
cfgp.parameter        = 'powspctrm';
cfgp.maskparameter    = 'mask';
cfgp.maskstyle        = 'outline';
cfgp.zlim             = clim;           % reuse your colour limits
cfgp.colormap         = color_map;
cfgp.channel          = 'all';          % or a ROI like {'O1','Oz','O2','POz',...}
cfgp.figure           = 'gcf';

close all
figure
set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
ft_singleplotTFR(cfgp, ga_diff);
cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]'); ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4); % 0–2 s × 8–14 Hz
set(gca, 'FontSize', fontSize);
title('Sternberg TFR — DIFF (HIGH - LOW), significant clusters outlined')
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_DIFF_SPL_trials_STATS.png');

%% Grand-average Scan Path Length time-course (LOW vs HIGH SPL)
grand_low = nanmean(scan_low, 1);
grand_high = nanmean(scan_high, 1);
sem_low = nanstd(scan_low, [], 1) ./ sqrt(sum(isfinite(scan_low), 1));
sem_high = nanstd(scan_high, [], 1) ./ sqrt(sum(isfinite(scan_high), 1));
t_plot = time_series(2:end);

figure;
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200]); hold on
shadedErrorBar(t_plot, grand_low, sem_low, 'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
shadedErrorBar(t_plot, grand_high, sem_high, 'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);
xlabel('Time [s]')
ylabel('Scan Path Length [px]')
title('Sternberg — Scan Path Length over time (LOW vs HIGH SPL trials)')
xlim([time_series(1) time_series(end)])
box on
set(gca, 'FontSize', 25)
legend({'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northwest')

% Paired t-tests (HIGH vs LOW) with FDR correction using mafdr
N = size(scan_low, 1);
T = size(scan_low, 2);

pvals = nan(1, T);
for t = 1:T
    [~, pvals(t)] = ttest(scan_high(:,t), scan_low(:,t));
end

% FDR correction (Benjamini–Hochberg)
qvals = mafdr(pvals, 'BHFDR', true);
sig_mask = qvals < 0.05;

% Add significance bar
ax = gca;
yl = ylim(ax);
yr = yl(2) - yl(1);
ybar = yl(1) + 0.03*yr;

d_sig   = diff([0, sig_mask, 0]);
on_sig  = find(d_sig == 1);
off_sig = find(d_sig == -1) - 1;

hold on
for k = 1:numel(on_sig)
    x0 = t_plot(on_sig(k));
    x1 = t_plot(off_sig(k));
    plot([x0 x1], [ybar ybar], 'k-', 'LineWidth', 12)
end

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_scanPathLength_LOWvsHIGH_SPL.png')

%% SPL and Alpha Power over Time
alpha_band = [8 14]; % Hz
cfg = [];
cfg.frequency = alpha_band;
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
gatfr_low_alpha = ft_selectdata(cfg, gatfr_low);
gatfr_high_alpha = ft_selectdata(cfg, gatfr_high);

%%
fontSize = 30;
close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200])

% Top subplot — Scan Path Length
subplot(3,1,1:2)
hold on
shadedErrorBar(t_plot, grand_low, sem_low, ...
    'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
shadedErrorBar(t_plot, grand_high, sem_high, ...
    'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);

% Mark significant intervals
for k = 1:numel(on_sig)
    x0 = t_plot(on_sig(k));
    x1 = t_plot(off_sig(k));
    plot([x0 x1], [ybar ybar], 'k-', 'LineWidth', 12)
end

ylabel('Scan Path Length [px]')
title('Sternberg — Scan Path Length and Alpha Power over time (LOW vs HIGH SPL trials)')
xlim([time_series(1) time_series(end)])
set(gca, 'FontSize', fontSize)
box on
legend({'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northwest', 'FontSize', fontSize*0.6)
set(gca, 'XColor', 'none')  % remove x-axis

% Bottom subplot — Alpha Power
subplot(3,1,3)
hold on
yyaxis right   % activate right axis

plot(gatfr_low_alpha.time, squeeze(gatfr_low_alpha.powspctrm), ...
    'LineWidth', 5, 'Color', colors(1,:), 'LineStyle', '--');
plot(gatfr_high_alpha.time, squeeze(gatfr_high_alpha.powspctrm), ...
    'LineWidth', 5, 'Color', colors(3,:), 'LineStyle', '--');

set(gca, 'YColor', [0 0 0])          % colour of right y-axis
ylabel('Alpha Power [\muV^2/Hz]')
xlabel('Time [s]')
xlim([-.5 2])
set(gca, 'FontSize', fontSize)
legend({'Alpha Power of LOW SPL Trials','Alpha Power of HIGH SPL Trials'}, 'Location','northwest', 'FontSize', fontSize*0.6)

% hide the left y-axis entirely
yyaxis left
set(gca, 'YColor', 'none', 'YTick', [])

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_SPLOverTime_AlphaPower_Subplots.png')


%% Done
disp('Completed LOW/HIGH SPL median-split TFRs and SPL time-course plots.')