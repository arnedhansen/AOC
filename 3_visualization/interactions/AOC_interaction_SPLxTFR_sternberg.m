%% AOC Sternberg — Median-split by within-subject Scan Path Length (SPL) to compare TFRs & SPL time-course

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
%subjects = subjects(1:10)

% aesthetics
color_map = cbrewer('seq', 'Reds', 64);
colors = color_def('AOC');
fontSize = 30;

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

% Full-resolution holders for stats (common full grid defined below)
fs_full     = 500;                                  % matches your gaze pipeline
t_full      = -0.5:1/fs_full:2;                     % sample grid incl. endpoints
t_plot_full = t_full(2:end);                        % step series aligns to t_full(2:end)
Tf          = numel(t_plot_full);                   % full-resolution length
scan_low_full  = nan(length(subjects), Tf);         % per-subject LOW-SPL full series
scan_high_full = nan(length(subjects), Tf);         % per-subject HIGH-SPL full series

datTS_ID    = [];
datTS_Trial = [];
datTS_Cond  = [];
datTS_tidx  = [];
datTS_t     = [];
datTS_SPL   = [];

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

    % Exclude NaNs
    good_idx = isfinite(spl) & isfinite(trlN);

    % Compute subject-specific z-scores and remove trials beyond 3 SD
    z_spl = (spl - nanmean(spl)) ./ nanstd(spl);
    good_idx = good_idx & abs(z_spl) <= 3;

    if ~any(good_idx)
        warning('No valid SPL for subject %s. Skipping subject.', subjects{subj})
        continue
    end

    % Random tie-breaking median split (within-subject, by SPL)
    rng(subjID, 'twister');
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
            error('ERROR IN CHANNEL SELECTION')
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
        load(datapath_gaze, 'ScanPathSeries', 'ScanPathSeriesBins', 'ScanPathSeriesT', 'trialinfo')
    catch
        warning('Missing gaze series for subject %subj, skipping gaze part.', subjects{subj})
    end

    if ~isempty(ScanPathSeriesBins)
        if size(trialinfo,2) < 2
            warning('Unexpected gaze trialinfo shape for subject %subj. Skipping gaze part.', subjects{subj})
        else
            % Interpolate each trial to the common grid of step times (time_series(2:end))
            subj_trials = nan(numel(ScanPathSeriesBins), T);
            ScanPathSeriesTBins = linspace(-0.5,2,50);
            for trl = 1:numel(ScanPathSeriesBins)
                srl = ScanPathSeriesBins{trl};
                tt  = ScanPathSeriesTBins;
                if isempty(srl) || isempty(tt) || numel(tt) ~= numel(srl)
                    continue
                end
                try
                    subj_trials(trl,:) = interp1(tt, srl, time_series(2:end), 'linear', NaN);
                catch
                    % leave as NaN
                end
            end

            gazeTrials = trialinfo(2,:);
            lowMask  = ismember(gazeTrials,  lowTrials);
            highMask = ismember(gazeTrials, highTrials);

            if any(lowMask)
                scan_low(subj,:)  = nanmean(subj_trials(lowMask,:), 1);
            end
            if any(highMask)
                scan_high(subj,:) = nanmean(subj_trials(highMask,:), 1);
            end
            % Interpolate each trial also to the common FULL grid of step times (t_full(2:end))
            subj_trials_full = nan(numel(ScanPathSeries), Tf);

            for trl = 1:numel(ScanPathSeries)
                srl_full = ScanPathSeries{trl}; % full-resolution step-length series for this trial
                tt_full  = ScanPathSeriesT{trl}; % its time vector (same length as srl_full+1 samples)

                if isempty(srl_full) || isempty(tt_full)
                    continue
                end

                % step series align to time vector
                try
                    subj_trials_full(trl, :) = interp1(tt_full, srl_full, t_plot_full, 'linear', NaN);
                catch
                    % leave as NaN
                end
            end

            % Use the same lowMask/highMask you already computed
            if any(lowMask)
                scan_low_full(subj, :)  = nanmean(subj_trials_full(lowMask, :), 1);
            end
            if any(highMask)
                scan_high_full(subj, :) = nanmean(subj_trials_full(highMask, :), 1);
            end
        end
    end

    % Collect LOW trials into long table
    if any(lowMask)
        trl_idx = find(lowMask);
        for r = 1:numel(trl_idx)
            tr = trl_idx(r);
            y  = subj_trials_full(tr, :);                % 1 x Tf
            good = isfinite(y);
            if any(good)
                n  = nnz(good);
                datTS_ID    = [datTS_ID;   repmat(subjID, n, 1)];
                datTS_Trial = [datTS_Trial; repmat(gazeTrials(tr), n, 1)];
                datTS_Cond  = [datTS_Cond;  repmat({'LOW'}, n, 1)];
                datTS_tidx  = [datTS_tidx;  find(good)'];
                datTS_t     = [datTS_t;     t_plot_full(good)'];
                datTS_SPL   = [datTS_SPL;   y(good)'];
            end
        end
    end

    % Collect HIGH trials into long table
    if any(highMask)
        trl_idx = find(highMask);
        for r = 1:numel(trl_idx)
            tr = trl_idx(r);
            y  = subj_trials_full(tr, :);
            good = isfinite(y);
            if any(good)
                n  = nnz(good);
                datTS_ID    = [datTS_ID;   repmat(subjID, n, 1)];
                datTS_Trial = [datTS_Trial; repmat(gazeTrials(tr), n, 1)];
                datTS_Cond  = [datTS_Cond;  repmat({'HIGH'}, n, 1)];
                datTS_tidx  = [datTS_tidx;  find(good)'];
                datTS_t     = [datTS_t;     t_plot_full(good)'];
                datTS_SPL   = [datTS_SPL;   y(good)'];
            end
        end
    end

end
datTS           = table;
datTS.ID        = datTS_ID;
datTS.Trial     = datTS_Trial;
datTS.Condition = categorical(datTS_Cond, {'LOW','HIGH'});  % set LOW as reference
datTS.t_index   = datTS_tidx;                               % 1..Tf
datTS.Time      = datTS_t;                                   % seconds
datTS.SPL       = datTS_SPL;

%% Grand-average TFRs (LOW vs HIGH SPL)
low_tfr_subs = low_tfr_subs(~cellfun(@isempty, low_tfr_subs));
high_tfr_subs = high_tfr_subs(~cellfun(@isempty, high_tfr_subs));

if isempty(low_tfr_subs) || isempty(high_tfr_subs)
    error('No subject TFRs available for grand average (LOW/HIGH SPL).')
end

gatfr_low = ft_freqgrandaverage([], low_tfr_subs{:});
gatfr_high = ft_freqgrandaverage([], high_tfr_subs{:});

%% Save VARIABLES
save '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_interaction_SPLxTFR_sternberg' ...
     gatfr_low gatfr_high scan_low scan_high scan_low_full scan_high_full time_series datTS t_plot_full

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
ft_singleplotTFR(cfg, gatfr_diff);
color_mapBlRd = customcolormap_preset('red-white-blue');
colormap(color_mapBlRd);
clim = ([-.5 .5]);
set(gca, 'CLim', clim);
cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4);
set(gca, 'FontSize', fontSize);
title('Sternberg TFR — DIFF (HIGH SPL - LOW SPL trials)');
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_DIFF_SPL_trials.png');

%% TFR DIFF STATS (CBPT)
% N = numel(low_tfr_subs);
% assert(N == numel(high_tfr_subs), 'LOW/HIGH cell arrays must be same length.');
% 
% % Neighbours (EEG): derive from electrode positions
% cfg_n               = [];
% cfg_n.method        = 'distance';          % safe default if you have .elec
% cfg_n.neighbourdist = 0.035;            % ~3.5 cm for ANT 128
% cfg_n.elec          = low_tfr_subs{1}.elec;
% neighbours          = ft_prepare_neighbours(cfg_n);
% 
% % Design (within-subject)
% design = zeros(2, 2*N);
% design(1, :) = [1:N,           1:N          ]; % subject id
% design(2, :) = [ones(1, N),    2*ones(1, N) ]; % condition code
% cfgs               = [];
% cfgs.parameter     = 'powspctrm';
% cfgs.method        = 'montecarlo';
% cfgs.statistic     = 'ft_statfun_depsamplesT';
% cfgs.correctm      = 'cluster';
% cfgs.clusteralpha  = 0.05;
% cfgs.clusterstatistic = 'maxsum';
% cfgs.minnbchan     = 2;                    % require ≥2 neighbouring chans in a cluster
% cfgs.neighbours    = neighbours;
% cfgs.tail          = 0;                    % two-sided
% cfgs.clustertail   = 0;
% cfgs.alpha         = 0.025;                % report threshold (two-sided)
% cfgs.numrandomization = 2000;              % increase if you want more precision
% cfgs.design        = design;
% cfgs.uvar          = 1;                    % unit of observation = subject
% cfgs.ivar          = 2;                    % independent variable = condition
% cfgs.channel       = 'all';                % or e.g., occipital ROI
% cfgs.latency       = [ 0  2 ];           % optionally restrict time window
% cfgs.frequency     = [ 4 30 ];            % optionally restrict frequency range
% 
% % Run stats across the full TFR (space × freq × time)
% stat_tfr = ft_freqstatistics(cfgs, low_tfr_subs{:}, high_tfr_subs{:});
% 
% % Grand averages for plotting (and DIFF)
% cfgGA          = [];
% cfgGA.parameter= 'powspctrm';
% ga_low         = ft_freqgrandaverage(cfgGA,  low_tfr_subs{:});
% ga_high        = ft_freqgrandaverage(cfgGA, high_tfr_subs{:});
% ga_diff        = ga_high; % copy metadata
% ga_diff.powspctrm = ga_high.powspctrm - ga_low.powspctrm;
% ga_diff.powspctrm = ga_diff.powspctrm(1:24, 1:27, find(ga_diff.time ==0):find(ga_diff.time ==2));
% 
% % Attach mask from stats (significant clusters)
% ga_diff.mask   = stat_tfr.mask; % logical mask time×freq×chan re-ordered internally by FT
% 
% % Plot masked DIFF TFR (average over chosen channels)
% cfgp                  = [];
% cfgp.parameter        = 'powspctrm';
% cfgp.maskparameter    = 'mask';
% cfgp.maskstyle        = 'outline';
% cfgp.zlim             = clim;           % reuse your colour limits
% cfgp.colormap         = color_map;
% cfgp.channel          = 'all';          % or a ROI like {'O1','Oz','O2','POz',...}
% cfgp.figure           = 'gcf';
% 
% close all
% figure
% set(gcf, 'Position', [100, 200, 2000, 1200], 'Color', 'w');
% ft_singleplotTFR(cfgp, ga_diff);
% cb = colorbar; ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', fontSize);
% xlabel('Time [s]'); ylabel('Frequency [Hz]');
% rectangle('Position', [0, 8, 2, 6], 'EdgeColor', 'k', 'LineWidth', 4); % 0–2 s × 8–14 Hz
% set(gca, 'FontSize', fontSize);
% title('Sternberg TFR — DIFF (HIGH - LOW), significant clusters outlined')
% saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_TFR_DIFF_SPL_trials_STATS.png');

%% Scan Path Length time-course (LOW vs HIGH SPL) TTEST
close all
clc
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

% Paired t-tests (HIGH vs LOW) with FDR correction using mafdr
N  = size(scan_low_full, 1);         % subjects
Tf = size(scan_low_full, 2);         % full-resolution time points

pvals = nan(1, Tf);
for t = 1:Tf
    [~, pvals(t)] = ttest(scan_high_full(:, t), scan_low_full(:, t));
end

% FDR correction (Benjamini–Hochberg) on full-res p-values
qvals    = mafdr(pvals, 'BHFDR', true);
sig_mask = qvals < 0.05;

% Convert the full-res significance mask into contiguous intervals (using t_plot_full)
d_sig   = diff([0, sig_mask, 0]);
on_sig  = find(d_sig == 1);
off_sig = find(d_sig == -1) - 1;

% y-position for the bar stays the same; only the x-locations change to full-res time
hold on
signBarHight = 0.15;
for k = 1:numel(on_sig)
    x0 = t_plot_full(on_sig(k));
    x1 = t_plot_full(off_sig(k));
    plot([x0 x1], [signBarHight signBarHight], 'k-', 'LineWidth', 12)
end
legend({'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northwest')

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_scanPathLength_LOWvsHIGH_SPL.png')

%% Scan Path Length time-course (LOW vs HIGH SPL) LMM
close all

% Optional: light trimming of extreme SPL at trial-level (robustness)
keep = isfinite(datTS.SPL) & abs(normalize(datTS.SPL, 'zscore')) <= 4;
datTS = datTS(keep, :);

% Per-timepoint LMMs: SPL ~ Condition + (1|ID)   [trial-level, full-res]
Tf = numel(t_plot_full);
p_lmm  = nan(1, Tf);
t_lmm  = nan(1, Tf);
b_lmm  = nan(1, Tf);  % fixed effect estimate for Condition HIGH vs LOW

for ti = 1:Tf
    idx = datTS.t_index == ti;
    if ~any(idx), continue; end

    tbl = datTS(idx, {'SPL','Condition','ID'});
    % Require at least some data in both levels to fit a contrast
    haveLOW  = any(tbl.Condition == 'LOW');
    haveHIGH = any(tbl.Condition == 'HIGH');
    if ~(haveLOW && haveHIGH), continue; end

    % Mixed model at this timepoint
    lme = fitlme(tbl, 'SPL ~ Condition + (1 + Condition | ID)');
    % estimate the population-level effect of Condition (HIGH vs LOW)
    % and allow each participant to have their own intercept and their own Condition slope

    % Extract the Condition effect (HIGH vs LOW)
    coefTab = lme.Coefficients;
    row = strcmp(coefTab.Name, 'Condition_HIGH'); 
    if any(row)
        b_lmm(ti) = coefTab.Estimate(row);
        t_lmm(ti) = coefTab.tStat(row);
        p_lmm(ti) = coefTab.pValue(row);
    end
end

% Multiple-comparisons control across time (BH-FDR)
q_lmm  = mafdr(p_lmm, 'BHFDR', true);
sig_lmm = q_lmm < 0.05;

% Find contiguous significant segments in full-res time
d_sig   = diff([0, sig_lmm, 0]);
on_sig  = find(d_sig == 1);
off_sig = find(d_sig == -1) - 1;

% Plot grand average of binned data
close all
clc
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

% Plot significance indication
hold on
signBarHight = 0.15;
for k = 1:numel(on_sig)
    x0 = t_plot_full(on_sig(k));
    x1 = t_plot_full(off_sig(k));
    plot([x0 x1], [signBarHight signBarHight], 'k-', 'LineWidth', 12)
end

%% SPL and Alpha Power over Time
alpha_band = [8 14]; % Hz
cfg = [];
cfg.frequency = alpha_band;
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
gatfr_low_alpha = ft_selectdata(cfg, gatfr_low);
gatfr_high_alpha = ft_selectdata(cfg, gatfr_high);
fontSize = 25

close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 2000 1200])

% Top subplot — Scan Path Length
subplot(4,1,1:2)
hold on
shadedErrorBar(t_plot, grand_low, sem_low, ...
    'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
shadedErrorBar(t_plot, grand_high, sem_high, ...
    'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);

% Mark significant intervals
for k = 1:numel(on_sig)
    x0 = t_plot_full(on_sig(k));
    x1 = t_plot_full(off_sig(k));
    %plot([x0 x1], [signBarHight signBarHight], 'k-', 'LineWidth', 12)
end

ylabel('Scan Path Length [px]')
title('Sternberg — Scan Path Length and Alpha Power over time (LOW vs HIGH SPL trials)')
xlim([-.5 2])
set(gca, 'FontSize', fontSize)
box on
legend({'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northwest', 'FontSize', fontSize*0.6)
set(gca, 'XColor', 'none')  % remove x-axis

% Bottom subplot — Alpha Power
subplot(4,1,3:4)
hold on
yyaxis right   % activate right axis
gatfr_diffs = gatfr_high_alpha;
gatfr_diffs.powspctrm = ((gatfr_high_alpha.powspctrm - gatfr_low_alpha.powspctrm) ./ (gatfr_high_alpha.powspctrm + gatfr_low_alpha.powspctrm)) * 100 ;
plot(gatfr_diffs.time, squeeze(gatfr_diffs.powspctrm), ...
    'LineWidth', 3, 'Color', 'k')%, 'LineStyle', '--');
% plot(gatfr_low_alpha.time, squeeze(gatfr_low_alpha.powspctrm), ...
%     'LineWidth', 5, 'Color', colors(1,:), 'LineStyle', '--');
% plot(gatfr_high_alpha.time, squeeze(gatfr_high_alpha.powspctrm), ...
%     'LineWidth', 5, 'Color', colors(3,:), 'LineStyle', '--');

set(gca, 'YColor', [0 0 0])          % colour of right y-axis
% ylabel('Alpha Power [\muV^2/Hz]')
xlabel('Time [s]')
xlim([-.5 2])
ylabel('Relative Alpha Power Difference [%]')
ylim([-5.5 5.5])
yticks([-5 -2.5 0 2.5 5])
yline(0, 'LineWidth', 1, 'LineStyle', '--')
set(gca, 'FontSize', fontSize)
legend({'Alpha Power Difference (HIGH-LOW SPL Trials)'}, 'Location','northeast', 'FontSize', fontSize*0.6)

% hide the left y-axis entirely
yyaxis left
set(gca, 'YColor', 'none', 'YTick', [])

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_SPLOverTime_AlphaPower_Subplots.png')
