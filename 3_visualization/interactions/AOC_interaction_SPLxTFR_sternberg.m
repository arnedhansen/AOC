%% AOC Interaction — Sternberg: SPL × TFR
% Median-splits trials by within-subject SPL, compares TFR and SPL time-course between low and high SPL. Loads merged_data_sternberg_trials. Saves figures.
%
% Key outputs:
%   TFR and SPL figures (low vs high SPL)

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
%subjects = subjects(1:10)

% Figure config
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
    % figure;
    % set(gcf, 'Color', 'w', 'Position', [0 0 1512 982]);
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
    % title(sprintf('Subject %s — within-subject median split (LOW=%d, HIGH=%d)', subjects{subj}, nHalf, numel(spl_sorted)-nHalf))
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
%load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_interaction_SPLxTFR_sternberg');

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
set(gcf, 'Color', 'w', 'Position', [0 0 1512 982]);
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
set(gcf, 'Color', 'w', 'Position', [0 0 1512 982]);
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
set(gcf, 'Color', 'w', 'Position', [0 0 1512 982]);
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

%% SPL and Alpha Power over Time
alpha_band = [8 14];
cfg = [];
cfg.frequency   = alpha_band;
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
fontSize = 25;

% Subject alpha time-courses
alpha_low_mat  = [];
alpha_high_mat = [];
t_alpha = [];
for s = 1:numel(low_tfr_subs)
    if isempty(low_tfr_subs{s}) || isempty(high_tfr_subs{s})
        continue
    end
    la = ft_selectdata(cfg, low_tfr_subs{s});   % 1 x time
    ha = ft_selectdata(cfg, high_tfr_subs{s});  % 1 x time

    if isempty(t_alpha)
        t_alpha = la.time;
    end

    alpha_low_mat(end+1, :)  = squeeze(la.powspctrm);
    alpha_high_mat(end+1, :) = squeeze(ha.powspctrm);
end

% Paired t-tests across subjects at each timepoint (HIGH vs LOW)
[~, p_alpha, ~, stats_alpha] = ttest(alpha_high_mat, alpha_low_mat);
tvals = stats_alpha.tstat;

% Cohen's dz for paired designs
diff_alpha = alpha_high_mat - alpha_low_mat;
dvals = nanmean(diff_alpha, 1) ./ nanstd(diff_alpha, 0, 1);

% FDR across time
toi = t_alpha >= 0 & t_alpha <= 2;
p_sub = p_alpha(toi & isfinite(p_alpha));
q_sub = mafdr(p_sub, 'BHFDR', true);
sig_alpha = false(size(p_alpha));
sig_alpha(toi) = q_sub < 0.05;

% Prepare SPL traces (LOW vs HIGH, binned grid)
grand_low  = nanmean(scan_low,  1);
grand_high = nanmean(scan_high, 1);
sem_low    = nanstd(scan_low,  [], 1) ./ sqrt(sum(isfinite(scan_low),  1));
sem_high   = nanstd(scan_high, [], 1) ./ sqrt(sum(isfinite(scan_high), 1));
t_plot     = time_series(2:end);

% Plot
close all
figure
set(gcf, 'Color', 'w', 'Position', [0 0 1512 900])

% Top subplot: Scan Path Length
subplot(4,1,1:3)
hold on
sebLow = shadedErrorBar(t_plot, grand_low,  sem_low,  'lineProps', {'-','Color',colors(1,:),'LineWidth',2.5}, 'transparent', true);
sebHigh = shadedErrorBar(t_plot, grand_high, sem_high, 'lineProps', {'-','Color',colors(3,:),'LineWidth',2.5}, 'transparent', true);
ylabel('Scan Path Length [px]')
title('Sternberg SPLxAlpha (LOW vs HIGH SPL trials)')
xlim([-.5 2])
set(gca, 'FontSize', fontSize)
box on
yline(0, '--')
legend([sebLow.mainLine, sebHigh.mainLine], {'LOW SPL ± SEM','HIGH SPL ± SEM'}, 'Location','northeast', 'FontSize', fontSize*0.6)
set(gca, 'XColor', 'none')

% Bottom subplot: Alpha
subplot(4,1,4)
hold on
y = dvals;
ylab = 'Cohen''s d_z';
h = plot(t_alpha, y, 'k-', 'LineWidth', 3);
yline(0, '--', 'LineWidth', 1)
xlim([-.5 2])
set(gca, 'FontSize', fontSize)
xlabel('Time [s]')
ylabel(ylab)

% Shade FDR-significant intervals
ylim = [-max(abs(dvals))*1.15 max(abs(dvals))*1.15];
yticks([-.5 -.25 0 .25 .5]);
d_sig   = diff([0, sig_alpha, 0]);
on_sig  = find(d_sig == 1);
off_sig = find(d_sig == -1) - 1;

for k = 1:numel(on_sig)
    x0 = t_alpha(on_sig(k));
    x1 = t_alpha(off_sig(k));
    patch([x0 x1 x1 x0], [ylim(1) ylim(1) ylim(2) ylim(2)], ...
        [0 0 0], 'FaceAlpha', 0.12, 'EdgeColor', 'none')
end
uistack(h, 'top')

saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/AOC_sternberg_SPLOxTFR_Subplots.png')