%% AOC Split Alpha Amp/Red Prep — Trial-Level Median Splits (Common Baseline)
% Recompute keeptrials TFR from dataEEG_TFR_*, apply one subject-pooled
% baseline in [-1.5 -0.5] s (dB), reduce occipital ROI x IAF band (IAF-4, IAF+2)
% x split latency, then median-split ALL trials within each subject (conditions pooled).
% Missing IAF for a condition falls back to fixed [8 14] Hz.
%
% Sternberg split metric: ERSD late window [1 2] s -> Low Alpha vs High Alpha
% N-back split metric:    ERSD full window [0 2] s -> Low Alpha vs High Alpha
%
% Saves trial indices + scalars for MS / GazeDev split scripts, plus sanity
% ERSD time-course figures (toggle single-subject figures below).
%
% Note: 2_feature_extraction/eeg/trials/AOC_eeg_fex_*_trials.m uses
% ft_freqbaseline with keeptrials (single-trial baseline) and is NOT used here.

%% Setup
startup
[subjects, paths, colors] = setup('AOC');
feat_dir = paths.features;
stats_dir = paths.splits_stats;
fig_dir = fullfile(paths.figures, 'splits', 'SplitERSERD', 'Prep');
if ~isfolder(stats_dir), mkdir(stats_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end

SAVE_SINGLE_SUBJECT_FIGS = false;  % true: save per-subject ERSD group time courses

fig_pos = [0 0 1512 982];
fontSize = 40;
baseline_window = [-1.5 -0.5];
alphaRange = [8 14];  % peak search window for IAF; ERSD uses (IAF-4, IAF+2)
plot_latency = [-0.5 2];

tasks(1).tag = 'sternberg';
tasks(1).data_fname = 'dataEEG_TFR_sternberg.mat';
tasks(1).data_var = 'dataTFR';
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_out = [2 4 6];
tasks(1).toi = -1.5:0.05:3;
tasks(1).winIAF = [1 2];
tasks(1).split_latency = [1 2];
tasks(1).ersd_var = 'ERSD_late';
tasks(1).group_lbl_low = 'Low Alpha';
tasks(1).group_lbl_high = 'High Alpha';

tasks(2).tag = 'nback';
tasks(2).data_fname = 'dataEEG_TFR_nback.mat';
tasks(2).data_var = 'dataTFR';
tasks(2).cond_codes = [21 22 23];
tasks(2).cond_out = [1 2 3];
tasks(2).toi = -1.5:0.05:2.25;
tasks(2).winIAF = [0 2];
tasks(2).split_latency = [0 2];
tasks(2).ersd_var = 'ERSD_full';
tasks(2).group_lbl_low = 'Low Alpha';
tasks(2).group_lbl_high = 'High Alpha';

fprintf('\n=== AOC Split AlphaAmpRed Prep (trial-level, common baseline) ===\n');
fprintf('SAVE_SINGLE_SUBJECT_FIGS = %d\n', SAVE_SINGLE_SUBJECT_FIGS);
fprintf('Output: %s\n', stats_dir);
fprintf('Figures: %s\n', fig_dir);

split_alpha_amp_red = struct();
split_alpha_amp_red.meta.created = datestr(now, 31);
split_alpha_amp_red.meta.script = 'AOC_splits_AlphaAmpRed_Prep.m';
split_alpha_amp_red.meta.baseline_window = baseline_window;
split_alpha_amp_red.meta.baseline_type = 'subject_pooled_db';
split_alpha_amp_red.meta.baseline_note = ...
    'Mean raw power over all trials and time in baseline_window (chan x freq), then 10*log10(pow/baseline) per trial';
split_alpha_amp_red.meta.alpha_band = 'IAF-4 to IAF+2 (condition-wise; fallback [8 14] Hz)';
split_alpha_amp_red.meta.alphaRange_peak = alphaRange;
split_alpha_amp_red.meta.roi = 'occipital O/I channels';
split_alpha_amp_red.meta.split_rule = 'within-subject median across all trials (conditions pooled)';

for ti = 1:numel(tasks)
    tk = tasks(ti);
    task_tag = tk.tag;
    fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));

    fig_dir_subj = fullfile(fig_dir, 'single_subjects', task_tag);
    if SAVE_SINGLE_SUBJECT_FIGS && ~isfolder(fig_dir_subj)
        mkdir(fig_dir_subj);
    end

    task_out = struct();
    task_out.tag = task_tag;
    task_out.split_label = 'splitMedian_trial';
    task_out.ersd_var = tk.ersd_var;
    task_out.split_latency = tk.split_latency;
    task_out.group_lbl_low = tk.group_lbl_low;
    task_out.group_lbl_high = tk.group_lbl_high;
    task_out.subjects = struct( ...
        'ID', {}, 'Trial', {}, 'Condition', {}, 'ERSD', {}, ...
        'median_thr', {}, 'idx_low', {}, 'idx_high', {}, ...
        'trial_ids_low', {}, 'trial_ids_high', {}, ...
        'ersd_tc_low', {}, 'ersd_tc_high', {}, 'time', {}, ...
        'tfr_low', {}, 'tfr_high', {}, 'n_low', {}, 'n_high', {});

    ga_low = [];
    ga_high = [];
    ga_time = [];
    n_ok = 0;
    n_fail = 0;

    for s = 1:numel(subjects)
        sid_str = subjects{s};
        sid = str2double(sid_str);
        clc
        fprintf('[PREP %s] Subject %d / %d (%s)\n', upper(task_tag), s, numel(subjects), sid_str);

        eeg_dir = fullfile(feat_dir, sid_str, 'eeg');
        data_file = fullfile(eeg_dir, tk.data_fname);
        if ~isfile(data_file)
            fprintf('  Missing %s\n', data_file);
            n_fail = n_fail + 1;
            continue
        end

        try
            S = load(data_file, tk.data_var);
            dataTFR = S.(tk.data_var);

            condCodes = dataTFR.trialinfo(:, 1);
            trialIDs = dataTFR.trialinfo(:, 2);
            condSet = unique(condCodes(:))';
            if ~all(ismember(condSet, tk.cond_codes))
                error('Unexpected condition codes: %s', mat2str(condSet));
            end

            % Keeptrials TFR (raw power)
            cfg = [];
            cfg.method = 'mtmconvol';
            cfg.output = 'pow';
            cfg.taper = 'hanning';
            cfg.foi = 2:2:40;
            cfg.t_ftimwin = ones(size(cfg.foi)) * 0.5;
            cfg.toi = tk.toi;
            cfg.pad = 'nextpow2';
            cfg.keeptrials = 'yes';
            tf = ft_freqanalysis(cfg, dataTFR);

            % Subject-pooled common baseline (NOT single-trial)
            [tf_bl, ~] = apply_subject_pooled_db_baseline(tf, baseline_window);

            chUse = occ_channels_from_labels(tf_bl.label);
            chIdx = find(ismember(tf_bl.label, chUse));
            if isempty(chIdx)
                chIdx = 1:numel(tf_bl.label);
            end
            tMaskSplit = tf_bl.time >= tk.split_latency(1) & tf_bl.time <= tk.split_latency(2);
            tMaskPlot = tf_bl.time >= plot_latency(1) & tf_bl.time <= plot_latency(2);
            time_plot = tf_bl.time(tMaskPlot);

            % Condition-wise IAF -> (IAF-4, IAF+2) band for each trial
            iaf_by_cond = nan(numel(tk.cond_codes), 1);
            for c = 1:numel(tk.cond_codes)
                trlIdx = find(condCodes == tk.cond_codes(c));
                [iaf_by_cond(c), ~] = iaf_from_retention_mtmfft(dataTFR, trlIdx, tk.winIAF, chUse, alphaRange);
            end

            nTrials = size(tf_bl.powspctrm, 1);
            ersd = nan(nTrials, 1);
            tc_all = nan(nTrials, numel(time_plot));
            for tr = 1:nTrials
                c = find(tk.cond_codes == condCodes(tr), 1);
                if isempty(c)
                    continue
                end
                band = ersd_alpha_band(iaf_by_cond(c), alphaRange);
                fMask = tf_bl.freq >= band(1) & tf_bl.freq <= band(2);
                if ~any(fMask)
                    continue
                end
                x = tf_bl.powspctrm(tr, chIdx, fMask, tMaskSplit);
                ersd(tr) = mean(x(:), 'omitnan');
                x = squeeze(mean(mean(tf_bl.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                tc_all(tr, :) = x(:)';
            end

            valid = isfinite(ersd);
            if sum(valid) < 4
                error('Too few finite trials for median split (%d).', sum(valid));
            end
            thr = median(ersd(valid), 'omitnan');
            idx_low = valid & (ersd < thr);
            idx_high = valid & (ersd >= thr);
            % If ties pile on the threshold, keep balanced by assigning equals
            % already to high; ensure both groups nonempty.
            if ~any(idx_low) || ~any(idx_high)
                error('Median split produced empty group (n_low=%d, n_high=%d).', sum(idx_low), sum(idx_high));
            end

            tc_low = mean(tc_all(idx_low, :), 1, 'omitnan');
            tc_high = mean(tc_all(idx_high, :), 1, 'omitnan');

            % Within-group trial-averaged TFR (for GazeDev EEG panels)
            tfr_low = average_trials_freq(tf_bl, idx_low);
            tfr_high = average_trials_freq(tf_bl, idx_high);

            condOut = nan(nTrials, 1);
            for c = 1:numel(tk.cond_codes)
                condOut(condCodes == tk.cond_codes(c)) = tk.cond_out(c);
            end

            entry = struct();
            entry.ID = sid;
            entry.Trial = trialIDs(:);
            entry.Condition = condOut(:);
            entry.ERSD = ersd(:);
            entry.median_thr = thr;
            entry.idx_low = idx_low(:);
            entry.idx_high = idx_high(:);
            entry.trial_ids_low = trialIDs(idx_low);
            entry.trial_ids_high = trialIDs(idx_high);
            entry.ersd_tc_low = tc_low(:)';
            entry.ersd_tc_high = tc_high(:)';
            entry.time = time_plot(:)';
            entry.tfr_low = tfr_low;
            entry.tfr_high = tfr_high;
            entry.n_low = sum(idx_low);
            entry.n_high = sum(idx_high);
            task_out.subjects(end + 1) = entry; %#ok<AGROW>

            if isempty(ga_time)
                ga_time = time_plot(:)';
                ga_low = nan(numel(subjects), numel(ga_time));
                ga_high = nan(numel(subjects), numel(ga_time));
            end
            % Interpolate if time vector length mismatches (should not)
            if numel(time_plot) == numel(ga_time) && all(abs(time_plot(:) - ga_time) < 1e-9)
                ga_low(s, :) = tc_low;
                ga_high(s, :) = tc_high;
            else
                ga_low(s, :) = interp1(time_plot(:), tc_low(:), ga_time, 'linear', NaN);
                ga_high(s, :) = interp1(time_plot(:), tc_high(:), ga_time, 'linear', NaN);
            end
            n_ok = n_ok + 1;

            if SAVE_SINGLE_SUBJECT_FIGS
                plot_ersd_group_timecourse(time_plot, tc_low, [], tc_high, [], ...
                    colors, tk.group_lbl_low, tk.group_lbl_high, fontSize, fig_pos, ...
                    sprintf('%s | subj %s (n_low=%d, n_high=%d)', task_tag, sid_str, sum(idx_low), sum(idx_high)), ...
                    fullfile(fig_dir_subj, sprintf('AOC_splitERSERD_prep_%s_subj%s_timecourse.png', task_tag, sid_str)));
            end

            % Per-subject lightweight save (indices + scalars; no large TFR)
            subj_split = rmfield(entry, {'tfr_low', 'tfr_high'});
            save(fullfile(eeg_dir, sprintf('alphaAmpRed_trialSplit_%s.mat', task_tag)), 'subj_split');

        catch ME
            fprintf('  FAILED %s: %s\n', sid_str, ME.message);
            n_fail = n_fail + 1;
        end
    end

    task_out.n_subjects_ok = n_ok;
    task_out.n_subjects_fail = n_fail;
    task_out.ga_time = ga_time;
    task_out.ga_low = ga_low;
    task_out.ga_high = ga_high;
    split_alpha_amp_red.(task_tag) = task_out;

    % Grand-average sanity figure (mean +/- SEM across subjects)
    if n_ok >= 2 && ~isempty(ga_time)
        mL = mean(ga_low, 1, 'omitnan');
        mH = mean(ga_high, 1, 'omitnan');
        nL = sum(isfinite(ga_low), 1);
        nH = sum(isfinite(ga_high), 1);
        sL = std(ga_low, 0, 1, 'omitnan') ./ max(sqrt(nL), 1);
        sH = std(ga_high, 0, 1, 'omitnan') ./ max(sqrt(nH), 1);
        plot_ersd_group_timecourse(ga_time, mL, sL, mH, sH, ...
            colors, tk.group_lbl_low, tk.group_lbl_high, fontSize, fig_pos, ...
            sprintf('%s grand average (n=%d subjects)', task_tag, n_ok), ...
            fullfile(fig_dir, sprintf('AOC_splitERSERD_prep_%s_GA_timecourse.png', task_tag)));
    end

    fprintf('Task %s done: ok=%d fail=%d\n', task_tag, n_ok, n_fail);
end

out_file = fullfile(stats_dir, 'AOC_splitAlphaAmpRed_trialSplits.mat');
save(out_file, 'split_alpha_amp_red', '-v7.3');
fprintf('\nSaved trial splits to: %s\n', out_file);
fprintf('Done.\n');

%% ========================= Local Functions =========================
function [tf_bl, bl_ref] = apply_subject_pooled_db_baseline(tf, bl_win)
% APPLY_SUBJECT_POOLED_DB_BASELINE
% Mean raw power over ALL trials and baseline time (chan x freq), then
% convert every trial to dB: 10*log10(pow / baseline).
P = tf.powspctrm; % rpt x chan x freq x time
if ndims(P) ~= 4
    error('Expected keeptrials powspctrm (rpt x chan x freq x time), got ndims=%d.', ndims(P));
end
tMask = tf.time >= bl_win(1) & tf.time <= bl_win(2);
if ~any(tMask)
    error('No time samples in baseline window [%.2f %.2f].', bl_win(1), bl_win(2));
end
% Average over time then trials -> chan x freq
bl_ref = squeeze(mean(mean(P(:, :, :, tMask), 4, 'omitnan'), 1, 'omitnan'));
if isvector(bl_ref)
    bl_ref = reshape(bl_ref, [numel(tf.label), numel(tf.freq)]);
end
bl_ref(~isfinite(bl_ref) | bl_ref <= 0) = NaN;

bl4 = reshape(bl_ref, [1 size(bl_ref, 1) size(bl_ref, 2) 1]);
P_db = 10 * log10(bsxfun(@rdivide, P, bl4));

tf_bl = tf;
tf_bl.powspctrm = P_db;
end

function tfr_avg = average_trials_freq(tf_bl, tr_mask)
tfr_avg = tf_bl;
P = tf_bl.powspctrm(tr_mask, :, :, :);
tfr_avg.powspctrm = squeeze(mean(P, 1, 'omitnan')); % chan x freq x time
if isfield(tfr_avg, 'trialinfo')
    tfr_avg.trialinfo = tf_bl.trialinfo(tr_mask, :);
end
tfr_avg.dimord = 'chan_freq_time';
% Drop rpt-related cumtapcnt if present
if isfield(tfr_avg, 'cumtapcnt')
    tfr_avg = rmfield(tfr_avg, 'cumtapcnt');
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end + 1} = lab; %#ok<AGROW>
    end
end
if isempty(ch)
    ch = labels;
end
end

function [IAF, powerIAF] = iaf_from_retention_mtmfft(dataTFR, trialinds, winSec, chLabs, alphaRange)
IAF = NaN;
powerIAF = NaN;
if isempty(trialinds) || isempty(chLabs)
    return
end
try
    cfg_sel = [];
    cfg_sel.latency = winSec;
    cfg_sel.trials = trialinds(:);
    cfg_sel.channel = chLabs(:);
    dw = ft_selectdata(cfg_sel, dataTFR);
    if isempty(dw.trial)
        return
    end
    cfgf = [];
    cfgf.method = 'mtmfft';
    cfgf.output = 'pow';
    cfgf.taper = 'dpss';
    cfgf.tapsmofrq = 2;
    cfgf.foilim = [6 18];
    cfgf.pad = 'nextpow2';
    cfgf.keeptrials = 'no';
    fr = ft_freqanalysis(cfgf, dw);
    ps = mean(fr.powspctrm, 1);
    [IAF, powerIAF] = iaf_peak_rules(fr.freq(:), ps(:), alphaRange);
catch
    IAF = NaN;
    powerIAF = NaN;
end
end

function [IAF, powerIAF] = iaf_peak_rules(freq, spec, alphaRange)
freq = freq(:);
spec = spec(:);
IAF = NaN;
powerIAF = NaN;
amask = freq >= alphaRange(1) & freq <= alphaRange(2);
alphaFreqs = freq(amask);
alphaSpec = spec(amask);
if numel(alphaSpec) < 3
    return
end
[pks, locs] = findpeaks(alphaSpec);
if isempty(pks)
    return
end
[~, ind] = max(pks);
IAF = alphaFreqs(locs(ind));
bandIdx = freq > (IAF - 4) & freq < (IAF + 2);
if any(bandIdx)
    powerIAF = mean(spec(bandIdx));
end
if locs(ind) == 1 || locs(ind) == numel(alphaFreqs)
    powerIAF = NaN;
end
end

function band = ersd_alpha_band(IAF, alphaRange)
% ERSD band: (IAF-4, IAF+2) when IAF is valid; otherwise fixed alphaRange ([8 14]).
if ~isfinite(IAF)
    band = alphaRange;
else
    band = [IAF - 4, IAF + 2];
end
if any(~isfinite(band)) || band(1) >= band(2)
    band = alphaRange;
end
end

function plot_ersd_group_timecourse(t, yLow, eLow, yHigh, eHigh, colors, lblLow, lblHigh, fsz, fig_pos, ttl, out_file)
figure('Position', fig_pos, 'Color', 'w');
hold on
legendFontSize = fsz * 0.666;
if isempty(eLow), eLow = zeros(size(yLow)); end
if isempty(eHigh), eHigh = zeros(size(yHigh)); end
ebL = shadedErrorBar(t, yLow, eLow, 'lineProps', {'-', 'Color', colors(1, :)});
ebH = shadedErrorBar(t, yHigh, eHigh, 'lineProps', {'-', 'Color', colors(3, :)});
set(ebL.mainLine, 'LineWidth', 5, 'Color', colors(1, :));
set(ebH.mainLine, 'LineWidth', 5, 'Color', colors(3, :));
set(ebL.patch, 'FaceColor', colors(1, :), 'FaceAlpha', 0.25);
set(ebH.patch, 'FaceColor', colors(3, :), 'FaceAlpha', 0.25);
set(ebL.edge(1), 'Color', 'none'); set(ebL.edge(2), 'Color', 'none');
set(ebH.edge(1), 'Color', 'none'); set(ebH.edge(2), 'Color', 'none');
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
yline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
set(gca, 'FontSize', fsz);
xlabel('Time [s]');
ylabel('Power Change [dB]');
xlim([-0.5 2]);
title(ttl, 'Interpreter', 'none', 'FontSize', fsz * 0.7);
leg1 = patch(NaN, NaN, colors(1, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(1, :), 'LineWidth', 1.5);
leg2 = patch(NaN, NaN, colors(3, :), 'FaceAlpha', 0.25, 'EdgeColor', colors(3, :), 'LineWidth', 1.5);
legend([leg1, leg2], {[' ' lblLow], [' ' lblHigh]}, ...
    'Location', 'southeast', 'FontSize', legendFontSize, 'Box', 'off');
box off
drawnow; pause(0.05);
saveas(gcf, out_file);
close(gcf);
end
