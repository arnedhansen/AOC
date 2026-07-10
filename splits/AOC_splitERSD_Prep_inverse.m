%% AOC Split Behavior Inverse Prep — Trial-Level Median Splits (Common Baseline)
% Within-subject median split on gaze deviation or microsaccade rate (conditions
% pooled), then compute occipital alpha ERSD time courses and TFRs per group.
%
% Sternberg outcome summary: ERSD_late [1 2] s
% N-back outcome summary:    ERSD_full [0 2] s
%
% Saves:
%   AOC_splitGazeDev_inverse_trialSplits.mat
%   AOC_splitMS_inverse_trialSplits.mat

%% Setup
startup
[subjects, paths, colors] = setup('AOC');
feat_dir = paths.features;
stats_dir = paths.splits_stats;
fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/splits/SplitERSERD/inverse';
if ~isfolder(stats_dir), mkdir(stats_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end

fig_pos = [0 0 1512 982];
fontSize = 40;
baseline_window = [-1.5 -0.5];
alpha_band = [8 14];
plot_latency = [-0.5 2];

tasks(1).tag = 'sternberg';
tasks(1).data_fname = 'dataEEG_TFR_sternberg.mat';
tasks(1).data_var = 'dataTFR';
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_out = [2 4 6];
tasks(1).toi = -1.5:0.05:3;
tasks(1).ersd_latency = [1 2];
tasks(1).ersd_var = 'ERSD_late';
tasks(1).gaze_trials_file = 'AOC_gaze_matrix_sternberg_trials.mat';
tasks(1).gaze_trials_var = 'gaze_data_sternberg_trials';

tasks(2).tag = 'nback';
tasks(2).data_fname = 'dataEEG_TFR_nback.mat';
tasks(2).data_var = 'dataTFR';
tasks(2).cond_codes = [21 22 23];
tasks(2).cond_out = [1 2 3];
tasks(2).toi = -1.5:0.05:2.25;
tasks(2).ersd_latency = [0 2];
tasks(2).ersd_var = 'ERSD_full';
tasks(2).gaze_trials_file = 'AOC_gaze_matrix_nback_trials.mat';
tasks(2).gaze_trials_var = 'gaze_data_nback_trials';

split_defs(1).key = 'gaze';
split_defs(1).split_var = 'GazeDeviationFullBL';
split_defs(1).group_lbl_low = 'Low Gaze Dev';
split_defs(1).group_lbl_high = 'High Gaze Dev';
split_defs(1).out_struct = 'split_gaze_dev_inverse';
split_defs(1).out_file = 'AOC_splitGazeDev_inverse_trialSplits.mat';

split_defs(2).key = 'ms';
split_defs(2).split_var = 'MSRateFullBL';
split_defs(2).group_lbl_low = 'Low MS';
split_defs(2).group_lbl_high = 'High MS';
split_defs(2).out_struct = 'split_ms_inverse';
split_defs(2).out_file = 'AOC_splitMS_inverse_trialSplits.mat';

fprintf('\n=== AOC Split Behavior Inverse Prep (trial-level) ===\n');

for di = 1:numel(split_defs)
    sd = split_defs(di);
    split_out = struct();
    split_out.meta.created = datestr(now, 31);
    split_out.meta.script = 'AOC_splitBehavior_inverse_Prep.m';
    split_out.meta.split_var = sd.split_var;
    split_out.meta.baseline_window = baseline_window;
    split_out.meta.baseline_type = 'subject_pooled_db';
    split_out.meta.alpha_band = alpha_band;
    split_out.meta.roi = 'occipital O/I channels';
    split_out.meta.split_rule = 'within-subject median across all trials (conditions pooled)';

    for ti = 1:numel(tasks)
        tk = tasks(ti);
        task_tag = tk.tag;
        fprintf('\n\n========== %s | %s ==========\n', upper(task_tag), sd.split_var);

        gaze_mat_file = fullfile(feat_dir, tk.gaze_trials_file);
        if ~isfile(gaze_mat_file)
            error('Missing gaze trial matrix: %s', gaze_mat_file);
        end
        Gmat = load(gaze_mat_file, tk.gaze_trials_var);
        Tg_all = struct2table(Gmat.(tk.gaze_trials_var));

        task_out = struct();
        task_out.tag = task_tag;
        task_out.split_label = 'splitMedian_trial';
        task_out.split_var = sd.split_var;
        task_out.ersd_var = tk.ersd_var;
        task_out.ersd_latency = tk.ersd_latency;
        task_out.group_lbl_low = sd.group_lbl_low;
        task_out.group_lbl_high = sd.group_lbl_high;
        task_out.subjects = struct( ...
            'ID', {}, 'Trial', {}, 'Condition', {}, 'SplitMetric', {}, 'ERSD', {}, ...
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
            fprintf('[PREP inverse %s - %s] Subject %d / %d (%s)\n', sd.key, upper(task_tag), s, numel(subjects), sid_str);

            eeg_dir = fullfile(feat_dir, sid_str, 'eeg');
            data_file = fullfile(eeg_dir, tk.data_fname);
            if ~isfile(data_file)
                fprintf('  Missing %s\n', data_file);
                n_fail = n_fail + 1;
                continue
            end

            try
                rows = Tg_all(Tg_all.ID == sid, :);
                if isempty(rows) || ~ismember(sd.split_var, rows.Properties.VariableNames)
                    error('No gaze rows or missing %s for subject %s.', sd.split_var, sid_str);
                end

                S = load(data_file, tk.data_var);
                dataTFR = S.(tk.data_var);
                condCodes = dataTFR.trialinfo(:, 1);
                trialIDs = dataTFR.trialinfo(:, 2);
                nTrials = numel(trialIDs);

                split_metric = nan(nTrials, 1);
                for tr = 1:nTrials
                    ix = find(rows.Trial == trialIDs(tr), 1);
                    if ~isempty(ix)
                        split_metric(tr) = rows.(sd.split_var)(ix);
                    end
                end

                valid = isfinite(split_metric);
                if sum(valid) < 4
                    error('Too few finite trials for median split (%d).', sum(valid));
                end
                thr = median(split_metric(valid), 'omitnan');
                idx_low = valid & (split_metric < thr);
                idx_high = valid & (split_metric >= thr);
                if ~any(idx_low) || ~any(idx_high)
                    error('Median split produced empty group (n_low=%d, n_high=%d).', sum(idx_low), sum(idx_high));
                end

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

                [tf_bl, ~] = apply_subject_pooled_db_baseline(tf, baseline_window);
                chUse = occ_channels_from_labels(tf_bl.label);
                chIdx = find(ismember(tf_bl.label, chUse));
                if isempty(chIdx), chIdx = 1:numel(tf_bl.label); end
                fMask = tf_bl.freq >= alpha_band(1) & tf_bl.freq <= alpha_band(2);
                tMaskErsd = tf_bl.time >= tk.ersd_latency(1) & tf_bl.time <= tk.ersd_latency(2);
                tMaskPlot = tf_bl.time >= plot_latency(1) & tf_bl.time <= plot_latency(2);
                time_plot = tf_bl.time(tMaskPlot);

                ersd = nan(nTrials, 1);
                for tr = 1:nTrials
                    x = tf_bl.powspctrm(tr, chIdx, fMask, tMaskErsd);
                    ersd(tr) = mean(x(:), 'omitnan');
                end

                tc_all = nan(nTrials, numel(time_plot));
                for tr = 1:nTrials
                    x = squeeze(mean(mean(tf_bl.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                    tc_all(tr, :) = x(:)';
                end
                tc_low = mean(tc_all(idx_low, :), 1, 'omitnan');
                tc_high = mean(tc_all(idx_high, :), 1, 'omitnan');

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
                entry.SplitMetric = split_metric(:);
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
                if numel(time_plot) == numel(ga_time)
                    ga_low(s, :) = tc_low;
                    ga_high(s, :) = tc_high;
                else
                    ga_low(s, :) = interp1(time_plot(:), tc_low(:), ga_time, 'linear', NaN);
                    ga_high(s, :) = interp1(time_plot(:), tc_high(:), ga_time, 'linear', NaN);
                end
                n_ok = n_ok + 1;

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
        split_out.(task_tag) = task_out;

        if n_ok >= 2 && ~isempty(ga_time)
            mL = mean(ga_low, 1, 'omitnan');
            mH = mean(ga_high, 1, 'omitnan');
            nL = sum(isfinite(ga_low), 1);
            nH = sum(isfinite(ga_high), 1);
            sL = std(ga_low, 0, 1, 'omitnan') ./ max(sqrt(nL), 1);
            sH = std(ga_high, 0, 1, 'omitnan') ./ max(sqrt(nH), 1);
            plot_ersd_group_timecourse(ga_time, mL, sL, mH, sH, ...
                colors, sd.group_lbl_low, sd.group_lbl_high, fontSize, fig_pos, ...
                sprintf('%s %s GA ERSD (n=%d)', task_tag, sd.key, n_ok), ...
                fullfile(fig_dir, sprintf('AOC_split_inverse_prep_%s_%s_GA_ersd.png', task_tag, sd.key)));
        end
        fprintf('Task %s (%s) done: ok=%d fail=%d\n', task_tag, sd.key, n_ok, n_fail);
    end

    out_path = fullfile(stats_dir, sd.out_file);
    eval(sprintf('%s = split_out;', sd.out_struct)); %#ok<EVLC>
    save(out_path, sd.out_struct, '-v7.3');
    fprintf('Saved: %s\n', out_path);
end

fprintf('\nDone.\n');

%% ========================= Local Functions =========================
function [tf_bl, bl_ref] = apply_subject_pooled_db_baseline(tf, bl_win)
P = tf.powspctrm;
if ndims(P) ~= 4
    error('Expected keeptrials powspctrm (rpt x chan x freq x time), got ndims=%d.', ndims(P));
end
tMask = tf.time >= bl_win(1) & tf.time <= bl_win(2);
if ~any(tMask)
    error('No time samples in baseline window [%.2f %.2f].', bl_win(1), bl_win(2));
end
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
tfr_avg.powspctrm = squeeze(mean(P, 1, 'omitnan'));
if isfield(tfr_avg, 'trialinfo')
    tfr_avg.trialinfo = tf_bl.trialinfo(tr_mask, :);
end
tfr_avg.dimord = 'chan_freq_time';
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
if isempty(ch), ch = labels; end
end

function plot_ersd_group_timecourse(t, yLow, eLow, yHigh, eHigh, colors, lblLow, lblHigh, fsz, fig_pos, ttl, out_file)
figure('Position', fig_pos, 'Color', 'w');
hold on
legendFontSize = fsz * 0.666;
if isempty(eLow), eLow = zeros(size(yLow)); end
if isempty(eHigh), eHigh = zeros(size(yHigh)); end
ebL = shadedErrorBar(t, yLow, eLow, 'lineProps', {'-', 'Color', colors(1, :)});
ebH = shadedErrorBar(t, yHigh, eHigh, 'lineProps', {'-', 'Color', colors(3, :)});
set(ebL.mainLine, 'LineWidth', 2, 'Color', colors(1, :));
set(ebH.mainLine, 'LineWidth', 2, 'Color', colors(3, :));
set(ebL.patch, 'FaceColor', colors(1, :), 'FaceAlpha', 0.20);
set(ebH.patch, 'FaceColor', colors(3, :), 'FaceAlpha', 0.20);
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
