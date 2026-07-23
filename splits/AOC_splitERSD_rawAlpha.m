%% AOC Split by Raw Alpha Power — Full Pipeline (Prep + GazeDev + MS)
% Within-subject median split on RAW occipital IAF-band alpha power
% (no subject-pooled dB baseline), then GazeDev and MS analyses with CBPT.
%
% Sternberg split window: [1 2] s | N-back split window: [0 2] s
% Alpha band: IAF-4 to IAF+2 (fallback [8 14] Hz)
%
% Figures:  .../figures/splits/SplitERSERD/rawAlpha/
% Stats:    paths.splits_stats (distinct *rawAlpha* filenames)

%% Setup
startup
[subjects, paths, colors, headmodel] = setup('AOC');
addpath(paths.seb_path);

feat_dir = paths.features;
stats_dir = paths.splits_stats;
fig_dir_root = fullfile(paths.figures, 'splits', 'SplitERSERD', 'rawAlpha');
fig_dir_prep = fullfile(fig_dir_root, 'Prep');
if ~isfolder(stats_dir), mkdir(stats_dir); end
if ~isfolder(fig_dir_root), mkdir(fig_dir_root); end
if ~isfolder(fig_dir_prep), mkdir(fig_dir_prep); end

SAVE_SINGLE_SUBJECT_FIGS = false;
fig_pos = [0 0 1512 982];
fontSize = 40;
alphaRange = [8 14];
plot_latency = [-0.5 2];
tfr_winsor_cfg = struct('enable', true, 'prctile', [2 98]);
split_label = 'splitMedian_trial_rawAlpha';

tasks(1).tag = 'sternberg';
tasks(1).data_fname = 'dataEEG_TFR_sternberg.mat';
tasks(1).data_var = 'dataTFR';
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_out = [2 4 6];
tasks(1).toi = -1.5:0.05:3;
tasks(1).winIAF = [1 2];
tasks(1).split_latency = [1 2];
tasks(1).topo_latency = [1 2];
tasks(1).split_var = 'AlphaRaw_late';
tasks(1).et_fname = 'dataET_sternberg.mat';
tasks(1).et_fname_ms = 'dataET_sternberg';
tasks(1).gaze_trials_file = 'AOC_gaze_matrix_sternberg_trials.mat';
tasks(1).gaze_trials_var = 'gaze_data_sternberg_trials';
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
tasks(2).topo_latency = [0 2];
tasks(2).split_var = 'AlphaRaw_full';
tasks(2).et_fname = 'dataET_nback.mat';
tasks(2).et_fname_ms = 'dataET_nback';
tasks(2).gaze_trials_file = 'AOC_gaze_matrix_nback_trials.mat';
tasks(2).gaze_trials_var = 'gaze_data_nback_trials';
tasks(2).group_lbl_low = 'Low Alpha';
tasks(2).group_lbl_high = 'High Alpha';

fprintf('\n=== AOC Split RAW Alpha — Full Pipeline (Prep + GazeDev + MS) ===\n');
fprintf('Figures: %s\n', fig_dir_root);

%% ========================= PART 1: PREP (raw alpha splits) =========================
split_alpha_amp_red_raw = struct();
split_alpha_amp_red_raw.meta.created = datestr(now, 31);
split_alpha_amp_red_raw.meta.script = 'AOC_splitERSD_Prep_rawAlpha.m';
split_alpha_amp_red_raw.meta.baseline_type = 'none_raw_power';
split_alpha_amp_red_raw.meta.baseline_note = ...
    'No baseline: split uses raw TFR power (occipital ROI x IAF band)';
split_alpha_amp_red_raw.meta.alpha_band = 'IAF-4 to IAF+2 (condition-wise; fallback [8 14] Hz)';
split_alpha_amp_red_raw.meta.alphaRange_peak = alphaRange;
split_alpha_amp_red_raw.meta.roi = 'occipital O/I channels';
split_alpha_amp_red_raw.meta.split_rule = 'within-subject median across all trials (conditions pooled)';

for ti = 1:numel(tasks)
    tk = tasks(ti);
    task_tag = tk.tag;
    fprintf('\n\n========== PREP RAW ALPHA: %s ==========\n', upper(task_tag));

    fig_dir_subj = fullfile(fig_dir_prep, 'single_subjects', task_tag);
    if SAVE_SINGLE_SUBJECT_FIGS && ~isfolder(fig_dir_subj)
        mkdir(fig_dir_subj);
    end

    task_out = struct();
    task_out.tag = task_tag;
    task_out.split_label = split_label;
    task_out.split_var = tk.split_var;
    task_out.split_latency = tk.split_latency;
    task_out.group_lbl_low = tk.group_lbl_low;
    task_out.group_lbl_high = tk.group_lbl_high;
    task_out.subjects = struct( ...
        'ID', {}, 'Trial', {}, 'Condition', {}, 'AlphaRaw', {}, ...
        'median_thr', {}, 'idx_low', {}, 'idx_high', {}, ...
        'trial_ids_low', {}, 'trial_ids_high', {}, ...
        'alpha_tc_low', {}, 'alpha_tc_high', {}, 'time', {}, ...
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
        fprintf('[PREP rawAlpha %s] Subject %d / %d (%s)\n', upper(task_tag), s, numel(subjects), sid_str);

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

            chUse = occ_channels_from_labels(tf.label);
            chIdx = find(ismember(tf.label, chUse));
            if isempty(chIdx)
                chIdx = 1:numel(tf.label);
            end
            tMaskSplit = tf.time >= tk.split_latency(1) & tf.time <= tk.split_latency(2);
            tMaskPlot = tf.time >= plot_latency(1) & tf.time <= plot_latency(2);
            time_plot = tf.time(tMaskPlot);

            iaf_by_cond = nan(numel(tk.cond_codes), 1);
            for c = 1:numel(tk.cond_codes)
                trlIdx = find(condCodes == tk.cond_codes(c));
                [iaf_by_cond(c), ~] = iaf_from_retention_mtmfft(dataTFR, trlIdx, tk.winIAF, chUse, alphaRange);
            end

            nTrials = size(tf.powspctrm, 1);
            alpha_raw = nan(nTrials, 1);
            tc_all = nan(nTrials, numel(time_plot));
            for tr = 1:nTrials
                c = find(tk.cond_codes == condCodes(tr), 1);
                if isempty(c)
                    continue
                end
                band = ersd_alpha_band(iaf_by_cond(c), alphaRange);
                fMask = tf.freq >= band(1) & tf.freq <= band(2);
                if ~any(fMask)
                    continue
                end
                x = tf.powspctrm(tr, chIdx, fMask, tMaskSplit);
                alpha_raw(tr) = mean(x(:), 'omitnan');
                x = squeeze(mean(mean(tf.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                tc_all(tr, :) = x(:)';
            end

            valid = isfinite(alpha_raw);
            if sum(valid) < 4
                error('Too few finite trials for median split (%d).', sum(valid));
            end
            thr = median(alpha_raw(valid), 'omitnan');
            idx_low = valid & (alpha_raw < thr);
            idx_high = valid & (alpha_raw >= thr);
            if ~any(idx_low) || ~any(idx_high)
                error('Median split produced empty group (n_low=%d, n_high=%d).', sum(idx_low), sum(idx_high));
            end

            tc_low = mean(tc_all(idx_low, :), 1, 'omitnan');
            tc_high = mean(tc_all(idx_high, :), 1, 'omitnan');
            tfr_low = average_trials_freq(tf, idx_low);
            tfr_high = average_trials_freq(tf, idx_high);

            condOut = nan(nTrials, 1);
            for c = 1:numel(tk.cond_codes)
                condOut(condCodes == tk.cond_codes(c)) = tk.cond_out(c);
            end

            entry = struct();
            entry.ID = sid;
            entry.Trial = trialIDs(:);
            entry.Condition = condOut(:);
            entry.AlphaRaw = alpha_raw(:);
            entry.median_thr = thr;
            entry.idx_low = idx_low(:);
            entry.idx_high = idx_high(:);
            entry.trial_ids_low = trialIDs(idx_low);
            entry.trial_ids_high = trialIDs(idx_high);
            entry.alpha_tc_low = tc_low(:)';
            entry.alpha_tc_high = tc_high(:)';
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
            if numel(time_plot) == numel(ga_time) && all(abs(time_plot(:) - ga_time) < 1e-9)
                ga_low(s, :) = tc_low;
                ga_high(s, :) = tc_high;
            else
                ga_low(s, :) = interp1(time_plot(:), tc_low(:), ga_time, 'linear', NaN);
                ga_high(s, :) = interp1(time_plot(:), tc_high(:), ga_time, 'linear', NaN);
            end
            n_ok = n_ok + 1;

            if SAVE_SINGLE_SUBJECT_FIGS
                plot_rawalpha_group_timecourse(time_plot, tc_low, [], tc_high, [], ...
                    colors, tk.group_lbl_low, tk.group_lbl_high, fontSize, fig_pos, ...
                    sprintf('%s | subj %s (n_low=%d, n_high=%d)', task_tag, sid_str, sum(idx_low), sum(idx_high)), ...
                    fullfile(fig_dir_subj, sprintf('AOC_split_rawAlpha_prep_%s_subj%s_timecourse.png', task_tag, sid_str)));
            end

            subj_split = rmfield(entry, {'tfr_low', 'tfr_high'});
            save(fullfile(eeg_dir, sprintf('alphaAmpRed_trialSplit_rawAlpha_%s.mat', task_tag)), 'subj_split');

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
    split_alpha_amp_red_raw.(task_tag) = task_out;

    if n_ok >= 2 && ~isempty(ga_time)
        mL = mean(ga_low, 1, 'omitnan');
        mH = mean(ga_high, 1, 'omitnan');
        nL = sum(isfinite(ga_low), 1);
        nH = sum(isfinite(ga_high), 1);
        sL = std(ga_low, 0, 1, 'omitnan') ./ max(sqrt(nL), 1);
        sH = std(ga_high, 0, 1, 'omitnan') ./ max(sqrt(nH), 1);
        plot_rawalpha_group_timecourse(ga_time, mL, sL, mH, sH, ...
            colors, tk.group_lbl_low, tk.group_lbl_high, fontSize, fig_pos, ...
            sprintf('%s grand average raw alpha (n=%d subjects)', task_tag, n_ok), ...
            fullfile(fig_dir_prep, sprintf('AOC_split_rawAlpha_prep_%s_GA_timecourse.png', task_tag)));
    end
    fprintf('Prep %s done: ok=%d fail=%d\n', task_tag, n_ok, n_fail);
end

out_file = fullfile(stats_dir, 'AOC_splitAlphaAmpRed_trialSplits_rawAlpha.mat');
save(out_file, 'split_alpha_amp_red_raw', '-v7.3');
fprintf('\nSaved trial splits to: %s\n', out_file);

%% ========================= PART 2+3: GazeDev + MS per task =========================
for ti = 1:numel(tasks)
tk = tasks(ti);
task_tag = tk.tag;
if ~isfield(split_alpha_amp_red_raw, task_tag)
    warning('Skipping outcomes for %s: missing prep field.', task_tag);
    continue
end
task_split = split_alpha_amp_red_raw.(task_tag);
subj_splits = task_split.subjects;
nSubj = numel(subj_splits);
if nSubj < 2
    warning('Skipping outcomes for %s: fewer than 2 subjects.', task_tag);
    continue
end

uIDs = [subj_splits.ID]';
split_info_str = sprintf('Within-subject median on raw %s (conditions pooled, no baseline)', tk.split_var);
fprintf('\n\n========== OUTCOMES: %s ==========\n%s\n', upper(task_tag), split_info_str);

alpha_mean_low = nan(nSubj, 1);
alpha_mean_high = nan(nSubj, 1);
alpha_thr = nan(nSubj, 1);
n_low = nan(nSubj, 1);
n_high = nan(nSubj, 1);
for s = 1:nSubj
    alpha_mean_low(s) = mean(subj_splits(s).AlphaRaw(subj_splits(s).idx_low), 'omitnan');
    alpha_mean_high(s) = mean(subj_splits(s).AlphaRaw(subj_splits(s).idx_high), 'omitnan');
    alpha_thr(s) = subj_splits(s).median_thr;
    n_low(s) = subj_splits(s).n_low;
    n_high(s) = subj_splits(s).n_high;
end
fprintf('Subjects: %d | mean n_low=%.1f | mean n_high=%.1f\n', nSubj, mean(n_low), mean(n_high));

%% --- GazeDev ---
fig_prefix_g = sprintf('AOC_split_rawAlpha_%s_GazeDev', task_tag);

figure('Position', fig_pos, 'Color', 'w'); hold on
x = (1:nSubj)';
plot(x, alpha_thr, 'k-', 'LineWidth', 1.2);
h_low = scatter(x, alpha_mean_low, 70, colors(1, :), 'filled', 'MarkerFaceAlpha', 0.85);
h_high = scatter(x, alpha_mean_high, 70, colors(3, :), 'filled', 'MarkerFaceAlpha', 0.85);
xlabel('Participant (index)'); ylabel('Mean trial raw alpha');
title(sprintf('Trial-level median split (raw alpha) [%s]', task_tag), 'Interpreter', 'none');
legend([h_low, h_high], {tk.group_lbl_low, tk.group_lbl_high}, 'Location', 'best', 'FontSize', fontSize-4, 'Box', 'off');
set(gca, 'FontSize', fontSize); box off
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir_root, sprintf('%s_inclusion.png', fig_prefix_g)));
close(gcf);

tfr_low_all = {};
tfr_high_all = {};
metrics_Alpha = nan(nSubj, 2);
metrics_Dev = nan(nSubj, 2);
gaze_cfg = init_gaze_tc_cfg();
t_plot = gaze_cfg.t_vec;
Tf = gaze_cfg.n_samp;
idx_viable = (t_plot >= 0) & (t_plot <= 2);
dev_tc = nan(nSubj, 2, Tf);
missing_gaze = {};

gaze_mat_file = fullfile(feat_dir, tk.gaze_trials_file);
Gmat = load(gaze_mat_file, tk.gaze_trials_var);
Tg = struct2table(Gmat.(tk.gaze_trials_var));

fprintf('\n=== Aggregating GazeDev within raw-alpha trial-split groups ===\n');
for s = 1:nSubj
    clc; fprintf('[RAW ALPHA GAZEDEV - %s] Subject %d / %d\n', upper(task_tag), s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder), subj_folder = sid_str; end

    ids_low = subj_splits(s).trial_ids_low(:);
    ids_high = subj_splits(s).trial_ids_high(:);
    metrics_Alpha(s, 1) = alpha_mean_low(s);
    metrics_Alpha(s, 2) = alpha_mean_high(s);

    if isfield(subj_splits(s), 'tfr_low') && ~isempty(subj_splits(s).tfr_low)
        tfr_low_all{end+1} = subj_splits(s).tfr_low; %#ok<AGROW>
        tfr_high_all{end+1} = subj_splits(s).tfr_high; %#ok<AGROW>
    end

    rows = Tg(Tg.ID == sid, :);
    if ~isempty(rows) && ismember('GazeDeviationFullBL', rows.Properties.VariableNames)
        metrics_Dev(s, 1) = mean(rows.GazeDeviationFullBL(ismember(rows.Trial, ids_low)), 'omitnan');
        metrics_Dev(s, 2) = mean(rows.GazeDeviationFullBL(ismember(rows.Trial, ids_high)), 'omitnan');
    end

    et_file = fullfile(feat_dir, subj_folder, 'gaze', tk.et_fname);
    if ~isfile(et_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end
    try
        [tc_low, tc_high] = build_gaze_tc_pct_by_trial_ids(et_file, ids_low, ids_high, gaze_cfg);
        dev_tc(s, 1, :) = tc_low;
        dev_tc(s, 2, :) = tc_high;
    catch ME
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        fprintf('  Warning: GazeDev TC failed for %s (%s)\n', sid_str, ME.message);
    end
end

channels = {};
if ~isempty(tfr_low_all)
    channels = occ_channels_from_labels(tfr_low_all{1}.label);
elseif ~isempty(tfr_high_all)
    channels = occ_channels_from_labels(tfr_high_all{1}.label);
end

color_map_tfr = customcolormap_preset('red-white-blue');
if ~isempty(tfr_low_all) && ~isempty(tfr_high_all)
    fprintf('\n=== Plotting collapsed TFRs (raw power) ===\n');
    plot_group_tfrs_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, color_map_tfr, ...
        fig_dir_root, fig_prefix_g, fig_pos, fontSize, tfr_winsor_cfg, tk.group_lbl_low, tk.group_lbl_high);
    fprintf('\n=== Plotting collapsed topoplots ===\n');
    plot_group_topos_collapsed_paired(tfr_low_all, tfr_high_all, channels, headmodel, tk.topo_latency, ...
        fig_dir_root, fig_prefix_g, fig_pos, fontSize, tk.group_lbl_low, tk.group_lbl_high);
end

fprintf('\n=== Plotting GazeDev rainclouds ===\n');
plot_paired_raincloud_positive(metrics_Alpha(:,1), metrics_Alpha(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'Raw alpha power', fig_dir_root, sprintf('%s_raincloud_alpha.png', fig_prefix_g), fig_pos, fontSize);
plot_paired_raincloud(metrics_Dev(:,1), metrics_Dev(:,2), colors, tk.group_lbl_low, tk.group_lbl_high, ...
    'Gaze deviation', fig_dir_root, sprintf('%s_raincloud_dev.png', fig_prefix_g), fig_pos, fontSize);

close all
fontSizeTC = 40;
rng(123)
t_vec = t_plot;
tc_window_idx = idx_viable;
tc_complete_min_frac = 0.80;
keep_tc = true(nSubj, 1);
for g = 1:2
    Xg = reshape(dev_tc(:, g, :), nSubj, Tf);
    [Xg, keep_g] = preprocess_gaze_subject_tc(Xg, idx_viable, gaze_cfg);
    dev_tc(:, g, :) = reshape(Xg, [nSubj, 1, Tf]);
    frac = mean(isfinite(Xg(:, tc_window_idx)), 2);
    keep_tc = keep_tc & keep_g & (frac >= tc_complete_min_frac);
end
dev_tc(~keep_tc, :, :) = NaN;
n_tc = sum(keep_tc);
fprintf('Included subjects for GazeDev TC: %d / %d\n', n_tc, nSubj);

low_group_timecourses = reshape(dev_tc(keep_tc, 1, :), n_tc, Tf);
high_group_timecourses = reshape(dev_tc(keep_tc, 2, :), n_tc, Tf);

gaze_ylabel = sprintf('Gaze Deviation\nChange [%%]');
cbpt_report_file = fullfile(stats_dir, sprintf('AOC_split_rawAlpha_GazeDev_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_splitERSD_Prep_rawAlpha.m', ...
    'task_tag', task_tag, ...
    'split_label', split_label, ...
    'split_info', split_info_str, ...
    'ersd_var', tk.split_var, ...
    'group_lbl_low', tk.group_lbl_low, ...
    'group_lbl_high', tk.group_lbl_high, ...
    'n_low_split', nSubj, ...
    'n_high_split', nSubj, ...
    'n_low_tc', n_tc, ...
    'n_high_tc', n_tc, ...
    'metric', 'Gaze deviation change [% baseline] (paired within subject)'));

tc_base = sprintf('AOC_split_rawAlpha_%s_GazeDev_timecourse', task_tag);
report_tag = sprintf('%s_%s_GazeDev', task_tag, split_label);
plot_paired_timecourse_CBPT(low_group_timecourses, high_group_timecourses, colors, gaze_ylabel, ...
    tc_base, report_tag, ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
plot_difference_timecourse_CBPT(low_group_timecourses, high_group_timecourses, ...
    sprintf('%s_HighMinusLow', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec);
plot_paired_timecourse_individuals(low_group_timecourses, high_group_timecourses, colors, gaze_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_individuals', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high);

t_win = (t_vec >= 0) & (t_vec <= 2);
dev_sum_low = mean(low_group_timecourses(:, t_win), 2, 'omitnan');
dev_sum_high = mean(high_group_timecourses(:, t_win), 2, 'omitnan');
ids_tc = uIDs(keep_tc);
ID_col = []; Group_col = strings(0,1); Included_col = [];
Alpha_col = []; GazeDev_col = []; GazeSummary_col = [];
for s = 1:nSubj
    for g = 1:2
        ID_col(end+1,1) = uIDs(s); %#ok<AGROW>
        if g == 1
            Group_col(end+1,1) = string(tk.group_lbl_low); %#ok<AGROW>
            Alpha_col(end+1,1) = metrics_Alpha(s,1); %#ok<AGROW>
            GazeDev_col(end+1,1) = metrics_Dev(s,1); %#ok<AGROW>
        else
            Group_col(end+1,1) = string(tk.group_lbl_high); %#ok<AGROW>
            Alpha_col(end+1,1) = metrics_Alpha(s,2); %#ok<AGROW>
            GazeDev_col(end+1,1) = metrics_Dev(s,2); %#ok<AGROW>
        end
        incl = keep_tc(s);
        Included_col(end+1,1) = incl; %#ok<AGROW>
        if incl
            ix = find(ids_tc == uIDs(s), 1);
            if g == 1, GazeSummary_col(end+1,1) = dev_sum_low(ix); %#ok<AGROW>
            else, GazeSummary_col(end+1,1) = dev_sum_high(ix); %#ok<AGROW>
            end
        else
            GazeSummary_col(end+1,1) = NaN; %#ok<AGROW>
        end
    end
end
stats_tbl = table(ID_col, Group_col, Included_col, Alpha_col, GazeDev_col, GazeSummary_col, ...
    'VariableNames', {'ID','Group','Included', tk.split_var, 'GazeDeviationFullBL', 'GazeDev_pct_0_2s'});
csv_out = fullfile(stats_dir, sprintf('AOC_split_rawAlpha_GazeDev_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('GazeDev CSV: %s | missing gaze: %d\n', csv_out, numel(unique(missing_gaze)));

%% --- MS ---
fig_prefix_m = sprintf('AOC_split_rawAlpha_%s_MS', task_tag);
ms_cfg = init_ms_tc_cfg();
metrics_MS = nan(nSubj, 2);
ms_tc = nan(nSubj, 2, ms_cfg.n_samp);
missing_et = {};

fprintf('\n=== Aggregating MS within raw-alpha trial-split groups ===\n');
for s = 1:nSubj
    clc; fprintf('[RAW ALPHA MS - %s] Subject %d / %d\n', upper(task_tag), s, nSubj);
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder), subj_folder = sid_str; end

    ids_low = subj_splits(s).trial_ids_low(:);
    ids_high = subj_splits(s).trial_ids_high(:);

    rows = Tg(Tg.ID == sid, :);
    if ~isempty(rows) && ismember('MSRateFullBL', rows.Properties.VariableNames)
        metrics_MS(s, 1) = mean(rows.MSRateFullBL(ismember(rows.Trial, ids_low)), 'omitnan');
        metrics_MS(s, 2) = mean(rows.MSRateFullBL(ismember(rows.Trial, ids_high)), 'omitnan');
    end

    et_file = fullfile(feat_dir, subj_folder, 'gaze', [tk.et_fname_ms, '.mat']);
    if ~isfile(et_file)
        missing_et{end+1} = sid_str; %#ok<AGROW>
        continue
    end
    try
        [tc_low, tc_high] = build_ms_tc_pct_by_trial_ids(et_file, ids_low, ids_high, ms_cfg);
        ms_tc(s, 1, :) = tc_low;
        ms_tc(s, 2, :) = tc_high;
    catch ME
        missing_et{end+1} = sid_str; %#ok<AGROW>
        fprintf('  Warning: MS TC failed for %s (%s)\n', sid_str, ME.message);
    end
end

fprintf('\n=== Plotting MS rainclouds ===\n');
plot_paired_raincloud(metrics_MS(:, 1), metrics_MS(:, 2), colors, ...
    tk.group_lbl_low, tk.group_lbl_high, 'Microsaccade rate', ...
    fig_dir_root, sprintf('%s_raincloud_ms.png', fig_prefix_m), fig_pos, fontSize);

close all
fprintf('\n=== Preparing microsaccade time courses ===\n');
t_vec = ms_cfg.t_vec;
idx_viable = (t_vec >= 0) & (t_vec <= 2);
ms_ylabel = sprintf('Microsaccade\nRate Change [%%]');
rng(123)
tc_window_idx = idx_viable;
keep_tc = true(nSubj, 1);
for g = 1:2
    Xg = reshape(ms_tc(:, g, :), nSubj, ms_cfg.n_samp);
    [Xg, keep_g] = preprocess_ms_subject_tc(Xg, idx_viable, ms_cfg);
    ms_tc(:, g, :) = reshape(Xg, [nSubj, 1, ms_cfg.n_samp]);
    frac = mean(isfinite(Xg(:, tc_window_idx)), 2);
    keep_tc = keep_tc & keep_g & (frac >= tc_complete_min_frac);
end
ms_tc(~keep_tc, :, :) = NaN;
n_tc = sum(keep_tc);
fprintf('Included subjects for MS TC: %d / %d\n', n_tc, nSubj);

low_group_timecourses = reshape(ms_tc(keep_tc, 1, :), n_tc, ms_cfg.n_samp);
high_group_timecourses = reshape(ms_tc(keep_tc, 2, :), n_tc, ms_cfg.n_samp);

cbpt_report_file = fullfile(stats_dir, sprintf('AOC_split_rawAlpha_MS_%s_%s_CBPT_report.txt', task_tag, split_label));
init_cbpt_report_file(cbpt_report_file, struct( ...
    'script', 'AOC_splitERSD_Prep_rawAlpha.m', ...
    'task_tag', task_tag, ...
    'split_label', split_label, ...
    'split_info', split_info_str, ...
    'ersd_var', tk.split_var, ...
    'group_lbl_low', tk.group_lbl_low, ...
    'group_lbl_high', tk.group_lbl_high, ...
    'n_low_split', nSubj, ...
    'n_high_split', nSubj, ...
    'n_low_tc', n_tc, ...
    'n_high_tc', n_tc, ...
    'metric', 'Microsaccade rate change [% baseline] (paired within subject)'));

tc_base = sprintf('AOC_split_rawAlpha_%s_MS_timecourse', task_tag);
report_tag = sprintf('%s_%s_MS', task_tag, split_label);
plot_paired_timecourse_CBPT(low_group_timecourses, high_group_timecourses, colors, ms_ylabel, ...
    tc_base, report_tag, ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high, cbpt_report_file);
plot_difference_timecourse_CBPT(low_group_timecourses, high_group_timecourses, ...
    sprintf('%s_HighMinusLow', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec);
plot_paired_timecourse_individuals(low_group_timecourses, high_group_timecourses, colors, ms_ylabel, 'Collapsed over conditions', ...
    sprintf('%s_individuals', tc_base), ...
    fig_dir_root, fig_pos, fontSizeTC, t_vec, ...
    tk.group_lbl_low, tk.group_lbl_high);

t_win = (t_vec >= 0) & (t_vec <= 2);
ms_sum_low = mean(low_group_timecourses(:, t_win), 2, 'omitnan');
ms_sum_high = mean(high_group_timecourses(:, t_win), 2, 'omitnan');
ids_tc = uIDs(keep_tc);
ID_col = []; Group_col = strings(0, 1); Included_col = [];
Alpha_col = []; MS_col = []; MSSummary_col = [];
for s = 1:nSubj
    sid = uIDs(s);
    for g = 1:2
        ID_col(end+1, 1) = sid; %#ok<AGROW>
        if g == 1
            Group_col(end+1, 1) = string(tk.group_lbl_low); %#ok<AGROW>
            Alpha_col(end+1, 1) = alpha_mean_low(s); %#ok<AGROW>
            MS_col(end+1, 1) = metrics_MS(s, 1); %#ok<AGROW>
        else
            Group_col(end+1, 1) = string(tk.group_lbl_high); %#ok<AGROW>
            Alpha_col(end+1, 1) = alpha_mean_high(s); %#ok<AGROW>
            MS_col(end+1, 1) = metrics_MS(s, 2); %#ok<AGROW>
        end
        incl = keep_tc(s);
        Included_col(end+1, 1) = incl; %#ok<AGROW>
        if incl
            ix = find(ids_tc == sid, 1);
            if g == 1
                MSSummary_col(end+1, 1) = ms_sum_low(ix); %#ok<AGROW>
            else
                MSSummary_col(end+1, 1) = ms_sum_high(ix); %#ok<AGROW>
            end
        else
            MSSummary_col(end+1, 1) = NaN; %#ok<AGROW>
        end
    end
end
stats_tbl = table(ID_col, Group_col, Included_col, Alpha_col, MS_col, MSSummary_col, ...
    'VariableNames', {'ID', 'Group', 'Included', tk.split_var, 'MSRateFullBL', 'MS_pct_0_2s'});
csv_out = fullfile(stats_dir, sprintf('AOC_split_rawAlpha_MS_%s_%s_stats_input.csv', task_tag, split_label));
writetable(stats_tbl, csv_out);
fprintf('MS CSV: %s | missing ET: %d\n', csv_out, numel(unique(missing_et)));
end

fprintf('\n=== Raw-alpha full pipeline done ===\n');
fprintf('Figures: %s\n', fig_dir_root);

%% ========================= Local Functions =========================

function gaze_cfg = init_gaze_tc_cfg()
gaze_cfg.fsample = 500;
gaze_cfg.screenW = 800;
gaze_cfg.screenH = 600;
gaze_cfg.centreX = 400;
gaze_cfg.centreY = 300;
gaze_cfg.blink_win = 50;
gaze_cfg.bl_win = [-1.5 -0.5];
gaze_cfg.t_win = [-0.5 2];
t_full = gaze_cfg.t_win(1):1/gaze_cfg.fsample:gaze_cfg.t_win(2);
gaze_cfg.t_vec = t_full(2:end);
gaze_cfg.n_samp = numel(gaze_cfg.t_vec);
gaze_cfg.min_trials_per_group = 3;
gaze_cfg.outlier_k_iqr = 1.5;
gaze_cfg.max_interp_gap_sec = 0.50;
gaze_cfg.min_subject_coverage = 0.60;
gaze_cfg.smooth_sec = 0.025;
gaze_cfg.win_sm = max(1, round(gaze_cfg.smooth_sec * gaze_cfg.fsample));
end

function [tc_low, tc_high] = build_gaze_tc_pct_by_trial_ids(et_file, ids_low, ids_high, gaze_cfg)
[D, err_msg] = load_dataETlong_with_tmp_fallback(et_file);
if isempty(D) || ~isfield(D, 'dataETlong') || ~isfield(D.dataETlong, 'trial')
    error('Could not load ET from %s (%s).', et_file, err_msg);
end
et = D.dataETlong;
if ~isfield(et, 'trialinfo') || size(et.trialinfo, 2) < 2
    error('Missing trial IDs in %s.', et_file);
end

t_plot = gaze_cfg.t_vec;
dev_low = [];
dev_high = [];

for trl = 1:numel(et.trial)
    tid = et.trialinfo(trl, 2);
    is_low = ismember(tid, ids_low);
    is_high = ismember(tid, ids_high);
    if ~is_low && ~is_high
        continue
    end
    if trl > numel(et.time)
        continue
    end

    data = double(et.trial{trl});
    t = double(et.time{trl});
    if size(data, 1) < 2 || isempty(t) || numel(t) ~= size(data, 2)
        continue
    end

    data = data(1:min(3, size(data, 1)), :);
    data(2, :) = gaze_cfg.screenH - data(2, :);
    valid_idx = data(1, :) >= 0 & data(1, :) <= gaze_cfg.screenW & ...
        data(2, :) >= 0 & data(2, :) <= gaze_cfg.screenH;
    data = data(:, valid_idx);
    t = t(valid_idx);
    if isempty(t)
        continue
    end
    data = remove_blinks(data, gaze_cfg.blink_win);

    idx_base = (t >= gaze_cfg.bl_win(1)) & (t <= gaze_cfg.bl_win(2));
    if ~any(idx_base)
        continue
    end
    dev = sqrt((data(1, :) - gaze_cfg.centreX).^2 + (data(2, :) - gaze_cfg.centreY).^2);
    gd_base = mean(dev(idx_base), 'omitnan');
    if ~isfinite(gd_base) || gd_base <= 0
        continue
    end
    dev_pct = 100 * (dev - gd_base) ./ gd_base;
    tc_interp = interp1(t, dev_pct, t_plot, 'linear', NaN);
    if is_low
        dev_low(end+1, :) = tc_interp;
    end
    if is_high
        dev_high(end+1, :) = tc_interp;
    end
end

tc_low = nan(1, gaze_cfg.n_samp);
tc_high = nan(1, gaze_cfg.n_samp);
if size(dev_low, 1) >= gaze_cfg.min_trials_per_group
    tc_low = mean(dev_low, 1, 'omitnan');
end
if size(dev_high, 1) >= gaze_cfg.min_trials_per_group
    tc_high = mean(dev_high, 1, 'omitnan');
end
end

function [X, keep_subj] = preprocess_gaze_subject_tc(X, idx_viable, gaze_cfg)
X(~isfinite(X)) = NaN;
med_metric = median(X(:, idx_viable), 2, 'omitnan');
[X, keep_subj] = exclude_outlier_trajectories(X, med_metric, gaze_cfg.outlier_k_iqr);
max_interp_gap_smp = max(1, round(gaze_cfg.max_interp_gap_sec * gaze_cfg.fsample));
for s = 1:size(X, 1)
    X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
end
subj_cov = mean(isfinite(X(:, idx_viable)), 2);
keep_subj = keep_subj & (subj_cov >= gaze_cfg.min_subject_coverage);
X(~keep_subj, :) = NaN;
if gaze_cfg.win_sm > 1
    X = movmean(X, gaze_cfg.win_sm, 2, 'omitnan');
end
end

function [X_keep, keep_subj] = exclude_outlier_trajectories(X, metric, k_iqr)
med_m = median(metric, 'omitnan');
iqr_m = iqr(metric);
if ~isfinite(iqr_m) || iqr_m == 0
    low_m = -inf; high_m = inf;
else
    low_m = med_m - k_iqr * iqr_m;
    high_m = med_m + k_iqr * iqr_m;
end
keep_subj = isfinite(metric) & (metric >= low_m) & (metric <= high_m);
X_keep = X;
X_keep(~keep_subj, :) = NaN;
end

function x = fill_short_nan_gaps(x, max_gap_smp)
if isempty(x), return, end
valid = isfinite(x);
if all(~valid) || all(valid), return, end
n = numel(x); i = 1;
while i <= n
    if valid(i), i = i + 1; continue, end
    j = i;
    while j <= n && ~valid(j), j = j + 1; end
    gap_start = i; gap_end = j - 1; gap_len = gap_end - gap_start + 1;
    left = gap_start - 1; right = gap_end + 1;
    if left >= 1 && right <= n && valid(left) && valid(right) && gap_len <= max_gap_smp
        x(gap_start:gap_end) = interp1([left right], [x(left) x(right)], gap_start:gap_end);
    end
    i = j;
end
end

function [D, err_msg] = load_dataETlong_with_tmp_fallback(et_file)
D = [];
err_msg = '';
try
    D = load(et_file, 'dataETlong');
    return
catch
    err_msg = 'direct load failed';
end
tmp_file = '';
try
    tmp_file = fullfile(tempdir, sprintf('AOC_tmp_%s.mat', char(java.util.UUID.randomUUID)));
    copyfile(et_file, tmp_file, 'f');
    D = load(tmp_file, 'dataETlong');
catch
    err_msg = 'direct and tmp load failed';
    D = [];
end
if ~isempty(tmp_file) && isfile(tmp_file)
    delete(tmp_file);
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end+1} = lab;
    end
end
if isempty(ch), ch = labels; end
end

function subj_folder = resolve_subject_folder(subjects, sid)
subj_folder = '';
for i = 1:numel(subjects)
    sval = str2double(subjects{i});
    if isfinite(sval) && sval == sid
        subj_folder = subjects{i}; return
    end
end
end

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm), return, end
P = freq_in.powspctrm;
if size(P, 1) < 3, return, end
P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function A = mean_over_channels_tfr(T, ch_idx)
P = T.powspctrm;
if isempty(P)
    A = nan(numel(T.freq), numel(T.time)); return
end
Psize = size(P);
chan_dim = find(Psize == numel(T.label), 1, 'first');
if isempty(chan_dim), chan_dim = 2; end
switch ndims(P)
    case 4
        if chan_dim == 2
            A = squeeze(mean(P(:, ch_idx, :, :), 2, 'omitnan'));
        else
            A = squeeze(mean(P(ch_idx, :, :, :), 1, 'omitnan'));
        end
        if ndims(A) == 3, A = squeeze(mean(A, 1, 'omitnan')); end
    case 3
        if chan_dim == 1
            A = squeeze(mean(P(ch_idx, :, :), 1, 'omitnan'));
        else
            A = squeeze(mean(P(:, ch_idx, :), 2, 'omitnan'));
            if size(A, 1) == numel(T.time), A = A'; end
        end
    otherwise
        error('Unsupported TFR dimensionality.');
end
end

function plot_group_tfrs_collapsed_paired(tfr_low, tfr_high, channels, headmodel, color_map, fig_dir, fig_prefix, fig_pos, fsz, winsor_cfg, lblLow, lblHigh)
cfg = []; cfg.keepindividual = 'yes';
ga_low = ft_freqgrandaverage(cfg, tfr_low{:});
ga_high = ft_freqgrandaverage(cfg, tfr_high{:});
if winsor_cfg.enable
    ga_low = winsorize_freq_subjects(ga_low, winsor_cfg.prctile);
    ga_high = winsorize_freq_subjects(ga_high, winsor_cfg.prctile);
end
[~, ch_idx] = ismember(channels, ga_low.label);
ch_idx = ch_idx(ch_idx > 0);
freq_idx = ga_low.freq >= 5 & ga_low.freq <= 30;
time_idx = ga_low.time >= -0.5 & ga_low.time <= 2;
Ared = mean_over_channels_tfr(ga_low, ch_idx);
Aamp = mean_over_channels_tfr(ga_high, ch_idx);
mx = max([max(Ared(freq_idx, time_idx), [], 'all'), max(Aamp(freq_idx, time_idx), [], 'all')]);
if ~isfinite(mx) || mx <= 0, mx = 0.1; end
clim_abs = [0, 0.9*mx];

cfg = [];
cfg.channel = channels; cfg.colorbar = 'no'; cfg.zlim = clim_abs;
cfg.xlim = [-0.5 2]; cfg.ylim = [5 30]; cfg.layout = headmodel.layANThead;
fsz_tfr = round(fsz * 0.8);
figure('Position', [0 0 1512*0.666 982], 'Color', 'w');
ax = subplot(2,1,1); cfg.figure = ax; ft_singleplotTFR(cfg, ga_low);
colormap(ax, color_map); set(ax, 'CLim', clim_abs); cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr); ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr); rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.Label.String = 'Raw Power'; title(ax, lblLow, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');
ax = subplot(2,1,2); cfg.figure = ax; ft_singleplotTFR(cfg, ga_high);
colormap(ax, color_map); set(ax, 'CLim', clim_abs); cbar = colorbar(ax);
xlabel(ax, 'Time [s]', 'FontSize', fsz_tfr); ylabel(ax, 'Frequency [Hz]', 'FontSize', fsz_tfr);
set(ax, 'FontSize', fsz_tfr); rectangle('Position', [0 8 2 6], 'EdgeColor', 'k', 'LineWidth', 2);
cbar.Label.String = 'Raw Power'; title(ax, lblHigh, 'FontSize', fsz_tfr+4, 'Interpreter', 'none');
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_tfr_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function plot_group_topos_collapsed_paired(tfr_low, tfr_high, channels, headmodel, topo_latency, fig_dir, fig_prefix, fig_pos, fsz, lblLow, lblHigh)
cfg_ga = []; cfg_ga.keepindividual = 'yes';
ga_low = ft_freqgrandaverage(cfg_ga, tfr_low{:});
ga_high = ft_freqgrandaverage(cfg_ga, tfr_high{:});
cfg_sel = []; cfg_sel.frequency = [8 14]; cfg_sel.avgoverfreq = 'yes';
ga_low = ft_selectdata(cfg_sel, ga_low);
ga_high = ft_selectdata(cfg_sel, ga_high);
cfg = [];
cfg.layout = headmodel.layANThead;
allch = cfg.layout.label;
cfg.channel = allch(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.marker = 'off'; cfg.highlight = 'on'; cfg.highlightchannel = channels;
cfg.highlightsymbol = '.'; cfg.highlightsize = 10; cfg.figure = 'gcf';
cfg.gridscale = 300; cfg.comment = 'no'; cfg.xlim = topo_latency;
cfg.colormap = rdbu_cmap(64); cfg.zlim = 'maxabs';
figure('Position', fig_pos, 'Color', 'w');
ax = subplot(1,2,1); cfg.figure = ax; ft_topoplotER(cfg, ga_low); colorbar(ax);
set(ax, 'FontSize', fsz); title(ax, sprintf('%s - collapsed', lblLow), 'Interpreter', 'none');
ax = subplot(1,2,2); cfg.figure = ax; ft_topoplotER(cfg, ga_high); colorbar(ax);
set(ax, 'FontSize', fsz); title(ax, sprintf('%s - collapsed', lblHigh), 'Interpreter', 'none');
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, sprintf('%s_topo_collapsedConditions.png', fig_prefix)));
close(gcf);
end

function cmap = rdbu_cmap(n)
rdbu_11 = [33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; ...
    253 219 199; 244 165 130; 214 96 77; 178 24 43] / 255;
x = linspace(0, 1, size(rdbu_11, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, rdbu_11, xi, 'linear');
end

function plot_paired_raincloud(yLow, yHigh, colors, lblLow, lblHigh, ylab, fig_dir, out_name, fig_pos, fsz)
all_vals = [yLow(:); yHigh(:)]; all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals), ymax = 1; else, ymax = max(abs(all_vals))*1.15; if ymax<=0, ymax=1; end, end
figure('Position', fig_pos, 'Color', 'w'); hold on
draw_one_cloud(yLow, 1, colors(1,:), 0.35, 96, 0.45);
draw_one_cloud(yHigh, 2, colors(3,:), 0.35, 96, 0.45);
for i = 1:numel(yLow)
    if isfinite(yLow(i)) && isfinite(yHigh(i))
        plot([1 2], [yLow(i) yHigh(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end
end
yline(0, '--', 'Color', [0.4 0.4 0.4]);
ylim([-ymax ymax]);
set(gca, 'XTick', [1 2], 'XTickLabel', {lblLow, lblHigh}, 'FontSize', fsz-2);
ylabel(ylab, 'Interpreter', 'none');
title('Within-subject trial split (conditions collapsed)', 'Interpreter', 'none');
box off; pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, out_name)); close(gcf);
end

function draw_one_cloud(yvals, xpos, col, box_w, dot_size, dot_alpha)
y = yvals(isfinite(yvals));
if numel(y) < 3, return, end
[f, xi] = ksdensity(y, 'NumPoints', 120);
if max(f) > 0, f = f / max(f) * 0.35; else, f = zeros(size(f)); end
x_den = xpos - 0.08; x_box = xpos + 0.03;
fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25); q3 = prctile(y, 75); med = median(y);
p5 = prctile(y, 5); p95 = prctile(y, 95);
plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [x_box-box_w/2, q1, box_w, max(q3-q1, eps)], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = (box_w/2) * 2 * (rand(numel(y),1)-0.5);
scatter(x_box + jit, y, dot_size, col, 'filled', 'MarkerFaceAlpha', dot_alpha, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end

function plot_paired_timecourse_individuals(low_group_timecourses, high_group_timecourses, colors, ylab, title_tag, out_name, fig_dir, fig_pos, fsz, t_vec, lblLow, lblHigh)
t_plot = t_vec(:);
low_curves = low_group_timecourses;
high_curves = high_group_timecourses;
low_curves(~isfinite(low_curves)) = NaN;
high_curves(~isfinite(high_curves)) = NaN;
colR_light = colors(1,:)*0.35+0.65; colA_light = colors(3,:)*0.35+0.65;
figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; plot(t_plot, low_curves', 'Color', colR_light, 'LineWidth', 0.8);
plot(t_plot, mean(low_curves,1,'omitnan'), 'Color', colors(1,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblLow, size(low_curves,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
nexttile; hold on; plot(t_plot, high_curves', 'Color', colA_light, 'LineWidth', 0.8);
plot(t_plot, mean(high_curves,1,'omitnan'), 'Color', colors(3,:), 'LineWidth', 2.5);
xline(0,'--k'); xlim([-0.5 2]); xlabel('Time [s]'); ylabel(ylab);
title(sprintf('%s (n=%d) - %s', lblHigh, size(high_curves,1), title_tag), 'Interpreter','none');
set(gca,'FontSize',fsz-6); box off
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function plot_paired_timecourse_CBPT(low_group_timecourses, high_group_timecourses, colors, ylab, out_name, report_tag, fig_dir, fig_pos, fsz, t_vec, lblLow, lblHigh, cbpt_report_file)
t_plot = t_vec(:)';
num_subjects_in = size(low_group_timecourses, 1);
[low_cb, high_cb, t_cb, n_pairs, dt_cb, drop_info] = prepare_paired_cbpt_poststim(low_group_timecourses, high_group_timecourses, t_plot);
cohens_d_cb = paired_cohens_dz_curve(low_cb, high_cb);

figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([2 1]); hold on
low_group_mean = mean(low_group_timecourses, 1, 'omitnan');
high_group_mean = mean(high_group_timecourses, 1, 'omitnan');
low_group_sem = std(low_group_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(low_group_timecourses), 1)), 1);
high_group_sem = std(high_group_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(high_group_timecourses), 1)), 1);
low_group_plot = shadedErrorBar(t_plot, low_group_mean, low_group_sem, 'lineProps', {'-'}, 'transparent', true);
high_group_plot = shadedErrorBar(t_plot, high_group_mean, high_group_sem, 'lineProps', {'-'}, 'transparent', true);
set(low_group_plot.mainLine, 'Color', colors(1,:), 'LineWidth', 5);
set(high_group_plot.mainLine, 'Color', colors(3,:), 'LineWidth', 5);
set(low_group_plot.patch, 'FaceColor', colors(1,:), 'FaceAlpha', 0.25);
set(high_group_plot.patch, 'FaceColor', colors(3,:), 'FaceAlpha', 0.25);
set(low_group_plot.edge(1), 'Color', 'none'); set(low_group_plot.edge(2), 'Color', 'none');
set(high_group_plot.edge(1), 'Color', 'none'); set(high_group_plot.edge(2), 'Color', 'none');
xline(0, '--k'); ylabel(ylab); xlim([-0.5 2]);
ylim(ylim_from_mean_sem(low_group_mean, low_group_sem, high_group_mean, high_group_sem));
box off; set(gca, 'FontSize', fsz-4);
leg_p1 = patch(NaN, NaN, colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(1,:), 'LineWidth', 1.5);
leg_p2 = patch(NaN, NaN, colors(3,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(3,:), 'LineWidth', 1.5);
legend([leg_p1 leg_p2], {[' ' lblLow], [' ' lblHigh]}, 'Location', 'best', 'FontSize', fsz*0.75, 'Box', 'off');

nexttile; hold on
n_perm = 10000; alpha_cbpt = 0.05; tail_cbpt = 'twotail';
[clusters, tvals_cl, thr] = ft_cluster_permutation_1d_paired(low_cb, high_cb, n_perm, alpha_cbpt, tail_cbpt, t_cb);
log_cbpt_report(cbpt_report_file, build_cbpt_report_lines(struct( ...
    'tag', report_tag, 'modality', 'gaze', 'nR', n_pairs, 'n_input', num_subjects_in, ...
    'lbl_low', lblLow, 'lbl_high', lblHigh, 'n_perm', n_perm, 'alpha', alpha_cbpt, 'tail', tail_cbpt, ...
    'nT_cb', numel(t_cb), 'clusters', clusters, ...
    'tvals', tvals_cl, 'thr', thr, 't_plot', t_cb, 'dt_cb', dt_cb, 'drop_info', drop_info)));
ylims = ylim_from_effect_curve(cohens_d_cb);
shade_sig_clusters(gca, clusters, t_cb, dt_cb, ylims, alpha_cbpt);
ylim(ylims);
plot(t_cb, cohens_d_cb, 'k-', 'LineWidth', 3.5);
yline(0, '--'); xline(0, '--k');
xlabel('Time [s]'); ylabel('Cohen''s d');
xlim([-0.5 2]); box off; set(gca, 'FontSize', fsz-4);
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function plot_difference_timecourse_CBPT(low_group_timecourses, high_group_timecourses, out_name, fig_dir, fig_pos, fsz, t_vec)
t_plot = t_vec(:)';
difference_timecourses = high_group_timecourses - low_group_timecourses;

figure('Position', fig_pos, 'Color', 'w');
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([2 1]); hold on
difference_mean = mean(difference_timecourses, 1, 'omitnan');
difference_sem = std(difference_timecourses, 0, 1, 'omitnan') ./ max(sqrt(sum(isfinite(difference_timecourses), 1)), 1);
difference_ci = 1.96 * difference_sem;
difference_plot = shadedErrorBar(t_plot, difference_mean, difference_ci, 'lineProps', {'-'}, 'transparent', true);
set(difference_plot.mainLine, 'Color', [0.1 0.1 0.1], 'LineWidth', 5);
set(difference_plot.patch, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.25);
set(difference_plot.edge(1), 'Color', 'none'); set(difference_plot.edge(2), 'Color', 'none');
yline(0, '--');
xline(0, '--k');
xlim([-0.5 2]);
ylim(ylim_from_single_mean_err(difference_mean, difference_ci));
ylabel('High - Low Change [%]');
box off; set(gca, 'FontSize', fsz-4);

nexttile; hold on
n_perm = 10000; alpha_cbpt = 0.05; tail_cbpt = 'twotail';
[low_cb, high_cb, t_cb, ~, dt_cb] = prepare_paired_cbpt_poststim(low_group_timecourses, high_group_timecourses, t_plot);
cohens_d_cb = paired_cohens_dz_curve(low_cb, high_cb);
[clusters, ~, ~] = ft_cluster_permutation_1d_paired(low_cb, high_cb, n_perm, alpha_cbpt, tail_cbpt, t_cb);
ylims = ylim_from_effect_curve(cohens_d_cb);
shade_sig_clusters(gca, clusters, t_cb, dt_cb, ylims, alpha_cbpt);
ylim(ylims);
plot(t_cb, cohens_d_cb, 'k-', 'LineWidth', 3.5);
yline(0, '--');
xline(0, '--k');
xlabel('Time [s]');
ylabel('Cohen''s d');
xlim([-0.5 2]);
box off; set(gca, 'FontSize', fsz-4);
finalize_tc_tiledlayout(tl);
pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, [out_name, '.png']));
close(gcf);
end

function [low_cb, high_cb, t_cb, n_pairs, dt_cb, drop_info] = prepare_paired_cbpt_poststim(low, high, t_vec)
min_pairs = 3;
n_subjects_input = size(low, 1);
t_vec = t_vec(:)';
post = t_vec >= 0 & t_vec <= 2;
low_post = low(:, post);
high_post = high(:, post);
t_post = t_vec(post);
n_bins_post_stim = numel(t_post);
complete = all(isfinite(low_post) & isfinite(high_post), 2);
n_subjects_dropped = sum(~complete);
low_cb = low_post(complete, :);
high_cb = high_post(complete, :);
n_pairs = size(low_cb, 1);
ok_t = sum(isfinite(low_cb) & isfinite(high_cb), 1) >= min_pairs;
n_bins_dropped = sum(~ok_t);
low_cb = low_cb(:, ok_t);
high_cb = high_cb(:, ok_t);
t_cb = t_post(ok_t);
if n_pairs < min_pairs
    error('Fewer than %d subjects with complete post-stimulus data for CBPT.', min_pairs);
end
if isempty(t_cb)
    error('No valid post-stimulus time points remain for CBPT.');
end
dt_cb = mean(diff(t_cb), 'omitnan');
drop_info = struct( ...
    'n_subjects_input', n_subjects_input, ...
    'n_subjects_cbpt', n_pairs, ...
    'n_subjects_dropped', n_subjects_dropped, ...
    'n_bins_post_stim', n_bins_post_stim, ...
    'n_bins_cbpt', numel(t_cb), ...
    'n_bins_dropped', n_bins_dropped, ...
    'min_pairs', min_pairs);
fprintf('CBPT exclusions: subjects %d -> %d (dropped %d); time points %d -> %d (dropped %d)\n', ...
    drop_info.n_subjects_input, drop_info.n_subjects_cbpt, drop_info.n_subjects_dropped, ...
    drop_info.n_bins_post_stim, drop_info.n_bins_cbpt, drop_info.n_bins_dropped);
end

function dz = paired_cohens_dz_curve(low, high, min_pairs)
if nargin < 3, min_pairs = 3; end
nT = size(low, 2);
dz = nan(1, nT);
for t = 1:nT
    paired_ok = isfinite(low(:, t)) & isfinite(high(:, t));
    if sum(paired_ok) < min_pairs, continue, end
    diff_t = high(paired_ok, t) - low(paired_ok, t);
    sd = std(diff_t);
    if sd > 0
        dz(t) = mean(diff_t) / sd;
    end
end
end

function ylims = ylim_from_effect_curve(d_curve)
valid_d = d_curve(isfinite(d_curve));
mx = max(abs(valid_d), [], 'omitnan');
if isempty(valid_d) || ~isfinite(mx) || mx == 0
    ylims = [-0.6 0.6];
else
    ylims = [-max(mx + 0.1, 0.6), max(mx + 0.1, 0.6)];
end
end

function shade_sig_clusters(ax, clusters, t_plot_cb, dt_cb, ylims, alpha_cbpt)
axes(ax); %#ok<LAXES>
sig = false(1, numel(t_plot_cb));
for k = 1:numel(clusters)
    if clusters(k).p < alpha_cbpt
        sig(clusters(k).idx) = true;
    end
end
run_start = [false, diff(sig) == 1];
run_end = [diff(sig) == -1, false];
if ~isempty(sig) && sig(1), run_start(1) = true; end
if ~isempty(sig) && sig(end), run_end(end) = true; end
starts = find(run_start);
ends = find(run_end);
patch_alpha = 0.4 * ~any(sig) + 0.25 * any(sig);
for k = 1:numel(starts)
    t1 = t_plot_cb(starts(k)) - dt_cb / 2;
    t2 = t_plot_cb(ends(k)) + dt_cb / 2;
    patch([t1 t2 t2 t1], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.5 0.5 0.5], ...
        'FaceAlpha', patch_alpha, 'EdgeColor', 'none');
end
end

function finalize_tc_tiledlayout(tl)
drawnow;
tl.OuterPosition = [0.06, 0.04, 0.90, 0.92];
end

function yl = ylim_from_single_mean_err(m, e)
lo = m - e;
hi = m + e;
if isempty(lo) || all(~isfinite(lo))
    yl = [-1, 1];
    return
end
span = max(max(hi, [], 'omitnan') - min(lo, [], 'omitnan'), eps);
pad = 0.10 * span;
yl = [min(lo, [], 'omitnan') - pad, max(hi, [], 'omitnan') + pad];
end

function yl = ylim_from_mean_sem(m1, s1, m2, s2)
lo = min([m1 - s1, m2 - s2], [], 'omitnan');
hi = max([m1 + s1, m2 + s2], [], 'omitnan');
if isempty(lo) || all(~isfinite(lo))
    yl = [-1, 1];
    return
end
span = max(max(hi, [], 'omitnan') - min(lo, [], 'omitnan'), eps);
pad = 0.10 * span;
yl = [min(lo, [], 'omitnan') - pad, max(hi, [], 'omitnan') + pad];
end

function init_cbpt_report_file(report_path, meta)
if isempty(report_path), return, end
if isfile(report_path), delete(report_path); end
lines = {};
lines{end+1} = '=== AOC Split ERS/ERD CBPT Report (PAIRED / trial-level) ===';
lines{end+1} = sprintf('Script: %s', meta.script);
lines{end+1} = sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
lines{end+1} = sprintf('Task: %s | Split: %s', meta.task_tag, meta.split_label);
lines{end+1} = sprintf('ERSD variable: %s', meta.ersd_var);
lines{end+1} = meta.split_info;
lines{end+1} = sprintf('Subjects contributing both groups: n=%d (TC n=%d)', meta.n_low_split, meta.n_low_tc);
lines{end+1} = sprintf('Outcome metric: %s', meta.metric);
lines{end+1} = 'CBPT method: FieldTrip ft_timelockstatistics (montecarlo, cluster maxsum, depsamplesT)';
lines{end+1} = 'Exclusions: no NaN imputation; subjects need complete post-stim time points; time points need >=3 paired subjects';
lines{end+1} = 'Defaults: n_perm=10000, clusteralpha=0.05, two-tailed, post-stimulus [0 2] s, no temporal binning';
lines{end+1} = sprintf('TC QC: subjects included for plotting if >=80%% finite samples in [0 2] s (n=%d of %d split subjects)', meta.n_low_tc, meta.n_low_split);
lines{end+1} = '';
append_lines_to_file(report_path, lines);
end

function lines = build_cbpt_report_lines(R)
lbl_lo = sanitize_label_for_fname(R.lbl_low);
lbl_hi = sanitize_label_for_fname(R.lbl_high);
nExtreme = sum(abs(R.tvals) > R.thr.tcrit & isfinite(R.tvals));
if isempty(R.clusters), maxClMass=0; maxClExtent=0;
else
    maxClMass = max([0, arrayfun(@(k) R.clusters(k).mass, 1:numel(R.clusters))]);
    maxClExtent = max([0, arrayfun(@(k) R.clusters(k).extent, 1:numel(R.clusters))]);
end
lines = {};
if isfield(R, 'n_input')
    n_line = sprintf('n=%d pairs (of %d input)', R.nR, R.n_input);
else
    n_line = sprintf('n=%d pairs', R.nR);
end
lines{end+1} = sprintf('  [%s] %s tcrit=%.2f |t|>tcrit at %d timepts; max mass=%.1f; max extent=%d', ...
    R.tag, n_line, R.thr.tcrit, nExtreme, maxClMass, maxClExtent);
lines{end+1} = sprintf('    Method: n_perm=%d, alpha=%.3f, two-tailed paired (%s vs %s), post-stim [0 2] s, no temporal binning', ...
    R.n_perm, R.alpha, lbl_lo, lbl_hi);
if isfield(R, 'drop_info') && ~isempty(R.drop_info)
    d = R.drop_info;
    lines{end+1} = sprintf('    Exclusions: subjects %d -> %d (dropped %d, %.1f%% kept)', ...
        d.n_subjects_input, d.n_subjects_cbpt, d.n_subjects_dropped, ...
        100 * d.n_subjects_cbpt / max(d.n_subjects_input, 1));
    lines{end+1} = sprintf('    Exclusions: time points %d -> %d (dropped %d, %.1f%% kept; min_pairs=%d)', ...
        d.n_bins_post_stim, d.n_bins_cbpt, d.n_bins_dropped, ...
        100 * d.n_bins_cbpt / max(d.n_bins_post_stim, 1), d.min_pairs);
end
if ~isempty(R.clusters)
    for k = 1:numel(R.clusters)
        idx = R.clusters(k).idx;
        t_start = R.t_plot(idx(1)) - R.dt_cb/2;
        t_end = R.t_plot(idx(end)) + R.dt_cb/2;
        status = 'n.s.'; if R.clusters(k).p < R.alpha, status = 'SIGNIFICANT'; end
        lines{end+1} = sprintf('    Cluster %d: [%.3f, %.3f] s; mass=%.1f; p=%.4f; %s', ...
            k, t_start, t_end, R.clusters(k).mass, R.clusters(k).p, status);
    end
else
    lines{end+1} = '    No clusters formed';
end
lines{end+1} = '';
end

function txt = sanitize_label_for_fname(label)
txt = lower(label);
txt = regexprep(txt, '[^a-z0-9]+', '_');
txt = regexprep(txt, '_+', '_');
txt = regexprep(txt, '^_|_$', '');
end

function log_cbpt_report(report_path, lines)
for i = 1:numel(lines), fprintf('%s\n', lines{i}); end
append_lines_to_file(report_path, lines);
end

function append_lines_to_file(file_path, lines)
if isempty(file_path) || isempty(lines), return, end
fid = fopen(file_path, 'a');
if fid < 0, warning('Could not open CBPT report: %s', file_path); return, end
cleanup = onCleanup(@() fclose(fid));
for i = 1:numel(lines), fprintf(fid, '%s\n', lines{i}); end
end

function [clusters, tvals, thr] = ft_cluster_permutation_1d_paired(low_group_timecourses, high_group_timecourses, nPerm, alpha, tail, t_plot_ds)
if nargin < 5, tail = 'twotail'; end
if nargin < 6, t_plot_ds = []; end
nS = size(low_group_timecourses, 1);
nT = size(low_group_timecourses, 2);
if size(high_group_timecourses, 1) ~= nS || size(high_group_timecourses, 2) ~= nT
    error('Low and high matrices must match for paired CBPT.');
end
if any(~isfinite(low_group_timecourses(:))) || any(~isfinite(high_group_timecourses(:)))
    error('CBPT input must not contain NaN. Use prepare_paired_cbpt_poststim first.');
end
df = nS - 1;
chan_label = 'metric';
if ~isempty(t_plot_ds) && numel(t_plot_ds) == nT
    time_vec = t_plot_ds(:)';
else
    time_vec = (0:nT-1)' / 500;
end
tl1 = struct('label', {{chan_label}}, 'time', time_vec, 'dimord', 'rpt_chan_time');
tl2 = tl1;
tl1.trial = reshape(low_group_timecourses, [nS, 1, nT]);
tl2.trial = reshape(high_group_timecourses, [nS, 1, nT]);
cfg_neigh = struct('label', chan_label, 'neighblabel', {{}});
cfg = struct();
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = cfg_neigh;
cfg.numrandomization = nPerm;
cfg.channel = chan_label;
cfg.latency = 'all';
if strcmpi(tail, 'onetail_pos')
    cfg.tail = 1; cfg.clustertail = 1; cfg.alpha = alpha;
elseif strcmpi(tail, 'onetail_neg')
    cfg.tail = -1; cfg.clustertail = -1; cfg.alpha = alpha;
else
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = alpha / 2;
end
cfg.design = [1:nS, 1:nS; ones(1, nS), 2*ones(1, nS)];
cfg.uvar = 1;
cfg.ivar = 2;
stat = ft_timelockstatistics(cfg, tl1, tl2);
tvals = stat.stat(1, :);
clusters = struct('idx', {}, 'mass', {}, 'extent', {}, 'p', {});
has_neg = isfield(stat, 'negclusters') && ~isempty(stat.negclusters);
has_pos = isfield(stat, 'posclusters') && ~isempty(stat.posclusters);
clist = [];
lmat = zeros(1, nT);
if has_pos
    clist = stat.posclusters(:);
    lmat = stat.posclusterslabelmat;
end
if has_neg
    lmat_neg = stat.negclusterslabelmat;
    npos = numel(clist);
    lmat(lmat_neg > 0) = npos + lmat_neg(lmat_neg > 0);
    if isempty(clist)
        clist = stat.negclusters(:);
    else
        clist = [clist; stat.negclusters(:)];
    end
end
for k = 1:numel(clist)
    idx = find(lmat(1, :) == k);
    if isempty(idx), continue, end
    idx = idx(1):idx(end);
    clusters(end+1).idx = idx; %#ok<AGROW>
    clusters(end).mass = sum(abs(tvals(idx)), 'omitnan');
    clusters(end).extent = numel(idx);
    clusters(end).p = clist(k).prob;
end
if strcmpi(tail, 'onetail_pos') || strcmpi(tail, 'onetail_neg')
    thr.tcrit = tinv(1 - alpha, df);
else
    thr.tcrit = tinv(1 - alpha/2, df);
end
end

function ms_cfg = init_ms_tc_cfg()
ms_cfg.fsample = 500;
ms_cfg.screenW = 800;
ms_cfg.screenH = 600;
ms_cfg.blink_win = 50;
ms_cfg.sigma_ms = 50;
sigma_samp = round(ms_cfg.sigma_ms / (1000 / ms_cfg.fsample));
kHalf = 3 * sigma_samp;
x_kern = -kHalf:kHalf;
ms_cfg.gKernel = exp(-x_kern.^2 / (2 * sigma_samp^2));
ms_cfg.gKernel = ms_cfg.gKernel / sum(ms_cfg.gKernel);
ms_cfg.t_comp = [-1.5 2.5];
ms_cfg.n_comp = round(diff(ms_cfg.t_comp) * ms_cfg.fsample) + 1;
ms_cfg.t_comp_vec = linspace(ms_cfg.t_comp(1), ms_cfg.t_comp(2), ms_cfg.n_comp);
ms_cfg.t_win = [-0.5 2];
[~, crop_start] = min(abs(ms_cfg.t_comp_vec - ms_cfg.t_win(1)));
[~, crop_end] = min(abs(ms_cfg.t_comp_vec - ms_cfg.t_win(2)));
ms_cfg.crop_idx = crop_start:crop_end;
ms_cfg.n_samp = numel(ms_cfg.crop_idx);
ms_cfg.t_vec = ms_cfg.t_comp_vec(ms_cfg.crop_idx);
ms_cfg.bl_win = [-1.5 -0.5];
ms_cfg.bl_idx_comp = ms_cfg.t_comp_vec >= ms_cfg.bl_win(1) & ms_cfg.t_comp_vec <= ms_cfg.bl_win(2);
ms_cfg.min_trials_per_group = 3;
ms_cfg.outlier_k_iqr = 1.5;
ms_cfg.max_interp_gap_sec = 0.50;
ms_cfg.min_subject_coverage = 0.60;
ms_cfg.smooth_sec = 0.025;
ms_cfg.win_sm = max(1, round(ms_cfg.smooth_sec * ms_cfg.fsample));
end

function [tc_low, tc_high] = build_ms_tc_pct_by_trial_ids(et_file, ids_low, ids_high, ms_cfg)
S = load(et_file);
if isfield(S, 'dataETlong')
    et = S.dataETlong;
elseif isfield(S, 'dataet')
    et = S.dataet;
else
    error('No ET struct in %s.', et_file);
end
fsample = ms_cfg.fsample;
n_comp = ms_cfg.n_comp;
spikes_low = [];
spikes_high = [];

for trl = 1:numel(et.trial)
    if size(et.trialinfo, 2) > 1
        tid = et.trialinfo(trl, 2);
    else
        continue
    end
    is_low = ismember(tid, ids_low);
    is_high = ismember(tid, ids_high);
    if ~is_low && ~is_high
        continue
    end

    raw = et.trial{trl};
    t = et.time{trl};
    raw = raw(1:min(3, size(raw, 1)), :);
    raw(2, :) = ms_cfg.screenH - raw(2, :);
    oob = raw(1,:) < 0 | raw(1,:) > ms_cfg.screenW | raw(2,:) < 0 | raw(2,:) > ms_cfg.screenH;
    raw(:, oob) = NaN;
    raw = remove_blinks(raw, ms_cfg.blink_win);

    idx_win = t >= ms_cfg.t_comp(1) & t <= ms_cfg.t_comp(2);
    gx = raw(1, idx_win);
    gy = raw(2, idx_win);
    if sum(isfinite(gx) & isfinite(gy)) < 50
        continue
    end
    velData = [gx; gy];
    [~, msDetails] = detect_microsaccades(fsample, velData, length(gx));
    spikeVec = zeros(1, length(gx));
    if ~isempty(msDetails.Onset)
        onsets = msDetails.Onset(msDetails.Onset >= 1 & msDetails.Onset <= length(gx));
        spikeVec(onsets) = 1;
    end
    if length(spikeVec) >= n_comp
        spikeVec = spikeVec(1:n_comp);
    else
        spikeVec(end+1:n_comp) = 0;
    end
    if is_low
        spikes_low(end+1, :) = spikeVec;
    end
    if is_high
        spikes_high(end+1, :) = spikeVec;
    end
end

tc_low = spikes_to_pct_tc(spikes_low, ms_cfg);
tc_high = spikes_to_pct_tc(spikes_high, ms_cfg);
end

function tc = spikes_to_pct_tc(spikes, ms_cfg)
tc = nan(1, ms_cfg.n_samp);
if size(spikes, 1) < ms_cfg.min_trials_per_group
    return
end
rate = mean(spikes, 1) * ms_cfg.fsample;
smoothed = conv(rate, ms_cfg.gKernel, 'same');
bl_mean = nanmean(smoothed(ms_cfg.bl_idx_comp));
if isfinite(bl_mean) && bl_mean > 0
    smoothed_pct = (smoothed - bl_mean) ./ bl_mean * 100;
else
    smoothed_pct = nan(size(smoothed));
end
tc = smoothed_pct(ms_cfg.crop_idx);
end

function [X, keep_subj] = preprocess_ms_subject_tc(X, idx_viable, ms_cfg)
X(~isfinite(X)) = NaN;
med_metric = median(X(:, idx_viable), 2, 'omitnan');
[X, keep_subj] = exclude_outlier_trajectories(X, med_metric, ms_cfg.outlier_k_iqr);
max_interp_gap_smp = max(1, round(ms_cfg.max_interp_gap_sec * ms_cfg.fsample));
for s = 1:size(X, 1)
    X(s, :) = fill_short_nan_gaps(X(s, :), max_interp_gap_smp);
end
subj_cov = mean(isfinite(X(:, idx_viable)), 2);
keep_subj = keep_subj & (subj_cov >= ms_cfg.min_subject_coverage);
X(~keep_subj, :) = NaN;
if ms_cfg.win_sm > 1
    X = movmean(X, ms_cfg.win_sm, 2, 'omitnan');
end
end


function tfr_avg = average_trials_freq(tf, tr_mask)
tfr_avg = tf;
P = tf.powspctrm(tr_mask, :, :, :);
tfr_avg.powspctrm = squeeze(mean(P, 1, 'omitnan')); % chan x freq x time
if isfield(tfr_avg, 'trialinfo')
    tfr_avg.trialinfo = tf.trialinfo(tr_mask, :);
end
tfr_avg.dimord = 'chan_freq_time';
if isfield(tfr_avg, 'cumtapcnt')
    tfr_avg = rmfield(tfr_avg, 'cumtapcnt');
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
% Alpha band: (IAF-4, IAF+2) when IAF is valid; otherwise fixed alphaRange ([8 14]).
if ~isfinite(IAF)
    band = alphaRange;
else
    band = [IAF - 4, IAF + 2];
end
if any(~isfinite(band)) || band(1) >= band(2)
    band = alphaRange;
end
end

function plot_rawalpha_group_timecourse(t, yLow, eLow, yHigh, eHigh, colors, lblLow, lblHigh, fsz, fig_pos, ttl, out_file)
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
set(gca, 'FontSize', fsz);
xlabel('Time [s]');
ylabel('Raw Alpha Power');
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


function plot_paired_raincloud_positive(yLow, yHigh, colors, lblLow, lblHigh, ylab, fig_dir, out_name, fig_pos, fsz)
all_vals = [yLow(:); yHigh(:)]; all_vals = all_vals(isfinite(all_vals));
figure('Position', fig_pos, 'Color', 'w'); hold on
draw_one_cloud(yLow, 1, colors(1,:), 0.35, 96, 0.45);
draw_one_cloud(yHigh, 2, colors(3,:), 0.35, 96, 0.45);
for i = 1:numel(yLow)
    if isfinite(yLow(i)) && isfinite(yHigh(i))
        plot([1 2], [yLow(i) yHigh(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
    end
end
if isempty(all_vals)
    ylim([0 1]);
else
    ymin = min(all_vals); ymax = max(all_vals);
    pad = max((ymax - ymin) * 0.10, eps);
    ylim([max(0, ymin - pad), ymax + pad]);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {lblLow, lblHigh}, 'FontSize', fsz-2);
ylabel(ylab, 'Interpreter', 'none');
title('Within-subject trial split (conditions collapsed)', 'Interpreter', 'none');
box off; pause(0.05); drawnow;
saveas(gcf, fullfile(fig_dir, out_name)); close(gcf);
end

