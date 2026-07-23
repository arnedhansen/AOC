%% AOC Split ERSD — Baseline Carryover Diagnostic (N-back + Sternberg)
% Tests whether trial-to-trial alpha alternation can arise from baselining:
%   high alpha on trial t elevates the baseline of trial t+1, so trial t+1
%   looks like strong ERD under single-trial dB baselining, then the pattern
%   flips again. Relevant for the late-timecourse "reversal" near the next
%   trial onset (especially N-back).
%
% For each task and subject this script:
%   1) Recomputes keeptrials TFR (same settings as AOC_splitERSD_Prep.m)
%   2) Extracts occipital IAF-band raw alpha in baseline / early / late / full
%   3) Builds ERSD under single-trial dB AND subject-pooled dB baselines
%   4) Pairs consecutive trials (global Trial ID differ by 1)
%   5) Saves scalars + multi-panel figures under figures/splits/.../diagnostics
%
% Toggle SUBJECT_LIMIT for a quick smoke test (Inf = all subjects).
% Set RELOAD_CACHED = true to skip TFR and only replot from a previous run.

%% Setup
startup
[subjects, paths, colors] = setup('AOC');
addpath(paths.seb_path);

feat_dir = paths.features;
stats_dir = fullfile(paths.splits_stats, 'diagnostics');
fig_dir = fullfile(paths.figures, 'splits', 'SplitERSERD', 'diagnostics', 'baselineCarryover');
if ~isfolder(stats_dir), mkdir(stats_dir); end
if ~isfolder(fig_dir), mkdir(fig_dir); end

SUBJECT_LIMIT = Inf;   % e.g. 3 for a quick test; Inf = all
RELOAD_CACHED = false; % true: load previous .mat and only remake figures
SAVE_INTERMEDIATE = true;

fig_pos = [0 0 1512 982];
fontSize = 28;
baseline_window = [-1.5 -0.5];
alphaRange = [8 14];
plot_latency = [-0.5 2];
edge_window = [1.5 2];  % late epoch edge used for "reversal" zoom

tasks(1).tag = 'sternberg';
tasks(1).data_fname = 'dataEEG_TFR_sternberg.mat';
tasks(1).data_var = 'dataTFR';
tasks(1).cond_codes = [22 24 26];
tasks(1).cond_out = [2 4 6];
tasks(1).toi = -1.5:0.05:3;
tasks(1).winIAF = [1 2];
tasks(1).split_latency = [1 2];
tasks(1).ersd_var = 'ERSD_late';

tasks(2).tag = 'nback';
tasks(2).data_fname = 'dataEEG_TFR_nback.mat';
tasks(2).data_var = 'dataTFR';
tasks(2).cond_codes = [21 22 23];
tasks(2).cond_out = [1 2 3];
tasks(2).toi = -1.5:0.05:2.25;
tasks(2).winIAF = [0 2];
tasks(2).split_latency = [0 2];
tasks(2).ersd_var = 'ERSD_full';

cache_file = fullfile(stats_dir, 'AOC_splitERSD_baselineCarryover_cache.mat');

fprintf('\n=== AOC Baseline Carryover Diagnostic ===\n');
fprintf('SUBJECT_LIMIT = %g | RELOAD_CACHED = %d\n', SUBJECT_LIMIT, RELOAD_CACHED);
fprintf('Figures: %s\n', fig_dir);
fprintf('Stats:    %s\n', stats_dir);

%% Extract or reload
if RELOAD_CACHED
    if ~isfile(cache_file)
        error('RELOAD_CACHED=true but missing cache:\n  %s', cache_file);
    end
    C = load(cache_file, 'carryover');
    carryover = C.carryover;
    fprintf('Loaded cache: %s\n', cache_file);
else
    carryover = struct();
    carryover.meta.created = datestr(now, 31);
    carryover.meta.script = 'AOC_splitERSD_baselineCarryover_diagnostic.m';
    carryover.meta.baseline_window = baseline_window;
    carryover.meta.plot_latency = plot_latency;
    carryover.meta.edge_window = edge_window;
    carryover.meta.note = [ ...
        'Compares raw alpha, single-trial dB ERSD, and subject-pooled dB ERSD ', ...
        'for lag-1 consecutive trials (Trial ID differ by 1).'];

    nSubjUse = min(numel(subjects), SUBJECT_LIMIT);
    for ti = 1:numel(tasks)
        tk = tasks(ti);
        task_tag = tk.tag;
        fprintf('\n\n========== TASK: %s ==========\n', upper(task_tag));

        trial_rows = repmat(struct( ...
            'ID', [], 'Trial', [], 'Condition', [], ...
            'raw_bl', [], 'raw_early', [], 'raw_late', [], 'raw_full', [], 'raw_split', [], ...
            'ersd_st_early', [], 'ersd_st_late', [], 'ersd_st_full', [], 'ersd_st_split', [], ...
            'ersd_pool_early', [], 'ersd_pool_late', [], 'ersd_pool_full', [], 'ersd_pool_split', [], ...
            'is_high_st', [], 'is_high_pool', [], 'is_high_raw', []), 0, 1);

        pair_rows = repmat(struct( ...
            'ID', [], 'Trial_t', [], 'Trial_tp1', [], 'Condition_t', [], 'Condition_tp1', [], ...
            'raw_split_t', [], 'raw_bl_tp1', [], 'raw_split_tp1', [], ...
            'ersd_st_split_t', [], 'ersd_st_split_tp1', [], ...
            'ersd_pool_split_t', [], 'ersd_pool_split_tp1', [], ...
            'is_high_st_t', [], 'is_high_st_tp1', [], ...
            'is_high_pool_t', [], 'is_high_pool_tp1', [], ...
            'is_high_raw_t', [], 'is_high_raw_tp1', []), 0, 1);

        % Subject-level next-trial timecourses conditioned on previous High/Low
        ga_time = [];
        tc_next_st_after_low = [];
        tc_next_st_after_high = [];
        tc_next_pool_after_low = [];
        tc_next_pool_after_high = [];
        tc_next_raw_after_low = [];
        tc_next_raw_after_high = [];
        n_ok = 0;
        n_fail = 0;

        for s = 1:nSubjUse
            sid_str = subjects{s};
            sid = str2double(sid_str);
            clc
            fprintf('[CARRYOVER %s] Subject %d / %d (%s)\n', upper(task_tag), s, nSubjUse, sid_str);

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

                [tf_pool, ~] = apply_subject_pooled_db_baseline(tf, baseline_window);
                tf_st = apply_single_trial_db_baseline(tf, baseline_window);

                chUse = occ_channels_from_labels(tf.label);
                chIdx = find(ismember(tf.label, chUse));
                if isempty(chIdx)
                    chIdx = 1:numel(tf.label);
                end

                tMaskBl = tf.time >= baseline_window(1) & tf.time <= baseline_window(2);
                tMaskEarly = tf.time >= 0 & tf.time <= 1;
                tMaskLate = tf.time >= 1 & tf.time <= 2;
                tMaskFull = tf.time >= 0 & tf.time <= 2;
                tMaskSplit = tf.time >= tk.split_latency(1) & tf.time <= tk.split_latency(2);
                tMaskPlot = tf.time >= plot_latency(1) & tf.time <= plot_latency(2);
                time_plot = tf.time(tMaskPlot);

                iaf_by_cond = nan(numel(tk.cond_codes), 1);
                for c = 1:numel(tk.cond_codes)
                    trlIdx = find(condCodes == tk.cond_codes(c));
                    [iaf_by_cond(c), ~] = iaf_from_retention_mtmfft(dataTFR, trlIdx, tk.winIAF, chUse, alphaRange);
                end

                nTrials = size(tf.powspctrm, 1);
                raw_bl = nan(nTrials, 1);
                raw_early = nan(nTrials, 1);
                raw_late = nan(nTrials, 1);
                raw_full = nan(nTrials, 1);
                raw_split = nan(nTrials, 1);
                ersd_st_early = nan(nTrials, 1);
                ersd_st_late = nan(nTrials, 1);
                ersd_st_full = nan(nTrials, 1);
                ersd_st_split = nan(nTrials, 1);
                ersd_pool_early = nan(nTrials, 1);
                ersd_pool_late = nan(nTrials, 1);
                ersd_pool_full = nan(nTrials, 1);
                ersd_pool_split = nan(nTrials, 1);
                tc_raw = nan(nTrials, numel(time_plot));
                tc_st = nan(nTrials, numel(time_plot));
                tc_pool = nan(nTrials, numel(time_plot));
                condOut = nan(nTrials, 1);

                for tr = 1:nTrials
                    c = find(tk.cond_codes == condCodes(tr), 1);
                    if isempty(c)
                        continue
                    end
                    condOut(tr) = tk.cond_out(c);
                    band = ersd_alpha_band(iaf_by_cond(c), alphaRange);
                    fMask = tf.freq >= band(1) & tf.freq <= band(2);
                    if ~any(fMask)
                        continue
                    end

                    raw_bl(tr) = mean_roi_band(tf.powspctrm(tr, :, :, :), chIdx, fMask, tMaskBl);
                    raw_early(tr) = mean_roi_band(tf.powspctrm(tr, :, :, :), chIdx, fMask, tMaskEarly);
                    raw_late(tr) = mean_roi_band(tf.powspctrm(tr, :, :, :), chIdx, fMask, tMaskLate);
                    raw_full(tr) = mean_roi_band(tf.powspctrm(tr, :, :, :), chIdx, fMask, tMaskFull);
                    raw_split(tr) = mean_roi_band(tf.powspctrm(tr, :, :, :), chIdx, fMask, tMaskSplit);

                    ersd_st_early(tr) = mean_roi_band(tf_st.powspctrm(tr, :, :, :), chIdx, fMask, tMaskEarly);
                    ersd_st_late(tr) = mean_roi_band(tf_st.powspctrm(tr, :, :, :), chIdx, fMask, tMaskLate);
                    ersd_st_full(tr) = mean_roi_band(tf_st.powspctrm(tr, :, :, :), chIdx, fMask, tMaskFull);
                    ersd_st_split(tr) = mean_roi_band(tf_st.powspctrm(tr, :, :, :), chIdx, fMask, tMaskSplit);

                    ersd_pool_early(tr) = mean_roi_band(tf_pool.powspctrm(tr, :, :, :), chIdx, fMask, tMaskEarly);
                    ersd_pool_late(tr) = mean_roi_band(tf_pool.powspctrm(tr, :, :, :), chIdx, fMask, tMaskLate);
                    ersd_pool_full(tr) = mean_roi_band(tf_pool.powspctrm(tr, :, :, :), chIdx, fMask, tMaskFull);
                    ersd_pool_split(tr) = mean_roi_band(tf_pool.powspctrm(tr, :, :, :), chIdx, fMask, tMaskSplit);

                    tc_raw(tr, :) = squeeze(mean(mean(tf.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                    tc_st(tr, :) = squeeze(mean(mean(tf_st.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                    tc_pool(tr, :) = squeeze(mean(mean(tf_pool.powspctrm(tr, chIdx, fMask, tMaskPlot), 2, 'omitnan'), 3, 'omitnan'));
                end

                valid = isfinite(raw_split) & isfinite(ersd_st_split) & isfinite(ersd_pool_split);
                if sum(valid) < 4
                    error('Too few finite trials (%d).', sum(valid));
                end

                thr_raw = median(raw_split(valid), 'omitnan');
                thr_st = median(ersd_st_split(valid), 'omitnan');
                thr_pool = median(ersd_pool_split(valid), 'omitnan');
                is_high_raw = valid & (raw_split >= thr_raw);
                is_high_st = valid & (ersd_st_split >= thr_st);
                is_high_pool = valid & (ersd_pool_split >= thr_pool);
                is_low_raw = valid & (raw_split < thr_raw);
                is_low_st = valid & (ersd_st_split < thr_st);
                is_low_pool = valid & (ersd_pool_split < thr_pool);

                for tr = 1:nTrials
                    if ~valid(tr)
                        continue
                    end
                    row = struct();
                    row.ID = sid;
                    row.Trial = trialIDs(tr);
                    row.Condition = condOut(tr);
                    row.raw_bl = raw_bl(tr);
                    row.raw_early = raw_early(tr);
                    row.raw_late = raw_late(tr);
                    row.raw_full = raw_full(tr);
                    row.raw_split = raw_split(tr);
                    row.ersd_st_early = ersd_st_early(tr);
                    row.ersd_st_late = ersd_st_late(tr);
                    row.ersd_st_full = ersd_st_full(tr);
                    row.ersd_st_split = ersd_st_split(tr);
                    row.ersd_pool_early = ersd_pool_early(tr);
                    row.ersd_pool_late = ersd_pool_late(tr);
                    row.ersd_pool_full = ersd_pool_full(tr);
                    row.ersd_pool_split = ersd_pool_split(tr);
                    row.is_high_st = double(is_high_st(tr));
                    row.is_high_pool = double(is_high_pool(tr));
                    row.is_high_raw = double(is_high_raw(tr));
                    trial_rows(end + 1, 1) = row; %#ok<AGROW>
                end

                % Consecutive pairs by global Trial ID
                [trial_sorted, ord] = sort(trialIDs(:));
                dTrial = diff(trial_sorted);
                pair_idx = find(dTrial == 1);
                for p = 1:numel(pair_idx)
                    i0 = ord(pair_idx(p));
                    i1 = ord(pair_idx(p) + 1);
                    if ~(valid(i0) && valid(i1))
                        continue
                    end
                    prow = struct();
                    prow.ID = sid;
                    prow.Trial_t = trialIDs(i0);
                    prow.Trial_tp1 = trialIDs(i1);
                    prow.Condition_t = condOut(i0);
                    prow.Condition_tp1 = condOut(i1);
                    prow.raw_split_t = raw_split(i0);
                    prow.raw_bl_tp1 = raw_bl(i1);
                    prow.raw_split_tp1 = raw_split(i1);
                    prow.ersd_st_split_t = ersd_st_split(i0);
                    prow.ersd_st_split_tp1 = ersd_st_split(i1);
                    prow.ersd_pool_split_t = ersd_pool_split(i0);
                    prow.ersd_pool_split_tp1 = ersd_pool_split(i1);
                    prow.is_high_st_t = double(is_high_st(i0));
                    prow.is_high_st_tp1 = double(is_high_st(i1));
                    prow.is_high_pool_t = double(is_high_pool(i0));
                    prow.is_high_pool_tp1 = double(is_high_pool(i1));
                    prow.is_high_raw_t = double(is_high_raw(i0));
                    prow.is_high_raw_tp1 = double(is_high_raw(i1));
                    pair_rows(end + 1, 1) = prow; %#ok<AGROW>
                end

                % Next-trial timecourses conditioned on previous High/Low (split metric)
                if isempty(ga_time)
                    ga_time = time_plot(:)';
                    nAlloc = nSubjUse;
                    nT = numel(ga_time);
                    tc_next_st_after_low = nan(nAlloc, nT);
                    tc_next_st_after_high = nan(nAlloc, nT);
                    tc_next_pool_after_low = nan(nAlloc, nT);
                    tc_next_pool_after_high = nan(nAlloc, nT);
                    tc_next_raw_after_low = nan(nAlloc, nT);
                    tc_next_raw_after_high = nan(nAlloc, nT);
                end

                next_st_L = [];
                next_st_H = [];
                next_pool_L = [];
                next_pool_H = [];
                next_raw_L = [];
                next_raw_H = [];
                for p = 1:numel(pair_idx)
                    i0 = ord(pair_idx(p));
                    i1 = ord(pair_idx(p) + 1);
                    if ~(valid(i0) && valid(i1))
                        continue
                    end
                    if is_low_st(i0)
                        next_st_L = [next_st_L; tc_st(i1, :)]; %#ok<AGROW>
                    elseif is_high_st(i0)
                        next_st_H = [next_st_H; tc_st(i1, :)]; %#ok<AGROW>
                    end
                    if is_low_pool(i0)
                        next_pool_L = [next_pool_L; tc_pool(i1, :)]; %#ok<AGROW>
                    elseif is_high_pool(i0)
                        next_pool_H = [next_pool_H; tc_pool(i1, :)]; %#ok<AGROW>
                    end
                    if is_low_raw(i0)
                        next_raw_L = [next_raw_L; tc_raw(i1, :)]; %#ok<AGROW>
                    elseif is_high_raw(i0)
                        next_raw_H = [next_raw_H; tc_raw(i1, :)]; %#ok<AGROW>
                    end
                end

                if ~isempty(next_st_L)
                    tc_next_st_after_low(s, :) = mean(next_st_L, 1, 'omitnan');
                end
                if ~isempty(next_st_H)
                    tc_next_st_after_high(s, :) = mean(next_st_H, 1, 'omitnan');
                end
                if ~isempty(next_pool_L)
                    tc_next_pool_after_low(s, :) = mean(next_pool_L, 1, 'omitnan');
                end
                if ~isempty(next_pool_H)
                    tc_next_pool_after_high(s, :) = mean(next_pool_H, 1, 'omitnan');
                end
                if ~isempty(next_raw_L)
                    tc_next_raw_after_low(s, :) = mean(next_raw_L, 1, 'omitnan');
                end
                if ~isempty(next_raw_H)
                    tc_next_raw_after_high(s, :) = mean(next_raw_H, 1, 'omitnan');
                end

                n_ok = n_ok + 1;
            catch ME
                fprintf('  FAILED %s: %s\n', sid_str, ME.message);
                n_fail = n_fail + 1;
            end
        end

        task_out = struct();
        task_out.tag = task_tag;
        task_out.ersd_var = tk.ersd_var;
        task_out.split_latency = tk.split_latency;
        task_out.n_subjects_ok = n_ok;
        task_out.n_subjects_fail = n_fail;
        task_out.trials = trial_rows;
        task_out.pairs = pair_rows;
        task_out.ga_time = ga_time;
        task_out.tc_next_st_after_low = tc_next_st_after_low;
        task_out.tc_next_st_after_high = tc_next_st_after_high;
        task_out.tc_next_pool_after_low = tc_next_pool_after_low;
        task_out.tc_next_pool_after_high = tc_next_pool_after_high;
        task_out.tc_next_raw_after_low = tc_next_raw_after_low;
        task_out.tc_next_raw_after_high = tc_next_raw_after_high;
        task_out.summary = summarize_carryover(trial_rows, pair_rows, edge_window, ga_time, ...
            tc_next_st_after_low, tc_next_st_after_high, ...
            tc_next_pool_after_low, tc_next_pool_after_high, ...
            tc_next_raw_after_low, tc_next_raw_after_high);
        carryover.(task_tag) = task_out;
        fprintf('Task %s done: ok=%d fail=%d | n_trials=%d n_pairs=%d\n', ...
            task_tag, n_ok, n_fail, numel(trial_rows), numel(pair_rows));
    end

    if SAVE_INTERMEDIATE
        save(cache_file, 'carryover', '-v7.3');
        fprintf('\nSaved cache: %s\n', cache_file);
    end
end

%% Figures + exported summary
summary_rows = {};
for ti = 1:numel(tasks)
    task_tag = tasks(ti).tag;
    if ~isfield(carryover, task_tag)
        warning('Missing task field: %s', task_tag);
        continue
    end
    T = carryover.(task_tag);
    if isempty(T.pairs)
        warning('No consecutive pairs for %s', task_tag);
        continue
    end
    col_low = colors(1, :);
    col_high = colors(3, :);
    fig_prefix = sprintf('AOC_baselineCarryover_%s', task_tag);

    P = struct2table(T.pairs);
    Tr = struct2table(T.trials);
    tvec = T.ga_time(:);

    % --- Fig 1: mechanism scatters (3 panels) ---
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on
    scatter(P.raw_split_t, P.raw_bl_tp1, 18, [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.25);
    r1 = corr(P.raw_split_t, P.raw_bl_tp1, 'rows', 'complete');
    lsline_safe(P.raw_split_t, P.raw_bl_tp1);
    xlabel(sprintf('Raw alpha split window trial t'));
    ylabel('Raw alpha baseline trial t+1');
    title(sprintf('Carry into baseline (r=%.3f)', r1), 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.55); box off

    nexttile; hold on
    scatter(P.raw_bl_tp1, P.ersd_st_split_tp1, 18, [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.25);
    r2 = corr(P.raw_bl_tp1, P.ersd_st_split_tp1, 'rows', 'complete');
    lsline_safe(P.raw_bl_tp1, P.ersd_st_split_tp1);
    xlabel('Raw alpha baseline trial t+1');
    ylabel(sprintf('Single-trial %s trial t+1', T.ersd_var));
    title(sprintf('Baseline -> single-trial ERSD (r=%.3f)', r2), 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.55); box off

    nexttile; hold on
    scatter(P.ersd_st_split_t, P.ersd_st_split_tp1, 18, [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.25);
    r3 = corr(P.ersd_st_split_t, P.ersd_st_split_tp1, 'rows', 'complete');
    lsline_safe(P.ersd_st_split_t, P.ersd_st_split_tp1);
    xlabel(sprintf('Single-trial %s trial t', T.ersd_var));
    ylabel(sprintf('Single-trial %s trial t+1', T.ersd_var));
    title(sprintf('Lag-1 single-trial ERSD (r=%.3f)', r3), 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.55); box off

    sgtitle(sprintf('%s | mechanism / lag-1 scalars', task_tag), 'Interpreter', 'none', 'FontSize', fontSize * 0.7);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_01_mechanism_scatters.png', fig_prefix)));
    close(gcf);

    % --- Fig 2: lag-1 under three baselining modes ---
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_lag1_panel(P.raw_split_t, P.raw_split_tp1, 'Raw alpha', fontSize);
    plot_lag1_panel(P.ersd_st_split_t, P.ersd_st_split_tp1, 'Single-trial dB ERSD', fontSize);
    plot_lag1_panel(P.ersd_pool_split_t, P.ersd_pool_split_tp1, 'Pooled dB ERSD', fontSize);
    sgtitle(sprintf('%s | lag-1 by baseline mode', task_tag), 'Interpreter', 'none', 'FontSize', fontSize * 0.7);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_02_lag1_by_baseline.png', fig_prefix)));
    close(gcf);

    % --- Fig 3: transition matrices (High/Low) ---
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_transition_heatmap(P.is_high_raw_t, P.is_high_raw_tp1, 'Raw split', fontSize);
    plot_transition_heatmap(P.is_high_st_t, P.is_high_st_tp1, 'Single-trial ERSD', fontSize);
    plot_transition_heatmap(P.is_high_pool_t, P.is_high_pool_tp1, 'Pooled ERSD', fontSize);
    sgtitle(sprintf('%s | High/Low transition probabilities', task_tag), 'Interpreter', 'none', 'FontSize', fontSize * 0.7);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_03_transition_matrices.png', fig_prefix)));
    close(gcf);

    % --- Fig 4: next-trial timecourses after High vs Low ---
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_conditioned_tc(tvec, T.tc_next_st_after_low, T.tc_next_st_after_high, ...
        col_low, col_high, 'Next trial | single-trial dB', 'Power change [dB]', fontSize, edge_window);
    plot_conditioned_tc(tvec, T.tc_next_pool_after_low, T.tc_next_pool_after_high, ...
        col_low, col_high, 'Next trial | pooled dB', 'Power change [dB]', fontSize, edge_window);
    plot_conditioned_tc(tvec, T.tc_next_raw_after_low, T.tc_next_raw_after_high, ...
        col_low, col_high, 'Next trial | raw power', 'Raw power [a.u.]', fontSize, edge_window);
    sgtitle(sprintf('%s | next-trial TC after previous Low/High (%s)', task_tag, T.ersd_var), ...
        'Interpreter', 'none', 'FontSize', fontSize * 0.65);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_04_nextTrial_timecourses.png', fig_prefix)));
    close(gcf);

    % --- Fig 5: High-minus-Low next-trial difference (edge focus) ---
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_diff_tc(tvec, T.tc_next_st_after_high, T.tc_next_st_after_low, ...
        'Single-trial dB: afterHigh - afterLow', fontSize, edge_window);
    plot_diff_tc(tvec, T.tc_next_pool_after_high, T.tc_next_pool_after_low, ...
        'Pooled dB: afterHigh - afterLow', fontSize, edge_window);
    plot_diff_tc(tvec, T.tc_next_raw_after_high, T.tc_next_raw_after_low, ...
        'Raw: afterHigh - afterLow', fontSize, edge_window);
    sgtitle(sprintf('%s | next-trial difference (carryover signature)', task_tag), ...
        'Interpreter', 'none', 'FontSize', fontSize * 0.7);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_05_nextTrial_diff_edge.png', fig_prefix)));
    close(gcf);

    % --- Fig 6: example subject trial sequence ---
    uIDs = unique(Tr.ID);
    nEx = min(3, numel(uIDs));
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(nEx, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for ex = 1:nEx
        nexttile; hold on
        sid = uIDs(ex);
        rows = Tr(Tr.ID == sid, :);
        rows = sortrows(rows, 'Trial');
        yyaxis left
        plot(rows.Trial, rows.raw_split, '-o', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'MarkerSize', 4);
        ylabel('Raw split');
        yyaxis right
        plot(rows.Trial, rows.ersd_st_split, '-o', 'Color', col_high, 'LineWidth', 1.2, 'MarkerSize', 4);
        ylabel('Single-trial ERSD');
        xlabel('Trial ID');
        title(sprintf('Subject %d', sid), 'Interpreter', 'none');
        set(gca, 'FontSize', fontSize * 0.5); box off
    end
    sgtitle(sprintf('%s | example trial sequences (raw vs single-trial ERSD)', task_tag), ...
        'Interpreter', 'none', 'FontSize', fontSize * 0.65);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_06_example_sequences.png', fig_prefix)));
    close(gcf);

    % --- Fig 7: subject-wise lag-1 correlations ---
    figure('Position', fig_pos, 'Color', 'w');
    hold on
    [r_subj_raw, r_subj_st, r_subj_pool, ids_r] = subjectwise_lag1(P);
    x = (1:numel(ids_r))';
    plot(x, r_subj_raw, '-o', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.2 0.2]);
    plot(x, r_subj_st, '-o', 'Color', col_high, 'LineWidth', 1.5, 'MarkerFaceColor', col_high);
    plot(x, r_subj_pool, '-o', 'Color', col_low, 'LineWidth', 1.5, 'MarkerFaceColor', col_low);
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    xlabel('Participant (index)');
    ylabel('Lag-1 Pearson r');
    legend({'Raw', 'Single-trial ERSD', 'Pooled ERSD'}, 'Location', 'best', 'Box', 'off', 'FontSize', fontSize * 0.55);
    title(sprintf('%s | subject-wise lag-1 correlations', task_tag), 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.6); box off
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, sprintf('%s_07_subjectwise_lag1_r.png', fig_prefix)));
    close(gcf);

    % Collect summary
    Ssum = T.summary;
    Ssum.task = task_tag;
    summary_rows{end + 1} = Ssum; %#ok<AGROW>

    fprintf('\n--- %s summary ---\n', upper(task_tag));
    fprintf('  n_pairs=%d\n', Ssum.n_pairs);
    fprintf('  r(raw_t, bl_t+1)=%.3f | r(bl_t+1, ERSD_st_t+1)=%.3f\n', Ssum.r_raw_t_vs_bl_tp1, Ssum.r_bl_tp1_vs_ersd_st_tp1);
    fprintf('  lag1 r raw/st/pool = %.3f / %.3f / %.3f\n', Ssum.r_lag1_raw, Ssum.r_lag1_st, Ssum.r_lag1_pool);
    fprintf('  P(High->Low) raw/st/pool = %.3f / %.3f / %.3f\n', ...
        Ssum.p_high_to_low_raw, Ssum.p_high_to_low_st, Ssum.p_high_to_low_pool);
    fprintf('  edge mean (afterHigh-afterLow) st/pool/raw = %.3f / %.3f / %.3f\n', ...
        Ssum.edge_diff_st, Ssum.edge_diff_pool, Ssum.edge_diff_raw);
end

% --- Fig 8: N-back vs Sternberg comparison bars ---
if numel(summary_rows) >= 2
    figure('Position', fig_pos, 'Color', 'w');
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile
    vals = [summary_rows{1}.r_lag1_st, summary_rows{2}.r_lag1_st; ...
            summary_rows{1}.r_lag1_pool, summary_rows{2}.r_lag1_pool; ...
            summary_rows{1}.r_lag1_raw, summary_rows{2}.r_lag1_raw];
    b = bar(vals);
    b(1).FaceColor = colors(1, :);
    b(2).FaceColor = colors(3, :);
    set(gca, 'XTickLabel', {'ST ERSD', 'Pooled', 'Raw'});
    ylabel('Lag-1 r');
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    legend({summary_rows{1}.task, summary_rows{2}.task}, 'Box', 'off', 'FontSize', fontSize * 0.5);
    title('Lag-1 autocorrelation', 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.55); box off

    nexttile
    vals = [summary_rows{1}.p_high_to_low_st, summary_rows{2}.p_high_to_low_st; ...
            summary_rows{1}.p_high_to_low_pool, summary_rows{2}.p_high_to_low_pool; ...
            summary_rows{1}.p_high_to_low_raw, summary_rows{2}.p_high_to_low_raw];
    b = bar(vals);
    b(1).FaceColor = colors(1, :);
    b(2).FaceColor = colors(3, :);
    set(gca, 'XTickLabel', {'ST ERSD', 'Pooled', 'Raw'});
    ylabel('P(High -> Low)');
    yline(0.5, '--', 'Color', [0.6 0.6 0.6]);
    title('Alternation after High', 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.55); box off

    nexttile
    vals = [summary_rows{1}.r_raw_t_vs_bl_tp1, summary_rows{2}.r_raw_t_vs_bl_tp1; ...
            summary_rows{1}.r_bl_tp1_vs_ersd_st_tp1, summary_rows{2}.r_bl_tp1_vs_ersd_st_tp1; ...
            summary_rows{1}.edge_diff_st, summary_rows{2}.edge_diff_st];
    b = bar(vals);
    b(1).FaceColor = colors(1, :);
    b(2).FaceColor = colors(3, :);
    set(gca, 'XTickLabel', {'raw_t->bl_t+1', 'bl->ERSD_st', 'edge diff ST'});
    ylabel('Effect size');
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    title('Mechanism + edge signature', 'Interpreter', 'none');
    set(gca, 'FontSize', fontSize * 0.5); box off

    sgtitle('N-back vs Sternberg | baseline carryover diagnostics', 'Interpreter', 'none', 'FontSize', fontSize * 0.7);
    drawnow; pause(0.05);
    saveas(gcf, fullfile(fig_dir, 'AOC_baselineCarryover_08_task_comparison.png'));
    close(gcf);
end

% Export summary + pair CSVs
if ~isempty(summary_rows)
    SumT = struct2table([summary_rows{:}]);
    writetable(SumT, fullfile(stats_dir, 'AOC_splitERSD_baselineCarryover_summary.csv'));
    fprintf('\nWrote summary CSV.\n');
end
for ti = 1:numel(tasks)
    task_tag = tasks(ti).tag;
    if ~isfield(carryover, task_tag) || isempty(carryover.(task_tag).pairs)
        continue
    end
    writetable(struct2table(carryover.(task_tag).pairs), ...
        fullfile(stats_dir, sprintf('AOC_splitERSD_baselineCarryover_pairs_%s.csv', task_tag)));
    writetable(struct2table(carryover.(task_tag).trials), ...
        fullfile(stats_dir, sprintf('AOC_splitERSD_baselineCarryover_trials_%s.csv', task_tag)));
end

fprintf('\nDone. Figures in:\n  %s\n', fig_dir);

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

function tf_bl = apply_single_trial_db_baseline(tf, bl_win)
P = tf.powspctrm;
tMask = tf.time >= bl_win(1) & tf.time <= bl_win(2);
if ~any(tMask)
    error('No time samples in baseline window [%.2f %.2f].', bl_win(1), bl_win(2));
end
bl = mean(P(:, :, :, tMask), 4, 'omitnan'); % rpt x chan x freq
bl(~isfinite(bl) | bl <= 0) = NaN;
bl4 = reshape(bl, [size(bl, 1), size(bl, 2), size(bl, 3), 1]);
P_db = 10 * log10(bsxfun(@rdivide, P, bl4));
tf_bl = tf;
tf_bl.powspctrm = P_db;
end

function m = mean_roi_band(P1, chIdx, fMask, tMask)
% P1 is 1 x chan x freq x time (single trial slice already indexed as 4D)
x = P1(1, chIdx, fMask, tMask);
m = mean(x(:), 'omitnan');
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
if ~isfinite(IAF)
    band = alphaRange;
else
    band = [IAF - 4, IAF + 2];
end
if any(~isfinite(band)) || band(1) >= band(2)
    band = alphaRange;
end
end

function S = summarize_carryover(trial_rows, pair_rows, edge_window, ga_time, ...
    tc_st_L, tc_st_H, tc_pool_L, tc_pool_H, tc_raw_L, tc_raw_H)
S = struct();
S.n_trials = numel(trial_rows);
S.n_pairs = numel(pair_rows);
if S.n_pairs < 2
    S.r_raw_t_vs_bl_tp1 = NaN;
    S.r_bl_tp1_vs_ersd_st_tp1 = NaN;
    S.r_lag1_raw = NaN;
    S.r_lag1_st = NaN;
    S.r_lag1_pool = NaN;
    S.p_high_to_low_raw = NaN;
    S.p_high_to_low_st = NaN;
    S.p_high_to_low_pool = NaN;
    S.edge_diff_st = NaN;
    S.edge_diff_pool = NaN;
    S.edge_diff_raw = NaN;
    return
end
P = struct2table(pair_rows);
S.r_raw_t_vs_bl_tp1 = corr(P.raw_split_t, P.raw_bl_tp1, 'rows', 'complete');
S.r_bl_tp1_vs_ersd_st_tp1 = corr(P.raw_bl_tp1, P.ersd_st_split_tp1, 'rows', 'complete');
S.r_lag1_raw = corr(P.raw_split_t, P.raw_split_tp1, 'rows', 'complete');
S.r_lag1_st = corr(P.ersd_st_split_t, P.ersd_st_split_tp1, 'rows', 'complete');
S.r_lag1_pool = corr(P.ersd_pool_split_t, P.ersd_pool_split_tp1, 'rows', 'complete');
S.p_high_to_low_raw = transition_p_high_to_low(P.is_high_raw_t, P.is_high_raw_tp1);
S.p_high_to_low_st = transition_p_high_to_low(P.is_high_st_t, P.is_high_st_tp1);
S.p_high_to_low_pool = transition_p_high_to_low(P.is_high_pool_t, P.is_high_pool_tp1);
S.edge_diff_st = edge_mean_diff(ga_time, tc_st_H, tc_st_L, edge_window);
S.edge_diff_pool = edge_mean_diff(ga_time, tc_pool_H, tc_pool_L, edge_window);
S.edge_diff_raw = edge_mean_diff(ga_time, tc_raw_H, tc_raw_L, edge_window);
end

function p = transition_p_high_to_low(high_t, high_tp1)
high_t = high_t(:);
high_tp1 = high_tp1(:);
ok = isfinite(high_t) & isfinite(high_tp1) & (high_t == 1);
if ~any(ok)
    p = NaN;
    return
end
p = mean(high_tp1(ok) == 0);
end

function d = edge_mean_diff(tvec, tc_after_high, tc_after_low, edge_window)
if isempty(tvec) || isempty(tc_after_high) || isempty(tc_after_low)
    d = NaN;
    return
end
mH = mean(tc_after_high, 1, 'omitnan');
mL = mean(tc_after_low, 1, 'omitnan');
diffv = mH - mL;
mask = tvec(:)' >= edge_window(1) & tvec(:)' <= edge_window(2);
if ~any(mask)
    d = NaN;
    return
end
d = mean(diffv(mask), 'omitnan');
end

function lsline_safe(x, y)
ok = isfinite(x) & isfinite(y);
if sum(ok) < 3
    return
end
p = polyfit(x(ok), y(ok), 1);
xx = linspace(min(x(ok)), max(x(ok)), 50);
plot(xx, polyval(p, xx), 'r-', 'LineWidth', 2);
end

function plot_lag1_panel(x, y, ttl, fontSize)
nexttile; hold on
scatter(x, y, 18, [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.25);
r = corr(x, y, 'rows', 'complete');
lsline_safe(x, y);
xlabel('Trial t');
ylabel('Trial t+1');
title(sprintf('%s (r=%.3f)', ttl, r), 'Interpreter', 'none');
set(gca, 'FontSize', fontSize * 0.55); box off
end

function plot_transition_heatmap(high_t, high_tp1, ttl, fontSize)
nexttile
M = zeros(2, 2);
ok = isfinite(high_t) & isfinite(high_tp1);
% rows: Low/High at t; cols: Low/High at t+1
for a = 0:1
    for b = 0:1
        M(a + 1, b + 1) = sum(ok & high_t == a & high_tp1 == b);
    end
end
rowSum = sum(M, 2);
Pmat = M ./ max(rowSum, 1);
imagesc(Pmat);
colormap(gca, parula);
caxis([0 1]);
colorbar;
set(gca, 'XTick', [1 2], 'XTickLabel', {'Low t+1', 'High t+1'}, ...
    'YTick', [1 2], 'YTickLabel', {'Low t', 'High t'});
title(sprintf('%s', ttl), 'Interpreter', 'none');
for i = 1:2
    for j = 1:2
        text(j, i, sprintf('%.2f\n(n=%d)', Pmat(i, j), M(i, j)), ...
            'HorizontalAlignment', 'center', 'FontSize', fontSize * 0.4, 'Color', 'w');
    end
end
set(gca, 'FontSize', fontSize * 0.5);
end

function plot_conditioned_tc(tvec, tc_after_low, tc_after_high, col_low, col_high, ttl, ylab, fontSize, edge_window)
nexttile; hold on
if isempty(tvec)
    title(ttl, 'Interpreter', 'none');
    return
end
mL = mean(tc_after_low, 1, 'omitnan');
mH = mean(tc_after_high, 1, 'omitnan');
nL = sum(isfinite(tc_after_low), 1);
nH = sum(isfinite(tc_after_high), 1);
eL = std(tc_after_low, 0, 1, 'omitnan') ./ max(sqrt(nL), 1);
eH = std(tc_after_high, 0, 1, 'omitnan') ./ max(sqrt(nH), 1);
ebL = shadedErrorBar(tvec, mL, eL, 'lineProps', {'-', 'Color', col_low});
ebH = shadedErrorBar(tvec, mH, eH, 'lineProps', {'-', 'Color', col_high});
set(ebL.mainLine, 'LineWidth', 2); set(ebH.mainLine, 'LineWidth', 2);
set(ebL.patch, 'FaceColor', col_low, 'FaceAlpha', 0.20);
set(ebH.patch, 'FaceColor', col_high, 'FaceAlpha', 0.20);
set(ebL.edge(1), 'Color', 'none'); set(ebL.edge(2), 'Color', 'none');
set(ebH.edge(1), 'Color', 'none'); set(ebH.edge(2), 'Color', 'none');
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
xline(edge_window(1), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlim([min(tvec) max(tvec)]);
xlabel('Time [s]');
ylabel(ylab);
title(ttl, 'Interpreter', 'none');
leg1 = patch(NaN, NaN, col_low, 'FaceAlpha', 0.25, 'EdgeColor', col_low);
leg2 = patch(NaN, NaN, col_high, 'FaceAlpha', 0.25, 'EdgeColor', col_high);
legend([leg1, leg2], {' after Low t', ' after High t'}, 'Location', 'best', ...
    'FontSize', fontSize * 0.4, 'Box', 'off');
set(gca, 'FontSize', fontSize * 0.5); box off
end

function plot_diff_tc(tvec, tc_after_high, tc_after_low, ttl, fontSize, edge_window)
nexttile; hold on
if isempty(tvec)
    title(ttl, 'Interpreter', 'none');
    return
end
diff_subj = tc_after_high - tc_after_low;
m = mean(diff_subj, 1, 'omitnan');
n = sum(isfinite(diff_subj), 1);
e = std(diff_subj, 0, 1, 'omitnan') ./ max(sqrt(n), 1);
eb = shadedErrorBar(tvec, m, e, 'lineProps', {'-', 'Color', [0.1 0.1 0.1]});
set(eb.mainLine, 'LineWidth', 2);
set(eb.patch, 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.20);
set(eb.edge(1), 'Color', 'none'); set(eb.edge(2), 'Color', 'none');
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
xline(0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
yl = ylim;
patch([edge_window(1) edge_window(2) edge_window(2) edge_window(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.9 0.7 0.2], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
xlim([min(tvec) max(tvec)]);
xlabel('Time [s]');
ylabel('Difference');
title(ttl, 'Interpreter', 'none');
set(gca, 'FontSize', fontSize * 0.5); box off
end

function [r_raw, r_st, r_pool, ids] = subjectwise_lag1(P)
ids = unique(P.ID);
r_raw = nan(numel(ids), 1);
r_st = nan(numel(ids), 1);
r_pool = nan(numel(ids), 1);
for i = 1:numel(ids)
    rows = P(P.ID == ids(i), :);
    if height(rows) < 5
        continue
    end
    r_raw(i) = corr(rows.raw_split_t, rows.raw_split_tp1, 'rows', 'complete');
    r_st(i) = corr(rows.ersd_st_split_t, rows.ersd_st_split_tp1, 'rows', 'complete');
    r_pool(i) = corr(rows.ersd_pool_split_t, rows.ersd_pool_split_tp1, 'rows', 'complete');
end
end
