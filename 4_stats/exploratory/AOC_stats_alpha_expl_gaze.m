%% AOC Alpha x Gaze (Scan Path Length, Eye Velocity, Deviation, Microsaccades)
% Computes scan path length, mean eye velocity, gaze deviation, and microsaccade
% rate in early, late, and full retention windows from raw eye-tracking data.
% Baselined gaze metrics use the standard AOC baseline window [-1.5 -0.5] s
% (percent change). Alpha is taken from merged subject matrices as
% window-matched ERSD (occipital 8-14 Hz, dB).
%
% Statistical tests (command window only, no figures):
%   Linear mixed models per task, window, and gaze metric:
%       GazeZ ~ AlphaZ * Condition + (1|ID)
%   followed by Benjamini-Hochberg FDR correction across all AlphaZ main effects.
%
% Tasks: Sternberg, N-back

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
feat_dir = paths.features;

fprintf('\n============================================================\n');
fprintf('AOC Alpha x Gaze: Scan Path Length and Eye Velocity\n');
fprintf('Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('============================================================\n');

% Gaze processing parameters (aligned with AOC gaze trial FEX)
fsample = 500;
screenH = 600;
centreX = 400;
centreY = 300;
blink_win = 50;
min_valid_samples = 100;
bounds_x = [0 800];
bounds_y = [0 screenH];
polyOrd = 3;
velZthr = 4;

t_base  = [-1.5 -0.5];
t_early = [0 1];
t_late  = [1 2];
t_full  = [0 2];

windows = {'early', 'late', 'full'};
t_windows = {t_early, t_late, t_full};
alpha_vars = {'ERSD_early', 'ERSD_late', 'ERSD_full'};
spl_fields = {'SPL_earlyBL', 'SPL_lateBL', 'SPL_fullBL'};
vel_fields = {'Vel_earlyBL', 'Vel_lateBL', 'Vel_fullBL'};
dev_fields = {'Dev_earlyBL', 'Dev_lateBL', 'Dev_fullBL'};
ms_fields  = {'MS_earlyBL',  'MS_lateBL',  'MS_fullBL'};

tasks(1).tag = 'sternberg';
tasks(1).merged_file = 'AOC_merged_data_sternberg.mat';
tasks(1).merged_var = 'merged_data_sternberg';
tasks(1).et_file = 'dataET_sternberg';
tasks(1).cond_vals = [2 4 6];
tasks(1).cond_labels = {'WM2', 'WM4', 'WM6'};

tasks(2).tag = 'nback';
tasks(2).merged_file = 'AOC_merged_data_nback.mat';
tasks(2).merged_var = 'merged_data_nback';
tasks(2).et_file = 'dataET_nback';
tasks(2).cond_vals = [1 2 3];
tasks(2).cond_labels = {'1back', '2back', '3back'};

all_glmm_p = [];
all_glmm_labels = {};

for ti = 1:numel(tasks)
    tk = tasks(ti);
    fprintf('\n\n==================== TASK: %s ====================\n', upper(tk.tag));

    %% Load merged alpha (ERSD per window)
    merged_path = fullfile(feat_dir, tk.merged_file);
    if ~isfile(merged_path)
        warning('Missing merged file for %s: %s. Skipping task.', tk.tag, merged_path);
        continue
    end
    S = load(merged_path, tk.merged_var);
    Tmerged = struct2table(S.(tk.merged_var));
    req_alpha = intersect(alpha_vars, Tmerged.Properties.VariableNames);
    if numel(req_alpha) < numel(alpha_vars)
        warning('Task %s: missing ERSD window columns. Found: %s', tk.tag, strjoin(req_alpha, ', '));
    end

    %% Compute gaze metrics per subject
    nSubj = numel(subjects);
    nCond = numel(tk.cond_vals);
    gaze_tbl = table();

    for si = 1:nSubj
        sid = str2double(subjects{si});
        if ~isfinite(sid)
            continue
        end

        et_file = fullfile(feat_dir, subjects{si}, 'gaze', [tk.et_file, '.mat']);
        if ~isfile(et_file)
            fprintf('  Subject %s: missing %s.mat\n', subjects{si}, tk.et_file);
            continue
        end

        G = load(et_file);
        if isfield(G, 'dataETlong')
            dataETlong = G.dataETlong;
        elseif isfield(G, 'dataet')
            dataETlong = G.dataet;
        else
            fprintf('  Subject %s: dataETlong not found in %s\n', subjects{si}, [tk.et_file, '.mat']);
            continue
        end
        nTrials = size(dataETlong.trialinfo, 1);
        cond_code = dataETlong.trialinfo(:, 1) - 20;

        for ci = 1:nCond
            cval = tk.cond_vals(ci);
            trl_idx = find(cond_code == cval);
            if isempty(trl_idx)
                continue
            end

            spl_early_bl = []; spl_late_bl = []; spl_full_bl = [];
            vel_early_bl = []; vel_late_bl = []; vel_full_bl = [];
            dev_early_bl = []; dev_late_bl = []; dev_full_bl = [];
            ms_early_bl  = []; ms_late_bl  = []; ms_full_bl  = [];

            for k = 1:numel(trl_idx)
                trl = trl_idx(k);
                met = compute_trial_gaze_metrics(dataETlong, trl, fsample, screenH, ...
                    centreX, centreY, blink_win, min_valid_samples, bounds_x, bounds_y, ...
                    t_base, t_early, t_late, t_full, polyOrd, velZthr);
                if isempty(met)
                    continue
                end
                spl_early_bl(end+1) = met.spl_early_bl; %#ok<AGROW>
                spl_late_bl(end+1)  = met.spl_late_bl;  %#ok<AGROW>
                spl_full_bl(end+1)  = met.spl_full_bl;  %#ok<AGROW>
                vel_early_bl(end+1) = met.vel_early_bl; %#ok<AGROW>
                vel_late_bl(end+1)  = met.vel_late_bl;  %#ok<AGROW>
                vel_full_bl(end+1)  = met.vel_full_bl;  %#ok<AGROW>
                dev_early_bl(end+1) = met.dev_early_bl; %#ok<AGROW>
                dev_late_bl(end+1)  = met.dev_late_bl;  %#ok<AGROW>
                dev_full_bl(end+1)  = met.dev_full_bl;  %#ok<AGROW>
                ms_early_bl(end+1)  = met.ms_early_bl;  %#ok<AGROW>
                ms_late_bl(end+1)   = met.ms_late_bl;   %#ok<AGROW>
                ms_full_bl(end+1)   = met.ms_full_bl;   %#ok<AGROW>
            end

            row = table();
            row.ID = sid;
            row.Condition = cval;
            row.SPL_earlyBL = mean(spl_early_bl, 'omitnan');
            row.SPL_lateBL  = mean(spl_late_bl, 'omitnan');
            row.SPL_fullBL  = mean(spl_full_bl, 'omitnan');
            row.Vel_earlyBL = mean(vel_early_bl, 'omitnan');
            row.Vel_lateBL  = mean(vel_late_bl, 'omitnan');
            row.Vel_fullBL  = mean(vel_full_bl, 'omitnan');
            row.Dev_earlyBL = mean(dev_early_bl, 'omitnan');
            row.Dev_lateBL  = mean(dev_late_bl, 'omitnan');
            row.Dev_fullBL  = mean(dev_full_bl, 'omitnan');
            row.MS_earlyBL  = mean(ms_early_bl,  'omitnan');
            row.MS_lateBL   = mean(ms_late_bl,   'omitnan');
            row.MS_fullBL   = mean(ms_full_bl,   'omitnan');
            row.N_trials = numel(trl_idx);
            gaze_tbl = [gaze_tbl; row]; %#ok<AGROW>
        end
    end

    if isempty(gaze_tbl)
        warning('No gaze metrics computed for task %s.', tk.tag);
        continue
    end

    %% Merge gaze with alpha from merged table
    alpha_keep = Tmerged(:, [{'ID', 'Condition'}, req_alpha]);
    D = outerjoin(gaze_tbl, alpha_keep, 'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'left');

  % Sanity check against precomputed SPL in merged data (if available)
    spl_merged_cols = {'ScanPathLengthEarlyBL', 'ScanPathLengthLateBL', 'ScanPathLengthFullBL'};
    have_merged_spl = all(ismember(spl_merged_cols, Tmerged.Properties.VariableNames));
    if have_merged_spl
        chk = outerjoin(gaze_tbl(:, {'ID', 'Condition', 'SPL_earlyBL', 'SPL_lateBL', 'SPL_fullBL'}), ...
            Tmerged(:, {'ID', 'Condition', spl_merged_cols{:}}), ...
            'Keys', {'ID', 'Condition'}, 'MergeKeys', true, 'Type', 'inner');
        r_early = corr(chk.SPL_earlyBL, chk.ScanPathLengthEarlyBL, 'rows', 'complete');
        r_late  = corr(chk.SPL_lateBL,  chk.ScanPathLengthLateBL,  'rows', 'complete');
        r_full  = corr(chk.SPL_fullBL,  chk.ScanPathLengthFullBL,  'rows', 'complete');
        fprintf('\n--- SPL sanity check vs merged gaze matrix (Pearson r) ---\n');
        fprintf('  early: r = %.3f | late: r = %.3f | full: r = %.3f (n = %d rows)\n', ...
            r_early, r_late, r_full, height(chk));
    end

    fprintf('\n--- Gaze metric summary (baselined %% change) ---\n');
    fprintf('Rows (subject x condition): %d | Subjects: %d\n', ...
        height(D), numel(unique(D.ID)));
    print_metric_summary(D, spl_fields, 'Scan path length');
    print_metric_summary(D, vel_fields, 'Mean eye velocity');

    %% GLMMs per window (window-matched alpha and gaze)
    fprintf('\n--- Linear mixed models (GazeZ ~ AlphaZ * Condition + (1|ID)) ---\n');
    for wi = 1:numel(windows)
        av = alpha_vars{wi};
        if ~ismember(av, D.Properties.VariableNames)
            continue
        end
        for gi = 1:4
            switch gi
                case 1
                    gvar = spl_fields{wi};
                    glab = 'SPL';
                case 2
                    gvar = vel_fields{wi};
                    glab = 'Velocity';
                case 3
                    gvar = dev_fields{wi};
                    glab = 'Deviation';
                case 4
                    gvar = ms_fields{wi};
                    glab = 'Microsaccades';
            end

            T = D(:, {'ID', 'Condition', av, gvar});
            T = rmmissing(T);
            if height(T) < 10 || numel(unique(T.ID)) < 5
                fprintf('  %s | %s: insufficient complete cases for GLMM (n = %d)\n', ...
                    windows{wi}, glab, height(T));
                continue
            end

            T.AlphaZ = within_subject_z(T.(av), T.ID);
            T.GazeZ  = within_subject_z(T.(gvar), T.ID);
            T.CondCat = categorical(T.Condition);
            T.SubjCat = categorical(T.ID);

            label = sprintf('%s_%s_%s', tk.tag, windows{wi}, glab);
            try
                formula_full = 'GazeZ ~ AlphaZ * CondCat + (1 | SubjCat)';
                formula_red  = 'GazeZ ~ AlphaZ + CondCat + (1 | SubjCat)';
                lme_full = fitlme(T, formula_full);
                lme_red  = fitlme(T, formula_red);
                lrt = compare(lme_red, lme_full);
                use_int = lrt.pValue(2) < 0.05;
                if use_int
                    lme_use = lme_full;
                    int_str = 'interaction retained';
                else
                    lme_use = lme_red;
                    int_str = 'interaction dropped';
                end

                coef_tbl = lme_use.Coefficients;
                a_idx = find(strcmp(lme_use.CoefficientNames, 'AlphaZ'), 1);
                beta_a = NaN; p_a = NaN;
                if ~isempty(a_idx)
                    beta_a = coef_tbl.Estimate(a_idx);
                    p_a = coef_tbl.pValue(a_idx);
                end

                fprintf('\n  [%s] %s window | %s | n = %d obs, %d subjects | %s\n', ...
                    label, windows{wi}, glab, height(T), numel(unique(T.ID)), int_str);
                fprintf('    LRT interaction: chi2(%.0f) = %.3f, p = %.4f\n', ...
                    lrt.DF(2), lrt.LRStat(2), lrt.pValue(2));
                fprintf('    beta(AlphaZ) = %+.4f, p = %.4f\n', beta_a, p_a);
                disp(lme_use);

                all_glmm_p(end+1) = p_a; %#ok<AGROW>
                all_glmm_labels{end+1} = label; %#ok<AGROW>
            catch ME
                fprintf('  [%s] GLMM failed: %s\n', label, ME.message);
            end
        end
    end

    %% Per-condition correlations
    fprintf('\n--- Per-condition Spearman correlations ---\n');
    for ci = 1:nCond
        fprintf('  Condition %s (%d):\n', tk.cond_labels{ci}, tk.cond_vals(ci));
        subc = D(D.Condition == tk.cond_vals(ci), :);
        for wi = 1:numel(windows)
            av = alpha_vars{wi};
            if ~ismember(av, subc.Properties.VariableNames)
                continue
            end
            for gi = 1:2
                if gi == 1, gvar = spl_fields{wi}; glab = 'SPL';
                else,       gvar = vel_fields{wi}; glab = 'Velocity';
                end
                x = subc.(av); y = subc.(gvar);
                vld = isfinite(x) & isfinite(y);
                if sum(vld) < 5
                    continue
                end
                [rs, ps] = corr(x(vld), y(vld), 'type', 'Spearman');
                fprintf('    %s | %s: n = %d, rho = %+.3f, p = %.4f\n', ...
                    windows{wi}, glab, sum(vld), rs, ps);
            end
        end
    end
end

%% FDR across GLMM alpha main effects
if ~isempty(all_glmm_p)
    p_fdr = fdr_bh(all_glmm_p);
    fprintf('\n--- FDR correction (Benjamini-Hochberg) across GLMM alpha main effects ---\n');
    for i = 1:numel(all_glmm_p)
        sig = '';
        if p_fdr(i) < 0.05
            sig = ' *';
        end
        fprintf('  %s: p = %.4f, p_fdr = %.4f%s\n', all_glmm_labels{i}, all_glmm_p(i), p_fdr(i), sig);
    end
end

fprintf('\n=== Analysis complete ===\n');

%% ========================= Local functions =========================

function met = compute_trial_gaze_metrics(dataETlong, trl, fsample, screenH, ...
    centreX, centreY, blink_win, min_valid_samples, bounds_x, bounds_y, ...
    t_base, t_early, t_late, t_full, polyOrd, velZthr)

met = [];
raw_dat = dataETlong.trial{trl};
t = dataETlong.time{trl};
if isempty(raw_dat) || isempty(t)
    return
end

raw_dat = double(raw_dat(1:3, :));
t = double(t(:)');
raw_dat(2, :) = screenH - raw_dat(2, :);
inb = raw_dat(1, :) >= bounds_x(1) & raw_dat(1, :) <= bounds_x(2) & ...
      raw_dat(2, :) >= bounds_y(1) & raw_dat(2, :) <= bounds_y(2);
raw_dat(:, ~inb) = NaN;
raw_dat = remove_blinks(raw_dat, blink_win);

idx_base  = t >= t_base(1)  & t <= t_base(2);
idx_early = t >= t_early(1) & t <= t_early(2);
idx_late  = t >= t_late(1)  & t <= t_late(2);
idx_full  = t >= t_full(1)  & t <= t_full(2);

ok_base  = sum(all(isfinite(raw_dat(1:2, idx_base)),  1)) >= min_valid_samples;
ok_early = sum(all(isfinite(raw_dat(1:2, idx_early)), 1)) >= min_valid_samples;
ok_late  = sum(all(isfinite(raw_dat(1:2, idx_late)),  1)) >= min_valid_samples;
ok_full  = sum(all(isfinite(raw_dat(1:2, idx_full)),  1)) >= min_valid_samples;

% SPL and gaze deviation baseline
if ok_base
    xb = raw_dat(1, idx_base); yb = raw_dat(2, idx_base);
    spl_base = nansum(sqrt(diff(xb).^2 + diff(yb).^2));
    dxb = xb - centreX;
    dyb = yb - centreY;
    dev_base = mean(sqrt(dxb.^2 + dyb.^2), 'omitnan');
else
    spl_base = NaN;
    dev_base = NaN;
end

% Velocity on full available series (linear fill for SG filter)
x_all = raw_dat(1, :);
y_all = raw_dat(2, :);
if any(~isfinite(x_all))
    x_all = fillmissing(x_all, 'linear', 'EndValues', 'nearest');
end
if any(~isfinite(y_all))
    y_all = fillmissing(y_all, 'linear', 'EndValues', 'nearest');
end
[vx, vy] = compute_velocity_sg(x_all, y_all, fsample, polyOrd);
[vx, vy] = clean_velocity_components(vx, vy, velZthr);
vel_all = hypot(vx, vy);

if ok_base
    vel_base = mean(vel_all(idx_base), 'omitnan');
else
    vel_base = NaN;
end

% Microsaccade baseline (events/s)
if ok_base
    vel_base_xy = [xb; yb];
    T_base = sum(all(isfinite(vel_base_xy), 1)) / fsample;
    if T_base > 0
        [~, msb] = detect_microsaccades(fsample, vel_base_xy, size(vel_base_xy, 2));
        ms_rate_base = numel(msb.Onset) / T_base;
        if ~isfinite(ms_rate_base) || ms_rate_base <= 0
            ms_rate_base = NaN;
        end
    else
        ms_rate_base = NaN;
    end
else
    ms_rate_base = NaN;
end

% Window metrics
met.spl_early_bl = window_bl_metric(raw_dat(1, idx_early), raw_dat(2, idx_early), ...
    ok_early, spl_base, 'spl');
met.spl_late_bl = window_bl_metric(raw_dat(1, idx_late), raw_dat(2, idx_late), ...
    ok_late, spl_base, 'spl');
met.spl_full_bl = window_bl_metric(raw_dat(1, idx_full), raw_dat(2, idx_full), ...
    ok_full, spl_base, 'spl');

met.vel_early_bl = window_bl_metric(vel_all(idx_early), [], ok_early, vel_base, 'vel');
met.vel_late_bl  = window_bl_metric(vel_all(idx_late),  [], ok_late,  vel_base, 'vel');
met.vel_full_bl  = window_bl_metric(vel_all(idx_full),  [], ok_full,  vel_base, 'vel');

met.dev_early_bl = window_bl_metric(raw_dat(1, idx_early), raw_dat(2, idx_early), ...
    ok_early, dev_base, 'dev');
met.dev_late_bl = window_bl_metric(raw_dat(1, idx_late), raw_dat(2, idx_late), ...
    ok_late, dev_base, 'dev');
met.dev_full_bl = window_bl_metric(raw_dat(1, idx_full), raw_dat(2, idx_full), ...
    ok_full, dev_base, 'dev');

met.ms_early_bl = window_bl_metric(raw_dat(:, idx_early), [], ok_early, ms_rate_base, 'ms', fsample);
met.ms_late_bl  = window_bl_metric(raw_dat(:, idx_late),  [], ok_late,  ms_rate_base, 'ms', fsample);
met.ms_full_bl  = window_bl_metric(raw_dat(:, idx_full),  [], ok_full,  ms_rate_base, 'ms', fsample);
end

function val_bl = window_bl_metric(x, y, ok_win, base_val, mode, fsample)
val_bl = NaN;
if ~ok_win || ~isfinite(base_val) || base_val <= 0
    return
end
switch mode
    case 'spl'
        if nargin < 2 || isempty(y)
            return
        end
        val = nansum(sqrt(diff(x).^2 + diff(y).^2));
    case 'vel'
        val = mean(x, 'omitnan');
    case 'dev'
        if nargin < 2 || isempty(y)
            return
        end
        dx = x - mean(x, 'omitnan');
        dy = y - mean(y, 'omitnan');
        val = mean(sqrt(dx.^2 + dy.^2), 'omitnan');
    case 'ms'
        if nargin < 6 || isempty(fsample)
            error('window_bl_metric:fsampleRequired', 'fsample required for ms mode.');
        end
        XY = x;
        if size(XY, 1) ~= 2
            XY = XY(1:2, :);
        end
        T_win = sum(all(isfinite(XY), 1)) / fsample;
        if T_win > 0
            [~, msw] = detect_microsaccades(fsample, XY, size(XY, 2));
            val = numel(msw.Onset) / T_win;
        else
            val = NaN;
        end
    otherwise
        return
end
if isfinite(val)
    val_bl = 100 * (val - base_val) / base_val;
end
end

function z = within_subject_z(vals, ids)
z = nan(size(vals));
u = unique(ids);
for i = 1:numel(u)
    idx = ids == u(i);
    v = vals(idx);
    mu = mean(v, 'omitnan');
    sd = std(v, 'omitnan');
    if sd > 0
        z(idx) = (v - mu) / sd;
    else
        z(idx) = 0;
    end
end
end

function print_metric_summary(D, fields, label)
fprintf('%s:\n', label);
for i = 1:numel(fields)
    v = D.(fields{i});
    v = v(isfinite(v));
    if isempty(v)
        fprintf('  %s: no finite values\n', fields{i});
    else
        fprintf('  %s: mean = %+.2f, SD = %.2f, median = %+.2f (n = %d)\n', ...
            fields{i}, mean(v), std(v), median(v), numel(v));
    end
end
end

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
X = double(X(:)');
Y = double(Y(:)');
if numel(X) ~= numel(Y)
    error('compute_velocity_sg:SizeMismatch', 'X and Y must have same length.');
end

Ts = 1 / fs;
L = numel(X);
framelen = min(21, L);
if mod(framelen, 2) == 0
    framelen = framelen - 1;
end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0
    minLegal = minLegal + 1;
end
if framelen < minLegal
    framelen = minLegal;
end
if framelen > L
    framelen = L - mod(L, 2) + 1;
end
useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1) / Ts) * G(:, 2)';
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end

function [vx, vy] = clean_velocity_components(vx, vy, velZthr)
halfwin = 11;
if numel(vx) > 2 * halfwin
    vx(1:halfwin) = NaN;
    vx(end - halfwin + 1:end) = NaN;
    vy(1:halfwin) = NaN;
    vy(end - halfwin + 1:end) = NaN;
end
zvx = (vx - nanmean(vx)) ./ (nanstd(vx) + eps);
zvy = (vy - nanmean(vy)) ./ (nanstd(vy) + eps);
bad = abs(zvx) > velZthr | abs(zvy) > velZthr;
vx(bad) = NaN;
vy(bad) = NaN;
end

function p_fdr = fdr_bh(p)
p_fdr = nan(size(p));
valid = isfinite(p);
if ~any(valid)
    return
end
pv = p(valid);
pv = pv(:);
m = numel(pv);
[ps, si] = sort(pv);
ranks = (1:m)';
adj = ps .* m ./ ranks;
adj = min(adj, 1);
for i = m - 1:-1:1
    adj(i) = min(adj(i), adj(i + 1));
end
out = nan(m, 1);
out(si) = adj;
p_fdr(valid) = out;
end
