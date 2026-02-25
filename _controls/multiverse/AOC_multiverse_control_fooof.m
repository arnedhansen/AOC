%% AOC Multiverse — FOOOF R² Quality Control
% Loads per-trial (and per-subject) FOOOF R² CSVs produced by the
% multiverse prep scripts and generates diagnostic scatterplots.
%
% Figures per task (4 × 2 = 8 figures total):
%   1. R² per subject — jittered scatterplot, faceted by latency, colored by electrode set
%   2. R² vs aperiodic exponent — faceted by latency
%   3. R² vs aperiodic offset — faceted by latency
%   4. R² distribution histogram — with threshold line at 0.90
%
% Input:  /Volumes/.../data/controls/multiverse/fooof_r2_{task}.csv
%         /Volumes/.../data/controls/multiverse/fooof_r2_{task}_subject.csv
% Output: /Volumes/.../figures/controls/multiverse/FOOOF/

%% Paths
base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
if ispc
    base_data = 'W:\Students\Arne\AOC';
end
csv_dir = fullfile(base_data, 'data', 'controls', 'multiverse');
multiverse_csv_dir = fullfile(base_data, 'data', 'multiverse');
fig_dir_base = fullfile(base_data, 'figures', 'controls', 'multiverse', 'FOOOF');
if ~isfolder(fig_dir_base), mkdir(fig_dir_base); end

%% Select method
method_selection = 'welch500_50'; % 'singleFFT' | 'welch500_50' | 'both'

%%
r2_thresh = 0.90;
tasks = {'sternberg', 'nback'};
levels = parse_levels('both');
methods = parse_methods(method_selection);

for mi = 1:length(methods)
    method = methods{mi};
    fig_dir = get_method_fig_dir(fig_dir_base, method);
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    fprintf('\n=== Method: %s | Levels: %s ===\n', method, strjoin(levels, ', '));

    for t = 1:length(tasks)
        task = tasks{t};
        for lv = 1:length(levels)
            level = levels{lv};
            if strcmp(level, 'trial')
                suffix = '';
            else
                suffix = '_subject';
            end
            csv_file = resolve_r2_csv(csv_dir, task, suffix, method);
            if isempty(csv_file)
                fprintf('Skipping %s %s-level (method=%s): no matching CSV found.\n', task, level, method);
                continue
            end

            T = readtable(csv_file);
            fprintf('Loaded %s %s-level (%s): %d FOOOF fits.\n', task, level, method, height(T));

            required_cols = {'r_squared', 'aperiodic_offset', 'aperiodic_exponent', ...
                             'electrodes', 'latency_ms', 'subjectID'};
            if ~all(ismember(required_cols, T.Properties.VariableNames))
                fprintf('  Missing required columns in %s; skipping.\n', csv_file);
                continue
            end

            [T, gaze_status] = add_gaze_from_multiverse(T, multiverse_csv_dir, task, level, method);
            fprintf('  %s\n', gaze_status);

            %% --- Exclude extreme aperiodic parameter values (IQR-based) ---
            iqr_factor = 3;
            n_before = height(T);

            valid_off = isfinite(T.aperiodic_offset);
            valid_exp = isfinite(T.aperiodic_exponent);
            if ~any(valid_off) || ~any(valid_exp)
                fprintf('  No finite aperiodic values in %s; skipping.\n', csv_file);
                continue
            end

            Q1_off = prctile(T.aperiodic_offset(valid_off), 25);
            Q3_off = prctile(T.aperiodic_offset(valid_off), 75);
            IQR_off = Q3_off - Q1_off;
            extreme_off = valid_off & ...
                (T.aperiodic_offset < Q1_off - iqr_factor * IQR_off | ...
                 T.aperiodic_offset > Q3_off + iqr_factor * IQR_off);

            Q1_exp = prctile(T.aperiodic_exponent(valid_exp), 25);
            Q3_exp = prctile(T.aperiodic_exponent(valid_exp), 75);
            IQR_exp = Q3_exp - Q1_exp;
            extreme_exp = valid_exp & ...
                (T.aperiodic_exponent < Q1_exp - iqr_factor * IQR_exp | ...
                 T.aperiodic_exponent > Q3_exp + iqr_factor * IQR_exp);

            extreme_mask = extreme_off | extreme_exp;
            n_extreme = sum(extreme_mask);
            fprintf('  Extreme aperiodic values (%.0f× IQR): %d / %d (%.1f%%)\n', ...
                iqr_factor, n_extreme, n_before, 100 * n_extreme / max(n_before, 1));
            fprintf('    Offset  bounds: [%.2f, %.2f]  (excluded %d)\n', ...
                Q1_off - iqr_factor * IQR_off, Q3_off + iqr_factor * IQR_off, sum(extreme_off));
            fprintf('    Exponent bounds: [%.2f, %.2f]  (excluded %d)\n', ...
                Q1_exp - iqr_factor * IQR_exp, Q3_exp + iqr_factor * IQR_exp, sum(extreme_exp));
            T(extreme_mask, :) = [];

            n_below = sum(T.r_squared < r2_thresh & isfinite(T.r_squared));
            n_below_filter = sum(T.r_squared < 0.60 & isfinite(T.r_squared));
            n_valid = sum(isfinite(T.r_squared));
            fprintf('  After exclusion: %d fits remaining.\n', height(T));
            fprintf('  R² < %.2f: %d / %d (%.1f%%)\n', r2_thresh, n_below, n_valid, ...
                100 * n_below / max(n_valid, 1));

            target_electrode = 'occipital';
            T = T(strcmp(T.electrodes, target_electrode), :);
            if isempty(T)
                fprintf('  No %s rows remaining after filtering; skipping plots.\n', target_electrode);
                continue
            end

            lat_labels = unique(T.latency_ms, 'stable');
            elec_labels = unique(T.electrodes, 'stable');
            subj_ids = unique(T.subjectID);
            n_subj = length(subj_ids);
            subj_map = containers.Map(subj_ids, 1:n_subj);

            colors_elec = [0.8 0.3 0.2];  % occipital = red

            %% --- Figure 1: R² per subject (faceted by latency; occipital only) ---
            fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
            task_title = task; task_title(1) = upper(task_title(1));
            sgtitle(sprintf('FOOOF R^2 per subject — %s (%s-level)', task_title, level), ...
                'FontSize', 16, 'FontWeight', 'bold');
            for il = 1:length(lat_labels)
                subplot(2, 2, il); hold on;
                mask = strcmp(T.latency_ms, lat_labels{il});
                for ie = 1:length(elec_labels)
                    sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                    if ~any(sub_mask), continue, end
                    subj_x = arrayfun(@(sid) subj_map(sid), T.subjectID(sub_mask));
                    jitter = (ie - 1.5) * 0.2 + 0.05 * randn(sum(sub_mask), 1);
                    scatter(subj_x + jitter, T.r_squared(sub_mask), 8, ...
                        colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.4, ...
                        'DisplayName', elec_labels{ie});
                end
                yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                    'Label', sprintf('R^2 = %.2f', r2_thresh), 'HandleVisibility', 'off');
                xticks(1:n_subj);
                xticklabels(arrayfun(@num2str, subj_ids, 'UniformOutput', false));
                xlabel('Subject'); ylabel('R^2');
                title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
                if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
                ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
                hold off;
            end
            fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_per_subject_%s%s.png', task, suffix));
            exportgraphics(fig1, fname, 'Resolution', 300);
            fprintf('  Saved: %s\n', fname);
            close(fig1);

            %% --- Figure 2: R² vs aperiodic exponent (faceted by latency) ---
            fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
            sgtitle(sprintf('FOOOF R^2 vs Aperiodic Exponent — %s (%s-level)', task_title, level), ...
                'FontSize', 16, 'FontWeight', 'bold');
            for il = 1:length(lat_labels)
                subplot(2, 2, il); hold on;
                mask = strcmp(T.latency_ms, lat_labels{il});
                for ie = 1:length(elec_labels)
                    sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                    if ~any(sub_mask), continue, end
                    scatter(T.aperiodic_exponent(sub_mask), T.r_squared(sub_mask), 10, ...
                        colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.35, ...
                        'DisplayName', elec_labels{ie});
                end
                yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
                xlabel('Aperiodic Exponent'); ylabel('R^2');
                title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
                if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
                ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
                corr_lines = build_corr_lines(T.aperiodic_exponent(mask), T.r_squared(mask), [0.90]);
                model_lines = build_exponent_model_lines(T(mask, :));
                xl = xlim;
                yl = ylim;
                text(xl(2) - 0.02 * range(xl), yl(1) + 0.02 * range(yl), ...
                    strjoin([corr_lines, model_lines], newline), ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 3, ...
                    'Interpreter', 'none');
                hold off;
            end
            fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_vs_exponent_%s%s.png', task, suffix));
            exportgraphics(fig2, fname, 'Resolution', 300);
            fprintf('  Saved: %s\n', fname);
            close(fig2);

            %% --- Figure 3: R² vs aperiodic offset (faceted by latency) ---
            fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
            sgtitle(sprintf('FOOOF R^2 vs Aperiodic Offset — %s (%s-level)', task_title, level), ...
                'FontSize', 16, 'FontWeight', 'bold');
            for il = 1:length(lat_labels)
                subplot(2, 2, il); hold on;
                mask = strcmp(T.latency_ms, lat_labels{il});
                for ie = 1:length(elec_labels)
                    sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                    if ~any(sub_mask), continue, end
                    scatter(T.aperiodic_offset(sub_mask), T.r_squared(sub_mask), 10, ...
                        colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.35, ...
                        'DisplayName', elec_labels{ie});
                end
                yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
                xlabel('Aperiodic Offset'); ylabel('R^2');
                title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
                if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
                ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
                corr_lines = build_corr_lines(T.aperiodic_offset(mask), T.r_squared(mask), [0.90]);
                xl = xlim;
                yl = ylim;
                text(xl(2) - 0.02 * range(xl), yl(1) + 0.02 * range(yl), ...
                    strjoin(corr_lines, newline), ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 3, ...
                    'Interpreter', 'none');
                hold off;
            end
            fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_vs_offset_%s%s.png', task, suffix));
            exportgraphics(fig3, fname, 'Resolution', 300);
            fprintf('  Saved: %s\n', fname);
            close(fig3);

            %% --- Figure 4: R² distribution histogram ---
            fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
            valid_r2 = T.r_squared(isfinite(T.r_squared));
            histogram(valid_r2, 50, 'FaceColor', [0.3 0.5 0.7], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.8);
            threshold_vals = [0.60, 0.70, 0.80, 0.85, 0.90];
            threshold_colors = [ ...
                0.80 0.50 0.00; ...
                0.55 0.55 0.55; ...
                0.35 0.35 0.35; ...
                0.75 0.25 0.25; ...
                0.90 0.00 0.00];
            for th = 1:numel(threshold_vals)
                xline(threshold_vals(th), '--', ...
                    sprintf('R^2 = %.2f', threshold_vals(th)), ...
                    'Color', threshold_colors(th, :), 'LineWidth', 1.4, ...
                    'LabelOrientation', 'horizontal', ...
                    'LabelVerticalAlignment', 'middle');
            end
            xlabel('FOOOF R^2', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold');
            pct_below = 100 * arrayfun(@(thr) sum(valid_r2 < thr), threshold_vals) / max(numel(valid_r2), 1);
            threshold_summary = sprintf('<0.60: %.1f%%; <0.70: %.1f%%; <0.80: %.1f%%; <0.85: %.1f%%; <0.90: %.1f%%', ...
                pct_below(1), pct_below(2), pct_below(3), pct_below(4), pct_below(5));
            title(sprintf('FOOOF R^2 Distribution — %s (%s-level)  [%s]', ...
                task_title, level, threshold_summary), ...
                'FontSize', 14, 'FontWeight', 'bold');
            set(gca, 'FontSize', 12);
            fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_distribution_%s%s.png', task, suffix));
            exportgraphics(fig4, fname, 'Resolution', 300);
            fprintf('  Saved: %s\n', fname);
            close(fig4);
        end
    end
end

fprintf('\n=== FOOOF R² QC DONE ===\n');

function levels = parse_levels(level_env)
if isempty(level_env)
    levels = {'trial', 'subject'};
    return
end

tokens = strsplit(lower(strtrim(level_env)), ',');
tokens = strtrim(tokens);
tokens = tokens(~cellfun(@isempty, tokens));
if isempty(tokens)
    levels = {'trial', 'subject'};
    return
end

if any(ismember(tokens, {'both', 'all'}))
    levels = {'trial', 'subject'};
    return
end

valid = {'trial', 'subject'};
levels = intersect(valid, tokens, 'stable');
if isempty(levels)
    levels = {'trial', 'subject'};
end
end

function methods = parse_methods(method_env)
if isempty(method_env)
    methods = {};
    return
end

tokens = strsplit(strtrim(method_env), ',');
tokens = strtrim(tokens);
tokens = tokens(~cellfun(@isempty, tokens));
if isempty(tokens)
    methods = {};
    return
end

if any(ismember(lower(tokens), {'both', 'all'}))
    methods = {'singleFFT', 'welch500_50'};
    return
end

methods = tokens;
end

function lines = build_corr_lines(x, y, y_thresholds)
lines = cell(1, numel(y_thresholds) + 1);
lines{1} = corr_line(x, y, 'all');
for i = 1:numel(y_thresholds)
    thr = y_thresholds(i);
    sub = y > thr;
    lines{i + 1} = corr_line(x(sub), y(sub), sprintf('R^2>%.2f', thr));
end
end

function line = corr_line(x, y, label)
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
n = numel(y);
if n < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    line = sprintf('%s: r=NA, p=NA, n=%d', label, n);
    return
end
[R, P] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
line = sprintf('%s: r=%.3f, p=%.3g, n=%d', label, R, P, n);
end

function lines = build_exponent_model_lines(Tsub)
lines = {};
gaze_col = detect_gaze_column(Tsub.Properties.VariableNames);
if isempty(gaze_col)
    lines = {'exp~gaze: unavailable (no gaze column in CSV)'};
    return
end

gaze_raw = Tsub.(gaze_col);
exp_raw = Tsub.aperiodic_exponent;
r2_raw = Tsub.r_squared;
[gaze_num, gaze_ok] = coerce_gaze_numeric(gaze_raw);

valid = isfinite(exp_raw) & isfinite(r2_raw) & gaze_ok;
if sum(valid) < 4
    lines = {sprintf('exp~gaze: unavailable (n=%d)', sum(valid))};
    return
end

tbl = table(exp_raw(valid), r2_raw(valid), gaze_num(valid), ...
    'VariableNames', {'exponent', 'r2', 'gaze'});

lm1 = fitlm(tbl, 'exponent ~ gaze');
lm2 = fitlm(tbl, 'exponent ~ r2');
lm3 = fitlm(tbl, 'exponent ~ gaze * r2');

lines = {
    model_line_from_lm(lm1, 'gaze', 'exp~gaze')
    model_line_from_lm(lm2, 'r2', 'exp~R^2')
    model_line_from_lm(lm3, 'gaze:r2', 'exp~gaze*R^2')
    };
end

function gaze_col = detect_gaze_column(var_names)
candidates = {'gaze_value', 'gaze', 'gaze_measure', 'baseline_gaze'};
gaze_col = '';
for i = 1:numel(candidates)
    idx = find(strcmpi(var_names, candidates{i}), 1);
    if ~isempty(idx)
        gaze_col = var_names{idx};
        return
    end
end
end

function [gaze_num, ok_mask] = coerce_gaze_numeric(gaze_raw)
if isnumeric(gaze_raw)
    gaze_num = double(gaze_raw);
    ok_mask = isfinite(gaze_num);
    return
end

if islogical(gaze_raw)
    gaze_num = double(gaze_raw);
    ok_mask = true(size(gaze_num));
    return
end

if iscellstr(gaze_raw) || isstring(gaze_raw) || iscategorical(gaze_raw)
    gaze_cat = categorical(gaze_raw);
    gaze_num = double(gaze_cat);
    ok_mask = ~isundefined(gaze_cat);
    return
end

gaze_num = nan(size(gaze_raw));
ok_mask = false(size(gaze_raw));
end

function line = model_line_from_lm(lm, term_name, label)
coef_tbl = lm.Coefficients;
match_idx = find(strcmp(coef_tbl.Properties.RowNames, term_name), 1);
if isempty(match_idx)
    line = sprintf('%s: term missing (n=%d, R^2=%.3f)', label, lm.NumObservations, lm.Rsquared.Ordinary);
    return
end
b = coef_tbl.Estimate(match_idx);
p = coef_tbl.pValue(match_idx);
line = sprintf('%s: b=%.3f, p=%.3g, R^2=%.3f, n=%d', ...
    label, b, p, lm.Rsquared.Ordinary, lm.NumObservations);
end

function [T_out, status_msg] = add_gaze_from_multiverse(T_in, multiverse_dir, task, level, method)
T_out = T_in;
if ismember('gaze_value', T_out.Properties.VariableNames)
    status_msg = 'Gaze merge: using existing gaze_value column.';
    return
end

gaze_file = resolve_multiverse_gaze_csv(multiverse_dir, task, level, method);
if isempty(gaze_file)
    status_msg = sprintf('Gaze merge: no source file found for %s/%s/%s.', task, level, method);
    return
end

if strcmp(level, 'trial')
    join_keys = {'subjectID', 'Condition', 'Trial', 'latency_ms', 'electrodes'};
else
    join_keys = {'subjectID', 'Condition', 'latency_ms', 'electrodes'};
end
join_keys = join_keys(ismember(join_keys, T_out.Properties.VariableNames));
if isempty(join_keys)
    status_msg = 'Gaze merge: no valid join keys in T; skipped.';
    return
end

try
    opts = detectImportOptions(gaze_file);
    needed = [join_keys, {'gaze_value'}];
    avail = intersect(needed, opts.VariableNames, 'stable');
    if ~ismember('gaze_value', avail)
        status_msg = sprintf('Gaze merge: gaze_value missing in %s.', gaze_file);
        return
    end
    opts.SelectedVariableNames = avail;
    G = readtable(gaze_file, opts);
catch ME
    status_msg = sprintf('Gaze merge: failed to read %s (%s).', gaze_file, ME.message);
    return
end

for k = 1:numel(join_keys)
    key = join_keys{k};
    if ismember(key, G.Properties.VariableNames)
        G.(key) = cast_to_like(G.(key), T_out.(key));
    end
end

G = G(:, intersect(needed, G.Properties.VariableNames, 'stable'));
group_keys = join_keys(ismember(join_keys, G.Properties.VariableNames));
if isempty(group_keys)
    status_msg = sprintf('Gaze merge: no join keys available in %s.', gaze_file);
    return
end

G = groupsummary(G, group_keys, 'mean', 'gaze_value');
if ismember('mean_gaze_value', G.Properties.VariableNames)
    G = renamevars(G, 'mean_gaze_value', 'gaze_value');
end
if ismember('GroupCount', G.Properties.VariableNames)
    G.GroupCount = [];
end

T_out = outerjoin(T_out, G, 'Keys', group_keys, 'MergeKeys', true, 'Type', 'left');
if ismember('gaze_value', T_out.Properties.VariableNames)
    n_gaze = sum(isfinite_if_numeric(T_out.gaze_value));
else
    n_gaze = 0;
end
status_msg = sprintf('Gaze merge: %d rows with finite gaze (source: %s).', n_gaze, gaze_file);
end

function out = cast_to_like(in, ref)
if isnumeric(ref)
    out = double(in);
    return
end

if iscell(ref)
    if isstring(in)
        out = cellstr(in);
    elseif iscategorical(in)
        out = cellstr(string(in));
    elseif ischar(in)
        out = cellstr(in);
    else
        out = in;
    end
    return
end

if isstring(ref)
    out = string(in);
    return
end

if iscategorical(ref)
    out = categorical(string(in));
    return
end

out = in;
end

function mask = isfinite_if_numeric(x)
if isnumeric(x)
    mask = isfinite(x);
else
    mask = false(size(x));
end
end

function gaze_file = resolve_multiverse_gaze_csv(multiverse_dir, task, level, method)
if strcmp(level, 'subject')
    candidates = {
        fullfile(multiverse_dir, sprintf('multiverse_%s_subject_%s.csv', task, method)), ...
        fullfile(multiverse_dir, sprintf('multiverse_%s_subject.csv', task))
    };
else
    if strcmpi(method, 'singleFFT')
        candidates = {
            fullfile(multiverse_dir, sprintf('multiverse_%s_singleFFT.csv', task)), ...
            fullfile(multiverse_dir, sprintf('multiverse_%s.csv', task))
        };
    else
        candidates = {
            fullfile(multiverse_dir, sprintf('multiverse_%s_%s.csv', task, method)), ...
            fullfile(multiverse_dir, sprintf('multiverse_%s.csv', task))
        };
    end
end

gaze_file = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        gaze_file = candidates{i};
        return
    end
end
end

function fig_dir = get_method_fig_dir(fig_dir_base, method)
if strcmpi(method, 'singleFFT')
    fig_dir = fig_dir_base;
elseif contains(lower(method), 'welch')
    fig_dir = fullfile(fig_dir_base, 'welch');
else
    fig_dir = fullfile(fig_dir_base, method);
end
end

function csv_file = resolve_r2_csv(csv_dir, task, suffix, method)
if strcmpi(method, 'singleFFT')
    candidates = {
        fullfile(csv_dir, ['fooof_r2_' task suffix '_singleFFT.csv']), ...
        fullfile(csv_dir, ['fooof_r2_' task suffix '.csv'])
    };
else
    candidates = {
        fullfile(csv_dir, ['fooof_r2_' task suffix '_' method '.csv']), ...
        fullfile(csv_dir, ['fooof_r2_' task suffix '.csv'])
    };
end

csv_file = '';
for i = 1:length(candidates)
    if isfile(candidates{i})
        csv_file = candidates{i};
        return
    end
end
end
