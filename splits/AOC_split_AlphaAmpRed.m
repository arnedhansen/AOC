%% AOC Split Alpha Amp/Red (Sternberg)
% Creates a CSV with split labels and per-subject/per-load metrics
% for downstream statistics in Python.

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
pathAOC = paths.features;

if ispc
    base_data = 'W:\Students\Arne\AOC';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
end

feat_dir = fullfile(base_data, 'data', 'features');
out_dir = fullfile(base_data, 'data', 'stats', 'splits');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

fprintf('\n=== AOC Split Alpha Amp/Red (Sternberg) ===\n');

%% Task definition (Sternberg only)
task_tag = 'sternberg';
merged_file = 'AOC_merged_data_sternberg.mat';
merged_var = 'merged_data_sternberg';
cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};
gaze_fname = 'gaze_series_sternberg_trials.mat';

alpha_zero_pct = 5;
alpha_ref_percentile = 95;

%% Load merged data
src_file = fullfile(feat_dir, merged_file);
if ~isfile(src_file)
    error('Missing merged file: %s', src_file);
end

S = load(src_file, merged_var);
if ~isfield(S, merged_var)
    error('Variable %s not found in %s', merged_var, src_file);
end
T = struct2table(S.(merged_var));

%% Define split at subject level from AlphaPower_FOOOF_bl
uIDs = unique(T.ID);
nSubj = numel(uIDs);
alpha_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    alpha_mean(i) = mean(T.AlphaPower_FOOOF_bl(T.ID == sid), 'omitnan');
end

split_valid = isfinite(alpha_mean);
vals = alpha_mean(split_valid);
if numel(vals) >= 4
    q1 = prctile(vals, 25);
    q3 = prctile(vals, 75);
    iqr_val = q3 - q1;
    if isfinite(iqr_val) && iqr_val > 0
        lo = q1 - 3 * iqr_val;
        hi = q3 + 3 * iqr_val;
        split_valid = split_valid & alpha_mean >= lo & alpha_mean <= hi;
    end
end

valid_alpha = alpha_mean(split_valid);
if isempty(valid_alpha)
    error('No valid subject-level alpha values for split.');
end
alpha_abs_ref = prctile(abs(valid_alpha), alpha_ref_percentile);
alpha_zero_margin = (alpha_zero_pct / 100) * alpha_abs_ref;

reduction_ids = uIDs(split_valid & (alpha_mean < -alpha_zero_margin));
amplification_ids = uIDs(split_valid & (alpha_mean > alpha_zero_margin));
zero_ids = uIDs(split_valid & (abs(alpha_mean) <= alpha_zero_margin));
invalid_ids = uIDs(~split_valid);

fprintf('\n=== Split Summary [%s] ===\n', task_tag);
fprintf('Subjects total: %d\n', nSubj);
fprintf('Cutoff (|alpha|): %.4f\n', alpha_zero_margin);
fprintf('Reduction: %d\n', numel(reduction_ids));
fprintf('Amplification: %d\n', numel(amplification_ids));
fprintf('Excluded near-zero: %d\n', numel(zero_ids));
fprintf('Excluded invalid: %d\n', numel(invalid_ids));

%% Preallocate metrics
metrics_Alpha = nan(nSubj, 3);
metrics_Dev = nan(nSubj, 3);

% Gaze-derived percent baseline deviation time-course
fs = 500;
t_full = -0.5:1/fs:3;
t_plot = t_full(2:end);
Tf = numel(t_plot);
dev_tc = nan(nSubj, 3, Tf);

missing_gaze = {};

%% Aggregate subject-level/load-level metrics
fprintf('\n=== Aggregating merged and gaze metrics ===\n');
for s = 1:nSubj
    sid = uIDs(s);
    sid_str = num2str(sid);
    subj_rows = T(T.ID == sid, :);

    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if any(cmask)
            metrics_Alpha(s, c) = mean(subj_rows.AlphaPower_FOOOF_bl(cmask), 'omitnan');
            metrics_Dev(s, c) = mean(subj_rows.GazeDeviationFullBL(cmask), 'omitnan');
        end
    end

    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end
    gaze_file = fullfile(pathAOC, subj_folder, 'gaze', gaze_fname);
    if ~isfile(gaze_file)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    G = load(gaze_file);
    if ~isfield(G, 'trialinfo') || ~isfield(G, 'gaze_x') || ~isfield(G, 'gaze_y')
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    conds = parse_trialinfo_conds(G.trialinfo);
    if isempty(conds)
        missing_gaze{end+1} = sid_str; %#ok<AGROW>
        continue
    end

    for c = 1:3
        tr_idx = find(conds == cond_codes(c));
        if isempty(tr_idx)
            continue
        end
        dev_mat = nan(numel(tr_idx), Tf);
        for k = 1:numel(tr_idx)
            tr = tr_idx(k);
            x = double(G.gaze_x{tr});
            y = double(G.gaze_y{tr});
            if isempty(x) || isempty(y) || numel(x) ~= numel(y)
                continue
            end
            if isfield(G, 'time') && numel(G.time) >= tr && ~isempty(G.time{tr})
                tt = double(G.time{tr});
                if numel(tt) ~= numel(x)
                    tt = linspace(tt(1), tt(end), numel(x));
                end
            else
                tt = linspace(-0.5, 3, numel(x));
            end
            dev = sqrt((x - 400).^2 + (y - 300).^2);
            try
                dev_mat(k, :) = interp1(tt, dev, t_plot, 'linear', NaN);
            catch
            end
        end
        dev_tc(s, c, :) = nanmean(dev_mat, 1);
    end
end

fprintf('Missing gaze files/fields: %d\n', numel(unique(missing_gaze)));

%% Baseline-corrected gaze deviation (%), then [0..2] s mean by load
t_vec = linspace(-0.5, 3, Tf);
bl_idx = t_vec >= -0.5 & t_vec <= -0.25;
task_idx = t_vec >= 0 & t_vec <= 2;

dev_bl = mean(dev_tc(:, :, bl_idx), 3, 'omitnan');
dev_bl_3d = repmat(dev_bl, [1, 1, Tf]);
dev_bl_3d(dev_bl_3d <= 0 | ~isfinite(dev_bl_3d)) = NaN;
dev_tc_pct = (dev_tc ./ dev_bl_3d - 1) * 100;
dev_tc_pct(~isfinite(dev_tc_pct)) = NaN;
dev_pct_by_load = squeeze(mean(dev_tc_pct(:, :, task_idx), 3, 'omitnan'));

%% Export long-format CSV for Python stats
split_label = repmat("Excluded", nSubj, 1);
split_label(ismember(uIDs, reduction_ids)) = "Reduction";
split_label(ismember(uIDs, amplification_ids)) = "Amplification";
split_label(ismember(uIDs, zero_ids)) = "ExcludedNearZero";
split_label(ismember(uIDs, invalid_ids)) = "ExcludedInvalid";
is_included = split_label == "Reduction" | split_label == "Amplification";

rows_n = nSubj * 3;
ID = nan(rows_n, 1);
LoadValue = nan(rows_n, 1);
LoadLabel = strings(rows_n, 1);
Group = strings(rows_n, 1);
Included = false(rows_n, 1);
AlphaPower_FOOOF_bl = nan(rows_n, 1);
GazeDeviationFullBL = nan(rows_n, 1);
GazeDev_pct_0_2s = nan(rows_n, 1);
AlphaSplitSubjectMean = nan(rows_n, 1);
AlphaZeroCutoff = nan(rows_n, 1);

r = 1;
for s = 1:nSubj
    for c = 1:3
        ID(r) = uIDs(s);
        LoadValue(r) = cond_vals(c);
        LoadLabel(r) = string(cond_labels{c});
        Group(r) = split_label(s);
        Included(r) = is_included(s);
        AlphaPower_FOOOF_bl(r) = metrics_Alpha(s, c);
        GazeDeviationFullBL(r) = metrics_Dev(s, c);
        GazeDev_pct_0_2s(r) = dev_pct_by_load(s, c);
        AlphaSplitSubjectMean(r) = alpha_mean(s);
        AlphaZeroCutoff(r) = alpha_zero_margin;
        r = r + 1;
    end
end

split_stats = table(ID, LoadValue, LoadLabel, Group, Included, ...
    AlphaPower_FOOOF_bl, GazeDeviationFullBL, GazeDev_pct_0_2s, ...
    AlphaSplitSubjectMean, AlphaZeroCutoff);

out_csv = fullfile(out_dir, 'AOC_splitAmpRed_sternberg_stats_input.csv');
writetable(split_stats, out_csv);

fprintf('\nSaved CSV: %s\n', out_csv);
fprintf('Rows: %d\n', height(split_stats));
fprintf('Included rows (Reduction/Amplification): %d\n', sum(split_stats.Included));

%% Local functions
function conds = parse_trialinfo_conds(trialinfo)
conds = [];
if isempty(trialinfo)
    return
end
if isvector(trialinfo)
    conds = trialinfo(:);
elseif size(trialinfo, 2) >= 1
    if size(trialinfo, 2) == 2
        conds = trialinfo(:, 1);
    elseif size(trialinfo, 1) == 2
        conds = trialinfo(1, :)';
    else
        conds = trialinfo(:, 1);
    end
end
end

function subj_folder = resolve_subject_folder(subjects, sid)
subj_folder = '';
for i = 1:numel(subjects)
    sval = str2double(subjects{i});
    if isfinite(sval) && sval == sid
        subj_folder = subjects{i};
        return
    end
end
end
