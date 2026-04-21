%% AOC Follow-up Input Table for Alpha Split (Sternberg)
% This script prepares a follow-up CSV for downstream statistics in Python
% (AOC_stats_splits_AmpRed.py). No inferential models are run here.
%
% The table contains subject-by-load rows with:
% - split group (Reduction/Amplification),
% - TOI gaze deviation (%) using baseline normalization,
% - specParam alpha, RT, accuracy, and trial-count proxy columns.
%
% TOI is derived from the permutation cluster while accounting for 100 ms
% time-bin uncertainty.

%% Setup
startup
[subjects, paths] = setup('AOC');
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

fprintf('\n=== AOC Split Follow-up Input Table (Sternberg) ===\n');
fprintf('Feature directory: %s\n', feat_dir);
fprintf('Output directory:  %s\n', out_dir);

%% Parameters
task_tag = 'sternberg';
merged_file = fullfile(feat_dir, 'AOC_merged_data_sternberg.mat');
merged_var = 'merged_data_sternberg';
gaze_fname = 'gaze_series_sternberg_trials.mat';
cond_vals = [2 4 6];
cond_codes = [22 24 26];
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

% Split settings (match split script defaults)
alpha_zero_pct = 5;
alpha_ref_percentile = 95;

% Cluster-derived window and binning uncertainty handling
cluster_window_sec = [0.452 1.352];
bin_width_sec = 0.100;
toi_mode = 'expanded_to_bin_edges'; % {'expanded_to_bin_edges','half_bin_margin'}

switch lower(toi_mode)
    case 'expanded_to_bin_edges'
        toi_sec = [ ...
            floor(cluster_window_sec(1) / bin_width_sec) * bin_width_sec, ...
            ceil(cluster_window_sec(2) / bin_width_sec) * bin_width_sec];
    case 'half_bin_margin'
        toi_sec = [ ...
            cluster_window_sec(1) - bin_width_sec/2, ...
            cluster_window_sec(2) + bin_width_sec/2];
    otherwise
        error('Unknown toi_mode: %s', toi_mode);
end

baseline_sec = [-0.5 -0.25];
fprintf('Cluster window: [%.3f %.3f] s\n', cluster_window_sec(1), cluster_window_sec(2));
fprintf('TOI mode: %s -> TOI [%.3f %.3f] s\n', toi_mode, toi_sec(1), toi_sec(2));

%% Load merged subject-level condition table
if ~isfile(merged_file)
    error('Missing merged file: %s', merged_file);
end

S = load(merged_file, merged_var);
if ~isfield(S, merged_var)
    error('Variable %s not found in %s', merged_var, merged_file);
end
T = struct2table(S.(merged_var));
if ~ismember('ID', T.Properties.VariableNames) || ~ismember('Condition', T.Properties.VariableNames)
    error('Merged table is missing required columns ID/Condition.');
end

%% Build alpha split groups
uIDs = unique(T.ID);
nSubj = numel(uIDs);
alpha_mean = nan(nSubj, 1);
for i = 1:nSubj
    sid = uIDs(i);
    mask = T.ID == sid;
    alpha_mean(i) = mean(T.AlphaPower_FOOOF_bl(mask), 'omitnan');
end

split_valid = isfinite(alpha_mean);
vals = alpha_mean(split_valid);
if numel(vals) >= 4
    q1_split = prctile(vals, 25);
    q3_split = prctile(vals, 75);
    iqr_split = q3_split - q1_split;
    if isfinite(iqr_split) && iqr_split > 0
        lo_split = q1_split - 3 * iqr_split;
        hi_split = q3_split + 3 * iqr_split;
        split_valid = split_valid & (alpha_mean >= lo_split) & (alpha_mean <= hi_split);
    end
end

valid_alpha_mean = alpha_mean(split_valid);
if isempty(valid_alpha_mean)
    error('No finite subject-level alpha values found for split.');
end
alpha_abs_ref = prctile(abs(valid_alpha_mean), alpha_ref_percentile);
alpha_zero_margin = (alpha_zero_pct / 100) * alpha_abs_ref;

reduction_ids = uIDs(split_valid & (alpha_mean < -alpha_zero_margin));
amplification_ids = uIDs(split_valid & (alpha_mean > alpha_zero_margin));
zero_ids = uIDs(split_valid & (abs(alpha_mean) <= alpha_zero_margin));

fprintf('\n=== Split Summary ===\n');
fprintf('Subjects total: %d\n', nSubj);
fprintf('Reduction: %d, Amplification: %d, Excluded zero-band: %d\n', ...
    numel(reduction_ids), numel(amplification_ids), numel(zero_ids));

%% Assemble modeling table (subject x load)
rows = [];

for s = 1:nSubj
    sid = uIDs(s);
    sid_str = num2str(sid);

    if ismember(sid, reduction_ids)
        grp_num = 1;
    elseif ismember(sid, amplification_ids)
        grp_num = 2;
    else
        continue
    end

    subj_folder = resolve_subject_folder(subjects, sid);
    if isempty(subj_folder)
        subj_folder = sid_str;
    end

    gaze_file = fullfile(feat_dir, subj_folder, 'gaze', gaze_fname);
    gaze_by_cond = nan(1, 3);
    if isfile(gaze_file)
        G = load(gaze_file);
        if isfield(G, 'trialinfo') && isfield(G, 'gaze_x') && isfield(G, 'gaze_y')
            conds = parse_trialinfo_conds(G.trialinfo);
            for c = 1:3
                tr_mask = conds == cond_codes(c);
                tr_idx = find(tr_mask);
                if isempty(tr_idx)
                    continue
                end
                tr_vals = nan(numel(tr_idx), 1);
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
                    bl_idx = tt >= baseline_sec(1) & tt <= baseline_sec(2);
                    toi_idx = tt >= toi_sec(1) & tt <= toi_sec(2);
                    if ~any(bl_idx) || ~any(toi_idx)
                        continue
                    end
                    bl = mean(dev(bl_idx), 'omitnan');
                    if ~isfinite(bl) || bl <= 0
                        continue
                    end
                    toi_mean = mean(dev(toi_idx), 'omitnan');
                    tr_vals(k) = ((toi_mean / bl) - 1) * 100;
                end
                gaze_by_cond(c) = mean(tr_vals, 'omitnan');
            end
        end
    end

    subj_rows = T(T.ID == sid, :);
    for c = 1:3
        cmask = subj_rows.Condition == cond_vals(c);
        if ~any(cmask)
            continue
        end

        alpha_val = mean(subj_rows.AlphaPower_FOOOF_bl(cmask), 'omitnan');
        rt_val = pull_first_available_mean(subj_rows(cmask, :), {'ReactionTime', 'RT'});
        acc_val = pull_first_available_mean(subj_rows(cmask, :), {'Accuracy', 'Acc'});
        ntr_val = pull_first_available_mean(subj_rows(cmask, :), {'NTrials', 'nTrials', 'NumTrials', 'TrialsN'});

        if ~isfinite(alpha_val) || ~isfinite(gaze_by_cond(c))
            continue
        end

        rows(end+1, :) = [sid, grp_num, cond_vals(c), gaze_by_cond(c), alpha_val, rt_val, acc_val, ntr_val]; %#ok<AGROW>
    end
end

if isempty(rows)
    error('No rows available for follow-up models.');
end

tbl = array2table(rows, 'VariableNames', ...
    {'Subject', 'GroupNum', 'LoadNum', 'GazeTOI_pct', 'AlphaPower_FOOOF_bl', 'ReactionTime', 'Accuracy', 'NTrials'});
tbl.Subject = categorical(tbl.Subject);
tbl.Group = categorical(tbl.GroupNum, [1 2], {'Reduction', 'Amplification'});
tbl.Load = categorical(tbl.LoadNum, cond_vals, cond_labels);

% z-standardize gaze for downstream model stability
tbl.GazeTOI_z = zscore(tbl.GazeTOI_pct);

fprintf('\n=== Modeling table summary ===\n');
fprintf('Rows: %d | Subjects: %d\n', height(tbl), numel(unique(tbl.Subject)));

%% Save final modeling table only (stats are run in Python)
tbl.Included = true(height(tbl), 1);
tbl.LoadValue = tbl.LoadNum;
tbl.LoadLabel = string(tbl.Load);
out_csv = fullfile(out_dir, 'AOC_splitAmpRed_followup_model_input.csv');
writetable(tbl, out_csv);
fprintf('\nSaved follow-up input table:\n%s\n', out_csv);

%% ========================= Local Functions =========================
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

function val = pull_first_available_mean(tbl_in, candidates)
val = NaN;
for i = 1:numel(candidates)
    nm = candidates{i};
    if ismember(nm, tbl_in.Properties.VariableNames)
        x = tbl_in.(nm);
        if iscell(x)
            continue
        end
        val = mean(double(x), 'omitnan');
        return
    end
end
end

