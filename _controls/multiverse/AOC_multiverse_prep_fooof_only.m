%% AOC Multiverse — Trial-Level FOOOF-Only Refresh (Sternberg & N-Back)
% Recomputes ONLY FOOOF-derived fields for trial-level multiverse CSVs:
%   - alpha (for fooof == 'FOOOFed' rows)
%   - aperiodic_offset / aperiodic_exponent (for fooof == 'FOOOFed' rows)
% Keeps all non-FOOOFed and gaze/baseline metadata from existing CSVs.
%
% Input base CSVs expected:
%   multiverse_sternberg.csv, multiverse_nback.csv
%
% Output per mode:
%   multiverse_sternberg_<mode>.csv, multiverse_nback_<mode>.csv
%   fooof_r2_sternberg_<mode>.csv, fooof_r2_nback_<mode>.csv

disp(upper('=== AOC MULTIVERSE FOOOF-ONLY PREP START ==='))

%% Setup
path_preproc = [];
subjects = {};
try
  if exist('startup', 'file')
    startup
    disp(upper('Startup run.'))
  end
  if exist('setup', 'file')
    [subjects, path_preproc, ~, ~] = setup('AOC');
    disp(upper(['Setup: ' num2str(length(subjects)) ' subjects, path_preproc = ' path_preproc]))
  end
catch ME
  disp(upper(['Setup failed: ' ME.message ' — using base_features for subject list.']))
end

%% Paths
if ispc
  base_data = 'W:\Students\Arne\AOC';
  base_features = fullfile(base_data, 'data', 'features');
else
  base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
  base_features = fullfile(base_data, 'data', 'features');
end
out_dir = fullfile(base_data, 'data', 'multiverse');
if ~isfolder(out_dir), mkdir(out_dir); end
r2_dir = fullfile(base_data, 'data', 'controls', 'multiverse');
if ~isfolder(r2_dir), mkdir(r2_dir); end
disp(upper(['Paths: base_data = ' base_data ', out_dir = ' out_dir ', r2_dir = ' r2_dir]))

if isempty(subjects)
  dirs = dir(base_features);
  dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
  subjects = {dirs.name};
  disp(upper(['Subject list from base_features: ' num2str(length(subjects)) ' folders.']))
end
if isempty(path_preproc)
  path_preproc = base_features;
  disp(upper('Using base_features as path_preproc.'))
end

%% Config
fooof_mode = 'welch500_50';   % 'singleFFT' | 'welch500_50' | 'BOTH'
fooof_mode_env = getenv('AOC_FOOOF_MODE');
if ~isempty(fooof_mode_env), fooof_mode = fooof_mode_env; end
welch_seg_len_sec = 0.5;
welch_overlap = 0.5;
valid_fooof_modes = {'singleFFT', 'welch500_50', 'BOTH'};
if ~ismember(fooof_mode, valid_fooof_modes)
  error('Invalid fooof_mode "%s". Valid options: %s', fooof_mode, strjoin(valid_fooof_modes, ', '));
end
if strcmpi(fooof_mode, 'BOTH')
  fooof_modes_to_run = {'singleFFT', 'welch500_50'};
else
  fooof_modes_to_run = {fooof_mode};
end
use_parfor = init_parallel_pool_science_cloud();
alphaRange = [8 14];
latency_opts = {'0_500ms', '0_1000ms', '0_2000ms', '1000_2000ms'};
lat_windows = {[0 0.5], [0 1], [0 2], [1 2]};
electrodes_opts = {'posterior', 'occipital'};

%% Resolve channel labels once
labels_master = [];
for s = 1:length(subjects)
  sid = subjects{s};
  eeg_dir = fullfile(path_preproc, sid, 'eeg');
  early_file = fullfile(eeg_dir, 'power_stern_early_trials.mat');
  if ~isfile(early_file)
    early_file = fullfile(base_features, sid, 'eeg', 'power_stern_early_trials.mat');
  end
  if isfile(early_file)
    load(early_file, 'powload2_early');
    labels_master = powload2_early.label;
    break
  end
end
if isempty(labels_master)
  error('Could not resolve channel labels from power_stern_early_trials.mat.')
end
[idx_post, idx_occ] = get_channel_indices(labels_master);
ch_label_sets = {labels_master(idx_post), labels_master(idx_occ)};

for im = 1:numel(fooof_modes_to_run)
  run_mode = fooof_modes_to_run{im};
  mode_tag = regexprep(run_mode, '[^A-Za-z0-9]+', '_');
  disp(upper(['=== RUN FOOOF MODE: ' run_mode ' ===']))
  for t = 1:2
    if t == 1
      task_name = 'sternberg';
      base_csv = fullfile(out_dir, 'multiverse_sternberg.csv');
      eeg_td_file = 'dataEEG_sternberg.mat';
      eeg_tfr_file = 'dataEEG_TFR_sternberg.mat';
      early_file = 'power_stern_early_trials.mat';
      full_file = 'power_stern_full_trials.mat';
      cond_vals = [2 4 6];
      cond_codes = [22 24 26];
    else
      task_name = 'nback';
      base_csv = fullfile(out_dir, 'multiverse_nback.csv');
      eeg_td_file = 'dataEEG_nback.mat';
      eeg_tfr_file = 'dataEEG_TFR_nback.mat';
      early_file = 'power_nback_early_trials.mat';
      full_file = 'power_nback_full_trials.mat';
      cond_vals = [1 2 3];
      cond_codes = [21 22 23];
    end
    if ~isfile(base_csv)
      error('Base multiverse CSV missing: %s', base_csv);
    end
    disp(upper(['--- TASK ' task_name ': loading base CSV ---']))
    tbl = readtable(base_csv);
    if ~ismember('fooof', tbl.Properties.VariableNames)
      error('Expected column "fooof" in %s', base_csv);
    end
    foo_mask = strcmp(tbl.fooof, 'FOOOFed');
    foo_tbl = tbl(foo_mask, :);
    if isempty(foo_tbl)
      warning('No FOOOFed rows found in %s', base_csv);
      continue
    end
    subj_ids = unique(foo_tbl.subjectID);
    updates_all = table();
    r2_all = table();
    cfg_fooof = []; cfg_fooof.method = 'mtmfft'; cfg_fooof.taper = 'hanning';
    cfg_fooof.foilim = [2 40]; cfg_fooof.pad = 5; cfg_fooof.output = 'fooof'; cfg_fooof.keeptrials = 'no';

    for si = 1:numel(subjects)
      sid = subjects{si};
      sid_num = str2double(sid);
      if isnan(sid_num) || ~ismember(sid_num, subj_ids), continue, end
      eeg_dir = fullfile(path_preproc, sid, 'eeg');
      if ~isfolder(eeg_dir), eeg_dir = fullfile(base_features, sid, 'eeg'); end
      if ~isfolder(eeg_dir), continue, end

      % Load time-domain EEG
      data_td = [];
      eeg_tfr_path = fullfile(eeg_dir, eeg_tfr_file);
      eeg_td_path = fullfile(eeg_dir, eeg_td_file);
      if isfile(eeg_tfr_path)
        load(eeg_tfr_path, 'dataTFR'); data_td = dataTFR;
      elseif isfile(eeg_td_path)
        load(eeg_td_path, 'dataEEG');
        if exist('dataEEG', 'var'), data_td = dataEEG; else, load(eeg_td_path, 'dataTFR'); data_td = dataTFR; end
      end
      if isempty(data_td), continue, end

      % Subject IAF from existing power files
      p_early = fullfile(eeg_dir, early_file);
      p_full = fullfile(eeg_dir, full_file);
      IAF_band = alphaRange;
      if isfile(p_early) && isfile(p_full)
        load(p_early); load(p_full);
        if strcmp(task_name, 'sternberg')
          pow_early = {powload2_early, powload4_early, powload6_early};
          pow_full = {powload2_full,  powload4_full,  powload6_full};
        else
          pow_early = {powload1_early, powload2_early, powload3_early};
          pow_full  = {powload1_full,  powload2_full,  powload3_full};
        end
        iaf_pow = pow_full{1}; if isempty(iaf_pow), iaf_pow = pow_early{1}; end
        if ~isempty(iaf_pow)
          IAF_band = get_IAF_band(iaf_pow, idx_occ, alphaRange);
        end
      end

      subj_rows = foo_tbl(foo_tbl.subjectID == sid_num, :);
      for ci = 1:numel(cond_vals)
        cval = cond_vals(ci);
        ccode = cond_codes(ci);
        cond_rows = subj_rows(subj_rows.Condition == cval, :);
        if isempty(cond_rows), continue, end

        td_inds = find(data_td.trialinfo(:,1) == ccode);
        if isempty(td_inds), continue, end
        td_trialnums = data_td.trialinfo(td_inds, 2);

        key_rows = unique(cond_rows(:, {'Trial','electrodes','latency_ms'}));
        n_key = height(key_rows);
        td_idx_vec = nan(n_key,1); ie_vec = nan(n_key,1); il_vec = nan(n_key,1);
        for k = 1:n_key
          lat_k = char(string(key_rows.latency_ms(k)));
          elec_k = char(string(key_rows.electrodes(k)));
          il = find(strcmp(latency_opts, lat_k), 1);
          ie = find(strcmp(electrodes_opts, elec_k), 1);
          m = find(td_trialnums == key_rows.Trial(k), 1);
          if isempty(il) || isempty(ie) || isempty(m), continue, end
          il_vec(k) = il; ie_vec(k) = ie; td_idx_vec(k) = td_inds(m);
        end

        ap_off = nan(n_key,1); ap_exp = nan(n_key,1); r2v = nan(n_key,1); erv = nan(n_key,1); nseg = nan(n_key,1);
        alpha_c = nan(n_key,1); alpha_i = nan(n_key,1);
        if use_parfor
          parfor k = 1:n_key
            if ~isfinite(td_idx_vec(k)), continue, end
            ie = ie_vec(k); il = il_vec(k);
            [f_osc, f_axis, ao, ae, rr, ee, ns] = run_fooof_from_raw(data_td, td_idx_vec(k), ch_label_sets{ie}, ...
              lat_windows{il}, cfg_fooof, run_mode, welch_seg_len_sec, welch_overlap);
            ap_off(k) = ao; ap_exp(k) = ae; r2v(k) = rr; erv(k) = ee; nseg(k) = ns;
            if isempty(f_osc), continue, end
            bi = (f_axis >= alphaRange(1) & f_axis <= alphaRange(2));
            if any(bi), alpha_c(k) = mean(f_osc(bi), 'omitnan'); end
            bi = (f_axis >= IAF_band(1) & f_axis <= IAF_band(2));
            if any(bi), alpha_i(k) = mean(f_osc(bi), 'omitnan'); end
          end
        else
          for k = 1:n_key
            if ~isfinite(td_idx_vec(k)), continue, end
            ie = ie_vec(k); il = il_vec(k);
            [f_osc, f_axis, ao, ae, rr, ee, ns] = run_fooof_from_raw(data_td, td_idx_vec(k), ch_label_sets{ie}, ...
              lat_windows{il}, cfg_fooof, run_mode, welch_seg_len_sec, welch_overlap);
            ap_off(k) = ao; ap_exp(k) = ae; r2v(k) = rr; erv(k) = ee; nseg(k) = ns;
            if isempty(f_osc), continue, end
            bi = (f_axis >= alphaRange(1) & f_axis <= alphaRange(2));
            if any(bi), alpha_c(k) = mean(f_osc(bi), 'omitnan'); end
            bi = (f_axis >= IAF_band(1) & f_axis <= IAF_band(2));
            if any(bi), alpha_i(k) = mean(f_osc(bi), 'omitnan'); end
          end
        end

        upd = table();
        upd.subjectID = sid_num .* ones(n_key,1);
        upd.Condition = cval .* ones(n_key,1);
        upd.Trial = key_rows.Trial;
        upd.electrodes = key_rows.electrodes;
        upd.latency_ms = key_rows.latency_ms;
        upd.alpha_canonical = alpha_c;
        upd.alpha_IAF = alpha_i;
        upd.aperiodic_offset_new = ap_off;
        upd.aperiodic_exponent_new = ap_exp;
        upd.fooof_mode_new = repmat({run_mode}, n_key, 1);
        updates_all = [updates_all; upd]; %#ok<AGROW>

        r2k = table();
        r2k.subjectID = sid_num .* ones(n_key,1);
        r2k.Condition = cval .* ones(n_key,1);
        r2k.Trial = key_rows.Trial;
        r2k.latency_ms = key_rows.latency_ms;
        r2k.electrodes = key_rows.electrodes;
        r2k.fooof_mode = repmat({run_mode}, n_key, 1);
        r2k.welch_n_segments = nseg;
        r2k.r_squared = r2v;
        r2k.fooof_error = erv;
        r2k.aperiodic_exponent = ap_exp;
        r2k.aperiodic_offset = ap_off;
        r2_all = [r2_all; r2k]; %#ok<AGROW>
      end
    end

    if isempty(updates_all)
      warning('No FOOOF updates generated for task %s (%s).', task_name, run_mode);
      continue
    end

    % Merge updates into existing FOOOFed rows
    foo_tbl.orig_idx = find(foo_mask);
    foo_tbl = join(foo_tbl, updates_all, ...
      'Keys', {'subjectID','Condition','Trial','electrodes','latency_ms'}, ...
      'Type', 'left', ...
      'RightVariables', {'alpha_canonical','alpha_IAF','aperiodic_offset_new','aperiodic_exponent_new','fooof_mode_new'});
    is_canon = strcmp(foo_tbl.alpha_type, 'canonical');
    new_alpha = nan(height(foo_tbl), 1);
    new_alpha(is_canon) = foo_tbl.alpha_canonical(is_canon);
    new_alpha(~is_canon) = foo_tbl.alpha_IAF(~is_canon);
    has_new_alpha = isfinite(new_alpha);
    foo_tbl.alpha(has_new_alpha) = new_alpha(has_new_alpha);
    has_new_off = isfinite(foo_tbl.aperiodic_offset_new);
    has_new_exp = isfinite(foo_tbl.aperiodic_exponent_new);
    foo_tbl.aperiodic_offset(has_new_off) = foo_tbl.aperiodic_offset_new(has_new_off);
    foo_tbl.aperiodic_exponent(has_new_exp) = foo_tbl.aperiodic_exponent_new(has_new_exp);

    drop_cols = {'alpha_canonical','alpha_IAF','aperiodic_offset_new','aperiodic_exponent_new','fooof_mode_new'};
    drop_cols = intersect(drop_cols, foo_tbl.Properties.VariableNames);
    foo_tbl(:, drop_cols) = [];
    foo_tbl = sortrows(foo_tbl, 'orig_idx');
    tbl_updated = tbl;
    tbl_updated(foo_tbl.orig_idx, :) = removevars(foo_tbl, {'orig_idx'});

    out_csv = fullfile(out_dir, ['multiverse_' task_name '_' mode_tag '.csv']);
    writetable(tbl_updated, out_csv);
    disp(upper(['Written: ' out_csv]))

    r2_path = fullfile(r2_dir, ['fooof_r2_' task_name '_' mode_tag '.csv']);
    writetable(r2_all, r2_path);
    disp(upper(['Written: ' r2_path]))
  end
end

disp(upper('=== AOC MULTIVERSE FOOOF-ONLY PREP DONE ==='))

%% ---------- local helpers ----------
function [idx_post, idx_occ] = get_channel_indices(labels)
  n = length(labels);
  idx_occ = []; idx_pari = [];
  for i = 1:n
    lb = labels{i};
    if contains(lb, 'O') || contains(lb, 'I'), idx_occ(end+1) = i; end
    if contains(lb, 'P') && ~contains(lb, 'F'), idx_pari(end+1) = i; end
  end
  idx_post = unique([idx_occ, idx_pari]);
  if isempty(idx_post), idx_post = idx_occ; end
end

function IAF_band = get_IAF_band(pow_full, chIdx, alphaRange)
  if isempty(chIdx), IAF_band = alphaRange; return, end
  spec = squeeze(mean(mean(pow_full.powspctrm(:, chIdx, :), 1), 2));
  f = pow_full.freq(:);
  aIdx = find(f >= alphaRange(1) & f <= alphaRange(2));
  alphaF = f(aIdx); alphaP = spec(aIdx);
  [pks, locs] = findpeaks(double(alphaP));
  if isempty(pks), IAF_band = alphaRange; return, end
  [~, im] = max(pks);
  IAF = alphaF(locs(im));
  if IAF <= alphaRange(1) || IAF >= alphaRange(2), IAF_band = alphaRange; return, end
  IAF_band = [max(f(1), IAF-4), min(f(end), IAF+2)];
end

function use_parfor = init_parallel_pool_science_cloud()
  use_parfor = false;
  try
    has_parallel = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
  catch
    has_parallel = false;
  end
  if ~has_parallel
    disp(upper('PARFOR DISABLED: PARALLEL TOOLBOX NOT AVAILABLE.'))
    return
  end
  try
    p = gcp('nocreate');
    if isempty(p), parpool('local'); end
    use_parfor = true;
    disp(upper('PARFOR ENABLED.'))
  catch ME
    disp(upper(['PARFOR DISABLED: ' ME.message]))
    use_parfor = false;
  end
end

function [f_osc, f_axis, ap_offset, ap_exponent, fooof_r2, fooof_err, fooof_n_segments] = ...
    run_fooof_from_raw(data_td, trial_idx, ch_labels, time_win, cfg_fooof, fooof_mode, welch_seg_len_sec, welch_overlap)
  f_osc = []; f_axis = []; ap_offset = NaN; ap_exponent = NaN; fooof_r2 = NaN; fooof_err = NaN; fooof_n_segments = NaN;
  if isempty(data_td) || isempty(ch_labels), return, end
  try
    cfg_sel = []; cfg_sel.trials = trial_idx;
    d = ft_selectdata(cfg_sel, data_td);
    cfg_lat = []; cfg_lat.latency = time_win;
    d = ft_selectdata(cfg_lat, d);
    cfg_ch = []; cfg_ch.channel = ch_labels; cfg_ch.avgoverchan = 'yes';
    d = ft_selectdata(cfg_ch, d);
    d.label = {'ROI'};
    if strcmp(fooof_mode, 'welch500_50')
      [d_fooof, fooof_n_segments] = build_welch_segments_from_roi(d, welch_seg_len_sec, welch_overlap);
    else
      d_fooof = d;
      fooof_n_segments = 1;
    end
    if ~exist('ft_freqanalysis_Arne_FOOOF', 'file'), return, end
    out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, d_fooof);
    f_axis = out.freq(:)';
    if iscell(out.fooofparams), rep = out.fooofparams{1}; else, rep = out.fooofparams; end
    if isfield(rep, 'aperiodic_params') && numel(rep.aperiodic_params) >= 2
      ap_offset = rep.aperiodic_params(1);
      ap_exponent = rep.aperiodic_params(2);
    end
    if isfield(rep, 'r_squared'), fooof_r2 = rep.r_squared; end
    if isfield(rep, 'error'), fooof_err = rep.error; end
    if isfield(rep, 'fooofed_spectrum') && ~isempty(rep.fooofed_spectrum)
      f_osc = rep.fooofed_spectrum(:)';
    elseif isfield(rep, 'peak_params') && ~isempty(rep.peak_params)
      pk = rep.peak_params;
      g = zeros(size(out.freq(:)'));
      for p = 1:size(pk,1)
        g = g + pk(p,2) .* exp(-(out.freq(:)' - pk(p,1)).^2 ./ (2*pk(p,3)^2));
      end
      f_osc = g;
    end
  catch
    f_osc = []; f_axis = [];
  end
end

function [d_seg, n_segments] = build_welch_segments_from_roi(d_roi, seg_len_sec, overlap_frac)
  d_seg = d_roi;
  n_segments = 0;
  if isempty(d_roi) || ~isfield(d_roi, 'trial') || isempty(d_roi.trial) || isempty(d_roi.trial{1})
    return
  end
  x = double(d_roi.trial{1});
  if size(x,1) > 1, x = mean(x, 1, 'omitnan'); end
  if isfield(d_roi, 'fsample') && ~isempty(d_roi.fsample)
    fs = double(d_roi.fsample);
  else
    if ~isfield(d_roi, 'time') || isempty(d_roi.time) || numel(d_roi.time{1}) < 2, return, end
    fs = 1 / median(diff(double(d_roi.time{1})));
  end
  seg_n = max(2, round(seg_len_sec * fs));
  hop_n = max(1, round(seg_n * (1 - overlap_frac)));
  n = numel(x);
  if n < seg_n, d_seg = d_roi; n_segments = 1; return, end
  starts = 1:hop_n:(n - seg_n + 1);
  n_segments = numel(starts);
  d_seg.trial = cell(1, n_segments);
  d_seg.time = cell(1, n_segments);
  d_seg.sampleinfo = nan(n_segments, 2);
  for k = 1:n_segments
    idx = starts(k):(starts(k) + seg_n - 1);
    d_seg.trial{k} = x(idx);
    if isfield(d_roi, 'time') && ~isempty(d_roi.time) && ~isempty(d_roi.time{1})
      t0 = d_roi.time{1}(idx(1));
      d_seg.time{k} = t0 + (0:seg_n-1) ./ fs;
    else
      d_seg.time{k} = (0:seg_n-1) ./ fs;
    end
    d_seg.sampleinfo(k,:) = [1 seg_n];
  end
  d_seg.label = {'ROI'};
  d_seg.fsample = fs;
end
