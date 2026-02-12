%% AOC Multiverse — Trial-Level Data Preparation (Sternberg & N-Back)
% One-shot build: loads/computes EEG and gaze per trial for all
% electrode/FOOOF/latency/alpha/gaze/baseline/freq-method combinations.
% Writes multiverse_sternberg.csv and multiverse_nback.csv.
%
% Decision grid (9 dimensions, 5120 universes per task):
%   Electrodes:     all, posterior, parietal, occipital (4)
%   1/f:            FOOOFed, non-FOOOFed (2)
%   Latency:        0-500, 0-1000, 0-2000, 1000-2000 ms (4)
%   Alpha band:     canonical 8-14 Hz, IAF (2)
%   Gaze measure:   fixations, SPL, velocity, microsaccades, BCEA (5)
%   EEG baseline:   raw, dB [-0.5 -0.25] s (2)
%   Gaze baseline:  raw, pct_change [-0.5 -0.25] s (2)
%   Freq method:    hanning, dpss (2)
%
% Model in R: alpha ~ gaze_value * Condition + (1|subjectID)

disp(upper('=== AOC MULTIVERSE PREP START ==='))

%% Setup: run startup and setup FIRST (startup may clear the workspace)
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

%% Paths (Science Cloud vs Mac) — defined AFTER startup/setup so clear all cannot wipe them
if ispc
    base_data = 'W:\Students\Arne\AOC';
    base_features = fullfile(base_data, 'data', 'features');
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
    base_features = fullfile(base_data, 'data', 'features');
end
out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir), out_dir = pwd; end
disp(upper(['Paths: base_data = ' base_data ', base_features = ' base_features ', out_dir = ' out_dir]))

if isempty(subjects)
    dirs = dir(base_features);
    dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {dirs.name};
    disp(upper(['Subject list from base_features: ' num2str(length(subjects)) ' folders.']))
end
if isempty(path_preproc)
    path_preproc = base_features;
    disp(upper('Using base_features as path_preproc (EEG and gaze under same root).'))
end

%% Decision options (9 dimensions)
electrodes_opts    = {'all', 'posterior', 'parietal', 'occipital'};
fooof_opts         = {'FOOOFed', 'nonFOOOFed'};
latency_opts       = {'0_500ms', '0_1000ms', '0_2000ms', '1000_2000ms'};
alpha_opts         = {'canonical', 'IAF'};
gaze_opts          = {'gaze_density', 'scan_path_length', 'gaze_velocity', 'microsaccades', 'BCEA'};
baseline_eeg_opts  = {'raw', 'dB'};
baseline_gaze_opts = {'raw', 'pct_change'};
freq_method_opts   = {'hanning', 'dpss'};
n_elec = 4; n_fooof = 2; n_lat = 4; n_alpha = 2; n_gaze = 5;
n_bl_eeg = 2; n_bl_gaze = 2; n_fm = 2;
n_universes = n_elec * n_fooof * n_lat * n_alpha * n_gaze * n_bl_eeg * n_bl_gaze * n_fm;
alphaRange = [8 14];
disp(upper(['Decision grid: ' num2str(n_universes) ' universes per task (9 dimensions).']))

%% Get channel labels and indices (from first subject with power file)
disp(upper('Resolving channel sets (occipital, parietal, posterior, all)...'))
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
    error('Could not load power_stern_early_trials from any subject.')
end
[idx_all, idx_post, idx_pari, idx_occ] = get_channel_indices(labels_master);
ch_sets = {idx_all, idx_post, idx_pari, idx_occ};
disp(upper(['Channels: all=' num2str(length(idx_all)) ' post=' num2str(length(idx_post)) ...
  ' pari=' num2str(length(idx_pari)) ' occ=' num2str(length(idx_occ))]))

%% Build Sternberg multiverse table
disp(upper('--- STERNBERG TASK ---'))
tbl_s = build_task_multiverse('sternberg', subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    baseline_eeg_opts, baseline_gaze_opts, freq_method_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_fm, n_universes, alphaRange);
disp(upper(['Sternberg table rows: ' num2str(height(tbl_s))]))
if height(tbl_s) == 0
  error('AOC_multiverse_prep:NoData', 'Sternberg table is empty.')
end
writetable(tbl_s, fullfile(out_dir, 'multiverse_sternberg.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_sternberg.csv')]))

%% Build N-back multiverse table
disp(upper('--- N-BACK TASK ---'))
tbl_n = build_task_multiverse('nback', subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    baseline_eeg_opts, baseline_gaze_opts, freq_method_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_fm, n_universes, alphaRange);
disp(upper(['N-back table rows: ' num2str(height(tbl_n))]))
if height(tbl_n) == 0
  error('AOC_multiverse_prep:NoData', 'N-back table is empty.')
end
writetable(tbl_n, fullfile(out_dir, 'multiverse_nback.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_nback.csv')]))

disp(upper('=== AOC MULTIVERSE PREP DONE ==='))

%% ========== LOCAL FUNCTIONS ==========

function [idx_all, idx_post, idx_pari, idx_occ] = get_channel_indices(labels)
  % Occipital: contains O or I; Parietal: contains P but NOT F; Posterior: union
  n = length(labels);
  idx_occ = []; idx_pari = [];
  for i = 1:n
    lb = labels{i};
    if contains(lb, 'O') || contains(lb, 'I'), idx_occ(end+1) = i; end
    if contains(lb, 'P') && ~contains(lb, 'F'), idx_pari(end+1) = i; end
  end
  idx_post = unique([idx_occ, idx_pari]);
  idx_all = 1:n;
  if isempty(idx_post), idx_post = idx_occ; end
  disp(['    Occipital (' num2str(length(idx_occ)) '): ' strjoin(labels(idx_occ), ', ')])
  disp(['    Parietal  (' num2str(length(idx_pari)) '): ' strjoin(labels(idx_pari), ', ')])
  disp(['    Posterior (' num2str(length(idx_post)) '): ' strjoin(labels(idx_post), ', ')])
end

function a = bandpower_trials(pow, chIdx, band)
  if isempty(chIdx), a = nan(size(pow.powspctrm, 1), 1); return, end
  f = pow.freq(:);
  bandIdx = f >= band(1) & f <= band(2);
  if sum(bandIdx) == 0, a = nan(size(pow.powspctrm, 1), 1); return, end
  x = pow.powspctrm(:, chIdx, bandIdx);
  a = squeeze(mean(mean(x, 2), 3));
  a = a(:);
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

function nfix = count_fixations_in_window(x, y, t, t_win, fsample, vel_thresh_px, min_dur_samp)
  if nargin < 7, min_dur_samp = 25; end
  if nargin < 6, vel_thresh_px = 30; end
  idx = t >= t_win(1) & t <= t_win(2);
  xw = x(idx); yw = y(idx); xw = xw(:)'; yw = yw(:)';
  valid = isfinite(xw) & isfinite(yw);
  if sum(valid) < min_dur_samp, nfix = 0; return, end
  vel = sqrt(diff(xw).^2 + diff(yw).^2) * fsample;
  low = [false, vel < vel_thresh_px];
  low(~[valid(1:end-1); valid(2:end)]) = false;
  runs = diff([0 low 0]);
  onsets = find(runs == 1); offsets = find(runs == -1);
  if length(onsets) ~= length(offsets), nfix = 0; return, end
  nfix = sum(offsets - onsets >= min_dur_samp);
end

function bcea = compute_bcea(xw, yw)
  k95 = chi2inv(0.95, 2) / 2;
  xw = double(xw(:)); yw = double(yw(:));
  if numel(xw) < 10, bcea = NaN; return, end
  sx = std(xw); sy = std(yw);
  if sx == 0 || sy == 0, bcea = NaN; return, end
  rho = corr(xw, yw);
  bcea = 2 * k95 * pi * sx * sy * sqrt(1 - rho^2);
end

function [gaze_raw, gaze_bl] = compute_gaze_one_window(x, y, t, tw, dur, fsample, t_base)
  gaze_raw = nan(1, 5); gaze_bl = nan(1, 5);
  idx = t >= tw(1) & t <= tw(2);
  xw = x(idx); yw = y(idx);
  validxy = isfinite(xw) & isfinite(yw);
  xw = xw(validxy); yw = yw(validxy);
  if length(xw) < 50, return, end
  nfix = count_fixations_in_window(x, y, t, tw, fsample);
  dx = diff(xw); dy = diff(yw);
  spl = sum(sqrt(dx.^2 + dy.^2), 'omitnan');
  vel = spl / dur;
  try
    [ms_rate, ms_det] = detect_microsaccades(fsample, [xw; yw], length(xw));
    if exist('ms_det', 'var') && isfield(ms_det, 'Onset')
      ms_cnt = numel(ms_det.Onset);
    else
      ms_cnt = ms_rate * dur;
    end
  catch, ms_cnt = NaN;
  end
  bcea_val = compute_bcea(xw, yw);
  gaze_raw = [nfix, spl, vel, ms_cnt, bcea_val];
  % Baseline
  idx_b = t >= t_base(1) & t <= t_base(2);
  xb = x(idx_b); yb = y(idx_b);
  validb = isfinite(xb) & isfinite(yb);
  xb = xb(validb); yb = yb(validb);
  dur_base = t_base(2) - t_base(1);
  if length(xb) < 20, return, end
  nfix_b = count_fixations_in_window(x, y, t, t_base, fsample);
  dx_b = diff(xb); dy_b = diff(yb);
  spl_b = sum(sqrt(dx_b.^2 + dy_b.^2), 'omitnan');
  vel_b = spl_b / dur_base;
  try
    [ms_rate_b, ms_det_b] = detect_microsaccades(fsample, [xb; yb], length(xb));
    if exist('ms_det_b', 'var') && isfield(ms_det_b, 'Onset')
      ms_cnt_b = numel(ms_det_b.Onset);
    else
      ms_cnt_b = ms_rate_b * dur_base;
    end
  catch, ms_cnt_b = NaN;
  end
  bcea_b = compute_bcea(xb, yb);
  base_vals = [nfix_b, spl_b, vel_b, ms_cnt_b, bcea_b];
  for g = 1:5
    if isfinite(base_vals(g)) && base_vals(g) ~= 0 && isfinite(gaze_raw(g))
      gaze_bl(g) = (gaze_raw(g) - base_vals(g)) / abs(base_vals(g)) * 100;
    end
  end
end

function pow = get_power_from_TFR(tfr_all, ind_cond, latency)
  % Extract power from TFR for a given latency window: average over time → rpt_chan_freq
  cfg = []; cfg.latency = latency; cfg.trials = ind_cond;
  sel = ft_selectdata(cfg, tfr_all);
  sel.powspctrm = mean(sel.powspctrm, 4);
  if isfield(sel, 'time'), sel = rmfield(sel, 'time'); end
  sel.dimord = 'rpt_chan_freq';
  pow = sel;
end

function pow = compute_pow_cond_window(data_td, cond_inds, time_win, cfg_freq)
  % Compute trial-level power from time-domain EEG for one condition and time window
  pow = [];
  if isempty(cond_inds), return, end
  try
    cfg_sel = []; cfg_sel.trials = cond_inds;
    d = ft_selectdata(cfg_sel, data_td);
    cfg_lat = []; cfg_lat.latency = time_win;
    d = ft_selectdata(cfg_lat, d);
    pow = ft_freqanalysis(cfg_freq, d);
  catch ME
    disp(upper(['    WARNING: compute_pow_cond_window failed: ' ME.message]))
    pow = [];
  end
end

function pow_bl = compute_db_baseline_spectra(pow_task, pow_base)
  % Spectrum-level dB baseline: 10*log10(task / mean_baseline)
  % Uses mean across trials for stable baseline estimate (short 250ms window)
  pow_bl = pow_task;
  base_mean = mean(pow_base.powspctrm, 1);  % 1 x chan x freq
  base_mean(base_mean <= 0) = NaN;
  pow_bl.powspctrm = 10 * log10(bsxfun(@rdivide, pow_task.powspctrm, base_mean));
end

function [f_osc, f_axis] = run_fooof_one_trial_full(pow_trial, chIdx, freq, cfg_fooof)
  % Run FOOOF on one trial's ROI-averaged spectrum; return full oscillatory spectrum
  f_osc = []; f_axis = freq(:)';
  if isempty(chIdx), return, end
  try
    spec = squeeze(mean(pow_trial.powspctrm(1, chIdx, :), 2));
    if all(~isfinite(spec)) || numel(spec) < 10, return, end
    fake = struct('powspctrm', reshape(spec, [1 1 numel(spec)]), 'freq', freq(:)', ...
      'label', {{'ROI'}}, 'dimord', 'rpt_chan_freq');
    if ~exist('ft_freqanalysis_Arne_FOOOF', 'file'), return, end
    out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, fake);
    f_axis = out.freq(:)';
    if iscell(out.fooofparams), rep = out.fooofparams{1}; else, rep = out.fooofparams; end
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
    f_osc = []; f_axis = freq(:)';
  end
end

%% ========== MAIN BUILD FUNCTION ==========

function tbl = build_task_multiverse(task_name, subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    baseline_eeg_opts, baseline_gaze_opts, freq_method_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_fm, n_universes, alphaRange)

  disp(upper(['Building multiverse table for task: ' task_name ' (' num2str(n_universes) ' universes)']))
  if strcmp(task_name, 'sternberg')
    cond_codes = [22 24 26]; cond_vals = [2 4 6];
    early_file = 'power_stern_early_trials.mat'; full_file = 'power_stern_full_trials.mat';
    tfr_file = 'tfr_stern_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_sternberg.mat';
    eeg_td_file = 'dataEEG_sternberg.mat';
    et_file = 'dataET_sternberg.mat'; et_var = 'dataETlong';
  else
    cond_codes = [21 22 23]; cond_vals = [1 2 3];
    early_file = 'power_nback_early_trials.mat'; full_file = 'power_nback_full_trials.mat';
    tfr_file = 'tfr_nback_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_nback.mat';
    eeg_td_file = 'dataEEG_nback.mat';
    et_file = 'dataET_nback.mat'; et_var = 'dataETlong';
  end

  % Freq analysis configs: {hanning, dpss}
  cfg_hann = []; cfg_hann.method = 'mtmfft'; cfg_hann.taper = 'hanning';
  cfg_hann.foilim = [2 40]; cfg_hann.pad = 5; cfg_hann.keeptrials = 'yes';
  cfg_dpss = []; cfg_dpss.method = 'mtmfft'; cfg_dpss.taper = 'dpss'; cfg_dpss.tapsmofrq = 2;
  cfg_dpss.foilim = [2 40]; cfg_dpss.pad = 5; cfg_dpss.keeptrials = 'yes';
  cfg_freq_all = {cfg_hann, cfg_dpss};

  % FOOOF configs per freq method
  cfg_fooof_hann = []; cfg_fooof_hann.method = 'mtmfft'; cfg_fooof_hann.taper = 'hanning';
  cfg_fooof_hann.foilim = [2 40]; cfg_fooof_hann.pad = 5; cfg_fooof_hann.output = 'fooof'; cfg_fooof_hann.keeptrials = 'no';
  cfg_fooof_dpss = []; cfg_fooof_dpss.method = 'mtmfft'; cfg_fooof_dpss.taper = 'dpss'; cfg_fooof_dpss.tapsmofrq = 2;
  cfg_fooof_dpss.foilim = [2 40]; cfg_fooof_dpss.pad = 5; cfg_fooof_dpss.output = 'fooof'; cfg_fooof_dpss.keeptrials = 'no';
  cfg_fooof_all = {cfg_fooof_hann, cfg_fooof_dpss};

  % Latency windows
  lat_windows = {[0 0.5], [0 1], [0 2], [1 2]};
  lat_durs    = [0.5, 1, 2, 1];
  t_base = [-0.5 -0.25];
  fsample = 500;

  % Preallocate output arrays
  task_cell = {}; u_id_cell = []; elec_cell = {}; fooof_cell = {}; lat_cell = {};
  alpha_type_cell = {}; gaze_meas_cell = {}; bl_eeg_cell = {}; bl_gaze_cell = {}; fm_cell = {};
  subjectID_cell = []; Trial_cell = []; Condition_cell = [];
  alpha_val_cell = []; gaze_val_cell = [];

  for s = 1:length(subjects)
    sid = subjects{s};
    sid_num = str2double(sid);
    if isnan(sid_num), continue; end
    disp(upper(['Subject ' sid ' (' num2str(s) '/' num2str(length(subjects)) ')']))

    eeg_dir = fullfile(path_preproc, sid, 'eeg');
    gaze_dir = fullfile(base_features, sid, 'gaze');
    if ~isfolder(eeg_dir), eeg_dir = fullfile(base_features, sid, 'eeg'); end
    if ~isfolder(gaze_dir), gaze_dir = fullfile(path_preproc, sid, 'gaze'); end

    %% ====== EEG: Load pre-computed power (Hanning, raw + _bl) ======
    disp(upper(['  Loading EEG: ' eeg_dir]))
    early_path = fullfile(eeg_dir, early_file);
    full_path  = fullfile(eeg_dir, full_file);
    if ~isfile(early_path) || ~isfile(full_path)
      disp(upper('  Skip: missing power files.')); continue
    end
    load(early_path); load(full_path);
    disp(upper('  Loaded power early (0-1 s) and full (0-2 s).'))

    if strcmp(task_name, 'sternberg')
      pow_early = {powload2_early, powload4_early, powload6_early};
      pow_full  = {powload2_full,  powload4_full,  powload6_full};
      if exist('powload2_early_bl','var')
        pow_early_bl = {powload2_early_bl, powload4_early_bl, powload6_early_bl};
      else
        disp(upper('  WARNING: _bl early not found.')); pow_early_bl = {[], [], []};
      end
      if exist('powload2_full_bl','var')
        pow_full_bl = {powload2_full_bl, powload4_full_bl, powload6_full_bl};
      else
        disp(upper('  WARNING: _bl full not found.')); pow_full_bl = {[], [], []};
      end
    else
      pow_early = {powload1_early, powload2_early, powload3_early};
      pow_full  = {powload1_full,  powload2_full,  powload3_full};
      if exist('powload1_early_bl','var')
        pow_early_bl = {powload1_early_bl, powload2_early_bl, powload3_early_bl};
      else
        pow_early_bl = {[], [], []};
      end
      if exist('powload1_full_bl','var')
        pow_full_bl = {powload1_full_bl, powload2_full_bl, powload3_full_bl};
      else
        pow_full_bl = {[], [], []};
      end
    end

    %% ====== Load TFR (Hanning, for arbitrary windows) ======
    tfr_path = fullfile(eeg_dir, tfr_file);
    eeg_tfr_path = fullfile(eeg_dir, eeg_tfr_file);
    eeg_td_path = fullfile(eeg_dir, eeg_td_file);

    tfr_loaded = []; tfr_bl_loaded = [];
    if isfile(tfr_path)
      disp(upper('  Loading precomputed TFR (Hanning).'))
      load(tfr_path, 'tfr_all');
      tfr_loaded = tfr_all;
      cfg_bl = []; cfg_bl.baseline = t_base; cfg_bl.baselinetype = 'db';
      tfr_bl_loaded = ft_freqbaseline(cfg_bl, tfr_loaded);
      disp(upper('  TFR loaded + dB-baselined copy created.'))
    end

    %% ====== Load time-domain EEG (for DPSS + Hanning fallback) ======
    data_td = [];
    if isfile(eeg_tfr_path)
      disp(upper('  Loading time-domain EEG (dataEEG_TFR).'))
      load(eeg_tfr_path, 'dataTFR');
      data_td = dataTFR;
    elseif isfile(eeg_td_path)
      disp(upper('  Loading time-domain EEG (dataEEG).'))
      load(eeg_td_path, 'dataEEG');
      if ~exist('dataEEG', 'var')
        load(eeg_td_path, 'dataTFR'); data_td = dataTFR;
      else
        data_td = dataEEG;
      end
    end
    if isempty(data_td) && isempty(tfr_loaded)
      disp(upper('  WARNING: No TFR or time-domain EEG. Hanning 0-500ms/1-2s and all DPSS will be NaN.'))
    end

    %% ====== Subject IAF (from Hanning full, occipital, raw) ======
    disp(upper('  Computing subject IAF from full-latency occipital power.'))
    IAF_band = get_IAF_band(pow_full{1}, ch_sets{4}, alphaRange);

    %% ====== Load gaze ======
    disp(upper(['  Loading gaze: ' gaze_dir]))
    et_path = fullfile(gaze_dir, et_file);
    if ~isfile(et_path)
      disp(upper('  Skip: missing ET file.')); continue
    end
    load(et_path);
    if ~exist(et_var, 'var'), et_var = 'dataET'; end
    if ~exist(et_var, 'var'), disp(upper('  Skip: ET variable not found.')); continue; end
    dataET = eval(et_var);
    disp(upper(['  Gaze loaded: ' num2str(size(dataET.trialinfo, 1)) ' trials.']))

    trial_num = dataET.trialinfo(:,2);
    cond_code = dataET.trialinfo(:,1);
    if strcmp(task_name, 'sternberg'), cond_code = cond_code - 20; end
    if strcmp(task_name, 'nback'),     cond_code = cond_code - 20; end

    %% ====== Loop over conditions ======
    for cond = 1:3
      cval = cond_vals(cond);
      trl_idx = find(cond_code == cval);
      if isempty(trl_idx), continue; end

      n_trials_cond = length(trl_idx);
      trials_list = trial_num(trl_idx);

      % EEG-ET trial alignment
      n_eeg_trials = size(pow_early{cond}.powspctrm, 1);
      if n_trials_cond ~= n_eeg_trials
        disp(upper(['  WARNING: ET=' num2str(n_trials_cond) ' vs EEG=' num2str(n_eeg_trials) ' for cond ' num2str(cval) '. Using min.']))
        n_trials_cond = min(n_trials_cond, n_eeg_trials);
        trl_idx = trl_idx(1:n_trials_cond);
        trials_list = trial_num(trl_idx);
      end
      disp(upper(['  Condition ' num2str(cval) ': ' num2str(n_trials_cond) ' trials.']))

      %% ====== Build power struct grid: pow_s{lat, fm}, pow_bl_s{lat, fm} ======
      % lat: 1=0-500ms, 2=0-1s, 3=0-2s, 4=1-2s
      % fm:  1=hanning,  2=dpss
      pow_s    = cell(4, 2);
      pow_bl_s = cell(4, 2);

      % --- Pre-computed Hanning (0-1s, 0-2s) ---
      pow_s{2, 1} = pow_early{cond};       % Hanning, 0-1s, raw
      pow_s{3, 1} = pow_full{cond};        % Hanning, 0-2s, raw
      if ~isempty(pow_early_bl{cond}), pow_bl_s{2, 1} = pow_early_bl{cond}; end  % Hanning, 0-1s, dB
      if ~isempty(pow_full_bl{cond}),  pow_bl_s{3, 1} = pow_full_bl{cond};  end  % Hanning, 0-2s, dB

      % --- TFR-based Hanning (0-500ms, 1-2s) ---
      if ~isempty(tfr_loaded)
        cond_inds_tfr = find(tfr_loaded.trialinfo(:,1) == cond_codes(cond));
        if isempty(pow_s{1, 1})
          pow_s{1, 1} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{1});
        end
        if isempty(pow_s{4, 1})
          pow_s{4, 1} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{4});
        end
        if ~isempty(tfr_bl_loaded)
          cond_inds_bl = find(tfr_bl_loaded.trialinfo(:,1) == cond_codes(cond));
          if isempty(pow_bl_s{1, 1})
            pow_bl_s{1, 1} = get_power_from_TFR(tfr_bl_loaded, cond_inds_bl, lat_windows{1});
          end
          if isempty(pow_bl_s{4, 1})
            pow_bl_s{4, 1} = get_power_from_TFR(tfr_bl_loaded, cond_inds_bl, lat_windows{4});
          end
        end
      end

      % --- Time-domain fill (all DPSS + remaining Hanning gaps) ---
      if ~isempty(data_td)
        cond_inds_td = find(data_td.trialinfo(:,1) == cond_codes(cond));
        for ifm = 1:2
          cfg_f = cfg_freq_all{ifm};
          for il = 1:4
            if isempty(pow_s{il, ifm})
              disp(upper(['    Computing power from time-domain: lat=' num2str(il) ' fm=' num2str(ifm)]))
              pow_s{il, ifm} = compute_pow_cond_window(data_td, cond_inds_td, lat_windows{il}, cfg_f);
            end
          end
          % Baseline spectrum for dB correction (mean across trials)
          pow_base_fm = compute_pow_cond_window(data_td, cond_inds_td, t_base, cfg_f);
          for il = 1:4
            if isempty(pow_bl_s{il, ifm}) && ~isempty(pow_s{il, ifm}) && ~isempty(pow_base_fm)
              pow_bl_s{il, ifm} = compute_db_baseline_spectra(pow_s{il, ifm}, pow_base_fm);
            end
          end
        end
      end

      %% ====== Alpha: cell array {lat, fm, bl, fo}(trial, elec_alpha_col) ======
      % fo: 1=FOOOFed, 2=nonFOOOFed
      n_cols = n_elec * n_alpha;
      alpha = cell(n_lat, n_fm, n_bl_eeg, n_fooof);
      for il = 1:n_lat; for ifm = 1:n_fm; for ibl = 1:n_bl_eeg; for ifo = 1:n_fooof
        alpha{il, ifm, ibl, ifo} = nan(n_trials_cond, n_cols);
      end; end; end; end

      disp(upper('  Computing alpha (raw + FOOOF × hanning/dpss × raw/dB) per trial.'))
      for il = 1:n_lat
        for ifm = 1:n_fm
          for ibl = 1:n_bl_eeg
            if ibl == 1, p = pow_s{il, ifm}; else, p = pow_bl_s{il, ifm}; end
            if isempty(p) || ~isfield(p, 'powspctrm'), continue; end

            % nonFOOOFed: vectorized bandpower
            for ie = 1:n_elec
              chIdx = ch_sets{ie};
              for ia = 1:n_alpha
                col = (ie-1)*n_alpha + ia;
                if ia == 1, band = alphaRange; else, band = IAF_band; end
                alpha{il, ifm, ibl, 2}(:, col) = bandpower_trials(p, chIdx, band);
              end
            end

            % FOOOFed: run once per (trial, electrode_set), extract both alpha bands
            cfg_fo = cfg_fooof_all{ifm};
            for ie = 1:n_elec
              chIdx = ch_sets{ie};
              for tr = 1:min(n_trials_cond, size(p.powspctrm, 1))
                pow_tr = ft_selectdata(struct('trials', tr), p);
                [f_osc, f_axis] = run_fooof_one_trial_full(pow_tr, chIdx, p.freq, cfg_fo);
                if isempty(f_osc), continue, end
                for ia = 1:n_alpha
                  col = (ie-1)*n_alpha + ia;
                  if ia == 1, band = alphaRange; else, band = IAF_band; end
                  bandIdx = f_axis >= band(1) & f_axis <= band(2);
                  if sum(bandIdx) > 0
                    alpha{il, ifm, ibl, 1}(tr, col) = mean(f_osc(bandIdx), 'omitnan');
                  end
                end
              end
            end
          end
        end
      end

      %% ====== Gaze: 5 metrics × 4 windows × raw + pct_change ======
      disp(upper('  Computing gaze (fix, SPL, vel, MS, BCEA) × 4 windows × raw + pct_change.'))
      gaze_raw_all = cell(n_lat, 1); gaze_bl_all = cell(n_lat, 1);
      for il = 1:n_lat
        gaze_raw_all{il} = nan(n_trials_cond, 5);
        gaze_bl_all{il}  = nan(n_trials_cond, 5);
      end
      for tr = 1:n_trials_cond
        trl_glob = trl_idx(tr);
        raw_et = dataET.trial{trl_glob};
        t = dataET.time{trl_glob};
        if size(raw_et,1) < 2, continue; end
        x = raw_et(1,:); y = raw_et(2,:);
        if length(t) ~= length(x), continue; end
        for il = 1:n_lat
          [gaze_raw_all{il}(tr,:), gaze_bl_all{il}(tr,:)] = ...
            compute_gaze_one_window(x, y, t, lat_windows{il}, lat_durs(il), fsample, t_base);
        end
      end

      %% ====== Append rows: trial × universe ======
      disp(upper(['  Appending ' num2str(n_trials_cond * n_universes) ' rows for subject ' sid ' cond ' num2str(cval) '.']))
      for tr = 1:n_trials_cond
        for u = 1:n_universes
          [ie, ifo, il, ia, ig, ibe, ibg, ifm] = ind2sub( ...
            [n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_fm], u);
          task_cell{end+1, 1}       = task_name;
          u_id_cell(end+1, 1)       = u;
          elec_cell{end+1, 1}       = electrodes_opts{ie};
          fooof_cell{end+1, 1}      = fooof_opts{ifo};
          lat_cell{end+1, 1}        = latency_opts{il};
          alpha_type_cell{end+1, 1} = alpha_opts{ia};
          gaze_meas_cell{end+1, 1}  = gaze_opts{ig};
          bl_eeg_cell{end+1, 1}     = baseline_eeg_opts{ibe};
          bl_gaze_cell{end+1, 1}    = baseline_gaze_opts{ibg};
          fm_cell{end+1, 1}         = freq_method_opts{ifm};
          subjectID_cell(end+1, 1)  = sid_num;
          Trial_cell(end+1, 1)      = trials_list(tr);
          Condition_cell(end+1, 1)  = cval;

          % Alpha lookup
          col_alpha = (ie-1)*n_alpha + ia;
          av = alpha{il, ifm, ibe, ifo}(tr, col_alpha);
          alpha_val_cell(end+1, 1) = av;

          % Gaze lookup
          if ibg == 1
            gv = gaze_raw_all{il}(tr, ig);
          else
            gv = gaze_bl_all{il}(tr, ig);
          end
          gaze_val_cell(end+1, 1) = gv;
        end
      end
    end
  end

  tbl = table(task_cell, u_id_cell, elec_cell, fooof_cell, lat_cell, alpha_type_cell, gaze_meas_cell, ...
    bl_eeg_cell, bl_gaze_cell, fm_cell, subjectID_cell, Trial_cell, Condition_cell, alpha_val_cell, gaze_val_cell, ...
    'VariableNames', {'task', 'universe_id', 'electrodes', 'fooof', 'latency_ms', 'alpha_type', 'gaze_measure', ...
    'baseline_eeg', 'baseline_gaze', 'freq_method', 'subjectID', 'Trial', 'Condition', 'alpha', 'gaze_value'});
end
