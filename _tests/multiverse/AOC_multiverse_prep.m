%% AOC Multiverse — Trial-Level Data Preparation (Sternberg & N-Back)
% One-shot build: loads/computes EEG (0-0.5, 0-1, 0-2 s) and gaze per trial,
% all electrode/FOOOF/latency/alpha/gaze combinations. Writes multiverse_*_sternberg.csv and multiverse_*_nback.csv.
% Electrodes: occipital (O|I), parietal (P), posterior (PO|O|I), all. Bilateral.
% 1/f: raw vs FOOOF (model minus aperiodic). Alpha: canonical 8-14 Hz vs IAF.
% Gaze: density (fixations/trial), SPL normalized, velocity, microsaccades/trial.
% Model in R: alpha ~ gaze_value * Condition + (1|subjectID).

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

%% Decision options
electrodes_opts = {'all', 'posterior', 'parietal', 'occipital'};
fooof_opts     = {'FOOOFed', 'nonFOOOFed'};
latency_opts   = {'0_500ms', '0_1000ms', '0_2000ms'};
alpha_opts     = {'canonical', 'IAF'};
gaze_opts      = {'gaze_density', 'scan_path_length', 'gaze_velocity', 'microsaccades'};
n_elec = 4; n_fooof = 2; n_lat = 3; n_alpha = 2; n_gaze = 4;
n_universes = n_elec * n_fooof * n_lat * n_alpha * n_gaze;
alphaRange = [8 14];
disp(upper(['Decision grid: ' num2str(n_universes) ' universes per task.']))

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
    error('Could not load power_stern_early_trials from any subject. Run EEG trial-level extraction first.')
end
[idx_all, idx_post, idx_pari, idx_occ] = get_channel_indices(labels_master);
ch_sets = {idx_all, idx_post, idx_pari, idx_occ};
disp(upper(['Channels: all=' num2str(length(idx_all)) ' post=' num2str(length(idx_post)) ' pari=' num2str(length(idx_pari)) ' occ=' num2str(length(idx_occ))]))

%% Build Sternberg multiverse table
disp(upper('--- STERNBERG TASK ---'))
tbl_s = build_task_multiverse('sternberg', subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_universes, alphaRange);
disp(upper(['Sternberg table rows: ' num2str(height(tbl_s))]))
if height(tbl_s) == 0
  error('AOC_multiverse_prep:NoData', 'Sternberg table is empty. Ensure at least one subject has both power_stern_*_trials.mat and dataET_sternberg.mat under base_features/<subjectID>/eeg and .../gaze.')
end
writetable(tbl_s, fullfile(out_dir, 'multiverse_sternberg.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_sternberg.csv')]))

%% Build N-back multiverse table
disp(upper('--- N-BACK TASK ---'))
tbl_n = build_task_multiverse('nback', subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_universes, alphaRange);
disp(upper(['N-back table rows: ' num2str(height(tbl_n))]))
if height(tbl_n) == 0
  error('AOC_multiverse_prep:NoData', 'N-back table is empty. Ensure at least one subject has both power_nback_*_trials.mat and dataET_nback.mat under base_features/<subjectID>/eeg and .../gaze.')
end
writetable(tbl_n, fullfile(out_dir, 'multiverse_nback.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_nback.csv')]))

disp(upper('=== AOC MULTIVERSE PREP DONE ==='))

%% ========== LOCAL FUNCTIONS ==========

function [idx_all, idx_post, idx_pari, idx_occ] = get_channel_indices(labels)
  n = length(labels);
  idx_occ = []; idx_pari = []; idx_post = [];
  for i = 1:n
    lb = labels{i};
    if contains(lb, 'O') || contains(lb, 'I'), idx_occ(end+1) = i; end
    if contains(lb, 'P'), idx_pari(end+1) = i; end
    if contains(lb, 'PO') || contains(lb, 'O') || contains(lb, 'I'), idx_post(end+1) = i; end
  end
  idx_all = 1:n;
  if isempty(idx_post), idx_post = idx_occ; end
end

function a = bandpower_trials(pow, chIdx, band)
  if isempty(chIdx)
    a = nan(size(pow.powspctrm, 1), 1); return
  end
  f = pow.freq(:);
  bandIdx = f >= band(1) & f <= band(2);
  if sum(bandIdx) == 0
    a = nan(size(pow.powspctrm, 1), 1); return
  end
  x = pow.powspctrm(:, chIdx, bandIdx);
  a = squeeze(mean(mean(x, 2), 3));
  if iscolumn(a), a = a'; end
  a = a(:);
end

function IAF_band = get_IAF_band(pow_full, chIdx, alphaRange)
  if isempty(chIdx), IAF_band = alphaRange; return, end
  spec = squeeze(mean(mean(pow_full.powspctrm(:, chIdx, :), 1), 2));
  f = pow_full.freq(:);
  aIdx = f >= alphaRange(1) & f <= alphaRange(2);
  alphaP = spec(aIdx);
  [pks, locs] = findpeaks(alphaP);
  if isempty(pks), IAF_band = alphaRange; return, end
  [~, im] = max(pks);
  IAF = f(aIdx(locs(im)));
  if IAF <= alphaRange(1) || IAF >= alphaRange(2), IAF_band = alphaRange; return, end
  IAF_band = [max(f(1), IAF-4), min(f(end), IAF+2)];
end

function nfix = count_fixations_in_window(x, y, t, t_win, fsample, vel_thresh_px, min_dur_samp)
  if nargin < 7, min_dur_samp = 25; end
  if nargin < 6, vel_thresh_px = 30; end
  idx = t >= t_win(1) & t <= t_win(2);
  xw = x(idx); yw = y(idx);
  xw = xw(:)'; yw = yw(:)';
  valid = isfinite(xw) & isfinite(yw);
  if sum(valid) < min_dur_samp, nfix = 0; return, end
  vel = sqrt(diff(xw).^2 + diff(yw).^2) * fsample;
  low = [false, vel < vel_thresh_px];
  low(~[valid(1:end-1); valid(2:end)]) = false;
  runs = diff([0 low 0]);
  onsets = find(runs == 1);
  offsets = find(runs == -1);
  if length(onsets) ~= length(offsets), nfix = 0; return, end
  nfix = sum(offsets - onsets >= min_dur_samp);
end

function pow_050 = get_power_050_from_TFR(tfr_all, ind_cond, latency_050)
  cfg = []; cfg.latency = latency_050; cfg.trials = ind_cond;
  sel = ft_selectdata(cfg, tfr_all);
  sel.powspctrm = mean(sel.powspctrm, 4);
  if isfield(sel, 'time'), sel = rmfield(sel, 'time'); end
  sel.dimord = 'rpt_chan_freq';
  pow_050 = sel;
end

function [pow_050_2, pow_050_4, pow_050_6] = get_power_050_from_timedomain(dataEEG, cond_codes, t_050, cfg_mtmfft)
  % Take long trial data, select 0-0.5 s epoch, run mtmfft → trial-level spectrum (rpt_chan_freq).
  ind2 = find(dataEEG.trialinfo(:,1) == cond_codes(1));
  ind4 = find(dataEEG.trialinfo(:,1) == cond_codes(2));
  ind6 = find(dataEEG.trialinfo(:,1) == cond_codes(3));
  cfg_lat = []; cfg_lat.latency = t_050;
  cfg_t2 = []; cfg_t2.trials = ind2; d2 = ft_selectdata(cfg_t2, dataEEG); d2 = ft_selectdata(cfg_lat, d2);
  cfg_t4 = []; cfg_t4.trials = ind4; d4 = ft_selectdata(cfg_t4, dataEEG); d4 = ft_selectdata(cfg_lat, d4);
  cfg_t6 = []; cfg_t6.trials = ind6; d6 = ft_selectdata(cfg_t6, dataEEG); d6 = ft_selectdata(cfg_lat, d6);
  pow_050_2 = ft_freqanalysis(cfg_mtmfft, d2);
  pow_050_4 = ft_freqanalysis(cfg_mtmfft, d4);
  pow_050_6 = ft_freqanalysis(cfg_mtmfft, d6);
end

function alpha_fooof = run_fooof_one_trial(pow_trial, chIdx, freq, band, cfg_fooof)
  alpha_fooof = NaN;
  if isempty(chIdx), return, end
  try
    spec = squeeze(mean(pow_trial.powspctrm(1, chIdx, :), 2));
    if all(~isfinite(spec)) || numel(spec) < 10, return, end
    fake = struct('powspctrm', reshape(spec, [1 1 numel(spec)]), 'freq', freq(:)', ...
      'label', {{'ROI'}}, 'dimord', 'rpt_chan_freq');
    if exist('ft_freqanalysis_Arne_FOOOF', 'file')
      out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, fake);
      if iscell(out.fooofparams), rep = out.fooofparams{1}; else, rep = out.fooofparams; end
      if isfield(rep, 'fooofed_spectrum') && ~isempty(rep.fooofed_spectrum)
        f_osc = rep.fooofed_spectrum(:);
      elseif isfield(rep, 'peak_params') && ~isempty(rep.peak_params)
        ap = rep.aperiodic_params(:);
        if numel(ap) >= 2
          ap_fit = ap(1) - ap(2) .* log10(out.freq(:));
        else
          ap_fit = zeros(size(out.freq(:)));
        end
        pk = rep.peak_params;
        g = zeros(size(out.freq(:)));
        for p = 1:size(pk,1)
          g = g + pk(p,2) .* exp(-(out.freq(:) - pk(p,1)).^2 ./ (2*pk(p,3)^2));
        end
        f_osc = ap_fit + g - ap_fit;
      else
        return
      end
      bandIdx = out.freq(:) >= band(1) & out.freq(:) <= band(2);
      if sum(bandIdx) > 0
        alpha_fooof = mean(f_osc(bandIdx), 'omitnan');
      end
    end
  catch
    alpha_fooof = NaN;
  end
end

function tbl = build_task_multiverse(task_name, subjects, path_preproc, base_features, ...
    ch_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_universes, alphaRange)

  disp(upper(['Building multiverse table for task: ' task_name]))
  if strcmp(task_name, 'sternberg')
    cond_codes = [22 24 26]; cond_vals = [2 4 6];
    early_file = 'power_stern_early_trials.mat'; full_file = 'power_stern_full_trials.mat';
    tfr_file = 'tfr_stern_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_sternberg.mat';
    et_file = 'dataET_sternberg.mat'; et_var = 'dataETlong';
  else
    cond_codes = [21 22 23]; cond_vals = [1 2 3];
    early_file = 'power_nback_early_trials.mat'; full_file = 'power_nback_full_trials.mat';
    tfr_file = 'tfr_nback_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_nback.mat';
    et_file = 'dataET_nback.mat'; et_var = 'dataETlong';
  end

  cfg_fooof = [];
  cfg_fooof.method = 'mtmfft'; cfg_fooof.taper = 'hanning';
  cfg_fooof.foilim = [2 40]; cfg_fooof.pad = 5; cfg_fooof.output = 'fooof';
  cfg_fooof.keeptrials = 'no';

  cfg_mtmfft = [];
  cfg_mtmfft.method = 'mtmfft'; cfg_mtmfft.taper = 'hanning';
  cfg_mtmfft.foilim = [2 40]; cfg_mtmfft.pad = 5; cfg_mtmfft.keeptrials = 'yes';

  task_cell = {}; u_id_cell = []; elec_cell = {}; fooof_cell = {}; lat_cell = {};
  alpha_type_cell = {}; gaze_meas_cell = {}; subjectID_cell = []; Trial_cell = []; Condition_cell = [];
  alpha_val_cell = []; gaze_val_cell = [];

  fsample = 500;
  t_050 = [0 0.5]; t_100 = [0 1]; t_200 = [0 2];
  dur_050 = 0.5; dur_100 = 1; dur_200 = 2;

  for s = 1:length(subjects)
    sid = subjects{s};
    sid_num = str2double(sid);
    if isnan(sid_num), continue; end
    disp(upper(['Subject ' sid ' (' num2str(s) '/' num2str(length(subjects)) ')']))

    eeg_dir = fullfile(path_preproc, sid, 'eeg');
    gaze_dir = fullfile(base_features, sid, 'gaze');
    if ~isfolder(eeg_dir), eeg_dir = fullfile(base_features, sid, 'eeg'); end
    if ~isfolder(gaze_dir), gaze_dir = fullfile(path_preproc, sid, 'gaze'); end

    %% Load EEG: power early/full; TFR for 0-0.5 or compute from raw TFR
    disp(upper(['  Loading EEG: ' eeg_dir]))
    early_path = fullfile(eeg_dir, early_file);
    full_path  = fullfile(eeg_dir, full_file);
    if ~isfile(early_path) || ~isfile(full_path)
      disp(upper(['  Skip: missing power files.'])); continue
    end
    load(early_path);
    load(full_path);
    disp(upper('  Loaded power early (0-1 s) and full (0-2 s).'))
    if strcmp(task_name, 'sternberg')
      pow_early_2 = powload2_early; pow_early_4 = powload4_early; pow_early_6 = powload6_early;
      pow_full_2  = powload2_full;  pow_full_4  = powload4_full;  pow_full_6  = powload6_full;
    else
      pow_early_2 = powload1_early; pow_early_4 = powload2_early; pow_early_6 = powload3_early;
      pow_full_2  = powload1_full;  pow_full_4  = powload2_full;  pow_full_6  = powload3_full;
    end

    pow_050_2 = []; pow_050_4 = []; pow_050_6 = [];
    tfr_path = fullfile(eeg_dir, tfr_file);
    eeg_tfr_path = fullfile(eeg_dir, eeg_tfr_file);
    if strcmp(task_name, 'sternberg')
      eeg_td_file = 'dataEEG_sternberg.mat';
    else
      eeg_td_file = 'dataEEG_nback.mat';
    end
    eeg_td_path = fullfile(eeg_dir, eeg_td_file);
    if isfile(tfr_path)
      disp(upper('  Loading TFR for 0-0.5 s from precomputed file.'))
      load(tfr_path, 'tfr_all');
      ind2 = find(tfr_all.trialinfo(:,1) == cond_codes(1));
      ind4 = find(tfr_all.trialinfo(:,1) == cond_codes(2));
      ind6 = find(tfr_all.trialinfo(:,1) == cond_codes(3));
      pow_050_2 = get_power_050_from_TFR(tfr_all, ind2, t_050);
      pow_050_4 = get_power_050_from_TFR(tfr_all, ind4, t_050);
      pow_050_6 = get_power_050_from_TFR(tfr_all, ind6, t_050);
    else
      data_td = [];
      if isfile(eeg_tfr_path)
        disp(upper('  Loading time-domain EEG (dataEEG_TFR); computing 0-0.5 s power from epoch (mtmfft).'))
        load(eeg_tfr_path, 'dataTFR');
        data_td = dataTFR;
      elseif isfile(eeg_td_path)
        disp(upper('  Loading time-domain EEG (dataEEG); computing 0-0.5 s power from epoch (mtmfft).'))
        load(eeg_td_path, 'dataEEG');
        if ~exist('dataEEG', 'var'), load(eeg_td_path, 'dataTFR'); data_td = dataTFR; else, data_td = dataEEG; end
      end
      if ~isempty(data_td)
        [pow_050_2, pow_050_4, pow_050_6] = get_power_050_from_timedomain(data_td, cond_codes, t_050, cfg_mtmfft);
        disp(upper('  0-500 ms power computed from long EEG epoch.'))
      else
        disp(upper('  No TFR or time-domain EEG found; 0-500 ms alpha will be NaN. Run preprocessing (dataEEG_* or dataEEG_TFR_*).'))
      end
    end

    %% Subject IAF from full (occipital)
    disp(upper('  Computing subject IAF from full-latency occipital power.'))
    IAF_band = get_IAF_band(pow_full_2, ch_sets{4}, alphaRange);

    %% Load gaze
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

    nTrials = size(dataET.trialinfo, 1);
    trial_num = dataET.trialinfo(:,2);
    cond_code = dataET.trialinfo(:,1);
    if strcmp(task_name, 'sternberg'), cond_code = cond_code - 20; end
    if strcmp(task_name, 'nback'), cond_code = cond_code - 20; end

    for cond = 1:3
      cval = cond_vals(cond);
      trl_idx = find(cond_code == cval);
      if isempty(trl_idx), continue; end

      if cond == 1
        pe_100 = pow_early_2; pf_200 = pow_full_2;
        pe_050 = pow_050_2; pf_050 = pow_050_2;
      elseif cond == 2
        pe_100 = pow_early_4; pf_200 = pow_full_4;
        pe_050 = pow_050_4; pf_050 = pow_050_4;
      else
        pe_100 = pow_early_6; pf_200 = pow_full_6;
        pe_050 = pow_050_6; pf_050 = pow_050_6;
      end

      n_trials_cond = length(trl_idx);
      trials_list = trial_num(trl_idx);
      disp(upper(['  Condition ' num2str(cval) ': ' num2str(n_trials_cond) ' trials. Computing alpha (raw + FOOOF) per trial per electrode set.']))

      %% Alpha: raw and FOOOF per trial, per electrode set, per latency, per alpha type
      alpha_raw_100 = nan(n_trials_cond, n_elec * n_alpha);
      alpha_raw_200 = nan(n_trials_cond, n_elec * n_alpha);
      alpha_raw_050 = nan(n_trials_cond, n_elec * n_alpha);
      alpha_fooof_100 = nan(n_trials_cond, n_elec * n_alpha);
      alpha_fooof_200 = nan(n_trials_cond, n_elec * n_alpha);
      alpha_fooof_050 = nan(n_trials_cond, n_elec * n_alpha);

      for ie = 1:n_elec
        chIdx = ch_sets{ie};
        for ia = 1:n_alpha
          col = (ie-1)*n_alpha + ia;
          if ia == 1, band = alphaRange; else, band = IAF_band; end
          alpha_raw_100(:, col) = bandpower_trials(pe_100, chIdx, band);
          alpha_raw_200(:, col) = bandpower_trials(pf_200, chIdx, band);
          if ~isempty(pf_050) && isfield(pf_050, 'powspctrm')
            alpha_raw_050(:, col) = bandpower_trials(pf_050, chIdx, band);
          end
          for tr = 1:n_trials_cond
            trl_sel = trl_idx(tr);
            pow1 = ft_selectdata(struct('trials', tr), pe_100);
            alpha_fooof_100(tr, col) = run_fooof_one_trial(pow1, chIdx, pe_100.freq, band, cfg_fooof);
            pow2 = ft_selectdata(struct('trials', tr), pf_200);
            alpha_fooof_200(tr, col) = run_fooof_one_trial(pow2, chIdx, pf_200.freq, band, cfg_fooof);
            if ~isempty(pf_050) && isfield(pf_050, 'powspctrm') && size(pf_050.powspctrm,1) >= tr
              pow05 = ft_selectdata(struct('trials', tr), pf_050);
              alpha_fooof_050(tr, col) = run_fooof_one_trial(pow05, chIdx, pf_050.freq, band, cfg_fooof);
            end
          end
        end
      end

      disp(upper('  Computing gaze metrics per trial (0-0.5, 0-1, 0-2 s): fixations, SPL norm, velocity, microsaccades.'))
      %% Gaze per trial for 0-0.5, 0-1, 0-2
      fix_050 = nan(n_trials_cond, 1); fix_100 = nan(n_trials_cond, 1); fix_200 = nan(n_trials_cond, 1);
      spl_050 = nan(n_trials_cond, 1); spl_100 = nan(n_trials_cond, 1); spl_200 = nan(n_trials_cond, 1);
      vel_050 = nan(n_trials_cond, 1); vel_100 = nan(n_trials_cond, 1); vel_200 = nan(n_trials_cond, 1);
      ms_050  = nan(n_trials_cond, 1); ms_100  = nan(n_trials_cond, 1); ms_200  = nan(n_trials_cond, 1);

      for tr = 1:n_trials_cond
        trl_glob = trl_idx(tr);
        raw = dataET.trial{trl_glob};
        t = dataET.time{trl_glob};
        if size(raw,1) < 2, continue; end
        x = raw(1,:); y = raw(2,:);
        if length(t) ~= length(x), continue; end

        for iw = 1:3
          if iw == 1, tw = t_050; dur = dur_050; elseif iw == 2, tw = t_100; dur = dur_100; else, tw = t_200; dur = dur_200; end
          idx = t >= tw(1) & t <= tw(2);
          xw = x(idx); yw = y(idx);
          xw = xw(isfinite(xw) & isfinite(yw)); yw = yw(isfinite(xw) & isfinite(yw));
          if length(xw) < 50
            if iw == 1, fix_050(tr)=0; spl_050(tr)=NaN; vel_050(tr)=NaN; ms_050(tr)=NaN; end
            if iw == 2, fix_100(tr)=0; spl_100(tr)=NaN; vel_100(tr)=NaN; ms_100(tr)=NaN; end
            if iw == 3, fix_200(tr)=0; spl_200(tr)=NaN; vel_200(tr)=NaN; ms_200(tr)=NaN; end
            continue
          end
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
          catch
            ms_cnt = NaN;
          end
          if iw == 1, fix_050(tr)=nfix; spl_050(tr)=spl/dur; vel_050(tr)=vel; ms_050(tr)=ms_cnt/dur; end
          if iw == 2, fix_100(tr)=nfix; spl_100(tr)=spl/dur; vel_100(tr)=vel; ms_100(tr)=ms_cnt/dur; end
          if iw == 3, fix_200(tr)=nfix; spl_200(tr)=spl/dur; vel_200(tr)=vel; ms_200(tr)=ms_cnt/dur; end
        end
      end

      gaze_050 = [fix_050, spl_050, vel_050, ms_050];
      gaze_100 = [fix_100, spl_100, vel_100, ms_100];
      gaze_200 = [fix_200, spl_200, vel_200, ms_200];

      disp(upper(['  Appending ' num2str(n_trials_cond * n_universes) ' rows for subject ' sid ' condition ' num2str(cval) '.']))
      %% Append rows for this subject x condition
      for tr = 1:n_trials_cond
        for u = 1:n_universes
          [ie, ifo, il, ia, ig] = ind2sub([n_elec, n_fooof, n_lat, n_alpha, n_gaze], u);
          task_cell{end+1, 1} = task_name;
          u_id_cell(end+1, 1) = u;
          elec_cell{end+1, 1} = electrodes_opts{ie};
          fooof_cell{end+1, 1} = fooof_opts{ifo};
          lat_cell{end+1, 1} = latency_opts{il};
          alpha_type_cell{end+1, 1} = alpha_opts{ia};
          gaze_meas_cell{end+1, 1} = gaze_opts{ig};
          subjectID_cell(end+1, 1) = sid_num;
          Trial_cell(end+1, 1) = trials_list(tr);
          Condition_cell(end+1, 1) = cval;
          col_alpha = (ie-1)*n_alpha + ia;
          if il == 1
            if ifo == 1, alpha_val_cell(end+1, 1) = alpha_fooof_050(tr, col_alpha); else, alpha_val_cell(end+1, 1) = alpha_raw_050(tr, col_alpha); end
            gv = gaze_050(tr, ig);
          elseif il == 2
            if ifo == 1, alpha_val_cell(end+1, 1) = alpha_fooof_100(tr, col_alpha); else, alpha_val_cell(end+1, 1) = alpha_raw_100(tr, col_alpha); end
            gv = gaze_100(tr, ig);
          else
            if ifo == 1, alpha_val_cell(end+1, 1) = alpha_fooof_200(tr, col_alpha); else, alpha_val_cell(end+1, 1) = alpha_raw_200(tr, col_alpha); end
            gv = gaze_200(tr, ig);
          end
          gaze_val_cell(end+1, 1) = gv;
        end
      end
    end
  end

  tbl = table(task_cell, u_id_cell, elec_cell, fooof_cell, lat_cell, alpha_type_cell, gaze_meas_cell, ...
    subjectID_cell, Trial_cell, Condition_cell, alpha_val_cell, gaze_val_cell, ...
    'VariableNames', {'task', 'universe_id', 'electrodes', 'fooof', 'latency_ms', 'alpha_type', 'gaze_measure', ...
    'subjectID', 'Trial', 'Condition', 'alpha', 'gaze_value'});
end
