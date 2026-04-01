%% Reciprocal subject×load LMEs: Gaze ~ Load*Alpha and Alpha ~ Load*Gaze
%
% Long-format rows: one per subject × WM load (2/4/6). Alpha = occipital 8–14 Hz
% TFR power [dB] (FOOOF-baselined TFR), averaged 0.5–1.5 s. Gaze = mean trial-level
% gaze deviation [dB]: 10*log10(task/baseline), baseline [-0.5 -0.25] s, task [0.5 1.5] s.
%
% Load is coded as (Load-4)/2 so 0 = load 4, +1 = load 6 (same as AOC_split_AlphaLoads_500_1500ms).
% In each model the *moderator* is grand-mean centered across all rows (Alpha_c or Gaze_c) so
% the interaction is interpretable; DVs stay on the raw scale (Gaze, Alpha).
%
% Models (random slope for Load by subject, with intercept-only fallback):
%   (1) Gaze ~ Load * Alpha_c + (Load|Subject)
%   (2) Alpha ~ Load * Gaze_c + (Load|Subject)
%
% These are associational within-subject×load data; they do not identify causal direction.

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;

if ispc
    base_data = 'W:\Students\Arne\AOC';
else
    base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
end
feat_dir = fullfile(base_data, 'data', 'features');
fig_dir = fullfile(base_data, 'figures', 'splits', 'ReciprocalAlphaGazeLME');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

cond_vals = [2 4 6];
t_analysis_eeg = [0.5 1.5];
t_base_gaze = [-0.5 -0.25];
t_task_gaze = [0.5 1.5];

winsor_cfg = struct('enable', true, 'prctile', [2 98]);

log_dir = fullfile(base_data, 'data', 'controls', 'logs');
if ~isfolder(log_dir)
    mkdir(log_dir);
end
cmdlog_file = fullfile(log_dir, sprintf('AOC_alpha_gaze_reciprocal_LME_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
diary('off');
diary(cmdlog_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Command window log file: %s\n', cmdlog_file);

fprintf('\n=== Reciprocal LMEs: Gaze ~ Load*Alpha_c, Alpha ~ Load*Gaze_c ===\n');
fprintf('Figure/log directory: %s\n', fig_dir);

%% Load TFR per subject
cfg_bl = [];
cfg_bl.baseline = [-.5 -.25];
cfg_bl.baselinetype = 'absolute';
load2 = cell(1, numel(subjects));
load4 = cell(1, numel(subjects));
load6 = cell(1, numel(subjects));

for subj = 1:numel(subjects)
    fprintf('LOADING TFR Subject %d / %d\n', subj, numel(subjects));
    eeg_dir = fullfile(path, subjects{subj}, 'eeg');
    f = fullfile(eeg_dir, 'tfr_stern.mat');
    if ~isfile(f)
        error('Missing: %s', f);
    end
    datTFR = load(f);
    if isfield(datTFR, 'tfr2_fooof_bl')
        load2{subj} = datTFR.tfr2_fooof_bl;
        load4{subj} = datTFR.tfr4_fooof_bl;
        load6{subj} = datTFR.tfr6_fooof_bl;
    else
        load2{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr2_fooof);
        load4{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr4_fooof);
        load6{subj} = ft_freqbaseline(cfg_bl, datTFR.tfr6_fooof);
    end
end
disp('TFR LOADING FINISHED');

%% Occipital channels
occ_channels = {};
labels = load2{1}.label;
for i = 1:numel(labels)
    label = labels{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label; %#ok<AGROW>
    end
end
channels = occ_channels;

%% Grand average + winsorize
cfg = [];
cfg.keepindividual = 'yes';
ga2 = ft_freqgrandaverage(cfg, load2{:});
ga4 = ft_freqgrandaverage(cfg, load4{:});
ga6 = ft_freqgrandaverage(cfg, load6{:});
if winsor_cfg.enable
    ga2 = winsorize_freq_subjects(ga2, winsor_cfg.prctile);
    ga4 = winsorize_freq_subjects(ga4, winsor_cfg.prctile);
    ga6 = winsorize_freq_subjects(ga6, winsor_cfg.prctile);
end

cfg = [];
cfg.frequency = [8 14];
cfg.latency = t_analysis_eeg;
cfg.channel = channels;
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
val2 = ft_selectdata(cfg, ga2);
val4 = ft_selectdata(cfg, ga4);
val6 = ft_selectdata(cfg, ga6);
alpha2 = val2.powspctrm;
alpha4 = val4.powspctrm;
alpha6 = val6.powspctrm;

nSubj = numel(subjects);

%% Gaze dB per load (ET)
fprintf('\nComputing gaze deviation [dB] per load from dataET_sternberg...\n');
[DEV2, DEV4, DEV6] = sternberg_gaze_deviation_dB_by_load(subjects, path, t_base_gaze, t_task_gaze);

%% Long table: subject × load
Subject_col = [];
Load_raw = [];
Gaze_col = [];
Alpha_col = [];

idx_ok = find(all(isfinite([alpha2(:), alpha4(:), alpha6(:), DEV2(:), DEV4(:), DEV6(:)]), 2));
fprintf('\nSubjects with finite alpha (2/4/6) and gaze dB (2/4/6): %d / %d\n', numel(idx_ok), nSubj);

for k = 1:numel(idx_ok)
    ii = idx_ok(k);
    av = [alpha2(ii), alpha4(ii), alpha6(ii)];
    gv = [DEV2(ii), DEV4(ii), DEV6(ii)];
    for c = 1:3
        Subject_col(end+1, 1) = ii; %#ok<AGROW>
        Load_raw(end+1, 1) = cond_vals(c); %#ok<AGROW>
        Alpha_col(end+1, 1) = av(c); %#ok<AGROW>
        Gaze_col(end+1, 1) = gv(c); %#ok<AGROW>
    end
end

if numel(Subject_col) < 15
    error('Too few observations for LME (need at least ~15 rows).');
end

mu_alpha_grand = mean(Alpha_col);
mu_gaze_grand = mean(Gaze_col);
Alpha_c = Alpha_col - mu_alpha_grand;
Gaze_c = Gaze_col - mu_gaze_grand;

Load_ctr = (Load_raw - 4) / 2;

tbl = table(Subject_col, Load_ctr, Gaze_col, Alpha_col, Alpha_c, Gaze_c, ...
    'VariableNames', {'Subject', 'Load', 'Gaze', 'Alpha', 'Alpha_c', 'Gaze_c'});
tbl.Subject = categorical(tbl.Subject);

fprintf('\n--- Data summary ---\n');
fprintf('Rows (subject×load): %d; subjects: %d\n', height(tbl), numel(unique(Subject_col)));
fprintf('Grand means: Alpha=%.6g [dB], Gaze=%.6g [dB] (subtracted from moderator in each model).\n', ...
    mu_alpha_grand, mu_gaze_grand);
fprintf('Load coding: Load = (WM_load - 4) / 2  (0 = load 4; +1 = load 6; -1 = load 2).\n');
fprintf('EEG window [%.2f %.2f] s; gaze task [%.2f %.2f] s; gaze baseline [%.2f %.2f] s.\n', ...
    t_analysis_eeg(1), t_analysis_eeg(2), t_task_gaze(1), t_task_gaze(2), t_base_gaze(1), t_base_gaze(2));

%% Model 1: Gaze ~ Load * Alpha_c + (Load|Subject)
fprintf('\n========== MODEL 1: Gaze ~ Load * Alpha_c + (Load|Subject) ==========\n');
fprintf('DV: Gaze [dB] at each load. Moderator: Alpha_c (Alpha grand-mean centered).\n');
try
    [lme1, f1] = fitlme_with_random_slope_fallback(tbl, 'Gaze ~ Load * Alpha_c');
    fprintf('Fitted: %s\n', f1);
    disp(anova(lme1));
    fprintf('Fixed effects:\n');
    disp(lme1.Coefficients);
    try
        fprintf('\n95%% CI (fixed):\n');
        disp(coefCI(lme1));
    catch %#ok<CTCH>
    end
catch ME
    fprintf('Model 1 failed: %s\n', ME.message);
end

%% Model 2: Alpha ~ Load * Gaze_c + (Load|Subject)
fprintf('\n========== MODEL 2: Alpha ~ Load * Gaze_c + (Load|Subject) ==========\n');
fprintf('DV: Alpha [dB] at each load. Moderator: Gaze_c (Gaze grand-mean centered).\n');
try
    [lme2, f2] = fitlme_with_random_slope_fallback(tbl, 'Alpha ~ Load * Gaze_c');
    fprintf('Fitted: %s\n', f2);
    disp(anova(lme2));
    fprintf('Fixed effects:\n');
    disp(lme2.Coefficients);
    try
        fprintf('\n95%% CI (fixed):\n');
        disp(coefCI(lme2));
    catch %#ok<CTCH>
    end
catch ME
    fprintf('Model 2 failed: %s\n', ME.message);
end

fprintf('\n=== Done. Associational models; interpret Load:moderator as moderation at fixed grand mean of the other variable. ===\n');
fprintf('Log: %s\n', cmdlog_file);

%% --- Local functions ---

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm)
    return
end
P = freq_in.powspctrm;
if size(P, 1) < 3
    return
end
P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function [lme, used_formula] = fitlme_with_random_slope_fallback(tbl_in, fixed_formula)
formula_rs = sprintf('%s + (Load|Subject)', fixed_formula);
formula_ri = sprintf('%s + (1|Subject)', fixed_formula);
try
    lme = fitlme(tbl_in, formula_rs);
    used_formula = formula_rs;
catch
    lme = fitlme(tbl_in, formula_ri);
    used_formula = formula_ri;
end
end

function [DEV2, DEV4, DEV6] = sternberg_gaze_deviation_dB_by_load(subjects, base_path, t_base, t_task)
% Per subject: mean across trials of 10*log10(gd_task/gd_baseline) for loads 2/4/6.
screenW = 800; screenH = 600;
centreX = 400; centreY = 300;
blink_win = 50;
min_valid_samples = 100;
bounds_x = [0 screenW];
bounds_y = [0 screenH];

nSubj = numel(subjects);
DEV2 = nan(nSubj, 1);
DEV4 = nan(nSubj, 1);
DEV6 = nan(nSubj, 1);

for si = 1:nSubj
    fprintf('Gaze ET Subject %d / %d\n', si, nSubj);
    et_file = fullfile(base_path, subjects{si}, 'gaze', 'dataET_sternberg.mat');
    if ~isfile(et_file)
        et_file = fullfile(base_path, subjects{si}, 'gaze', 'dataET_sternberg');
    end
    if ~isfile(et_file)
        warning('Missing ET file for subject %s', subjects{si});
        continue
    end
    tmp = load(et_file);
    dataETlong = select_dataetlong(tmp, subjects{si});
    nTrials = size(dataETlong.trialinfo, 1);
    cond_code = dataETlong.trialinfo(:, 1) - 20;

    gd_by_cond = cell(3, 1);
    for c = 1:3
        gd_by_cond{c} = [];
    end

    for trl = 1:nTrials
        raw_dat = dataETlong.trial{trl};
        t = dataETlong.time{trl};
        if isempty(raw_dat) || isempty(t)
            continue
        end
        raw_dat = raw_dat(1:3, :);
        raw_dat(2, :) = screenH - raw_dat(2, :);

        inb = raw_dat(1, :) >= bounds_x(1) & raw_dat(1, :) <= bounds_x(2) & ...
            raw_dat(2, :) >= bounds_y(1) & raw_dat(2, :) <= bounds_y(2);
        raw_dat(:, ~inb) = NaN;

        raw_dat = remove_blinks(raw_dat, blink_win);

        idx_b = t >= t_base(1) & t <= t_base(2);
        idx_t = t >= t_task(1) & t <= t_task(2);
        dat_b = raw_dat(:, idx_b);
        dat_t = raw_dat(:, idx_t);

        ok_b = sum(all(isfinite(dat_b(1:2, :)), 1)) >= min_valid_samples;
        ok_t = sum(all(isfinite(dat_t(1:2, :)), 1)) >= min_valid_samples;
        if ~ok_b || ~ok_t
            continue
        end

        dx_b = dat_b(1, :) - centreX;
        dy_b = dat_b(2, :) - centreY;
        gaze_dev_b = nanmean(sqrt(dx_b.^2 + dy_b.^2));

        dx_t = dat_t(1, :) - centreX;
        dy_t = dat_t(2, :) - centreY;
        gaze_dev_t = nanmean(sqrt(dx_t.^2 + dy_t.^2));

        if ~(isfinite(gaze_dev_b) && gaze_dev_b > 0 && isfinite(gaze_dev_t) && gaze_dev_t > 0)
            continue
        end
        gd_dB = 10 * log10(gaze_dev_t / gaze_dev_b);
        if ~isfinite(gd_dB)
            continue
        end

        cc = cond_code(trl);
        if cc == 2
            gd_by_cond{1}(end+1, 1) = gd_dB;
        elseif cc == 4
            gd_by_cond{2}(end+1, 1) = gd_dB;
        elseif cc == 6
            gd_by_cond{3}(end+1, 1) = gd_dB;
        end
    end

    if ~isempty(gd_by_cond{1})
        DEV2(si) = mean(gd_by_cond{1}, 'omitnan');
    end
    if ~isempty(gd_by_cond{2})
        DEV4(si) = mean(gd_by_cond{2}, 'omitnan');
    end
    if ~isempty(gd_by_cond{3})
        DEV6(si) = mean(gd_by_cond{3}, 'omitnan');
    end
end
end

function dataETlong = select_dataetlong(tmp, subj_label)
if isfield(tmp, 'dataETlong')
    dataETlong = tmp.dataETlong;
elseif isfield(tmp, 'dataet')
    dataETlong = tmp.dataet;
elseif isfield(tmp, 'dataET')
    dataETlong = tmp.dataET;
else
    error('No ET data structure found in dataET_sternberg for subject %s.', subj_label);
end
if ~isfield(dataETlong, 'trial') || ~isfield(dataETlong, 'trialinfo') || ~isfield(dataETlong, 'time')
    error('ET structure for subject %s is missing trial/trialinfo/time fields.', subj_label);
end
end
