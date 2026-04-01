%% AOC Split Sternberg Alpha Loads — 500–1500 ms window (minimal: coupling + LME)
% Same logic as AOC_split_AlphaLoads.m (alpha–gaze slope coupling figure), but:
%   - Occipital alpha (FOOOF TFR, baseline [-0.5 -0.25] s) averaged over 0.5–1.5 s post-stimulus
%   - Gaze deviation recomputed from dataET_sternberg: baseline [-0.5 -0.25] s, task [0.5 1.5] s,
%     baselined in dB: 10*log10(gaze_dev_task / gaze_dev_baseline) (same convention as AlphaAmpRed)
%
% Outputs: inclusion histogram + alpha–gaze coupling figure + LME (Dev ~ Load * AlphaSlope_c).
% Figure directory: .../figures/splits/SplitAlphaLoads500_1500/

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
fig_dir = fullfile(base_data, 'figures', 'splits', 'SplitAlphaLoads500_1500');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end
cond_vals = [2 4 6];

t_analysis_eeg = [0.5 1.5]; % s post-stimulus (alpha TFR average)
t_base_gaze = [-0.5 -0.25];
t_task_gaze = [0.5 1.5];

log_dir = fullfile(base_data, 'data', 'controls', 'logs');
if ~isfolder(log_dir)
    mkdir(log_dir);
end
cmdlog_file = fullfile(log_dir, sprintf('AOC_split_AlphaLoads_500_1500ms_commandwindow_%s.log', datestr(now,'yyyymmdd_HHMMSS')));
diary('off');
diary(cmdlog_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Command window log file: %s\n', cmdlog_file);

%% Figure setup
fig_pos = [0 0 1512 982];
fontSize = 15;
winsor_cfg = struct();
winsor_cfg.enable = true;
winsor_cfg.prctile = [2 98];
split_cfg = struct();
split_cfg.use_threshold_override = true;
split_cfg.threshold = 0.015;

%% Loop over subjects - load EEG TFR (specParam, baselined)
cfg_bl = [];
cfg_bl.baseline     = [-.5 -.25];
cfg_bl.baselinetype = 'absolute';
for subj = 1:length(subjects)
    fprintf('LOADING Subject %d / %d\n', subj, length(subjects))
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
disp('LOADING FINISHED')

%% Occipital channels
occ_channels = {};
labels = load2{1, 1}.label;
for i = 1:length(labels)
    label = labels{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Grand average
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

%% Alpha power (0.5–1.5 s)
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

%% Subject-level slopes across load
nSubj = length(subjects);
slope = nan(nSubj, 1);
for subj = 1:nSubj
    y = [alpha2(subj), alpha4(subj), alpha6(subj)];
    Xmat = [ones(3,1), [2;4;6]];
    b = Xmat \ y';
    slope(subj,1) = b(2);
end

%% Split
nSubj = numel(slope);
if split_cfg.use_threshold_override
    thr = split_cfg.threshold;
    if ~isscalar(thr) || ~isfinite(thr) || thr < 0
        error('split_cfg.threshold must be a finite, non-negative scalar.');
    end
    idx_inc = slope > thr;
    idx_dec = slope < -thr;
    idx_flat = ~(idx_inc | idx_dec);
else
    k = floor(nSubj/3);
    if k < 1
        error('Not enough subjects for a zero-centered tail split.');
    end
    pos_idx = find(slope > 0);
    [~, ord_pos] = sort(slope(pos_idx), 'descend');
    pos_sel = pos_idx(ord_pos(1:min(k, numel(pos_idx))));
    neg_idx = find(slope < 0);
    [~, ord_neg] = sort(slope(neg_idx), 'ascend');
    neg_sel = neg_idx(ord_neg(1:min(k, numel(neg_idx))));
    k_eff = min(numel(pos_sel), numel(neg_sel));
    if k_eff < k
        warning('Only %d participants per tail possible (requested %d).', k_eff, k);
    end
    pos_sel = pos_sel(1:k_eff);
    neg_sel = neg_sel(1:k_eff);
    idx_inc = false(nSubj,1);
    idx_dec = false(nSubj,1);
    idx_flat = true(nSubj,1);
    idx_inc(pos_sel) = true;
    idx_dec(neg_sel) = true;
    idx_flat(idx_inc | idx_dec) = false;
end

n_inc = sum(idx_inc);
n_dec = sum(idx_dec);
n_f = sum(idx_flat);

if any(idx_dec)
    t1 = max(slope(idx_dec));
else
    t1 = NaN;
end
if any(idx_inc)
    t2 = min(slope(idx_inc));
else
    t2 = NaN;
end

%% Inclusion figure
close all
figure('Position', fig_pos, 'Color', 'w');
hold on
% bin width from tail cutoffs
bin_w = (t2 - t1) / 8;
if ~isfinite(bin_w) || bin_w <= 0
    all_s = slope(isfinite(slope));
    if isempty(all_s)
        bin_w = 1;
    else
        bin_w = (max(all_s) - min(all_s)) / 12;
        if ~isfinite(bin_w) || bin_w <= 0
            bin_w = 1e-3;
        end
    end
end
min_s = min(slope);
max_s = max(slope);
k_left = ceil((t1 - min_s) / bin_w);
k_right = ceil((max_s - t1) / bin_w);
bin_edges = t1 + (-k_left:k_right) * bin_w;
histogram(slope(idx_inc), 'BinEdges', bin_edges, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6);
histogram(slope(idx_dec),  'BinEdges', bin_edges, 'FaceColor', [0 0 0.8], 'FaceAlpha', 0.6);
histogram(slope(idx_flat),   'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);
xline(t1, 'k--', 'LineWidth', 2);
xline(t2, 'k--', 'LineWidth', 2);
xlabel('Alpha power slope')
ylabel('Participants')
title(sprintf('Alpha power slope (0.5–1.5 s), WM load 2/4/6'), 'FontSize', 18)
legend({sprintf('Increase (N=%d)', n_inc), sprintf('Decrease (N=%d)', n_dec), ...
    sprintf('intermediate (N=%d)', n_f), 'cutoffs'}, 'Box', 'off')
box on
set(gca, 'FontSize', 15)
drawnow
saveas(gcf, fullfile(fig_dir, 'AOC_split_AlphaLoads_500_1500ms_inclusion.png'));

fprintf('\nSplit summary (alpha TFR 0.5–1.5 s):\n');
fprintf('Increase: %d, Decrease: %d, Intermediate: %d\n', n_inc, n_dec, n_f);

%% Gaze: mean dB deviation per load (recomputed; not merged table)
fprintf('\nComputing gaze deviation (dB) per load from dataET_sternberg...\n');
[DEV2, DEV4, DEV6] = sternberg_gaze_deviation_dB_by_load(subjects, path, t_base_gaze, t_task_gaze);

%% Alpha–gaze coupling (subject-level)
disp(upper('Alpha–gaze coupling over WM load (subject-level models)...'))

gaze_dev_slope = nan(nSubj, 1);
for ii = 1:nSubj
    ydev = [DEV2(ii), DEV4(ii), DEV6(ii)];
    if all(isfinite(ydev))
        Xmat = [ones(3, 1), [2; 4; 6]];
        b = Xmat \ ydev(:);
        gaze_dev_slope(ii) = b(2);
    end
end

mask_couple = isfinite(slope(:)) & isfinite(gaze_dev_slope(:));
mask_couple_fig = mask_couple & (idx_inc | idx_dec);
n_couple = sum(mask_couple);
n_couple_fig = sum(mask_couple_fig);
mI = mask_couple_fig & idx_inc;
mD = mask_couple_fig & idx_dec;
n_inc_fig = sum(mI);
n_dec_fig = sum(mD);
r_p_inc = NaN; p_p_inc = NaN; r_s_inc = NaN; p_s_inc = NaN; p_pearson_inv_inc = NaN;
r_p_dec = NaN; p_p_dec = NaN; r_s_dec = NaN; p_s_dec = NaN; p_pearson_inv_dec = NaN;

fprintf('\n========== SUBJECT-LEVEL ALPHA–GAZE COUPLING (WM LOAD 2–6) ==========\n');
fprintf(['Gaze: 10*log10(dev_task/dev_baseline), task [%.2f %.2f] s, baseline [%.2f %.2f] s; ', ...
    'slope = OLS of mean gaze dB on load.\n'], t_task_gaze(1), t_task_gaze(2), t_base_gaze(1), t_base_gaze(2));
fprintf('Alpha: occipital 8–14 Hz, TFR baseline [-0.5 -0.25] s, average [%.2f %.2f] s.\n', ...
    t_analysis_eeg(1), t_analysis_eeg(2));

mask_gaze_fin = isfinite(gaze_dev_slope);
n_gaze_fin = sum(mask_gaze_fin);
fprintf('\n--- Descriptive statistics ---\n');
fprintf('N with finite gaze deviation slope: %d\n', n_gaze_fin);
if n_gaze_fin > 0
    gs_all = gaze_dev_slope(mask_gaze_fin);
    fprintf('  Gaze dev slope (dB per item): min=%.5g, max=%.5g, mean=%.5g, SD=%.5g, median=%.5g\n', ...
        min(gs_all), max(gs_all), mean(gs_all), std(gs_all), median(gs_all));
end
fprintf('N with finite alpha and gaze slopes (any): %d; increase/decrease (figure): %d\n', ...
    n_couple, n_couple_fig);
if n_couple_fig > 0
    sa_d = slope(mask_couple_fig);
    sg_d = gaze_dev_slope(mask_couple_fig);
    fprintf('  Alpha slope (coupling N): min=%.5g, max=%.5g, mean=%.5g, SD=%.5g\n', ...
        min(sa_d), max(sa_d), mean(sa_d), std(sa_d));
    fprintf('  Gaze dev slope (coupling N): min=%.5g, max=%.5g, mean=%.5g, SD=%.5g\n', ...
        min(sg_d), max(sg_d), mean(sg_d), std(sg_d));
end

fprintf('\n--- Within-group correlation (alpha slope vs gaze dev slope), per group ---\n');
if n_inc_fig >= 3
    sa_i = slope(mI); sg_i = gaze_dev_slope(mI);
    [r_p_inc, p_p_inc] = corr(sa_i, sg_i, 'Type', 'Pearson', 'rows', 'complete');
    [r_s_inc, p_s_inc] = corr(sa_i, sg_i, 'Type', 'Spearman', 'rows', 'complete');
    df_i = n_inc_fig - 2;
    if df_i > 0 && abs(r_p_inc) < 1 - eps
        t_ri = r_p_inc * sqrt(df_i / max(eps, 1 - r_p_inc^2));
        p_pearson_inv_inc = tcdf(t_ri, df_i);
    end
    fprintf('Alpha increase (N=%d): Pearson r=%.5f, p=%.6g; H1 inv. coupling (rho<0): p=%.6g; Spearman rho=%.5f, p=%.6g\n', ...
        n_inc_fig, r_p_inc, p_p_inc, p_pearson_inv_inc, r_s_inc, p_s_inc);
else
    fprintf('Alpha increase: N=%d (skipped correlation; need N>=3)\n', n_inc_fig);
end
if n_dec_fig >= 3
    sa_dv = slope(mD); sg_dv = gaze_dev_slope(mD);
    [r_p_dec, p_p_dec] = corr(sa_dv, sg_dv, 'Type', 'Pearson', 'rows', 'complete');
    [r_s_dec, p_s_dec] = corr(sa_dv, sg_dv, 'Type', 'Spearman', 'rows', 'complete');
    df_d = n_dec_fig - 2;
    if df_d > 0 && abs(r_p_dec) < 1 - eps
        t_rd = r_p_dec * sqrt(df_d / max(eps, 1 - r_p_dec^2));
        p_pearson_inv_dec = tcdf(t_rd, df_d);
    end
    fprintf('Alpha decrease (N=%d): Pearson r=%.5f, p=%.6g; H1 inv. coupling (rho<0): p=%.6g; Spearman rho=%.5f, p=%.6g\n', ...
        n_dec_fig, r_p_dec, p_p_dec, p_pearson_inv_dec, r_s_dec, p_s_dec);
else
    fprintf('Alpha decrease: N=%d (skipped correlation; need N>=3)\n', n_dec_fig);
end
if n_inc_fig >= 2
    p_ols_inc = polyfit(slope(mI), gaze_dev_slope(mI), 1);
    fprintf('OLS, alpha increase: gaze_dev_slope = %.5g + %.5g * alpha_slope\n', p_ols_inc(2), p_ols_inc(1));
end
if n_dec_fig >= 2
    p_ols_dec = polyfit(slope(mD), gaze_dev_slope(mD), 1);
    fprintf('OLS, alpha decrease: gaze_dev_slope = %.5g + %.5g * alpha_slope\n', p_ols_dec(2), p_ols_dec(1));
end

col_inc = [0.8 0 0];
col_dec = [0 0 0.8];
col_flat = [0.5 0.5 0.5];
mu_gs = nan(1, 2);
sem_gs = nan(1, 2);
n_gs = zeros(1, 2);
for gi = 1:2
    if gi == 1
        gmask = mask_couple_fig & idx_inc;
    else
        gmask = mask_couple_fig & idx_dec;
    end
    vals = gaze_dev_slope(gmask);
    vals = vals(isfinite(vals));
    n_gs(gi) = numel(vals);
    if n_gs(gi) >= 1
        mu_gs(gi) = mean(vals);
        if n_gs(gi) >= 2
            sem_gs(gi) = std(vals) / sqrt(n_gs(gi));
        else
            sem_gs(gi) = 0;
        end
    end
end

fprintf('\n--- Building alpha–gaze coupling figures ---\n');
figure('Position', fig_pos, 'Color', 'w');
axS = subplot(2, 2, [1 2]);
hold(axS, 'on');
h1 = scatter(axS, NaN, NaN, 40, col_inc, 'filled', 'MarkerEdgeColor', col_inc * 0.35);
h2 = scatter(axS, NaN, NaN, 40, col_dec, 'filled', 'MarkerEdgeColor', col_dec * 0.35);
idx_plot = find(mask_couple_fig);
for ii = idx_plot(:)'
    if idx_inc(ii)
        scatter(axS, slope(ii), gaze_dev_slope(ii), 36, col_inc, 'filled', 'MarkerFaceAlpha', 0.78, ...
            'MarkerEdgeColor', col_inc * 0.35, 'LineWidth', 0.45);
    elseif idx_dec(ii)
        scatter(axS, slope(ii), gaze_dev_slope(ii), 36, col_dec, 'filled', 'MarkerFaceAlpha', 0.78, ...
            'MarkerEdgeColor', col_dec * 0.35, 'LineWidth', 0.45);
    end
end
if sum(mI) >= 2
    sa_i = slope(mI);
    p_inc = polyfit(sa_i, gaze_dev_slope(mI), 1);
    xl_i = [min(sa_i), max(sa_i)];
    pad_i = 0.05 * max(abs(diff(xl_i)), eps);
    xl_i = [xl_i(1) - pad_i, xl_i(2) + pad_i];
    plot(axS, xl_i, polyval(p_inc, xl_i), '-', 'Color', col_inc * 0.55, 'LineWidth', 1.65);
end
if sum(mD) >= 2
    sa_d = slope(mD);
    p_dec = polyfit(sa_d, gaze_dev_slope(mD), 1);
    xl_d = [min(sa_d), max(sa_d)];
    pad_d = 0.05 * max(abs(diff(xl_d)), eps);
    xl_d = [xl_d(1) - pad_d, xl_d(2) + pad_d];
    plot(axS, xl_d, polyval(p_dec, xl_d), '-', 'Color', col_dec * 0.55, 'LineWidth', 1.65);
end
xline(axS, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
yline(axS, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
hold(axS, 'off');
sa_sc = slope(mask_couple_fig);
sg_sc = gaze_dev_slope(mask_couple_fig);
sx_sc = max(abs(sa_sc(isfinite(sa_sc))));
sy_sc = max(abs(sg_sc(isfinite(sg_sc))));
if ~isfinite(sx_sc) || sx_sc <= 0, sx_sc = eps; end
if ~isfinite(sy_sc) || sy_sc <= 0, sy_sc = eps; end
xlim(axS, [-1 1] * sx_sc * 1.05);
ylim(axS, [-1 1] * sy_sc * 1.05);
xlabel(axS, 'Alpha power slope (a.u. per item)', 'FontSize', fontSize - 1);
ylabel(axS, 'Gaze deviation slope (dB per item)', 'FontSize', fontSize - 1);
title(axS, ' ', 'FontSize', fontSize, 'Interpreter', 'tex');
set(axS, 'FontSize', fontSize - 2);
legend(axS, [h1, h2], {'Alpha increase', 'Alpha decrease'}, ...
    'Location', 'best', 'Box', 'off', 'Interpreter', 'tex', 'FontSize', fontSize - 3);
txt = {};
if n_inc_fig >= 3
    txt{end+1} = sprintf(['Alpha increase (N=%d): Pearson \\it r\\rm = %.3f, \\it p\\rm = %.4g; ', ...
        'H_1: \\rho < 0: \\it p\\rm = %.4g; Spearman \\rho = %.3f, \\it p\\rm = %.4g'], ...
        n_inc_fig, r_p_inc, p_p_inc, p_pearson_inv_inc, r_s_inc, p_s_inc);
else
    txt{end+1} = sprintf('Alpha increase: N = %d (correlation requires N \\geq 3)', n_inc_fig);
end
if n_dec_fig >= 3
    txt{end+1} = sprintf(['Alpha decrease (N=%d): Pearson \\it r\\rm = %.3f, \\it p\\rm = %.4g; ', ...
        'H_1: \\rho < 0: \\it p\\rm = %.4g; Spearman \\rho = %.3f, \\it p\\rm = %.4g'], ...
        n_dec_fig, r_p_dec, p_p_dec, p_pearson_inv_dec, r_s_dec, p_s_dec);
else
    txt{end+1} = sprintf('Alpha decrease: N = %d (correlation requires N \\geq 3)', n_dec_fig);
end
text(axS, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', fontSize - 4, 'Interpreter', 'tex');

fprintf('\n--- Gaze deviation slope by group (dB per item; same sample as scatter) ---\n');
fprintf('  group                n    mean(gaze slope)       SEM\n');
tail_lbl = {'alpha increase', 'alpha decrease'};
for gi = 1:2
    if n_gs(gi) < 1
        fprintf('  %-20s   --          --              --\n', tail_lbl{gi});
    else
        fprintf('  %-20s  %3d  %18.5g  %12.5g\n', tail_lbl{gi}, n_gs(gi), mu_gs(gi), sem_gs(gi));
    end
end

axAlpha = subplot(2, 2, 3);
hold(axAlpha, 'on');
if numel(bin_edges) >= 2 && all(isfinite(bin_edges(:)))
    histogram(axAlpha, slope(idx_inc), 'BinEdges', bin_edges, 'FaceColor', col_inc, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_dec), 'BinEdges', bin_edges, 'FaceColor', col_dec, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_flat), 'BinEdges', bin_edges, 'FaceColor', col_flat, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    lim_a = max(abs([bin_edges(1), bin_edges(end)]));
else
    histogram(axAlpha, slope(idx_inc), 'FaceColor', col_inc, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_dec), 'FaceColor', col_dec, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    histogram(axAlpha, slope(idx_flat), 'FaceColor', col_flat, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    all_s = slope(isfinite(slope));
    lim_a = max(abs([min(all_s), max(all_s)]));
end
xline(axAlpha, 0, 'k-', 'LineWidth', 1);
if isfinite(t1), xline(axAlpha, t1, 'k--', 'LineWidth', 1.5); end
if isfinite(t2), xline(axAlpha, t2, 'k--', 'LineWidth', 1.5); end
hold(axAlpha, 'off');
if isfinite(lim_a) && lim_a > 0
    xlim(axAlpha, [-lim_a lim_a]);
end
xlabel(axAlpha, 'Alpha power slope', 'FontSize', fontSize - 1);
ylabel(axAlpha, 'Participants', 'FontSize', fontSize - 1);
title(axAlpha, 'Alpha power slope across WM load', 'FontSize', fontSize - 1, 'Interpreter', 'none');
set(axAlpha, 'FontSize', fontSize - 2);

axGaze = subplot(2, 2, 4);
idx_g = find(mask_couple_fig);
gv_raw = gaze_dev_slope(idx_g);
[~, ord] = sort(gv_raw);
idx_ord = idx_g(ord);
gv_sorted = gv_raw(ord);
n_b = numel(idx_ord);
C_b = zeros(n_b, 3);
for k = 1:n_b
    ii = idx_ord(k);
    if idx_inc(ii)
        C_b(k, :) = col_inc;
    else
        C_b(k, :) = col_dec;
    end
end
if n_b >= 1
    h_bar = ones(n_b, 1);
    b = bar(axGaze, gv_sorted, h_bar, 'BarWidth', 5, 'EdgeColor', 'none');
    b.FaceColor = 'flat';
    b.CData = C_b;
    xline(axGaze, 0, 'k-', 'LineWidth', 1);
    mg = max(abs(gv_sorted(isfinite(gv_sorted))));
    if isfinite(mg) && mg > 0
        xlim(axGaze, [-mg mg] * 1.05);
    end
    ylim(axGaze, [0 1.12]);
    set(axGaze, 'YTick', []);
    hold(axGaze, 'on');
    plot(axGaze, NaN, NaN, 's', 'MarkerFaceColor', col_inc, 'MarkerEdgeColor', col_inc * 0.35, 'MarkerSize', 9);
    plot(axGaze, NaN, NaN, 's', 'MarkerFaceColor', col_dec, 'MarkerEdgeColor', col_dec * 0.35, 'MarkerSize', 9);
    legend(axGaze, {'Alpha increase', 'Alpha decrease'}, 'Location', 'northeast', ...
        'Box', 'off', 'FontSize', fontSize - 4, 'Interpreter', 'none');
    hold(axGaze, 'off');
end
xlabel(axGaze, 'Gaze deviation slope (dB per item)', 'FontSize', fontSize - 1);
ylabel(axGaze, '', 'FontSize', fontSize - 1);
title(axGaze, 'Gaze deviation slopes', 'FontSize', fontSize - 1, 'Interpreter', 'none');
set(axGaze, 'FontSize', fontSize - 2);

sgtitle({'Subject-level Alpha–gaze coupling'; ...
    'EEG: 0.5–1.5 s; Gaze: 10·log_{10}(task/baseline), task 0.5–1.5 s'}, ...
    'FontSize', fontSize + 1, 'FontWeight', 'bold', 'Interpreter', 'tex');
drawnow;
fig_coupling_path = fullfile(fig_dir, 'AOC_split_AlphaLoads_500_1500ms_alpha_gaze_slope_coupling.png');
saveas(gcf, fig_coupling_path);
fprintf('Saved figure: %s\n', fig_coupling_path);

%% LME
Subject_ag = [];
Load_ag = [];
Dev_ag = [];
Alpha_ag = [];
idx_ag = find(all(isfinite([DEV2, DEV4, DEV6]), 2));
mu_alpha = mean(slope(idx_ag));
for k = 1:numel(idx_ag)
    ii = idx_ag(k);
    dev_vals = [DEV2(ii), DEV4(ii), DEV6(ii)];
    for c = 1:3
        Subject_ag = [Subject_ag; ii];
        Load_ag = [Load_ag; cond_vals(c)];
        Dev_ag = [Dev_ag; dev_vals(c)];
        Alpha_ag = [Alpha_ag; slope(ii) - mu_alpha];
    end
end
tbl_ag = table(Subject_ag, Load_ag, Dev_ag, Alpha_ag, ...
    'VariableNames', {'Subject', 'Load', 'Dev', 'AlphaSlope_c'});
tbl_ag.Subject = categorical(tbl_ag.Subject);
tbl_ag.Load = (tbl_ag.Load - 4) / 2;

disp('--- LME: Dev [dB] ~ Load * AlphaSlope (between-subject) + (Load|Subject) ---')
fprintf('Load coded as (Load-4)/2 (0 = load 4; +1 = load 6).\n');
fprintf('Dev = mean trial-level gaze deviation in dB (10*log10(task/baseline)) per load.\n');
try
    n_subj_ag = numel(categories(categorical(tbl_ag.Subject)));
    n_obs_ag = height(tbl_ag);
    fprintf('LME data: %d subjects with complete gaze at loads 2/4/6, %d rows total.\n', ...
        n_subj_ag, n_obs_ag);
    fprintf('Grand mean alpha slope (subtracted as AlphaSlope_c): %.6g\n', mu_alpha);

    [lme_ag, lme_ag_formula] = fitlme_with_random_slope_fallback(tbl_ag, 'Dev ~ Load * AlphaSlope_c');
    fprintf('Fitted formula: %s\n', lme_ag_formula);
    fprintf('\nANOVA (marginal tests):\n');
    disp(anova(lme_ag))

    fprintf('Fixed-effects coefficients:\n');
    disp(lme_ag.Coefficients)

    coef_ag = lme_ag.Coefficients;
    idx_int = contains(coef_ag.Name, ':') & (contains(coef_ag.Name, 'Load') & contains(coef_ag.Name, 'AlphaSlope'));
    ix_int = find(idx_int, 1, 'first');
    if ~isempty(ix_int)
        row = coef_ag(ix_int, :);
        fprintf('Primary test (Load:AlphaSlope_c): inverse coupling if estimate < 0.\n');
        fprintf('  Estimate=%.6g, SE=%.6g, t=%.4g, DF=%.4g, p=%.6g\n', ...
            row.Estimate(1), row.SE(1), row.tStat(1), row.DF(1), row.pValue(1));
    else
        fprintf('Note: Load:AlphaSlope interaction row not found in coefficient table (check names).\n');
    end
    try
        fprintf('\n95%% CI for fixed effects (coefCI):\n');
        disp(coefCI(lme_ag));
    catch %#ok<CTCH>
    end
    fprintf('\n========== END SUBJECT-LEVEL ALPHA–GAZE COUPLING OUTPUT ==========\n\n');
end

fprintf('Done. Figures in: %s\n', fig_dir);

%% --- Local functions ---

function freq_out = winsorize_freq_subjects(freq_in, prct_bounds)
freq_out = freq_in;
if ~isfield(freq_in, 'powspctrm') || isempty(freq_in.powspctrm)
    return
end

P = freq_in.powspctrm;
if size(P, 1) < 3
    % Too few subjects for stable percentile clipping.
    return
end

P_size = size(P);
P_2d = reshape(P, P_size(1), []);
lo = prctile(P_2d, prct_bounds(1), 1);
hi = prctile(P_2d, prct_bounds(2), 1);
P_2d = min(max(P_2d, lo), hi);
freq_out.powspctrm = reshape(P_2d, P_size);
end

function [DEV2, DEV4, DEV6] = sternberg_gaze_deviation_dB_by_load(subjects, base_path, t_base, t_task)
% Per subject: mean across trials of 10*log10(gd_task/gd_baseline) for loads 2/4/6.
% Preprocessing aligned with AOC_gaze_fex_sternberg_trials.m (screen mask, Y flip, blinks).

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
    clc
    fprintf('Computing gaze deviation for Subject %d', si)
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
