%% AOC Power Spectrum — Baseline Window (FOOOFed + baselined)
% Builds power spectra from the baseline interval [-1.5 -0.5] s only
% using time-resolved FOOOF+BL TFR files, for both Sternberg and N-back.

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;
baseline_window = [-1.5 -0.5];
freq_range = [3 30];

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');

% Skip when TFR feature extraction was run with RUN_FOOOF=false (no FOOOF fields in MAT).
chk = fullfile(path, subjects{1}, 'eeg', 'tfr_stern.mat');
Schk = load(chk, 'tfr2_fooof');
if ~isfield(Schk, 'tfr2_fooof')
    warning('AOC_eeg_powspctrm_baseline_window_fooof: tfr2_fooof not in %s (FOOOF disabled?). Skipping script.', chk);
    return
end

%% Channels
datapath = fullfile(path, subjects{1}, 'eeg');
cd(datapath);
D0 = load(fullfile(datapath, 'tfr_stern.mat'), 'tfr2_fooof');
ref_tfr = D0.tfr2_fooof;

occ_channels = {};
for i = 1:length(ref_tfr.label)
    label = ref_tfr.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Sternberg: load baseline-window power from tfr_stern
clc
disp('LOADING STERNBERG BASELINE-WINDOW POWERSPECTRA...')
for subj = 1:length(subjects)
    clc; fprintf('[VIZ POWSPCTRM BL - STERNBERG] Loading baseline-window FOOOF power for Subject %d / %d \n', subj, length(subjects))
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath);
    stern_file = fullfile(datapath, 'tfr_stern.mat');
    D = load(stern_file, 'tfr2_fooof', 'tfr4_fooof', 'tfr6_fooof');
    powl2{subj} = tfr_to_pow_baseline(D.tfr2_fooof, baseline_window, freq_range);
    powl4{subj} = tfr_to_pow_baseline(D.tfr4_fooof, baseline_window, freq_range);
    powl6{subj} = tfr_to_pow_baseline(D.tfr6_fooof, baseline_window, freq_range);
end

[powl2, powl4, powl6, keep_idx_stern] = exclude_empty_subjects(powl2, powl4, powl6, subjects, 'Sternberg');
% powl* are trimmed to complete subjects only; grand average over that list.
gapow2 = ft_freqgrandaverage_nanrobust([], powl2{:});
gapow4 = ft_freqgrandaverage_nanrobust([], powl4{:});
gapow6 = ft_freqgrandaverage_nanrobust([], powl6{:});

%% Sternberg plot
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.linewidth = 3;
hold on;
yline(0, '--')

elecs = ismember(gapow2.label, cfg.channel);
freqs = gapow2.freq;

pow2_e = gapow2.powspctrm(elecs, :);
pow4_e = gapow4.powspctrm(elecs, :);
pow6_e = gapow6.powspctrm(elecs, :);

m2 = mean(pow2_e, 1, 'omitnan');
n2 = sum(isfinite(pow2_e), 1);
se2 = std(pow2_e, 0, 1, 'omitnan') ./ sqrt(n2);
m4 = mean(pow4_e, 1, 'omitnan');
n4 = sum(isfinite(pow4_e), 1);
se4 = std(pow4_e, 0, 1, 'omitnan') ./ sqrt(n4);
m6 = mean(pow6_e, 1, 'omitnan');
n6 = sum(isfinite(pow6_e), 1);
se6 = std(pow6_e, 0, 1, 'omitnan') ./ sqrt(n6);

se2(n2 < 2) = NaN;
se4(n4 < 2) = NaN;
se6(n6 < 2) = NaN;

eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(1,:)});
eb4 = shadedErrorBar(freqs, m4, se4, 'lineProps', {'-','Color',colors(2,:)});
eb6 = shadedErrorBar(freqs, m6, se6, 'lineProps', {'-','Color',colors(3,:)});
ebs = [eb2, eb4, eb6];
for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(k,:));
    set(ebs(k).patch, 'FaceColor', colors(k,:), 'FaceAlpha', 0.20);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end

set(gca, 'FontSize', 20);
box off
xlim([5 20]);
ylim([0 .65]);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
leg_p2 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p6 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], {'WM load 2', 'WM load 4', 'WM load 6'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');

saveas(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_sternberg_fooof_baselineWindow.png'));

%% N-back: load baseline-window power from tfr_nback
clc
disp('LOADING N-BACK BASELINE-WINDOW POWERSPECTRA...')
for subj = 1:length(subjects)
    clc; fprintf('[VIZ POWSPCTRM BL - NBACK] Loading baseline-window FOOOF power for Subject %d / %d \n', subj, length(subjects))
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath);
    nback_file = fullfile(datapath, 'tfr_nback.mat');
    dat = load(nback_file, 'tfr1_fooof', 'tfr2_fooof', 'tfr3_fooof');
    powl1{subj} = tfr_to_pow_baseline(dat.tfr1_fooof, baseline_window, freq_range);
    powl2_nb{subj} = tfr_to_pow_baseline(dat.tfr2_fooof, baseline_window, freq_range);
    powl3{subj} = tfr_to_pow_baseline(dat.tfr3_fooof, baseline_window, freq_range);
end

[powl1, powl2_nb, powl3, keep_idx_nback] = exclude_empty_subjects(powl1, powl2_nb, powl3, subjects, 'N-back');
% powl* are trimmed to complete subjects only; grand average over that list.
gapow1 = ft_freqgrandaverage_nanrobust([], powl1{:});
gapow2_nb = ft_freqgrandaverage_nanrobust([], powl2_nb{:});
gapow3 = ft_freqgrandaverage_nanrobust([], powl3{:});

%% N-back plot
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.linewidth = 3;
hold on;
yline(0, '--')

elecs = ismember(gapow1.label, cfg.channel);
freqs = gapow1.freq;

pow1_e = gapow1.powspctrm(elecs, :);
pow2_e = gapow2_nb.powspctrm(elecs, :);
pow3_e = gapow3.powspctrm(elecs, :);

m1 = mean(pow1_e, 1, 'omitnan');
n1 = sum(isfinite(pow1_e), 1);
se1 = std(pow1_e, 0, 1, 'omitnan') ./ sqrt(n1);
m2 = mean(pow2_e, 1, 'omitnan');
n2 = sum(isfinite(pow2_e), 1);
se2 = std(pow2_e, 0, 1, 'omitnan') ./ sqrt(n2);
m3 = mean(pow3_e, 1, 'omitnan');
n3 = sum(isfinite(pow3_e), 1);
se3 = std(pow3_e, 0, 1, 'omitnan') ./ sqrt(n3);

se1(n1 < 2) = NaN;
se2(n2 < 2) = NaN;
se3(n3 < 2) = NaN;

eb1 = shadedErrorBar(freqs, m1, se1, 'lineProps', {'-','Color',colors(1,:)});
eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(2,:)});
eb3 = shadedErrorBar(freqs, m3, se3, 'lineProps', {'-','Color',colors(3,:)});
ebs = [eb1, eb2, eb3];
for k = 1:numel(ebs)
    set(ebs(k).mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(k,:));
    set(ebs(k).patch, 'FaceColor', colors(k,:), 'FaceAlpha', 0.20);
    set(ebs(k).edge(1), 'Color', 'none');
    set(ebs(k).edge(2), 'Color', 'none');
end

set(gca, 'FontSize', 20);
box off
xlim([5 20]);
ylim([0 .65]);
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p3 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1, leg_p2, leg_p3], {'1-back', '2-back', '3-back'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');

saveas(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_nback_fooof_baselineWindow.png'));

fprintf('Sternberg included N = %d subjects\n', numel(keep_idx_stern));
fprintf('N-back included N = %d subjects\n', numel(keep_idx_nback));

function S_pow = tfr_to_pow_baseline(S_tfr, baseline_window, freq_range)
cfg = [];
cfg.latency = baseline_window;
cfg.frequency = freq_range;
S_sel = ft_selectdata(cfg, S_tfr);
S_pow = S_sel;
S_pow.powspctrm = mean(S_sel.powspctrm, 3, 'omitnan');
S_pow.dimord = 'chan_freq';
if isfield(S_pow, 'time')
    S_pow = rmfield(S_pow, 'time');
end
end

function [a, b, c, keep_idx] = exclude_empty_subjects(a, b, c, subjects, label)
is_empty_a = cellfun(@isempty, a);
is_empty_b = cellfun(@isempty, b);
is_empty_c = cellfun(@isempty, c);
exclude_mask = is_empty_a | is_empty_b | is_empty_c;
exclude_idx = find(exclude_mask);
keep_idx = find(~exclude_mask);

fprintf('%s excluded subjects (any empty condition): %d\n', label, numel(exclude_idx));
if ~isempty(exclude_idx)
    fprintf('%s excluded subject indices: %s\n', label, mat2str(exclude_idx));
    fprintf('%s excluded subject IDs: %s\n', label, strjoin(subjects(exclude_idx), ', '));
end
fprintf('%s included subjects for grand average: %d\n', label, numel(keep_idx));

if isempty(keep_idx)
    error('%s: No valid subjects remain after excluding empty condition data.', label);
end

a = a(keep_idx);
b = b(keep_idx);
c = c(keep_idx);
end
