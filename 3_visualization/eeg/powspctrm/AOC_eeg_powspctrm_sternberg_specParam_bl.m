%% AOC Power Spectrum — Sternberg (specParam + baseline, late window [1 2]s)
% Loads pow*_fooof_bl_late from power_stern_fooof_TFR.mat (AOC_eeg_fex_sternberg_TFR.m),
% grand-averages across WM load 2/4/6, plots spectra, and saves figures.
%
% For non-specParam Hanning spectra use AOC_eeg_powspctrm_sternberg.m (raw) or
% AOC_eeg_powspctrm_sternberg_bl.m (dB baselined), both from power_stern_windows.mat.

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end

%% Define channels (reference: first subject)
datapath = fullfile(path, subjects{1}, 'eeg');
cd(datapath);
ref_file = fullfile(datapath, 'power_stern_fooof_TFR.mat');
tmp = load(ref_file, 'pow2_fooof_bl_late');
pow2_fooof_bl_late = tmp.pow2_fooof_bl_late;

occ_channels = {};
for i = 1:length(pow2_fooof_bl_late.label)
    label = pow2_fooof_bl_late.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
clc
disp('LOADING SPECPARAM+BL LATE DATA...')
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath)
    clc
    disp('LOADING SPECPARAM+BL LATE DATA...')
    disp(subj)
    stern_power_file = fullfile(datapath, 'power_stern_fooof_TFR.mat');
    D = load(stern_power_file, 'pow2_fooof_bl_late', 'pow4_fooof_bl_late', 'pow6_fooof_bl_late');
    powl2{subj} = D.pow2_fooof_bl_late;
    powl4{subj} = D.pow4_fooof_bl_late;
    powl6{subj} = D.pow6_fooof_bl_late;
end

% Grand-average per condition
is_empty2 = cellfun(@isempty, powl2);
is_empty4 = cellfun(@isempty, powl4);
is_empty6 = cellfun(@isempty, powl6);

miss2 = find(is_empty2);
miss4 = find(is_empty4);
miss6 = find(is_empty6);
if ~isempty(miss2)
    fprintf('WM load 2 missing (%d): %s\n', numel(miss2), strjoin(subjects(miss2), ', '));
end
if ~isempty(miss4)
    fprintf('WM load 4 missing (%d): %s\n', numel(miss4), strjoin(subjects(miss4), ', '));
end
if ~isempty(miss6)
    fprintf('WM load 6 missing (%d): %s\n', numel(miss6), strjoin(subjects(miss6), ', '));
end

powl2_valid = powl2(~is_empty2);
powl4_valid = powl4(~is_empty4);
powl6_valid = powl6(~is_empty6);

fprintf('Grand average N: WM load 2 = %d, load 4 = %d, load 6 = %d\n', ...
    numel(powl2_valid), numel(powl4_valid), numel(powl6_valid));

if isempty(powl2_valid) || isempty(powl4_valid) || isempty(powl6_valid)
    error('At least one WM load has no subject data; cannot plot all three conditions.');
end

gapow2 = ft_freqgrandaverage_nanrobust([], powl2_valid{:});
gapow4 = ft_freqgrandaverage_nanrobust([], powl4_valid{:});
gapow6 = ft_freqgrandaverage_nanrobust([], powl6_valid{:});

%% Plot grand-average power spectrum
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel   = channels;
cfg.figure    = 'gcf';
cfg.linewidth = 3;
hold on;
yline(0, '--')

elecs = ismember(gapow2.label, cfg.channel);
freqs = gapow2.freq;

pow2_e = gapow2.powspctrm(elecs,:);
pow4_e = gapow4.powspctrm(elecs,:);
pow6_e = gapow6.powspctrm(elecs,:);

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
ylim([-.16 .16])
ylabel('Power [dB]');
xlabel('Frequency [Hz]');
leg_p2 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p4 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p6 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p2, leg_p4, leg_p6], {' WM load 2', ' WM load 4', ' WM load 6'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('');

saveas(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_sternberg_fooof_bl_late.png'));
