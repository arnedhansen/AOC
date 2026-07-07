%% AOC Power Spectrum N Back specParam Baselined
% Plot grand average specParam baselined power spectra across N back conditions in the full window.
% Output: `AOC_powspctrm_nback_fooof_bl_full.png`.

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
path = paths.features;

fig_dir_pow = fullfile(paths.figures, 'eeg', 'powspctrm');
if ~isfolder(fig_dir_pow), mkdir(fig_dir_pow); end

%% Define channels (reference: first subject)
datapath = fullfile(path, subjects{1}, 'eeg');
cd(datapath);
nback_power_file = fullfile(datapath, 'power_nback_fooof_TFR.mat');
if ~isfile(nback_power_file)
    error('No n-back specParam power file for subject 1 (expected power_nback_fooof_TFR.mat).');
end
D0 = load(nback_power_file, 'pow1_fooof_bl_full');
pow1_fooof_bl_full = D0.pow1_fooof_bl_full;
occ_channels = {};
for i = 1:length(pow1_fooof_bl_full.label)
    label = pow1_fooof_bl_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
clc
disp('LOADING SPECPARAM+BL FULL DATA...')
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj}, 'eeg');
    cd(datapath)
    clc
    disp('LOADING SPECPARAM+BL FULL DATA...')
    fprintf('[VIZ POWSPCTRM NBACK SPECPARAM] Loading spectra for Subject %d/%d (%s)\n', subj, length(subjects), subjects{subj});
    nback_power_file = fullfile(datapath, 'power_nback_fooof_TFR.mat');
    if ~isfile(nback_power_file)
        powl1{subj} = [];
        powl2{subj} = [];
        powl3{subj} = [];
        continue
    end
    D = load(nback_power_file, 'pow1_fooof_bl_full', 'pow2_fooof_bl_full', 'pow3_fooof_bl_full');
    powl1{subj} = D.pow1_fooof_bl_full;
    powl2{subj} = D.pow2_fooof_bl_full;
    powl3{subj} = D.pow3_fooof_bl_full;
end

% Exclude subjects with any empty condition before grand-average.
is_empty1 = cellfun(@isempty, powl1);
is_empty2 = cellfun(@isempty, powl2);
is_empty3 = cellfun(@isempty, powl3);
exclude_mask = is_empty1 | is_empty2 | is_empty3;
exclude_idx = find(exclude_mask);
keep_idx = find(~exclude_mask);

fprintf('Excluded subjects (any empty condition): %d\n', numel(exclude_idx));
if ~isempty(exclude_idx)
    fprintf('Excluded subject indices: %s\n', mat2str(exclude_idx));
    fprintf('Excluded subject IDs: %s\n', strjoin(subjects(exclude_idx), ', '));
end
fprintf('Included subjects for grand average: %d\n', numel(keep_idx));

if isempty(keep_idx)
    error('No valid subjects remain after excluding empty condition data.');
end

powl1 = powl1(keep_idx);
powl2 = powl2(keep_idx);
powl3 = powl3(keep_idx);

gapow1 = ft_freqgrandaverage_nanrobust([], powl1{:});
gapow2 = ft_freqgrandaverage_nanrobust([], powl2{:});
gapow3 = ft_freqgrandaverage_nanrobust([], powl3{:});

%% Plot grand-average power spectrum
close all
figure('Position', [0 0 1512*0.4 982], 'Color', 'w');
cfg = [];
cfg.channel = channels;
cfg.figure = 'gcf';
cfg.linewidth = 3;
hold on;
yline(0, '--')

channels_seb = ismember(gapow1.label, cfg.channel);
pow1_e = gapow1.powspctrm(channels_seb, :);
pow2_e = gapow2.powspctrm(channels_seb, :);
pow3_e = gapow3.powspctrm(channels_seb, :);

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

freqs = gapow1.freq;
eb1 = shadedErrorBar(freqs, m1, se1, 'lineProps', {'-','Color',colors(1,:)});
eb2 = shadedErrorBar(freqs, m2, se2, 'lineProps', {'-','Color',colors(2,:)});
eb3 = shadedErrorBar(freqs, m3, se3, 'lineProps', {'-','Color',colors(3,:)});

set(eb1.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(1, :));
set(eb2.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(2, :));
set(eb3.mainLine, 'LineWidth', cfg.linewidth, 'Color', colors(3, :));
set(eb1.patch, 'FaceColor', colors(1, :), 'FaceAlpha', 0.20);
set(eb2.patch, 'FaceColor', colors(2, :), 'FaceAlpha', 0.20);
set(eb3.patch, 'FaceColor', colors(3, :), 'FaceAlpha', 0.20);
set(eb1.edge(1), 'Color', 'none');
set(eb1.edge(2), 'Color', 'none');
set(eb2.edge(1), 'Color', 'none');
set(eb2.edge(2), 'Color', 'none');
set(eb3.edge(1), 'Color', 'none');
set(eb3.edge(2), 'Color', 'none');

set(gca, 'Fontsize', 20);
box off
xlim([5 20]);
ylim([-.16 .16])
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
leg_p1 = patch(NaN, NaN, colors(1,:), 'EdgeColor', 'none');
leg_p2 = patch(NaN, NaN, colors(2,:), 'EdgeColor', 'none');
leg_p3 = patch(NaN, NaN, colors(3,:), 'EdgeColor', 'none');
legend([leg_p1, leg_p2, leg_p3], {' 1-back', ' 2-back', ' 3-back'}, ...
    'FontName', 'Arial', 'FontSize', 20, 'Box', 'off');
title('')

saveas(gcf, fullfile(fig_dir_pow, 'AOC_powspctrm_nback_fooof_bl_full.png'));
