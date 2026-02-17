%% AOC EEG Power Spectrum — Baseline Effects
% Control script: compares grand-average power spectra during the baseline
% period [-0.5 -0.25]s vs the retention period [0 2]s for both tasks.
%
% Loads per-subject TFR data (raw and FOOOF'd), averages over occipital
% channels and time within each window, then grand-averages across subjects.
%
% Figure layout (two figures per task — raw and FOOOF):
%   Left panel:   Raw baseline [-0.5 -0.25]s power spectra per condition
%   Middle panel: Raw retention [0 2]s power spectra per condition
%   Right panel:  Baselined retention [0 2]s power spectra per condition
%
% Data sources (per subject):
%   Raw + FOOOF TFR  → tfr_nback.mat   (tfr1/2/3, tfr1/2/3_fooof)
%                     → tfr_stern.mat   (tfr2/4/6, tfr2/4/6_fooof)
%   Baselined TFR    → tfr_nback.mat   (tfr1/2/3_bl, tfr1/2/3_fooof_bl)
%                     → tfr_stern.mat   (tfr2/4/6_bl, tfr2/4/6_fooof_bl)
%
% Key outputs:
%   AOC_eeg_powspctrm_baseline_effects_nback_raw.png
%   AOC_eeg_powspctrm_baseline_effects_nback_fooof.png
%   AOC_eeg_powspctrm_baseline_effects_sternberg_raw.png
%   AOC_eeg_powspctrm_baseline_effects_sternberg_fooof.png

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

%% Define occipital channels (from first subject)
datapath = strcat(path, subjects{1}, filesep, 'eeg');
cd(datapath);
load('tfr_nback.mat', 'tfr1');
occ_channels = {};
for i = 1:length(tfr1.label)
    lab = tfr1.label{i};
    if contains(lab, 'O') || contains(lab, 'I')
        occ_channels{end+1} = lab;
    end
end
channels = occ_channels;
clear tfr1

%% Task configurations
tasks = struct( ...
    'name',        {'nback',                                  'sternberg'}, ...
    'tfr_file',    {'tfr_nback.mat',                          'tfr_stern.mat'}, ...
    'raw_vars',    {{'tfr1','tfr2','tfr3'},                   {'tfr2','tfr4','tfr6'}}, ...
    'fooof_vars',  {{'tfr1_fooof','tfr2_fooof','tfr3_fooof'}, {'tfr2_fooof','tfr4_fooof','tfr6_fooof'}}, ...
    'raw_bl_vars', {{'tfr1_bl','tfr2_bl','tfr3_bl'},          {'tfr2_bl','tfr4_bl','tfr6_bl'}}, ...
    'fooof_bl_vars', {{'tfr1_fooof_bl','tfr2_fooof_bl','tfr3_fooof_bl'}, {'tfr2_fooof_bl','tfr4_fooof_bl','tfr6_fooof_bl'}}, ...
    'cond_labels', {{'1-back','2-back','3-back'},             {'WM load 2','WM load 4','WM load 6'}}, ...
    'title',       {'N-Back',                                 'Sternberg'});

blWin  = [-0.5 -0.25];   % Baseline window [s]
retWin = [0 2];           % Retention window [s]

% Save path
savePath = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/controls/';
if ~exist(savePath, 'dir'), mkdir(savePath); end

%% Loop over tasks
for t = 1:numel(tasks)
    task   = tasks(t);
    nConds = numel(task.raw_vars);
    nSubj  = length(subjects);
    init   = false;

    % ================================================================
    %  Collect per-subject spectra
    % ================================================================
    for subj = 1:nSubj
        clc
        fprintf('%s — Loading TFR for Subject %s (%d/%d)\n', ...
            upper(task.name), subjects{subj}, subj, nSubj);
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);

        % Load all TFR variables (raw + FOOOF + baselined)
        S = load(task.tfr_file);

        % First subject: initialise arrays
        % Raw and FOOOF may have different frequency grids
        if ~init
            refRaw = S.(task.raw_vars{1});
            refFF  = S.(task.fooof_vars{1});

            freqsRaw = refRaw.freq;
            freqsFF  = refFF.freq;
            nFreqRaw = numel(freqsRaw);
            nFreqFF  = numel(freqsFF);
            chIdx    = find(ismember(refRaw.label, channels));

            % Raw: baseline period + retention period (absolute)
            bl_raw   = nan(nSubj, nConds, nFreqRaw);
            ret_raw  = nan(nSubj, nConds, nFreqRaw);
            % Raw: retention period (baselined, dB)
            retbl_raw = nan(nSubj, nConds, nFreqRaw);

            % FOOOF: baseline period + retention period (absolute)
            bl_ff    = nan(nSubj, nConds, nFreqFF);
            ret_ff   = nan(nSubj, nConds, nFreqFF);
            % FOOOF: retention period (baselined, absolute)
            retbl_ff = nan(nSubj, nConds, nFreqFF);

            init = true;
        end

        for c = 1:nConds
            raw   = S.(task.raw_vars{c});
            ff    = S.(task.fooof_vars{c});
            rawBL = S.(task.raw_bl_vars{c});
            ffBL  = S.(task.fooof_bl_vars{c});

            % --- Raw: mean over occ channels × time window ---
            tBL_r  = raw.time >= blWin(1)  & raw.time <= blWin(2);
            tRet_r = raw.time >= retWin(1) & raw.time <= retWin(2);
            bl_raw(subj, c, :)  = mean(mean(raw.powspctrm(chIdx, :, tBL_r),  3, 'omitnan'), 1, 'omitnan');
            ret_raw(subj, c, :) = mean(mean(raw.powspctrm(chIdx, :, tRet_r), 3, 'omitnan'), 1, 'omitnan');

            % --- Raw baselined: retention only ---
            tRetBL_r = rawBL.time >= retWin(1) & rawBL.time <= retWin(2);
            retbl_raw(subj, c, :) = mean(mean(rawBL.powspctrm(chIdx, :, tRetBL_r), 3, 'omitnan'), 1, 'omitnan');

            % --- FOOOF: same (may have different freq grid) ---
            chIdxFF = find(ismember(ff.label, channels));
            tBL_f   = ff.time >= blWin(1)  & ff.time <= blWin(2);
            tRet_f  = ff.time >= retWin(1) & ff.time <= retWin(2);
            bl_ff(subj, c, :)  = mean(mean(ff.powspctrm(chIdxFF, :, tBL_f),  3, 'omitnan'), 1, 'omitnan');
            ret_ff(subj, c, :) = mean(mean(ff.powspctrm(chIdxFF, :, tRet_f), 3, 'omitnan'), 1, 'omitnan');

            % --- FOOOF baselined: retention only ---
            chIdxFFBL = find(ismember(ffBL.label, channels));
            tRetBL_f  = ffBL.time >= retWin(1) & ffBL.time <= retWin(2);
            retbl_ff(subj, c, :) = mean(mean(ffBL.powspctrm(chIdxFFBL, :, tRetBL_f), 3, 'omitnan'), 1, 'omitnan');
        end
    end

    % ================================================================
    %  Grand averages (mean ± SEM across subjects)
    % ================================================================
    nValid = sum(~isnan(bl_raw(:, 1, 1)));

    bl_raw_m    = nan(nConds, nFreqRaw);  bl_raw_s    = nan(nConds, nFreqRaw);
    ret_raw_m   = nan(nConds, nFreqRaw);  ret_raw_s   = nan(nConds, nFreqRaw);
    retbl_raw_m = nan(nConds, nFreqRaw);  retbl_raw_s = nan(nConds, nFreqRaw);
    bl_ff_m     = nan(nConds, nFreqFF);   bl_ff_s     = nan(nConds, nFreqFF);
    ret_ff_m    = nan(nConds, nFreqFF);   ret_ff_s    = nan(nConds, nFreqFF);
    retbl_ff_m  = nan(nConds, nFreqFF);   retbl_ff_s  = nan(nConds, nFreqFF);
    for c = 1:nConds
        bl_raw_m(c, :)    = mean(squeeze(bl_raw(:, c, :)),    1, 'omitnan');
        bl_raw_s(c, :)    = std(squeeze(bl_raw(:, c, :)),     0, 1, 'omitnan') / sqrt(nValid);
        ret_raw_m(c, :)   = mean(squeeze(ret_raw(:, c, :)),   1, 'omitnan');
        ret_raw_s(c, :)   = std(squeeze(ret_raw(:, c, :)),    0, 1, 'omitnan') / sqrt(nValid);
        retbl_raw_m(c, :) = mean(squeeze(retbl_raw(:, c, :)), 1, 'omitnan');
        retbl_raw_s(c, :) = std(squeeze(retbl_raw(:, c, :)),  0, 1, 'omitnan') / sqrt(nValid);
        bl_ff_m(c, :)     = mean(squeeze(bl_ff(:, c, :)),     1, 'omitnan');
        bl_ff_s(c, :)     = std(squeeze(bl_ff(:, c, :)),      0, 1, 'omitnan') / sqrt(nValid);
        ret_ff_m(c, :)    = mean(squeeze(ret_ff(:, c, :)),    1, 'omitnan');
        ret_ff_s(c, :)    = std(squeeze(ret_ff(:, c, :)),     0, 1, 'omitnan') / sqrt(nValid);
        retbl_ff_m(c, :)  = mean(squeeze(retbl_ff(:, c, :)),  1, 'omitnan');
        retbl_ff_s(c, :)  = std(squeeze(retbl_ff(:, c, :)),   0, 1, 'omitnan') / sqrt(nValid);
    end

    % ================================================================
    %  FIGURE 1 — Raw Power Spectra (3 subplots)
    % ================================================================
    close all
    figure('Position', [0, 0, 1512, 982], 'Color', 'w');

    % ---- Left: Raw Baseline [-0.5 -0.25]s ----
    subplot(1, 3, 1); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsRaw, bl_raw_m(c, :), bl_raw_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title('Baseline [-0.5 -0.25]s', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    % ---- Middle: Raw Retention [0 2]s ----
    subplot(1, 3, 2); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsRaw, ret_raw_m(c, :), ret_raw_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title('Retention [0 2]s', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    % ---- Right: Baselined Retention [0 2]s (dB) ----
    subplot(1, 3, 3); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsRaw, retbl_raw_m(c, :), retbl_raw_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    title('Retention [0 2]s (baselined)', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    sgtitle(sprintf('%s — Raw Power Spectra (N = %d)', ...
        task.title, nValid), 'FontSize', 22, 'FontWeight', 'bold');

    fname = sprintf('AOC_eeg_powspctrm_baseline_effects_%s_raw.png', task.name);
    saveas(gcf, fullfile(savePath, fname));
    fprintf('Saved: %s\n', fullfile(savePath, fname));

    % ================================================================
    %  FIGURE 2 — FOOOF Power Spectra (3 subplots)
    % ================================================================
    figure('Position', [0, 0, 1512, 982], 'Color', 'w');

    % ---- Left: FOOOF Baseline [-0.5 -0.25]s ----
    subplot(1, 3, 1); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsFF, bl_ff_m(c, :), bl_ff_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [FOOOF log]');
    title('Baseline [-0.5 -0.25]s', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    % ---- Middle: FOOOF Retention [0 2]s ----
    subplot(1, 3, 2); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsFF, ret_ff_m(c, :), ret_ff_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [FOOOF log]');
    title('Retention [0 2]s', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    % ---- Right: FOOOF Baselined Retention [0 2]s ----
    subplot(1, 3, 3); hold on;
    hLines = gobjects(nConds, 1);
    for c = 1:nConds
        eb = shadedErrorBar(freqsFF, retbl_ff_m(c, :), retbl_ff_s(c, :), ...
            'lineProps', {'-', 'Color', colors(c, :), 'LineWidth', 2});
        set(eb.patch, 'FaceColor', colors(c, :), 'FaceAlpha', 0.25);
        set(eb.edge(1), 'Color', colors(c, :));
        set(eb.edge(2), 'Color', colors(c, :));
        hLines(c) = eb.mainLine;
    end
    set(gca, 'FontSize', 15);
    xlim([3 30]);
    xlabel('Frequency [Hz]');
    ylabel('Power [FOOOF log, baselined]');
    title('Retention [0 2]s (baselined)', 'FontSize', 20);
    legend(hLines, task.cond_labels, 'Location', 'northeast', 'FontSize', 12);
    box on; hold off;

    sgtitle(sprintf('%s — FOOOF Power Spectra (N = %d)', ...
        task.title, nValid), 'FontSize', 22, 'FontWeight', 'bold');

    fname = sprintf('AOC_eeg_powspctrm_baseline_effects_%s_fooof.png', task.name);
    saveas(gcf, fullfile(savePath, fname));
    fprintf('Saved: %s\n', fullfile(savePath, fname));
end
