%% AOC Interaction — Horizontal Gaze Bias Index × Posterior Alpha Lateralization
% Exploratory analysis testing the relationship between horizontal gaze
% position bias and posterior alpha lateralization at the trial level.
%
% Horizontal Gaze Bias Index (HGBI):
%   HGBI = (time_left − time_right) / (time_left + time_right)
%   where time_left  = number of gaze samples left of screen centre (x < 400)
%         time_right = number of gaze samples right of screen centre (x > 400)
%   Range: −1 (all left) to +1 (all right); 0 = symmetric.
%
% Posterior Alpha Lateralization (Stroganova et al., 2007):
%   ALI = (R_alpha − L_alpha) / (R_alpha + L_alpha)
%   on occipital channels (O, I), odd numbers = left, even = right.
%   Pre-computed per trial in lateralization_*_trials.mat or per condition
%   in lateralization_*.mat.
%
% Analyses (per task: Sternberg, N-back, resting):
%   1. Trial-level HGBI computation from L-GAZE-X during full retention [0 2] s
%   2. Within-subject Spearman correlations (Fisher-z, group t-test)
%   3. GLMM: HGBI ~ AlphaLatZ * Condition + (1 | ID)   [tasks only]
%   4. Pooled analysis (collapsed across conditions) + condition interaction
%
% Visualizations:
%   - Scatter plot (HGBI vs alpha lateralization) coloured by condition
%   - Topographic maps: alpha power split by gaze direction (left vs right)
%   - Time-resolved HGBI and alpha lateralization (per-trial sliding window)
%   - Raincloud of within-subject correlations
%
% Prerequisites:
%   - dataET_sternberg.mat / dataET_nback.mat / dataET_resting.mat
%   - dataEEG_TFR_sternberg.mat / dataEEG_TFR_nback.mat / dataEEG_resting.mat
%   - lateralization_sternberg.mat / lateralization_nback.mat (condition-level)
%   - power_*_full_trials.mat (for topographic split)

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');

set(0, 'DefaultFigureColor', 'w')
fontSize = 26;

fs       = 500;           % EEG & ET sampling rate (Hz)
anaWin   = [0 2];         % retention / analysis window (s)
screenCX = 400;           % screen centre X (pixels, 800×600)

% Task definitions
tasks        = {'sternberg', 'nback', 'resting'};
taskETFiles  = {'dataET_sternberg', 'dataET_nback', 'dataET_resting'};
taskEEGFiles = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback', 'dataEEG_resting'};
taskIAFFiles = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat', ''};

% Alpha band (fallback; overridden by IAF when available)
alphaRange = [8 14];

% Output directory
outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions/HGBIxAlphaLat';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% Define posterior ROI (occipital: O, I channels) and left/right split
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat', 'powload2_late');
allLabs = powload2_late.label;

% Occipital channels
occ_channels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I')
        occ_channels{end+1} = L; %#ok<SAGROW>
    end
end

% Split into left (odd numbers) and right (even numbers)
left_channels  = {};
right_channels = {};
for i = 1:numel(occ_channels)
    ch = occ_channels{i};
    numStr = regexp(ch, '\d+', 'match');
    if ~isempty(numStr)
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch; %#ok<SAGROW>
        else
            right_channels{end+1} = ch; %#ok<SAGROW>
        end
    end
end

% Full posterior ROI (for topographies)
roi_labels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I') || startsWith(L, 'P')
        roi_labels{end+1} = L; %#ok<SAGROW>
    end
end

fprintf('Occipital L (%d): %s\n', numel(left_channels), strjoin(left_channels, ', '));
fprintf('Occipital R (%d): %s\n', numel(right_channels), strjoin(right_channels, ', '));

%% Ridge stabilisation for lateralization denominator
epsP      = 1e-12;
ridgeFrac = 0.01;

%% Results container
results = struct();

%% ====================================================================
%%                         TASK LOOP
%% ====================================================================
for taskIdx = 1:numel(tasks)
    taskName = tasks{taskIdx};
    isResting = strcmp(taskName, 'resting');
    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName), repmat('=',1,60));

    %% Holders — trial-by-trial
    all_hgbi       = [];
    all_alphalat   = [];
    all_subj       = [];
    all_cond       = [];

    %% Holders — time-resolved (per-subject means)
    timeRes_hgbi     = {};
    timeRes_alphalat = {};
    timeRes_time     = {};

    %% Holders — topographic split
    topo_gazeLeft_allSubj  = [];
    topo_gazeRight_allSubj = [];
    topo_template          = [];   % FieldTrip freq struct template for topoplot

    %% ================================================================
    %%                       SUBJECT LOOP
    %% ================================================================
    for s = 1:numel(subjects)

        fprintf('Task: %s | Subject %s (%d/%d)\n', taskName, subjects{s}, s, numel(subjects));

        dpeeg  = fullfile(path, subjects{s}, 'eeg');
        dpgaze = fullfile(path, subjects{s}, 'gaze');

        %% Load gaze data
        etFile = fullfile(dpgaze, [taskETFiles{taskIdx} '.mat']);
        if ~exist(etFile, 'file')
            fprintf('  Missing ET data, skipping.\n');
            continue;
        end
        cd(dpgaze);

        if isResting
            try
                load(taskETFiles{taskIdx}, 'dataET');
            catch
                fprintf('  Load error (ET resting), skipping.\n'); continue;
            end
        else
            try
                load(taskETFiles{taskIdx}, 'dataETlong');
            catch
                fprintf('  Load error (ET), skipping.\n'); continue;
            end
        end

        %% Load EEG data
        eegFile = fullfile(dpeeg, [taskEEGFiles{taskIdx} '.mat']);
        if ~exist(eegFile, 'file')
            fprintf('  Missing EEG data, skipping.\n');
            continue;
        end
        cd(dpeeg);

        if isResting
            try
                load(taskEEGFiles{taskIdx}, 'dataEEG', 'segInfo');
            catch
                fprintf('  Load error (EEG resting), skipping.\n'); continue;
            end
        else
            try
                load(taskEEGFiles{taskIdx}, 'dataTFR');
            catch
                fprintf('  Load error (EEG), skipping.\n'); continue;
            end
        end

        %% Load IAF
        IAF = NaN;
        if ~isResting
            iafFile = fullfile(dpeeg, taskIAFFiles{taskIdx});
            if exist(iafFile, 'file')
                tmp = load(iafFile, 'IAF_subj');
                if isfinite(tmp.IAF_subj)
                    IAF = tmp.IAF_subj;
                end
            end
        else
            % For resting: try Sternberg IAF, then N-back IAF
            for ic = {'IAF_sternberg_subject.mat', 'IAF_nback_subject.mat'}
                iafPath = fullfile(dpeeg, ic{1});
                if exist(iafPath, 'file')
                    tmp = load(iafPath, 'IAF_subj');
                    if isfinite(tmp.IAF_subj)
                        IAF = tmp.IAF_subj;
                        break;
                    end
                end
            end
        end

        if isfinite(IAF)
            bandAlpha = [IAF-4, IAF+2];
        else
            bandAlpha = alphaRange;
        end

        %% ============================================================
        %%  TASK: Sternberg / N-back (trial-level)
        %% ============================================================
        if ~isResting
            % Select posterior ROI and compute alpha lateralization per trial
            nT = min(numel(dataTFR.trial), size(dataETlong.trialinfo, 1));

            % Per-trial storage
            subj_hgbi     = nan(nT, 1);
            subj_alphalat = nan(nT, 1);
            subj_cond     = nan(nT, 1);

            % For topographic split: accumulate power spectra
            topo_gazeLeft  = [];
            topo_gazeRight = [];

            % Frequency analysis for alpha lateralization (trial-level)
            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.output     = 'pow';
            cfg.taper      = 'hanning';
            cfg.foi        = 1:0.5:30;
            cfg.keeptrials = 'yes';

            % Select retention window
            cfgsel         = [];
            cfgsel.latency = anaWin;
            dataSel        = ft_selectdata(cfgsel, dataTFR);

            % Compute trial-level power spectra
            freq_trials = ft_freqanalysis(cfg, dataSel);

            % Get channel indices for lateralization
            left_idx  = find(ismember(freq_trials.label, left_channels));
            right_idx = find(ismember(freq_trials.label, right_channels));
            alpha_fidx = freq_trials.freq >= bandAlpha(1) & freq_trials.freq <= bandAlpha(2);

            for tr = 1:nT
                %% HGBI from gaze
                t_gaze = dataETlong.time{tr};
                idx_ana = t_gaze >= anaWin(1) & t_gaze <= anaWin(2);

                X = double(dataETlong.trial{tr}(1, idx_ana));  % L-GAZE-X

                % Blink removal (use area channel)
                if size(dataETlong.trial{tr}, 1) >= 3
                    A = double(dataETlong.trial{tr}(3, idx_ana));
                    blink = ~isfinite(A) | (A <= 0);
                    X(blink) = NaN;
                end

                % Remove NaN samples for HGBI
                X_valid = X(isfinite(X));
                if numel(X_valid) < 50  % minimum viable samples
                    continue;
                end

                time_left  = sum(X_valid < screenCX);
                time_right = sum(X_valid > screenCX);
                total      = time_left + time_right;

                if total > 0
                    subj_hgbi(tr) = (time_left - time_right) / total;
                end

                %% Alpha lateralization from power spectra
                if tr <= size(freq_trials.powspctrm, 1)
                    L_pow = mean(freq_trials.powspctrm(tr, left_idx, alpha_fidx), [2 3]);
                    R_pow = mean(freq_trials.powspctrm(tr, right_idx, alpha_fidx), [2 3]);
                    denom = R_pow + L_pow + ridgeFrac * mean([R_pow, L_pow]) + epsP;
                    subj_alphalat(tr) = (R_pow - L_pow) / denom;
                end

                %% Condition code
                subj_cond(tr) = dataETlong.trialinfo(tr, 1);

                %% Topographic split by gaze direction
                if isfinite(subj_hgbi(tr)) && tr <= size(freq_trials.powspctrm, 1)
                    trialPow = squeeze(freq_trials.powspctrm(tr, :, alpha_fidx));
                    trialPow_alpha = mean(trialPow, 2);  % [chan x 1]
                    if subj_hgbi(tr) < 0  % gaze left
                        topo_gazeLeft = [topo_gazeLeft, trialPow_alpha]; %#ok<AGROW>
                    elseif subj_hgbi(tr) > 0  % gaze right
                        topo_gazeRight = [topo_gazeRight, trialPow_alpha]; %#ok<AGROW>
                    end
                end
            end % trial loop

            %% Time-resolved HGBI and alpha lateralization (sliding window)
            winDur  = 0.5;  % 500 ms window
            winStep = 0.1;  % 100 ms step
            winSamp = round(winDur * fs);
            stepSamp = round(winStep * fs);

            % Use the first valid trial to get time vector length
            t_example = dataETlong.time{1};
            winCenters = (anaWin(1) + winDur/2):winStep:(anaWin(2) - winDur/2);
            nWins = numel(winCenters);

            hgbi_timeSeries = nan(nT, nWins);
            for tr = 1:nT
                t_gaze = dataETlong.time{tr};
                X = double(dataETlong.trial{tr}(1, :));
                if size(dataETlong.trial{tr}, 1) >= 3
                    A = double(dataETlong.trial{tr}(3, :));
                    blink = ~isfinite(A) | (A <= 0);
                    X(blink) = NaN;
                end

                for wi = 1:nWins
                    wStart = winCenters(wi) - winDur/2;
                    wEnd   = winCenters(wi) + winDur/2;
                    idx_w  = t_gaze >= wStart & t_gaze < wEnd;
                    Xw = X(idx_w);
                    Xw_valid = Xw(isfinite(Xw));
                    if numel(Xw_valid) < 20, continue; end
                    tL = sum(Xw_valid < screenCX);
                    tR = sum(Xw_valid > screenCX);
                    if (tL + tR) > 0
                        hgbi_timeSeries(tr, wi) = (tL - tR) / (tL + tR);
                    end
                end
            end

            % Store time-resolved (subject mean)
            timeRes_hgbi{s}     = nanmean(hgbi_timeSeries, 1);
            timeRes_time{s}     = winCenters;

            %% Store topographic data (subject-mean)
            if ~isempty(topo_gazeLeft)
                topo_gazeLeft_allSubj = cat(3, topo_gazeLeft_allSubj, mean(topo_gazeLeft, 2));
            end
            if ~isempty(topo_gazeRight)
                topo_gazeRight_allSubj = cat(3, topo_gazeRight_allSubj, mean(topo_gazeRight, 2));
            end

            %% Append subject data
            valid = isfinite(subj_hgbi) & isfinite(subj_alphalat);
            all_hgbi     = [all_hgbi;     subj_hgbi(valid)]; %#ok<AGROW>
            all_alphalat = [all_alphalat;  subj_alphalat(valid)]; %#ok<AGROW>
            all_subj     = [all_subj;      repmat(str2double(subjects{s}), sum(valid), 1)]; %#ok<AGROW>
            all_cond     = [all_cond;      subj_cond(valid)]; %#ok<AGROW>

            % Save template for topographic plotting (first subject)
            if isempty(topo_template) && exist('freq_trials', 'var')
                topo_template = freq_trials;
                topo_template = rmfield(topo_template, 'powspctrm');
                topo_template.dimord = 'chan_freq';
                topo_template.freq   = mean(bandAlpha);
            end

            clear dataTFR dataETlong freq_trials dataSel

        %% ============================================================
        %%  RESTING STATE (sliding-window approach)
        %% ============================================================
        else
            % Resting state: continuous data, no trials
            % Use sliding windows for HGBI and alpha lateralization

            try
                segRange = segInfo.fixcross;  % analyse fixation-cross period
            catch
                fprintf('  Missing segInfo, skipping.\n');
                clear dataEEG dataET; continue;
            end

            % ET channels
            etLabels = dataET.label;
            xIdx = find(strcmp(etLabels, 'L-GAZE-X'), 1);
            aIdx = find(strcmp(etLabels, 'L-AREA'), 1);
            if isempty(xIdx)
                fprintf('  L-GAZE-X not found, skipping.\n');
                clear dataEEG dataET; continue;
            end

            gazeX = double(dataET.trial{1}(xIdx, :));
            if ~isempty(aIdx)
                gazeA = double(dataET.trial{1}(aIdx, :));
            else
                gazeA = ones(1, numel(gazeX));
            end

            % Blink removal
            blinkMask = ~isfinite(gazeA) | (gazeA <= 0) | ~isfinite(gazeX);
            gazeX(blinkMask) = NaN;

            % EEG: extract posterior channels
            post_idx = ismember(dataEEG.label, occ_channels);
            left_eeg_idx  = ismember(dataEEG.label, left_channels);
            right_eeg_idx = ismember(dataEEG.label, right_channels);

            if sum(post_idx) == 0
                fprintf('  Posterior channels not found, skipping.\n');
                clear dataEEG dataET; continue;
            end

            eeg_data = dataEEG.trial{1};
            t_eeg = dataEEG.time{1};

            % Sliding window
            winDur_rest   = 2;    % 2-second windows
            winStep_rest  = 1;    % 1-second steps
            winSamp_rest  = round(winDur_rest * fs);
            stepSamp_rest = round(winStep_rest * fs);

            % Buffer from segment edges
            bufferSamp = 2 * fs;
            s1 = segRange(1) + bufferSamp;
            s2 = segRange(2) - bufferSamp;
            if s1 >= s2
                fprintf('  Segment too short, skipping.\n');
                clear dataEEG dataET; continue;
            end

            winStarts = s1:stepSamp_rest:(s2 - winSamp_rest + 1);
            nWin = numel(winStarts);

            if nWin < 10
                fprintf('  Too few windows (%d), skipping.\n', nWin);
                clear dataEEG dataET; continue;
            end

            win_hgbi     = nan(nWin, 1);
            win_alphalat = nan(nWin, 1);

            for w = 1:nWin
                idx = winStarts(w):(winStarts(w) + winSamp_rest - 1);

                % HGBI
                Xw = gazeX(idx);
                Xw_valid = Xw(isfinite(Xw));
                if numel(Xw_valid) >= 50
                    tL = sum(Xw_valid < screenCX);
                    tR = sum(Xw_valid > screenCX);
                    if (tL + tR) > 0
                        win_hgbi(w) = (tL - tR) / (tL + tR);
                    end
                end

                % Alpha lateralization from EEG (FFT on window)
                eeg_left  = mean(eeg_data(left_eeg_idx, idx), 1);
                eeg_right = mean(eeg_data(right_eeg_idx, idx), 1);

                if any(~isfinite(eeg_left)) || any(~isfinite(eeg_right))
                    continue;
                end

                % Power via FFT
                nfft = 2^nextpow2(numel(idx));
                f = (0:nfft/2) * (fs / nfft);
                alpha_fidx = f >= bandAlpha(1) & f <= bandAlpha(2);

                pL = abs(fft(eeg_left - mean(eeg_left), nfft)).^2 / numel(idx);
                pR = abs(fft(eeg_right - mean(eeg_right), nfft)).^2 / numel(idx);
                pL = pL(1:nfft/2+1);
                pR = pR(1:nfft/2+1);

                L_alpha = mean(pL(alpha_fidx));
                R_alpha = mean(pR(alpha_fidx));

                denom = R_alpha + L_alpha + ridgeFrac * mean([R_alpha, L_alpha]) + epsP;
                win_alphalat(w) = (R_alpha - L_alpha) / denom;
            end

            % Treat each window as a "trial"
            valid = isfinite(win_hgbi) & isfinite(win_alphalat);
            all_hgbi     = [all_hgbi;     win_hgbi(valid)]; %#ok<AGROW>
            all_alphalat = [all_alphalat;  win_alphalat(valid)]; %#ok<AGROW>
            all_subj     = [all_subj;      repmat(str2double(subjects{s}), sum(valid), 1)]; %#ok<AGROW>
            all_cond     = [all_cond;      zeros(sum(valid), 1)]; %#ok<AGROW>  % 0 = resting

            clear dataEEG dataET gazeX gazeA eeg_data

        end % task vs resting

    end % subject loop

    %% ================================================================
    %%    GROUP STATISTICS
    %% ================================================================
    fprintf('\n--- Group Statistics: %s ---\n', upper(taskName));

    if isempty(all_hgbi)
        fprintf('  No valid data. Skipping.\n');
        continue;
    end

    % Convert condition codes for tasks
    if ~isResting
        all_cond = all_cond - 20;  % Sternberg: 22/24/26→2/4/6; N-back: 21/22/23→1/2/3
    end

    %% Within-subject Spearman correlations (pooled across conditions)
    unique_ids = unique(all_subj);
    nS = numel(unique_ids);
    r_ws = nan(nS, 1);

    for si = 1:nS
        idx = all_subj == unique_ids(si);
        h = all_hgbi(idx);
        a = all_alphalat(idx);
        if sum(isfinite(h) & isfinite(a)) >= 5
            r_ws(si) = corr(h, a, 'Type', 'Spearman', 'Rows', 'complete');
        end
    end

    % Group-level: Fisher-z → one-sample t-test vs 0
    valid_r = isfinite(r_ws);
    z_ws = atanh(r_ws(valid_r));
    if numel(z_ws) >= 3
        [~, p_group, ~, stats_group] = ttest(z_ws);
        mean_r = tanh(mean(z_ws));
        se_z   = std(z_ws) / sqrt(numel(z_ws));
        ci_r   = tanh([mean(z_ws) - 1.96*se_z, mean(z_ws) + 1.96*se_z]);

        fprintf('  Within-subject r: mean = %.3f, 95%% CI = [%.3f, %.3f]\n', mean_r, ci_r(1), ci_r(2));
        fprintf('  t(%d) = %.2f, p = %.4f\n', stats_group.df, stats_group.tstat, p_group);
    else
        p_group = NaN; mean_r = NaN; ci_r = [NaN NaN];
        stats_group = struct('df', NaN, 'tstat', NaN);
        fprintf('  Too few subjects with valid correlations.\n');
    end

    %% GLMM: HGBI ~ AlphaLatZ * Condition + (1 | ID)
    if ~isResting && numel(unique(all_cond)) > 1
        T = table(all_hgbi, all_alphalat, categorical(all_cond), categorical(all_subj), ...
            'VariableNames', {'HGBI', 'AlphaLat', 'Condition', 'ID'});

        % Z-score alpha lateralization within subject
        T.AlphaLatZ = nan(height(T), 1);
        for si = 1:nS
            idx = T.ID == categorical(unique_ids(si));
            vals = T.AlphaLat(idx);
            mu = mean(vals, 'omitnan');
            sd = std(vals, 'omitnan');
            if sd > 0
                T.AlphaLatZ(idx) = (vals - mu) / sd;
            else
                T.AlphaLatZ(idx) = 0;
            end
        end

        try
            lme_full    = fitlme(T, 'HGBI ~ AlphaLatZ * Condition + (1 | ID)');
            lme_reduced = fitlme(T, 'HGBI ~ AlphaLatZ + Condition + (1 | ID)');
            lrt = compare(lme_reduced, lme_full);

            fprintf('\n  GLMM: HGBI ~ AlphaLatZ * Condition + (1|ID)\n');
            fprintf('  LRT for interaction: Chi2 = %.2f, p = %.4f\n', lrt.LRStat(2), lrt.pValue(2));

            if lrt.pValue(2) < 0.05
                fprintf('  → Significant interaction; reporting full model\n');
                lme_final = lme_full;
            else
                fprintf('  → No significant interaction; reporting reduced model\n');
                lme_final = lme_reduced;
            end
            disp(lme_final)
        catch ME
            fprintf('  GLMM failed: %s\n', ME.message);
            lme_final = [];
        end
    else
        % Resting or single condition: simple GLMM without condition
        T = table(all_hgbi, all_alphalat, categorical(all_subj), ...
            'VariableNames', {'HGBI', 'AlphaLat', 'ID'});

        T.AlphaLatZ = nan(height(T), 1);
        for si = 1:nS
            idx = T.ID == categorical(unique_ids(si));
            vals = T.AlphaLat(idx);
            mu = mean(vals, 'omitnan');
            sd = std(vals, 'omitnan');
            if sd > 0
                T.AlphaLatZ(idx) = (vals - mu) / sd;
            else
                T.AlphaLatZ(idx) = 0;
            end
        end

        try
            lme_final = fitlme(T, 'HGBI ~ AlphaLatZ + (1 | ID)');
            fprintf('\n  GLMM: HGBI ~ AlphaLatZ + (1|ID)\n');
            disp(lme_final)
        catch ME
            fprintf('  GLMM failed: %s\n', ME.message);
            lme_final = [];
        end
    end

    %% ================================================================
    %%    FIGURE 1: Scatter Plot (HGBI vs Alpha Lateralization by Condition)
    %% ================================================================
    close all
    figure('Position', [0 0 1000 800]); hold on

    % Remove between-subject means for display
    hgbi_ws     = all_hgbi;
    alphalat_ws = all_alphalat;
    for si = 1:nS
        idx = all_subj == unique_ids(si);
        hgbi_ws(idx)     = all_hgbi(idx)     - mean(all_hgbi(idx),     'omitnan');
        alphalat_ws(idx) = all_alphalat(idx) - mean(all_alphalat(idx), 'omitnan');
    end

    unique_conds = sort(unique(all_cond));
    nCond = numel(unique_conds);
    h_sc = gobjects(nCond, 1);

    for ci = 1:nCond
        idx = all_cond == unique_conds(ci);
        cIdx = min(ci, size(colors, 1));
        h_sc(ci) = scatter(alphalat_ws(idx), hgbi_ws(idx), 20, colors(cIdx,:), ...
            'filled', 'MarkerFaceAlpha', 0.3);
    end

    % Regression line
    vld = isfinite(alphalat_ws) & isfinite(hgbi_ws);
    p_fit = polyfit(alphalat_ws(vld), hgbi_ws(vld), 1);
    xl = xlim; xfit = linspace(xl(1), xl(2), 100);
    h_reg = plot(xfit, polyval(p_fit, xfit), 'k-', 'LineWidth', 3);

    xlabel('Alpha Lateralization (within-subject centred)')
    ylabel('HGBI (within-subject centred)')

    if strcmp(taskName, 'sternberg')
        condLabels = arrayfun(@(x) sprintf('Load %d', x), unique_conds, 'UniformOutput', false);
    elseif strcmp(taskName, 'nback')
        condLabels = arrayfun(@(x) sprintf('%d-back', x), unique_conds, 'UniformOutput', false);
    else
        condLabels = {'Resting'};
    end
    legend([h_sc; h_reg], [condLabels; {'Regression'}], 'Location', 'best', 'FontSize', fontSize - 8)

    title(sprintf('%s: HGBI vs Alpha Lateralization\nr = %.3f [%.3f, %.3f], p = %.4f', ...
        upper(taskName), mean_r, ci_r(1), ci_r(2), p_group), 'FontSize', fontSize - 4)
    set(gca, 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_HGBIxAlphaLat_scatter_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 2: Raincloud of Within-Subject Correlations
    %% ================================================================
    figure('Position', [0 0 700 700]); hold on
    vals = r_ws(valid_r);

    % Jittered points
    xj = 0.85 + 0.3 * rand(size(vals));
    scatter(xj, vals, 60, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.6);

    % Mean + 95% CI
    if ~isnan(mean_r)
        plot([1 1], ci_r, 'k-', 'LineWidth', 3)
        plot(1, mean_r, 's', 'MarkerSize', 14, 'MarkerFaceColor', colors(1,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2)
    end
    yline(0, 'k:', 'LineWidth', 1.5)

    xlim([0.5 1.5])
    ylabel('Within-subject Spearman r')
    title(sprintf('%s: Within-Subject Correlations\nHGBI \\times \\alpha Lateralization\nr = %.3f, p = %.4f', ...
        upper(taskName), mean_r, p_group), 'FontSize', fontSize - 6)
    set(gca, 'XTick', [], 'FontSize', fontSize - 4)
    saveas(gcf, fullfile(outdir, ['AOC_HGBIxAlphaLat_raincloud_' taskName '.png']));

    %% ================================================================
    %%    FIGURE 3: Topographic Maps — Alpha Power by Gaze Direction
    %% ================================================================
    if ~isResting && ~isempty(topo_gazeLeft_allSubj) && ~isempty(topo_gazeRight_allSubj) && ~isempty(topo_template)
        % Grand-average across subjects
        topo_L = mean(topo_gazeLeft_allSubj, 3);
        topo_R = mean(topo_gazeRight_allSubj, 3);
        topo_diff = topo_R - topo_L;

        % Use saved FieldTrip freq struct template for topoplot
        topo_struct = topo_template;

        figure('Position', [0 0 1500 500]);

        % Panel 1: Gaze Left trials
        subplot(1,3,1)
        topo_struct.powspctrm = topo_L;
        cfg_topo = [];
        cfg_topo.xlim     = mean(bandAlpha) * [1 1];
        cfg_topo.zlim     = 'maxabs';
        cfg_topo.layout   = 'EEG1005.lay';
        cfg_topo.colorbar = 'yes';
        cfg_topo.comment  = 'no';
        cfg_topo.style    = 'straight';
        ft_topoplotER(cfg_topo, topo_struct);
        title('Gaze Left Trials', 'FontSize', fontSize - 4)

        % Panel 2: Gaze Right trials
        subplot(1,3,2)
        topo_struct.powspctrm = topo_R;
        ft_topoplotER(cfg_topo, topo_struct);
        title('Gaze Right Trials', 'FontSize', fontSize - 4)

        % Panel 3: Difference (Right - Left)
        subplot(1,3,3)
        topo_struct.powspctrm = topo_diff;
        cfg_topo.zlim = 'maxabs';
        ft_topoplotER(cfg_topo, topo_struct);
        title('Difference (R - L gaze)', 'FontSize', fontSize - 4)

        sgtitle(sprintf('%s: Alpha Power Topography by Gaze Direction', upper(taskName)), 'FontSize', fontSize)
        saveas(gcf, fullfile(outdir, ['AOC_HGBIxAlphaLat_topo_' taskName '.png']));
    end

    %% ================================================================
    %%    FIGURE 4: Time-Resolved HGBI
    %% ================================================================
    if ~isResting
        % Collect all subject time series (same time axis)
        validTS = ~cellfun(@isempty, timeRes_hgbi);
        if sum(validTS) >= 3
            % Stack into matrix
            hgbi_mat = cell2mat(timeRes_hgbi(validTS)');
            tVec = timeRes_time{find(validTS, 1, 'first')};

            ga_hgbi  = mean(hgbi_mat, 1, 'omitnan');
            sem_hgbi = std(hgbi_mat, [], 1, 'omitnan') / sqrt(sum(validTS));

            figure('Position', [0 0 1200 600]); hold on

            fill([tVec, fliplr(tVec)], [ga_hgbi - sem_hgbi, fliplr(ga_hgbi + sem_hgbi)], ...
                colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            plot(tVec, ga_hgbi, '-', 'LineWidth', 3, 'Color', colors(1,:))
            yline(0, 'k:', 'LineWidth', 1.5)
            xlabel('Time (s)')
            ylabel('HGBI')
            title(sprintf('%s: Time-Resolved Horizontal Gaze Bias Index', upper(taskName)), 'FontSize', fontSize)
            set(gca, 'FontSize', fontSize - 4)
            saveas(gcf, fullfile(outdir, ['AOC_HGBIxAlphaLat_timeseries_' taskName '.png']));
        end
    end

    %% ================================================================
    %%    FIGURE 5: Condition Comparison (HGBI and Alpha Lat by condition)
    %% ================================================================
    if ~isResting && nCond > 1
        figure('Position', [0 0 1400 600]);

        % Panel A: HGBI by condition
        subplot(1,2,1); hold on
        for ci = 1:nCond
            cond_vals = nan(nS, 1);
            for si = 1:nS
                idx = all_subj == unique_ids(si) & all_cond == unique_conds(ci);
                if sum(idx) >= 3
                    cond_vals(si) = mean(all_hgbi(idx), 'omitnan');
                end
            end
            vld_c = isfinite(cond_vals);
            xj = ci + 0.15 * randn(sum(vld_c), 1);
            cIdx = min(ci, size(colors, 1));
            scatter(xj, cond_vals(vld_c), 60, colors(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.6);
            m_val = mean(cond_vals(vld_c));
            se_val = std(cond_vals(vld_c)) / sqrt(sum(vld_c));
            plot([ci-0.2, ci+0.2], [m_val m_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 3);
            plot([ci ci], [m_val-se_val m_val+se_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 2);
        end
        yline(0, 'k:', 'LineWidth', 1.5)
        set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', fontSize - 4)
        ylabel('HGBI')
        title('HGBI by Condition', 'FontSize', fontSize - 4)

        % Panel B: Alpha Lateralization by condition
        subplot(1,2,2); hold on
        for ci = 1:nCond
            cond_vals = nan(nS, 1);
            for si = 1:nS
                idx = all_subj == unique_ids(si) & all_cond == unique_conds(ci);
                if sum(idx) >= 3
                    cond_vals(si) = mean(all_alphalat(idx), 'omitnan');
                end
            end
            vld_c = isfinite(cond_vals);
            xj = ci + 0.15 * randn(sum(vld_c), 1);
            cIdx = min(ci, size(colors, 1));
            scatter(xj, cond_vals(vld_c), 60, colors(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.6);
            m_val = mean(cond_vals(vld_c));
            se_val = std(cond_vals(vld_c)) / sqrt(sum(vld_c));
            plot([ci-0.2, ci+0.2], [m_val m_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 3);
            plot([ci ci], [m_val-se_val m_val+se_val], '-', 'Color', colors(cIdx,:), 'LineWidth', 2);
        end
        yline(0, 'k:', 'LineWidth', 1.5)
        set(gca, 'XTick', 1:nCond, 'XTickLabel', condLabels, 'FontSize', fontSize - 4)
        ylabel('Alpha Lateralization')
        title('\alpha Lateralization by Condition', 'FontSize', fontSize - 4)

        sgtitle(sprintf('%s: Condition Comparison', upper(taskName)), 'FontSize', fontSize)
        saveas(gcf, fullfile(outdir, ['AOC_HGBIxAlphaLat_conditions_' taskName '.png']));
    end

    %% Store results
    results.(taskName).r_within  = r_ws;
    results.(taskName).mean_r    = mean_r;
    results.(taskName).ci_r      = ci_r;
    results.(taskName).p_group   = p_group;
    if exist('stats_group', 'var') && isstruct(stats_group)
        results.(taskName).t_group = stats_group.tstat;
        results.(taskName).df_group = stats_group.df;
    end
    if exist('lme_final', 'var') && ~isempty(lme_final)
        results.(taskName).lme = lme_final;
    end
    results.(taskName).nTrials   = numel(all_hgbi);
    results.(taskName).nSubjects = nS;

end % task loop

%% Save all results
save(fullfile(outdir, 'AOC_HGBIxAlphaLat_results.mat'), 'results');
fprintf('\nDone! Results and figures saved to:\n  %s\n', outdir);
