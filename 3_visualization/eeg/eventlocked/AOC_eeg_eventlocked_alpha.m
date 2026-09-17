%% AOC event-locked posterior alpha (microsaccades and saccades)
% Exploratory analysis (not part of the registered-report confirmatory models).
% For Sternberg and N-back, EEG is time-locked to microsaccade onsets (0.1-1.0 deg)
% and to saccade onsets (> 1.0 deg) during retention. Posterior alpha is quantified
% as Hilbert envelope power and as within-trial phase locking (WTPL; Karvat,
% Ofir and Landau, JOCN 2024). Real events are compared with within-trial jittered
% surrogates. Load is retained as a factor.
%
% Event detection uses detect_microsaccades (Engbert and Kliegl, 2003) on the
% full-trial gaze trace, then splits by 2D amplitude at the paradigm ppd.
% Hilbert envelope uses IAF-4/+2 Hz with fallback [8 14] Hz (same rule as ERSD).
% WTPL is computed at 8:1:14 Hz with 5-cycle Morlet wavelets (Karvat Box 1).
%
% Outputs:
%   figures/eeg/eventlocked/*.png
%   data/stats/eventlocked/AOC_eeg_eventlocked_alpha_results.mat
%   data/stats/eventlocked/AOC_eeg_eventlocked_alpha_subject_summary.csv

%% Setup
startup
[subjects, paths, colors, ~] = setup('AOC');
featDir = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'eventlocked');
statDir = fullfile(paths.stats, 'eventlocked');
if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(statDir), mkdir(statDir); end

fs = 500;
screenW = 800;
screenH = 600;
blinkWin = 50;          % samples; matches gaze feature extraction
minISI_s = 0.020;       % merge Engbert fragments closer than 20 ms
ampMS = [0.1 1.0];      % deg
ampSacMin = 1.0;        % deg
ppd = aoc_ppd();
periPlot_s = 0.75;
waveletCycles = 5;
foiWTPL = 8:1:14;       % canonical alpha for WTPL (Karvat 7-13 Hz)
pad_s = waveletCycles / min(foiWTPL);
periExtract_s = periPlot_s + pad_s;
retWin = [0.5 2.0];     % retention; edges trimmed to available epoch length
blWin = [-0.75 -0.50];  % s relative to event, far from peri-event dip
preWin = [-0.40 -0.05];
postWin = [0.05 0.40];
minEvents = 5;
nPerm = 1000;
alphaThr = 0.05;
surrExclude_s = 0.050;
fontSize = 28;
figPos = [0 0 1512 982];
alphaRange = [8 14];

extractSamp = round(periExtract_s * fs);
periTime = -periPlot_s:0.010:periPlot_s;
nPlot = numel(periTime);
nExtract = 2 * extractSamp + 1;
extractTime = (-extractSamp:extractSamp) / fs;
if numel(extractTime) ~= nExtract
    error('Extract time vector does not match the sampling grid.');
end

tasks(1).tag = 'sternberg';
tasks(1).eegFile = 'dataEEG_TFR_sternberg.mat';
tasks(1).etFile = 'dataET_sternberg.mat';
tasks(1).iafFile = 'IAF_sternberg.mat';
tasks(1).iafVars = {'IAF2', 'IAF4', 'IAF6'};
tasks(1).condCodes = [22 24 26];
tasks(1).condLoads = [2 4 6];
tasks(1).condLabels = {'Load 2', 'Load 4', 'Load 6'};

tasks(2).tag = 'nback';
tasks(2).eegFile = 'dataEEG_TFR_nback.mat';
tasks(2).etFile = 'dataET_nback.mat';
tasks(2).iafFile = 'IAF_nback.mat';
tasks(2).iafVars = {'IAF1', 'IAF2', 'IAF3'};
tasks(2).condCodes = [21 22 23];
tasks(2).condLoads = [1 2 3];
tasks(2).condLabels = {'1-back', '2-back', '3-back'};

eventTypes(1).tag = 'ms';
eventTypes(1).label = 'Microsaccade';
eventTypes(1).ampMin = ampMS(1);
eventTypes(1).ampMax = ampMS(2);

eventTypes(2).tag = 'saccade';
eventTypes(2).label = 'Saccade';
eventTypes(2).ampMin = ampSacMin;
eventTypes(2).ampMax = Inf;

occLabs = [];
results = struct();
summaryRows = [];

fprintf('\n=== AOC event-locked alpha (MS and saccades) ===\n');
fprintf('ppd = %.2f px/deg | peri plot = +/-%.0f ms | extract = +/-%.0f ms\n', ...
    ppd, periPlot_s * 1000, periExtract_s * 1000);

for ti = 1:numel(tasks)
    tk = tasks(ti);
    for ei = 1:numel(eventTypes)
        ev = eventTypes(ei);
        fprintf('\n---------- %s | %s ----------\n', upper(tk.tag), upper(ev.label));
        nSubj = numel(subjects);
        nLoad = numel(tk.condCodes);

        envReal = nan(nSubj, nPlot);
        envSurr = nan(nSubj, nPlot);
        wtplReal = nan(nSubj, nPlot);
        wtplSurr = nan(nSubj, nPlot);
        envRealLoad = nan(nSubj, nLoad, nPlot);
        envSurrLoad = nan(nSubj, nLoad, nPlot);
        wtplRealLoad = nan(nSubj, nLoad, nPlot);
        wtplSurrLoad = nan(nSubj, nLoad, nPlot);
        nEvSubj = zeros(nSubj, 1);
        nEvLoad = zeros(nSubj, nLoad);
        nSurrSubj = zeros(nSubj, 1);
        iafUsed = nan(nSubj, 1);
        detectWinUsed = nan(nSubj, 2);

        for s = 1:nSubj
            subj = subjects{s};
            try
            eegPath = fullfile(featDir, subj, 'eeg', tk.eegFile);
            etPath = fullfile(featDir, subj, 'gaze', tk.etFile);
            if ~isfile(eegPath) || ~isfile(etPath)
                fprintf('  %s: missing EEG or ET, skip\n', subj);
                continue
            end

            E = load(eegPath, 'dataTFR');
            dataEEG = E.dataTFR;
            G = load(etPath);
            if isfield(G, 'dataETlong')
                dataET = G.dataETlong;
            else
                fn = fieldnames(G);
                dataET = G.(fn{1});
            end

            if isempty(occLabs)
                occLabs = occipital_labels(dataEEG.label);
                fprintf('  Occipital ROI (%d): %s\n', numel(occLabs), strjoin(occLabs, ', '));
            end

            iafMap = load_iaf_map(fullfile(featDir, subj, 'eeg', tk.iafFile), ...
                tk.iafVars, tk.condCodes, alphaRange);
            iafMean = mean(iafMap.iaf(isfinite(iafMap.iaf)));
            if ~isfinite(iafMean)
                iafMean = 10;
                bandAlpha = alphaRange;
            else
                bandAlpha = iafMap.band;
                if any(~isfinite(bandAlpha))
                    bandAlpha = alphaRange;
                end
            end
            iafUsed(s) = iafMean;
            foi = foiWTPL;

            useLabs = intersect(occLabs, dataEEG.label, 'stable');
            if numel(useLabs) < 3
                fprintf('  %s: fewer than 3 occipital channels, skip\n', subj);
                continue
            end
            cfg = [];
            cfg.channel = useLabs;
            evalc('dataROI = ft_selectdata(cfg, dataEEG);');
            cfg = [];
            cfg.avgoverchan = 'yes';
            evalc('dataOcc = ft_selectdata(cfg, dataROI);');

            cfg = [];
            cfg.demean = 'yes';
            cfg.bpfilter = 'yes';
            cfg.bpfreq = bandAlpha;
            cfg.bpfilttype = 'but';
            cfg.bpfiltord = 4;
            cfg.hilbert = 'abs';
            evalc('envDat = ft_preprocessing(cfg, dataOcc);');

            [idxEEG, idxET] = match_trials(dataEEG.trialinfo, dataET.trialinfo);
            if numel(idxEEG) < 5
                fprintf('  %s: too few matched trials (%d), skip\n', subj, numel(idxEEG));
                continue
            end

            tEEG0 = dataOcc.time{idxEEG(1)};
            detectWin = [max(retWin(1), tEEG0(1) + periExtract_s), ...
                min(retWin(2), tEEG0(end) - periExtract_s)];
            detectWinUsed(s, :) = detectWin;
            if detectWin(2) - detectWin(1) < 0.3
                fprintf('  %s: detect window too short [%.2f %.2f], skip\n', ...
                    subj, detectWin(1), detectWin(2));
                continue
            end

            envEpochs = [];
            envSurrEpochs = [];
            rawEpochs = {};
            rawSurrEpochs = {};
            evLoad = [];
            surrLoad = [];

            rng(1000 * str2double(subj) + ei, 'twister');

            for k = 1:numel(idxEEG)
                trE = idxEEG(k);
                trG = idxET(k);
                tEEG = dataOcc.time{trE};
                tGaze = dataET.time{trG};
                cond = dataEEG.trialinfo(trE, 1);
                loadIdx = find(tk.condCodes == cond, 1);
                if isempty(loadIdx)
                    continue
                end

                gaze = prepare_gaze(dataET.trial{trG}, tGaze, screenW, screenH, blinkWin);
                if sum(~gaze.blink) < 100
                    continue
                end
                evOnsets = detect_amp_events(gaze, fs, ppd, ev.ampMin, ev.ampMax, minISI_s);
                if isempty(evOnsets)
                    tEv = [];
                else
                    tEv = gaze.time(evOnsets);
                    tEv = tEv(tEv >= detectWin(1) & tEv <= detectWin(2));
                end

                envFull = double(envDat.trial{trE}(1, :));
                rawFull = double(dataOcc.trial{trE}(1, :));
                envFull = fillmissing(envFull, 'linear', 'EndValues', 'nearest');
                rawFull = fillmissing(rawFull, 'linear', 'EndValues', 'nearest');

                nBefore = size(envEpochs, 1);
                for ie = 1:numel(tEv)
                    [envSnip, rawSnip] = cut_epoch(tEEG, envFull, rawFull, tEv(ie), extractSamp);
                    if isempty(envSnip), continue; end
                    envEpochs = [envEpochs; envSnip]; %#ok<AGROW>
                    rawEpochs{end + 1} = rawSnip; %#ok<AGROW>
                    evLoad(end + 1, 1) = loadIdx; %#ok<AGROW>
                end
                nAdded = size(envEpochs, 1) - nBefore;
                if nAdded < 1
                    continue
                end
                tSurr = jitter_times(tEv, detectWin, nAdded, surrExclude_s, tEEG, periExtract_s);
                for ie = 1:numel(tSurr)
                    [envSnip, rawSnip] = cut_epoch(tEEG, envFull, rawFull, tSurr(ie), extractSamp);
                    if isempty(envSnip), continue; end
                    envSurrEpochs = [envSurrEpochs; envSnip]; %#ok<AGROW>
                    rawSurrEpochs{end + 1} = rawSnip; %#ok<AGROW>
                    surrLoad(end + 1, 1) = loadIdx; %#ok<AGROW>
                end
            end

            nEvSubj(s) = size(envEpochs, 1);
            nSurrSubj(s) = size(envSurrEpochs, 1);
            for li = 1:nLoad
                nEvLoad(s, li) = sum(evLoad == li);
            end
            if nEvSubj(s) < minEvents || nSurrSubj(s) < minEvents
                fprintf('  %s: %d real / %d surr events, skip GA\n', ...
                    subj, nEvSubj(s), nSurrSubj(s));
                continue
            end

            envReal(s, :) = interp1(extractTime, mean(envEpochs, 1, 'omitnan'), periTime, 'linear');
            envSurr(s, :) = interp1(extractTime, mean(envSurrEpochs, 1, 'omitnan'), periTime, 'linear');
            for li = 1:nLoad
                idxL = evLoad == li;
                idxLs = surrLoad == li;
                if sum(idxL) >= 1
                    envRealLoad(s, li, :) = interp1(extractTime, mean(envEpochs(idxL, :), 1, 'omitnan'), periTime, 'linear');
                end
                if sum(idxLs) >= 1
                    envSurrLoad(s, li, :) = interp1(extractTime, mean(envSurrEpochs(idxLs, :), 1, 'omitnan'), periTime, 'linear');
                end
            end

            wtplR = wtpl_on_epochs(rawEpochs, extractTime, periTime, foi, waveletCycles, fs);
            wtplS = wtpl_on_epochs(rawSurrEpochs, extractTime, periTime, foi, waveletCycles, fs);
            if ~isempty(wtplR)
                wtplReal(s, :) = mean(wtplR, 1, 'omitnan');
            end
            if ~isempty(wtplS)
                wtplSurr(s, :) = mean(wtplS, 1, 'omitnan');
            end
            for li = 1:nLoad
                if ~isempty(wtplR) && any(evLoad == li)
                    wtplRealLoad(s, li, :) = mean(wtplR(evLoad == li, :), 1, 'omitnan');
                end
                if ~isempty(wtplS) && any(surrLoad == li)
                    wtplSurrLoad(s, li, :) = mean(wtplS(surrLoad == li, :), 1, 'omitnan');
                end
            end

            fprintf('  %s (%d/%d): %d %s events (loads %s), IAF = %.1f Hz\n', ...
                subj, s, nSubj, nEvSubj(s), ev.tag, mat2str(nEvLoad(s, :)), iafMean);
            catch ME
                fprintf('  %s: error (%s), skip\n', subj, ME.message);
            end
        end

        tag = sprintf('%s_%s', tk.tag, ev.tag);
        valid = nEvSubj >= minEvents & nSurrSubj >= minEvents & ...
            all(isfinite(envReal), 2) & all(isfinite(envSurr), 2) & ...
            all(isfinite(wtplReal), 2) & all(isfinite(wtplSurr), 2);
        nValid = sum(valid);
        if nValid < 1
            warning('No valid subjects for %s. Skipping stats and figures.', tag);
            continue
        end
        fprintf('  Valid subjects: %d / %d | events/subj: median %.0f (range %d-%d)\n', ...
            nValid, nSubj, median(nEvSubj(valid)), min(nEvSubj(valid)), max(nEvSubj(valid)));
        if nValid < 5
            warning('Too few valid subjects for %s. Skipping stats and figures.', tag);
            continue
        end

        envRealV = envReal(valid, :);
        envSurrV = envSurr(valid, :);
        wtplRealV = wtplReal(valid, :);
        wtplSurrV = wtplSurr(valid, :);
        envRealLoadV = envRealLoad(valid, :, :);
        envSurrLoadV = envSurrLoad(valid, :, :);
        wtplRealLoadV = wtplRealLoad(valid, :, :);
        wtplSurrLoadV = wtplSurrLoad(valid, :, :);

        [envRealBl, envSurrBl] = baseline_pct(envRealV, envSurrV, periTime, blWin);
        [wtplRealBl, wtplSurrBl] = baseline_sub(wtplRealV, wtplSurrV, periTime, blWin);
        envDiff = envRealBl - envSurrBl;
        wtplDiff = wtplRealBl - wtplSurrBl;

        [clusEnv, tEnv, thrEnv] = cluster_permutation_1d(envDiff, nPerm, alphaThr);
        [clusW, tW, thrW] = cluster_permutation_1d(wtplDiff, nPerm, alphaThr);
        print_clusters(sprintf('%s envelope', tag), clusEnv, periTime);
        print_clusters(sprintf('%s WTPL', tag), clusW, periTime);

        preIdx = periTime >= preWin(1) & periTime <= preWin(2);
        postIdx = periTime >= postWin(1) & periTime <= postWin(2);
        envPre = mean(envDiff(:, preIdx), 2, 'omitnan');
        envPost = mean(envDiff(:, postIdx), 2, 'omitnan');
        wtplPre = mean(wtplDiff(:, preIdx), 2, 'omitnan');
        wtplPost = mean(wtplDiff(:, postIdx), 2, 'omitnan');
        [~, pEnvPP, ~, stEnvPP] = ttest(envPost, envPre);
        [~, pWPP, ~, stWPP] = ttest(wtplPost, wtplPre);
        fprintf('  Envelope post vs pre: t(%d) = %.2f, p = %.4f\n', stEnvPP.df, stEnvPP.tstat, pEnvPP);
        fprintf('  WTPL post vs pre: t(%d) = %.2f, p = %.4f\n', stWPP.df, stWPP.tstat, pWPP);

        hi = nLoad;
        lo = 1;
        envHi = squeeze(envRealLoadV(:, hi, :));
        envLo = squeeze(envRealLoadV(:, lo, :));
        loadPair = all(isfinite(envHi), 2) & all(isfinite(envLo), 2);
        clusLoad = struct('idx', {}, 'mass', {}, 'sign', {}, 'p', {});
        if sum(loadPair) >= 5
            envHiBl = baseline_pct_single(envHi(loadPair, :), periTime, blWin);
            envLoBl = baseline_pct_single(envLo(loadPair, :), periTime, blWin);
            [clusLoad, ~, ~] = cluster_permutation_1d(envHiBl - envLoBl, nPerm, alphaThr);
            print_clusters(sprintf('%s envelope load %d vs %d', tag, tk.condLoads(hi), tk.condLoads(lo)), ...
                clusLoad, periTime);
        end

        plot_main_2x2(periTime, envRealBl, envSurrBl, envDiff, clusEnv, ...
            wtplRealBl, wtplSurrBl, wtplDiff, clusW, colors, ev.label, tk.tag, ...
            fontSize, figPos, fullfile(figDir, sprintf('AOC_eventlocked_%s_main.png', tag)));
        plot_by_load(periTime, envRealLoadV, envSurrLoadV, periTime, blWin, colors, ...
            tk.condLabels, ev.label, tk.tag, fontSize, figPos, ...
            fullfile(figDir, sprintf('AOC_eventlocked_%s_byLoad.png', tag)));
        plot_prepost(envPre, envPost, wtplPre, wtplPost, colors, ev.label, tk.tag, ...
            fontSize, figPos, fullfile(figDir, sprintf('AOC_eventlocked_%s_prepost.png', tag)));

        results.(tag).task = tk.tag;
        results.(tag).event = ev.tag;
        results.(tag).periTime = periTime;
        results.(tag).nValid = nValid;
        results.(tag).nEvSubj = nEvSubj;
        results.(tag).nEvLoad = nEvLoad;
        results.(tag).iafUsed = iafUsed;
        results.(tag).detectWin = detectWinUsed;
        results.(tag).envRealBl = envRealBl;
        results.(tag).envSurrBl = envSurrBl;
        results.(tag).envDiff = envDiff;
        results.(tag).wtplRealBl = wtplRealBl;
        results.(tag).wtplSurrBl = wtplSurrBl;
        results.(tag).wtplDiff = wtplDiff;
        results.(tag).envRealLoad = envRealLoadV;
        results.(tag).envSurrLoad = envSurrLoadV;
        results.(tag).wtplRealLoad = wtplRealLoadV;
        results.(tag).wtplSurrLoad = wtplSurrLoadV;
        results.(tag).condLoads = tk.condLoads;
        results.(tag).clustersEnv = clusEnv;
        results.(tag).clustersWTPL = clusW;
        results.(tag).thrEnv = thrEnv;
        results.(tag).thrWTPL = thrW;
        results.(tag).tEnv = tEnv;
        results.(tag).tWTPL = tW;
        results.(tag).envPre = envPre;
        results.(tag).envPost = envPost;
        results.(tag).wtplPre = wtplPre;
        results.(tag).wtplPost = wtplPost;
        results.(tag).pEnvPrePost = pEnvPP;
        results.(tag).pWTPLPrePost = pWPP;
        results.(tag).clustersLoad = clusLoad;
        results.(tag).ppd = ppd;
        results.(tag).foiBand = alphaRange;
        results.(tag).valid = valid;

        validIdx = find(valid);
        for si = 1:nValid
            row = table;
            row.ID = str2double(subjects{validIdx(si)});
            row.Task = string(tk.tag);
            row.Event = string(ev.tag);
            row.NEvents = nEvSubj(validIdx(si));
            row.IAF = iafUsed(validIdx(si));
            row.EnvPre = envPre(si);
            row.EnvPost = envPost(si);
            row.WTPLPre = wtplPre(si);
            row.WTPLPost = wtplPost(si);
            summaryRows = [summaryRows; row]; %#ok<AGROW>
        end
    end
end

results.ppd = ppd;
results.foiWTPL = foiWTPL;
results.periPlot_s = periPlot_s;
results.blWin = blWin;
results.preWin = preWin;
results.postWin = postWin;
results.ampMS = ampMS;
results.ampSacMin = ampSacMin;
save(fullfile(statDir, 'AOC_eeg_eventlocked_alpha_results.mat'), 'results');
if ~isempty(summaryRows)
    writetable(summaryRows, fullfile(statDir, 'AOC_eeg_eventlocked_alpha_subject_summary.csv'));
end
fprintf('\nSaved figures to %s\n', figDir);
fprintf('Saved results to %s\n', statDir);

%% Local functions
function ppd = aoc_ppd()
viewDist_mm = 680;
ppm = 3.6;
ecc_dva = 3;
mpd = (viewDist_mm / 2) * tan(deg2rad(2 * ecc_dva)) / ecc_dva;
ppd = ppm * mpd;
end

function labs = occipital_labels(allLabs)
labs = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, 'O') || contains(L, 'I')
        labs{end + 1} = L; %#ok<AGROW>
    end
end
end

function iafMap = load_iaf_map(iafFile, iafVars, condCodes, alphaRange)
iafMap.cond = condCodes(:);
iafMap.iaf = nan(numel(condCodes), 1);
iafMap.band = alphaRange;
if ~isfile(iafFile)
    return
end
S = load(iafFile);
for i = 1:numel(iafVars)
    if isfield(S, iafVars{i})
        iafMap.iaf(i) = S.(iafVars{i});
    end
end
iafOk = iafMap.iaf(isfinite(iafMap.iaf) & iafMap.iaf > alphaRange(1) & iafMap.iaf < alphaRange(2));
if ~isempty(iafOk)
    iafM = mean(iafOk);
    iafMap.band = [iafM - 4, iafM + 2];
end
end

function [idxEEG, idxET] = match_trials(infoEEG, infoET)
nE = size(infoEEG, 1);
nG = size(infoET, 1);
if nE < 1 || nG < 1
    idxEEG = [];
    idxET = [];
    return
end
if size(infoEEG, 2) >= 2 && size(infoET, 2) >= 2
    idE = infoEEG(:, 2);
    idG = infoET(:, 2);
    if all(isfinite(idE)) && all(isfinite(idG))
        [~, idxEEG, idxET] = intersect(idE, idG, 'stable');
        return
    end
end
n = min(nE, nG);
idxEEG = (1:n)';
idxET = (1:n)';
end

function gaze = prepare_gaze(raw, t, screenW, screenH, blinkWin)
raw = double(raw(1:min(3, size(raw, 1)), :));
X = raw(1, :);
Yorig = raw(2, :);
oob = X < 0 | X > screenW | Yorig < 0 | Yorig > screenH | ~isfinite(X) | ~isfinite(Yorig);
dat = raw;
dat(1, :) = X;
dat(2, :) = screenH - Yorig;
if size(dat, 1) >= 3
    dat(3, oob) = NaN;
end
dat(1, oob) = NaN;
dat(2, oob) = NaN;
dat = remove_blinks(dat, blinkWin);
blink = ~isfinite(dat(1, :)) | ~isfinite(dat(2, :));
Xi = dat(1, :);
Yi = dat(2, :);
if any(isfinite(Xi))
    Xi = fillmissing(Xi, 'linear', 'EndValues', 'nearest');
    Yi = fillmissing(Yi, 'linear', 'EndValues', 'nearest');
end
gaze.X = Xi;
gaze.Y = Yi;
gaze.blink = blink;
gaze.time = t;
end

function onsets = detect_amp_events(gaze, fs, ppd, ampMin, ampMax, minISI_s)
onsets = [];
if numel(gaze.X) < 50
    return
end
[~, det] = detect_microsaccades(fs, [gaze.X; gaze.Y], numel(gaze.X));
if isempty(det.Onset)
    return
end
on = det.Onset(:);
off = det.Offset(:);
minISI = max(1, round(minISI_s * fs));
keepOn = on(1);
keepOff = off(1);
onM = [];
offM = [];
for i = 2:numel(on)
    if on(i) - keepOff <= minISI
        keepOff = off(i);
    else
        onM(end + 1, 1) = keepOn; %#ok<AGROW>
        offM(end + 1, 1) = keepOff; %#ok<AGROW>
        keepOn = on(i);
        keepOff = off(i);
    end
end
onM(end + 1, 1) = keepOn;
offM(end + 1, 1) = keepOff;

nPos = numel(gaze.X);
keep = false(size(onM));
for i = 1:numel(onM)
    if onM(i) < 1 || offM(i) > nPos || offM(i) < onM(i)
        continue
    end
    if any(gaze.blink(onM(i):offM(i)))
        continue
    end
    amp = hypot(gaze.X(offM(i)) - gaze.X(onM(i)), gaze.Y(offM(i)) - gaze.Y(onM(i))) / ppd;
    keep(i) = isfinite(amp) && amp >= ampMin && amp <= ampMax;
end
onsets = onM(keep);
end

function [envSnip, rawSnip] = cut_epoch(tEEG, envFull, rawFull, t0, extractSamp)
envSnip = [];
rawSnip = [];
[~, i0] = min(abs(tEEG - t0));
lo = i0 - extractSamp;
hi = i0 + extractSamp;
if lo < 1 || hi > numel(envFull)
    return
end
envSnip = envFull(lo:hi);
rawSnip = rawFull(lo:hi);
if any(~isfinite(envSnip)) || any(~isfinite(rawSnip))
    envSnip = [];
    rawSnip = [];
end
end

function tSurr = jitter_times(tEv, detectWin, nNeed, exclude_s, tEEG, periExtract_s)
tMin = max(detectWin(1), tEEG(1) + periExtract_s);
tMax = min(detectWin(2), tEEG(end) - periExtract_s);
tSurr = nan(nNeed, 1);
if ~(tMax > tMin)
    return
end
cand = tMin + (tMax - tMin) * rand(nNeed * 20, 1);
keep = false(size(cand));
for i = 1:numel(cand)
    if all(abs(cand(i) - tEv(:)) > exclude_s)
        keep(i) = true;
    end
end
cand = cand(keep);
if numel(cand) >= nNeed
    tSurr = cand(1:nNeed);
else
    tSurr = tMin + (tMax - tMin) * rand(nNeed, 1);
end
end

function wtplBand = wtpl_on_epochs(rawEpochs, extractTime, periTime, foi, widthCyc, fs)
wtplBand = [];
nEv = numel(rawEpochs);
if nEv < 1
    return
end
dataEv = [];
dataEv.label = {'occ'};
dataEv.fsample = fs;
dataEv.trial = rawEpochs;
dataEv.time = repmat({extractTime}, 1, nEv);
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'fourier';
cfg.foi = foi;
cfg.width = widthCyc;
cfg.toi = periTime;
cfg.keeptrials = 'yes';
cfg.pad = 'nextpow2';
evalc('freq = ft_freqanalysis(cfg, dataEv);');
nRpt = size(freq.fourierspctrm, 1);
nFr = numel(freq.freq);
nTi = numel(freq.time);
F = reshape(freq.fourierspctrm(:, 1, :, :), [nRpt, nFr, nTi]);
WTPL = wtpl_vectorized(F, freq.freq, freq.time, [-1 1]);
wtplF = squeeze(mean(WTPL, 2, 'omitnan'));
if nRpt == 1
    wtplF = wtplF(:)';
end
if numel(freq.time) == numel(periTime) && max(abs(freq.time(:)' - periTime(:)')) < 1e-8
    wtplBand = wtplF;
else
    wtplBand = nan(nRpt, numel(periTime));
    for i = 1:nRpt
        wtplBand(i, :) = interp1(freq.time, wtplF(i, :), periTime, 'linear');
    end
end
end

function WTPL = wtpl_vectorized(fourier, freqVec, timeVec, lags)
% Within-trial phase lock (Karvat, Ofir and Landau, JOCN 2024, Box 1 / Eq. 2).
% fourier is rpt x freq x time. lags are in cycles, excluding 0.
nRpt = size(fourier, 1);
nFreq = size(fourier, 2);
nTime = size(fourier, 3);
dt = median(diff(timeVec));
lags = lags(lags ~= 0);
nLag = numel(lags);
WTPL = nan(nRpt, nFreq, nTime);
if nLag < 1
    return
end
for fi = 1:nFreq
    if ~(freqVec(fi) > 0)
        continue
    end
    phi = reshape(angle(fourier(:, fi, :)), [nRpt, nTime]);
    acc = zeros(nRpt, nTime);
    nUsed = zeros(1, nTime);
    for li = 1:nLag
        shiftSamp = round((lags(li) / freqVec(fi)) / dt);
        if shiftSamp == 0
            continue
        end
        src = (1:nTime) + shiftSamp;
        valid = src >= 1 & src <= nTime;
        dphi = nan(nRpt, nTime);
        dphi(:, valid) = phi(:, valid) - phi(:, src(valid));
        acc(:, valid) = acc(:, valid) + exp(1i * dphi(:, valid));
        nUsed(valid) = nUsed(valid) + 1;
    end
    mag = abs(acc);
    den = nUsed;
    den(den == 0) = NaN;
    WTPL(:, fi, :) = mag ./ den;
end
end

function [realBl, surrBl] = baseline_pct(realMat, surrMat, t, blWin)
bl = t >= blWin(1) & t <= blWin(2);
realBl = nan(size(realMat));
surrBl = nan(size(surrMat));
for i = 1:size(realMat, 1)
    bR = mean(realMat(i, bl), 'omitnan');
    bS = mean(surrMat(i, bl), 'omitnan');
    if bR > 0
        realBl(i, :) = (realMat(i, :) - bR) / bR * 100;
    end
    if bS > 0
        surrBl(i, :) = (surrMat(i, :) - bS) / bS * 100;
    end
end
end

function out = baseline_pct_single(mat, t, blWin)
bl = t >= blWin(1) & t <= blWin(2);
out = nan(size(mat));
for i = 1:size(mat, 1)
    b = mean(mat(i, bl), 'omitnan');
    if b > 0
        out(i, :) = (mat(i, :) - b) / b * 100;
    end
end
end

function [realBl, surrBl] = baseline_sub(realMat, surrMat, t, blWin)
bl = t >= blWin(1) & t <= blWin(2);
realBl = realMat - mean(realMat(:, bl), 2, 'omitnan');
surrBl = surrMat - mean(surrMat(:, bl), 2, 'omitnan');
end

function [clusters, tvals, thr] = cluster_permutation_1d(S, nPerm, alpha)
% One-sample cluster permutation vs 0 (Maris and Oostenveld, 2007), time only.
ns = size(S, 1);
nFin = sum(isfinite(S), 1);
se = std(S, [], 1, 'omitnan') ./ sqrt(max(nFin, 1));
m = mean(S, 1, 'omitnan');
tvals = m ./ se;
tvals(~isfinite(tvals)) = 0;
tcrit = tinv(1 - 0.5 * alpha, max(ns - 1, 1));
thr.tcrit = tcrit;
clusters = compute_clusters(tvals, tcrit);
maxMass = zeros(1, nPerm);
for p = 1:nPerm
    flips = (rand(ns, 1) > 0.5) * 2 - 1;
    Sprm = S .* flips;
    seP = std(Sprm, [], 1, 'omitnan') ./ sqrt(max(sum(isfinite(Sprm), 1), 1));
    tP = mean(Sprm, 1, 'omitnan') ./ seP;
    tP(~isfinite(tP)) = 0;
    clP = compute_clusters(tP, tcrit);
    if isempty(clP)
        maxMass(p) = 0;
    else
        maxMass(p) = max([clP.mass]);
    end
end
thr.mass = prctile(maxMass, 100 * (1 - alpha));
for k = 1:numel(clusters)
    clusters(k).p = max(mean(maxMass >= clusters(k).mass), 1 / nPerm);
end
end

function clusters = compute_clusters(tvals, tcrit)
above = abs(tvals) > tcrit;
clusters = struct('idx', {}, 'mass', {}, 'sign', {}, 'p', {});
if ~any(above)
    return
end
d = diff([0, above, 0]);
on = find(d == 1);
off = find(d == -1) - 1;
for c = 1:numel(on)
    idx = on(c):off(c);
    clusters(end + 1).idx = idx; %#ok<AGROW>
    clusters(end).mass = sum(abs(tvals(idx)));
    clusters(end).sign = sign(mean(tvals(idx)));
    clusters(end).p = NaN;
end
end

function print_clusters(name, clusters, t)
nSig = 0;
for k = 1:numel(clusters)
    if clusters(k).p < 0.05
        nSig = nSig + 1;
        fprintf('  %s cluster %d: [%+.0f, %+.0f] ms, p = %.4f\n', ...
            name, nSig, t(clusters(k).idx(1)) * 1000, t(clusters(k).idx(end)) * 1000, clusters(k).p);
    end
end
if nSig == 0
    fprintf('  %s: no significant clusters\n', name);
end
end

function plot_sem(t, M, col, ls, lw)
ga = mean(M, 1, 'omitnan');
sem = std(M, [], 1, 'omitnan') / sqrt(sum(all(isfinite(M), 2)));
fill([t, fliplr(t)], [ga - sem, fliplr(ga + sem)], col, ...
    'FaceAlpha', 0.22, 'EdgeColor', 'none');
plot(t, ga, ls, 'LineWidth', lw, 'Color', col);
end

function shade_clusters(clusters, t, y, col)
yl = ylim;
for k = 1:numel(clusters)
    if ~(clusters(k).p < 0.05)
        continue
    end
    idx = clusters(k).idx;
    xx = [t(idx(1)), t(idx(end)), t(idx(end)), t(idx(1))];
    yy = [yl(1), yl(1), yl(2), yl(2)];
    patch(xx, yy, col, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
end
if nargin >= 3 && ~isempty(y)
    % y unused; ylim already applied
end
end

function plot_main_2x2(t, envR, envS, envD, clusE, wtR, wtS, wtD, clusW, ...
    colors, evLab, taskTag, fsz, figPos, outFile)
tms = t * 1000;
grey = [0.55 0.55 0.55];
figure('Position', figPos, 'Color', 'w');

subplot(2, 2, 1); hold on
plot_sem(tms, envS, grey, '--', 2);
plot_sem(tms, envR, colors(1, :), '-', 2.5);
xline(0, 'k:'); yline(0, 'k:');
xlabel('Time from event onset [ms]');
ylabel('Alpha envelope [% of baseline]');
title(sprintf('%s %s: power', upper(taskTag), evLab));
legend({'Surrogate', 'Event'}, 'Location', 'best', 'Box', 'off');
set(gca, 'FontSize', fsz - 8); box off

subplot(2, 2, 2); hold on
plot_sem(tms, envD, colors(2, :), '-', 2.5);
xline(0, 'k:'); yline(0, 'k:');
shade_clusters(clusE, tms, [], colors(2, :));
xlabel('Time from event onset [ms]');
ylabel('Power difference (event minus surrogate) [%]');
title('Power: event minus surrogate');
set(gca, 'FontSize', fsz - 8); box off

subplot(2, 2, 3); hold on
plot_sem(tms, wtS, grey, '--', 2);
plot_sem(tms, wtR, colors(3, :), '-', 2.5);
xline(0, 'k:'); yline(0, 'k:');
xlabel('Time from event onset [ms]');
ylabel('WTPL (baseline subtracted)');
title(sprintf('%s %s: WTPL', upper(taskTag), evLab));
legend({'Surrogate', 'Event'}, 'Location', 'best', 'Box', 'off');
set(gca, 'FontSize', fsz - 8); box off

subplot(2, 2, 4); hold on
plot_sem(tms, wtD, colors(3, :), '-', 2.5);
xline(0, 'k:'); yline(0, 'k:');
shade_clusters(clusW, tms, [], colors(3, :));
xlabel('Time from event onset [ms]');
ylabel('WTPL difference (event minus surrogate)');
title('WTPL: event minus surrogate');
set(gca, 'FontSize', fsz - 8); box off

sgtitle(sprintf('Event-locked posterior alpha  |  %s  |  %s', upper(taskTag), evLab), ...
    'FontSize', fsz, 'FontWeight', 'bold');
saveas(gcf, outFile);
close(gcf);
end

function plot_by_load(t, envRealLoad, envSurrLoad, periTime, blWin, colors, ...
    condLabels, evLab, taskTag, fsz, figPos, outFile)
tms = t * 1000;
nLoad = size(envRealLoad, 2);
figure('Position', figPos, 'Color', 'w');
hold on
h = gobjects(nLoad, 1);
for li = 1:nLoad
    R = squeeze(envRealLoad(:, li, :));
    S = squeeze(envSurrLoad(:, li, :));
    good = all(isfinite(R), 2) & all(isfinite(S), 2);
    if sum(good) < 3
        continue
    end
    [Rb, Sb] = baseline_pct(R(good, :), S(good, :), periTime, blWin);
    D = Rb - Sb;
    ga = mean(D, 1, 'omitnan');
    sem = std(D, [], 1, 'omitnan') / sqrt(sum(good));
    fill([tms, fliplr(tms)], [ga - sem, fliplr(ga + sem)], colors(li, :), ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    h(li) = plot(tms, ga, '-', 'LineWidth', 2.5, 'Color', colors(li, :));
end
xline(0, 'k:'); yline(0, 'k:');
xlabel('Time from event onset [ms]');
ylabel('Power difference (event minus surrogate) [%]');
title(sprintf('%s %s: load', upper(taskTag), evLab));
ok = isgraphics(h);
if any(ok)
    legend(h(ok), condLabels(ok), 'Location', 'best', 'Box', 'off');
end
set(gca, 'FontSize', fsz - 6); box off
saveas(gcf, outFile);
close(gcf);
end

function plot_prepost(envPre, envPost, wtplPre, wtplPost, colors, evLab, taskTag, ...
    fsz, figPos, outFile)
figure('Position', figPos, 'Color', 'w');
subplot(1, 2, 1);
plot_paired_dots(envPre, envPost, colors(1, :), 'Alpha envelope difference [%]', fsz);
title(sprintf('%s %s: power', upper(taskTag), evLab));
subplot(1, 2, 2);
plot_paired_dots(wtplPre, wtplPost, colors(3, :), 'WTPL difference', fsz);
title(sprintf('%s %s: WTPL', upper(taskTag), evLab));
saveas(gcf, outFile);
close(gcf);
end

function plot_paired_dots(pre, post, col, ylab, fsz)
hold on
n = numel(pre);
for i = 1:n
    if isfinite(pre(i)) && isfinite(post(i))
        plot([1 2], [pre(i) post(i)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    end
end
scatter(ones(n, 1) + 0.04 * randn(n, 1), pre, 70, col, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(2 * ones(n, 1) + 0.04 * randn(n, 1), post, 70, col, 'filled', 'MarkerFaceAlpha', 0.7);
plot([1 2], [mean(pre, 'omitnan') mean(post, 'omitnan')], 'k-', 'LineWidth', 2.5);
yline(0, 'k:');
xlim([0.5 2.5]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre [-400 -50] ms', 'Post [50 400] ms'}, ...
    'FontSize', fsz - 8);
ylabel(ylab);
box off
end
