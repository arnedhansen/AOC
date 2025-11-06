%% AOC Sternberg — α ↔ oculomotor couplings
% 1) CCF curves (α × speed @500 Hz; α × scan-path length; α × microsaccade rate) + SEM
% 2) For each measure: per-lag t-curve with cluster-permutation shading
% 3) For each measure: window-level subject metrics with group inference

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');
% subjects  = subjects(1:10);

set(0,'DefaultFigureColor','w')
fontSize  = 26;

fs        = 500;                  % Hz (EEG & ET)
epochWin  = [-.5 2.5];            % padding for filtering
anaWin    = [0 2];                % retention window
maxLagS   = 1;                    % seconds (+/-)
maxLag_500= round(maxLagS * fs);  % ±lag at 500 Hz

binDur    = 0.05;                 % 50 ms bins for SPL & MS
fsb       = 1/binDur;             % 20 Hz
maxLag_20 = round(maxLagS * fsb); % ±lag at 20 Hz

velZthr   = 4;                    % |z|>4 → interpolate
polyOrd   = 3;                    % SG derivative order

doSSD     = true;                 % use SSD α component instead of ROI average
doPrewhite= true;                 % AR(1) prewhitening before xcorr
nPerm     = 2000;                 % cluster-permutation count
alphaThr  = 0.05;                 % cluster-level alpha

% Stats windows (ms), positive = gaze follows α
winPos    = [0   250];             % early positive window
winNeg    = [250 500];            % late negative window (often negative)

outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions';
if ~exist(outdir,'dir'); mkdir(outdir); end

%% Occipital ROI labels containing 'O' or 'I'
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat');
allLabs = powload2_late.label;
roi_labels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, {'O'}) || contains(L, {'I'})
        roi_labels{end+1} = L; %#ok<SAGROW>
    end
end

%% Holders
CCF_speed_allSubs = {};   % per-subject: [nTrials × nLags(500 Hz)]
CCF_spl_allSubs   = {};   % per-subject: [nTrials × nLags(20 Hz)]
CCF_ms_allSubs    = {};   % per-subject: [nTrials × nLags(20 Hz)]

lags_500   = -maxLag_500:maxLag_500;
lags_20    = -maxLag_20:maxLag_20;
lags_ms_500= (lags_500/fs)*1000;
lags_ms_20 = (lags_20/fsb)*1000;
xmax500 = max(lags_ms_500);
xmax20  = max(lags_ms_20);
tickstep = 100;

subj_ids = cell(numel(subjects),1);

%% Loop subjects
for s = 1:numel(subjects)
    clc; fprintf('Subject %s (%d/%d)\n', subjects{s}, s, numel(subjects));
    subj_ids{s} = subjects{s};

    dpeeg  = fullfile(path, subjects{s}, 'eeg');
    dpgaze = fullfile(path, subjects{s}, 'gaze');

    % EEG time series
    cd(dpeeg)
    load dataEEG_TFR_sternberg

    % Subject IAF
    IAF_file = fullfile(dpeeg, 'IAF_sternberg_subject.mat');
    if exist(IAF_file, 'file'), tmp = load(IAF_file, 'IAF_subj'); IAF = tmp.IAF_subj;
    else, IAF = NaN; end
    if isfinite(IAF), bandAlpha = [IAF-4, IAF+2]; else, bandAlpha = [8 14]; end
    flank1 = [bandAlpha(1)-4, bandAlpha(1)-2];   % lower noise flank
    flank2 = [bandAlpha(2)+2, bandAlpha(2)+4];   % upper noise flank

    % Build α source: SSD component or ROI mean
    if doSSD
        alphaSrc = ssd_alpha_component(dataTFR, roi_labels, bandAlpha, flank1, flank2);
    else
        cfg = []; cfg.channel = roi_labels; roi = ft_selectdata(cfg, dataTFR);
        cfg = []; cfg.avgoverchan = 'yes'; alphaSrc = ft_selectdata(cfg, roi);
    end
    % Hilbert envelope
    cfg = []; cfg.bpfilter='yes'; cfg.bpfreq = bandAlpha;
    cfg.demean='yes'; cfg.hilbert='complex';
    aH = ft_preprocessing(cfg, alphaSrc);
    aH = ft_selectdata(struct('latency',epochWin), aH);
    aH = ft_selectdata(struct('latency',anaWin),   aH);

    nT_env = numel(aH.trial);
    aenv_500 = cell(1,nT_env);
    aenv_20  = cell(1,nT_env);
    for k = 1:nT_env
        ae = abs(aH.trial{k}(1,:));
        % interpolate residual NaNs
        if any(isnan(ae)), ae = fillmissing(ae,'linear','EndValues','nearest'); end
        % prewhiten (optional)
        if doPrewhite, ae = ar1_prewhiten(ae); end
        aenv_500{k} = ae;
        aenv_20{k}  = bin_series(ae, fs, binDur, 'mean');
    end

    % Eye-tracking (left eye)
    cd(dpgaze)
    load dataET_sternberg dataETlong

    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y','L-AREA'};
    et  = ft_selectdata(cfg, dataETlong);
    et  = ft_selectdata(struct('latency',epochWin), et);
    et  = ft_selectdata(struct('latency',anaWin),   et);

    nT = min([numel(aenv_500), numel(et.trial), size(dataETlong.trialinfo,1)]);
    aenv_500 = aenv_500(1:nT);
    aenv_20  = aenv_20(1:nT);

    ccf_speed = nan(nT, numel(lags_500));
    ccf_spl   = nan(nT, numel(lags_20));
    ccf_ms    = nan(nT, numel(lags_20));

    for tr = 1:nT
        X = double(et.trial{tr}(1,:));
        Y = double(et.trial{tr}(2,:));
        A = double(et.trial{tr}(3,:));

        % invert Y (Cartesian up)
        Y = -Y;

        % Blink removal by AREA + linear interp with edge padding
        blink = ~isfinite(A) | (A<=0);
        if any(blink)
            X(blink)=NaN; Y(blink)=NaN;
            if isnan(X(1)),  i1=find(isfinite(X),1,'first'); if ~isempty(i1), X(1:i1-1)=X(i1); end, end
            if isnan(X(end)), i2=find(isfinite(X),1,'last');  if ~isempty(i2), X(i2+1:end)=X(i2); end, end
            if isnan(Y(1)),  i1=find(isfinite(Y),1,'first'); if ~isempty(i1), Y(1:i1-1)=Y(i1); end, end
            if isnan(Y(end)), i2=find(isfinite(Y),1,'last');  if ~isempty(i2), Y(i2+1:end)=Y(i2); end, end
            X = fillmissing(X,'linear');  Y = fillmissing(Y,'linear');
        end

        % --- Speed @ 500 Hz
        [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd);
        % outlier handle
        zvx = (vx - nanmean(vx)) / nanstd(vx); zvy = (vy - nanmean(vy)) / nanstd(vy);
        bad = abs(zvx)>velZthr | abs(zvy)>velZthr;
        if any(bad)
            vx(bad)=NaN; vy(bad)=NaN;
            if isnan(vx(1)),  i1=find(isfinite(vx),1,'first'); if ~isempty(i1), vx(1:i1-1)=vx(i1); end, end
            if isnan(vx(end)), i2=find(isfinite(vx),1,'last');  if ~isempty(i2), vx(i2+1:end)=vx(i2); end, end
            if isnan(vy(1)),  i1=find(isfinite(vy),1,'first'); if ~isempty(i1), vy(1:i1-1)=vy(i1); end, end
            if isnan(vy(end)), i2=find(isfinite(vy),1,'last');  if ~isempty(i2), vy(i2+1:end)=vy(i2); end, end
            vx = fillmissing(vx,'linear'); vy = fillmissing(vy,'linear');
        end
        spd = hypot(vx, vy);
        if doPrewhite, spd = ar1_prewhiten(spd); end

        ae = aenv_500{tr};
        aZ = zscore(ae); sZ = zscore(spd);
        ccf_speed(tr,:) = xcorr(aZ, sZ, maxLag_500, 'coeff');

        % --- Scan-path length (50 ms sums)
        step = sqrt(diff(X).^2 + diff(Y).^2);
        step = [step, step(end)];
        if doPrewhite, step = ar1_prewhiten(step); end
        spl_bins = bin_series(step, fs, binDur, 'sum');
        ae20 = aenv_20{tr};
        nb = min(numel(ae20), numel(spl_bins));
        aZb = zscore(ae20(1:nb));
        splZ= zscore(spl_bins(1:nb));
        ccf_spl(tr,1:numel(lags_20)) = xcorr(aZb, splZ, maxLag_20, 'coeff');

        % --- Microsaccade rate (events/s in 50 ms)
        ms = ms_detect_engbert(X, Y, fs, 6, 6e-3, 0.02);
        ms_on = zeros(size(X));
        if ~isempty(ms.onsets)
            on = ms.onsets; on = on(on>=1 & on<=numel(ms_on));
            ms_on(on) = 1;
        end
        ms_rate = bin_series(ms_on, fs, binDur, 'sum')/binDur; % events/s
        nb = min(numel(ae20), numel(ms_rate));
        aZb = zscore(ae20(1:nb));
        msZ = zscore(ms_rate(1:nb));
        ccf_ms(tr,1:numel(lags_20)) = xcorr(aZb, msZ, maxLag_20, 'coeff');
    end

    CCF_speed_allSubs{end+1} = ccf_speed;
    CCF_spl_allSubs{end+1}   = ccf_spl;
    CCF_ms_allSubs{end+1}    = ccf_ms;
end

%% Group curves (trial-pooled, then subject-balanced for stats)
[CCF_speed_mu,  CCF_speed_sem,  CCF_speed_subj]  = group_curve(CCF_speed_allSubs);
[CCF_spl_mu,    CCF_spl_sem,    CCF_spl_subj]    = group_curve(CCF_spl_allSubs);
[CCF_ms_mu,     CCF_ms_sem,     CCF_ms_subj]     = group_curve(CCF_ms_allSubs);

%% Figure 1: α × speed (curve)
close all
figure('Position',[0 0 1200 1000]);
fill([lags_ms_500, fliplr(lags_ms_500)], [CCF_speed_mu-CCF_speed_sem, fliplr(CCF_speed_mu+CCF_speed_sem)], 0.85*[1 1 1], 'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_500, CCF_speed_mu, 'LineWidth', 3, 'Color', colors(3,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)'), ylabel('Correlation (r)')
% xticks(-500:50:500), ylim(1.25*[-max(abs(CCF_speed_mu)) max(abs(CCF_speed_mu))])
xticks(-xmax500:tickstep:xmax500);
title('CCF: Alpha Envelope \times Gaze Speed','FontSize',fontSize), set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_alpha_speed_CCF.png'));

%% Figure 2: α × scan-path length (curve)
figure('Position',[0 0 1200 1000]);
fill([lags_ms_20, fliplr(lags_ms_20)], [CCF_spl_mu-CCF_spl_sem, fliplr(CCF_spl_mu+CCF_spl_sem)], 0.85*[1 1 1], 'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_20, CCF_spl_mu, 'LineWidth', 3, 'Color', colors(2,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)'), ylabel('Correlation (r)')
% xticks(-500:50:500), ylim(1.25*[-max(abs(CCF_spl_mu)) max(abs(CCF_spl_mu))])
xticks(-xmax20:tickstep:xmax20);
title('CCF: Alpha Envelope \times Scan-Path Length (50 ms bins)','FontSize',fontSize), set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_alpha_spl_CCF.png'));

%% Figure 3: α × microsaccade rate (curve)
figure('Position',[0 0 1200 1000]);
fill([lags_ms_20, fliplr(lags_ms_20)], [CCF_ms_mu-CCF_ms_sem, fliplr(CCF_ms_mu+CCF_ms_sem)], 0.85*[1 1 1], 'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_20, CCF_ms_mu, 'LineWidth', 3, 'Color', colors(1,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)'), ylabel('Correlation (r)')
% xticks(-500:50:500), ylim(1.25*[-max(abs(CCF_ms_mu)) max(abs(CCF_ms_mu))]);
xticks(-xmax20:tickstep:xmax20)
title('CCF: Alpha Envelope \times Microsaccade Rate (50 ms bins)','FontSize',fontSize), set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_alpha_ms_CCF.png'));

%% Cluster-permutation across lags (subject-balanced)
% Prepare subject means over trials (subjects × lags)
S_speed = subj_mean_over_trials(CCF_speed_allSubs);
S_spl   = subj_mean_over_trials(CCF_spl_allSubs);
S_ms    = subj_mean_over_trials(CCF_ms_allSubs);

% Run cluster permutation (one-sample vs 0)
[clusS, tS, thrS] = cluster_permutation_1d(S_speed, nPerm, alphaThr);
[clusP, tP, thrP] = cluster_permutation_1d(S_spl,   nPerm, alphaThr);
[clusM, tM, thrM] = cluster_permutation_1d(S_ms,    nPerm, alphaThr);

% Plot per-lag t and significant clusters
plot_tcurve_with_clusters(lags_ms_500, tS, clusS, thrS, 'Speed', colors(3,:), outdir, fontSize);
plot_tcurve_with_clusters(lags_ms_20,  tP, clusP, thrP, 'ScanPath', colors(2,:), outdir, fontSize);
plot_tcurve_with_clusters(lags_ms_20,  tM, clusM, thrM, 'Microsaccades', colors(1,:), outdir, fontSize);

%% Window-level subject stats (Fisher-z; group t; CI)
winMaskS_pos = lags_ms_500>=winPos(1) & lags_ms_500<=winPos(2);
winMaskS_neg = lags_ms_500>=winNeg(1) & lags_ms_500<=winNeg(2);
winMaskP_pos = lags_ms_20 >=winPos(1) & lags_ms_20 <=winPos(2);
winMaskP_neg = lags_ms_20 >=winNeg(1) & lags_ms_20 <=winNeg(2);

stats_speed_pos = window_stats(S_speed, winMaskS_pos);
stats_speed_neg = window_stats(S_speed, winMaskS_neg);
stats_spl_pos   = window_stats(S_spl,   winMaskP_pos);
stats_spl_neg   = window_stats(S_spl,   winMaskP_neg);
stats_ms_pos    = window_stats(S_ms,    winMaskP_pos);
stats_ms_neg    = window_stats(S_ms,    winMaskP_neg);

% Visualise subject metrics as rain/bar hybrid
plot_window_stats(stats_speed_pos, 'Speed (+ window)', colors(3,:), outdir, fontSize);
plot_window_stats(stats_speed_neg, 'Speed (late window)', colors(3,:), outdir, fontSize);
plot_window_stats(stats_spl_pos,   'ScanPath (+ window)', colors(2,:), outdir, fontSize);
plot_window_stats(stats_spl_neg,   'ScanPath (late window)', colors(2,:), outdir, fontSize);
plot_window_stats(stats_ms_pos,    'Microsaccades (+ window)', colors(1,:), outdir, fontSize);
plot_window_stats(stats_ms_neg,    'Microsaccades (late window)', colors(1,:), outdir, fontSize);

%% Save summary tables
save(fullfile(outdir,'AOC_alpha_windows_stats.mat'), ...
    'stats_speed_pos','stats_speed_neg', ...
    'stats_spl_pos','stats_spl_neg', ...
    'stats_ms_pos','stats_ms_neg', ...
    'winPos','winNeg','subjects','doSSD','doPrewhite','nPerm','alphaThr');

%% FUNCTIONS
%%

function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
Ts = 1/fs; L = numel(X);
framelen = min(21, L); if mod(framelen,2)==0, framelen=framelen-1; end
minLegal = polyOrd+3; if mod(minLegal,2)==0, minLegal=minLegal+1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L, framelen = L - mod(L,2) + 1; end
useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1)/(Ts^1)) * G(:,2)';    % 1st-derivative kernel
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X)*fs;
    vy = gradient(Y)*fs;
end
end

function out = bin_series(x, fs, binDur, how)
if iscolumn(x), x = x.'; end
binSamp = round(binDur * fs);
nb = floor(numel(x)/binSamp);
out = nan(1, nb);
for b = 1:nb
    i1 = (b-1)*binSamp + 1;
    i2 = b*binSamp;
    seg = x(i1:i2);
    switch lower(how)
        case 'mean', out(b) = mean(seg,'omitnan');
        case 'sum',  out(b) = nansum(seg);
        otherwise, error('how must be ''mean'' or ''sum''.');
    end
end
end

function ms = ms_detect_engbert(X, Y, fs, lambda, minDurSec, minISISec)
if nargin < 4 || isempty(lambda),    lambda = 6;      end
if nargin < 5 || isempty(minDurSec), minDurSec = 0.006; end
if nargin < 6 || isempty(minISISec), minISISec = 0.02;  end

[vx, vy] = compute_velocity_sg(X, Y, fs, 3);

% Robust median-based SD
sx = 1.4826 * median(abs(vx - median(vx,'omitnan')),'omitnan');
sy = 1.4826 * median(abs(vy - median(vy,'omitnan')),'omitnan');
tx = lambda * max(sx, eps); ty = lambda * max(sy, eps);

rad2 = (vx./tx).^2 + (vy./ty).^2;
cand = isfinite(rad2) & (rad2 > 1);

d = diff([0; cand(:); 0]);
on  = find(d==1);
off = find(d==-1) - 1;

minDur = max(1, round(minDurSec * fs));
keep = (off - on + 1) >= minDur;
on  = on(keep); off = off(keep);

% Merge close events
minISI = max(1, round(minISISec * fs));
if ~isempty(on)
    O = on(1); F = off(1);
    onM = []; offM = [];
    for i = 2:numel(on)
        if on(i) - F <= minISI
            F = off(i);
        else
            onM = [onM; O]; offM = [offM; F]; %#ok<AGROW>
            O = on(i); F = off(i);
        end
    end
    onM = [onM; O]; offM = [offM; F];
    on = onM; off = offM;
end

ms.onsets = on;
ms.offsets = off;
end

function y = ar1_prewhiten(x)
% Remove AR(1) autocorrelation via Yule–Walker estimate
x = x(:).';
if any(isnan(x))
    x = fillmissing(x,'linear','EndValues','nearest');
end
x = x - mean(x,'omitnan');
% lag-1 autocorr
r0 = sum(x.^2); r1 = sum(x(1:end-1).*x(2:end));
phi = 0;
if r0 > 0, phi = r1 / r0; end
phi = max(min(phi, 0.99), -0.99);
% inverse filter: y_t = x_t - phi*x_{t-1}
y = x; y(2:end) = x(2:end) - phi*x(1:end-1);
end

function comp = ssd_alpha_component(ft_data, roi_labels, bandAlpha, flank1, flank2)
% Compute a single SSD component from ROI channels and return it as FT data struct
cfg = []; cfg.channel = roi_labels;
roi = ft_selectdata(cfg, ft_data);

% Build long matrices across trials
X = []; % channels × time
for k = 1:numel(roi.trial)
    X = [X, roi.trial{k}]; %#ok<AGROW>
end

% Band-pass filter copies for covariance estimation
cfgB = []; cfgB.bpfilter = 'yes'; cfgB.demean = 'yes';
cfgB.bpfreq   = bandAlpha;
roiB = ft_preprocessing(cfgB, roi);
cfgN1 = cfgB; cfgN1.bpfreq = flank1;
cfgN2 = cfgB; cfgN2.bpfreq = flank2;
roiN1 = ft_preprocessing(cfgN1, roi);
roiN2 = ft_preprocessing(cfgN2, roi);

XB = []; XN1 = []; XN2 = [];
for k = 1:numel(roi.trial)
    XB  = [XB,  roiB.trial{k}];  %#ok<AGROW>
    XN1 = [XN1, roiN1.trial{k}]; %#ok<AGROW>
    XN2 = [XN2, roiN2.trial{k}]; %#ok<AGROW>
end

CB  = cov(XB.');             % band covariance
CN  = cov([XN1, XN2].');     % noise covariance (flanks)

% Generalised eigen-decomposition
[W, D] = eig(CB, CN);
[~, ix] = sort(diag(D), 'descend');
W = W(:, ix);

% Project trials onto the first SSD component
proj = cell(size(roi.trial));
w1 = W(:,1);
for k = 1:numel(roi.trial)
    proj{k} = w1.' * roi.trial{k};
end

% Return as FT data (1 channel)
comp = [];
comp.label = {'SSD-Alpha-1'};
comp.fsample = roi.fsample;
comp.trial = proj;
comp.time  = roi.time;
end

function [mu, sem, subjMeans] = group_curve(CCF_allSubs)
% trial-pooled mean/SEM for the display; plus subject-balanced means for stats
CCF = cat(1, CCF_allSubs{:});
CCF = CCF(any(isfinite(CCF),2), :);
mu  = nanmean(CCF,1);
nPL = sum(isfinite(CCF),1);
sem = nanstd(CCF,[],1) ./ max(1, sqrt(nPL));

% Subject-balanced: mean across trials per subject (rows = subjects)
nS = numel(CCF_allSubs);
subjMeans = nan(nS, size(CCF,2));
for s = 1:nS
    M = CCF_allSubs{s};
    if isempty(M), continue; end
    subjMeans(s,:) = mean(M, 1, 'omitnan');
end
end

function [clusters, tvals, thr] = cluster_permutation_1d(S, nPerm, alpha)
% S: subjects × lags (subject means over trials)
% One-sample t-test vs 0 at each lag, cluster by |t| > tcrit, sign preserved.
% Permutation: random sign-flip per subject (good null for mean=0).

[ns, L] = size(S);
% t-values
se = std(S,[],1,'omitnan') ./ sqrt(sum(isfinite(S),1));
m  = mean(S,1,'omitnan');
tvals = m ./ se;

% Threshold from t-distribution
tcrit = tinv(1 - 0.5*alpha, ns-1); % two-sided
thr.tcrit = tcrit;

% Observed clusters
clusters = compute_clusters(tvals, tcrit);

% Permutation null (max cluster mass)
maxMass = zeros(1,nPerm);
for p = 1:nPerm
    flips = (rand(ns,1) > 0.5)*2 - 1;
    Sprm  = S .* flips;
    seP = std(Sprm,[],1,'omitnan') ./ sqrt(sum(isfinite(Sprm),1));
    mP  = mean(Sprm,1,'omitnan');
    tP  = mP ./ seP;
    clP = compute_clusters(tP, tcrit);
    if isempty(clP), maxMass(p) = 0; else, maxMass(p) = max([clP.mass]); end
end
maxMass = sort(maxMass);
thr.mass = maxMass(round((1-alpha)*nPerm));

% Assign cluster p-values
for k = 1:numel(clusters)
    clusters(k).p = 1 - mean(maxMass <= clusters(k).mass);
end
end

function clusters = compute_clusters(tvals, tcrit)
L = numel(tvals);
above = abs(tvals) > tcrit;
clusters = struct('idx',{},'mass',{},'sign',{});
if ~any(above), return; end
d = diff([0, above, 0]);
on  = find(d==1);
off = find(d==-1) - 1;
for c = 1:numel(on)
    idx = on(c):off(c);
    sgn = sign(mean(tvals(idx)));
    clusters(end+1).idx  = idx; %#ok<AGROW>
    clusters(end).mass   = sum(abs(tvals(idx)));
    clusters(end).sign   = sgn;
end
end

function plot_tcurve_with_clusters(lags_ms, tvals, clusters, thr, name, col, outdir, fontSize)
figure('Position',[0 0 1200 600]);
plot(lags_ms, tvals, 'LineWidth', 2, 'Color', col); hold on
yline(thr.tcrit,'k--'); yline(-thr.tcrit,'k--'); xline(0,'k:');
xlabel('Lag [ms] (positive: gaze follows \alpha)')
ylabel('t-statistic (one-sample vs 0)')
title(['Per-lag t: ', name], 'FontSize', fontSize)
set(gca,'FontSize',fontSize)

% Shade significant clusters at cluster-mass level
for k = 1:numel(clusters)
    if clusters(k).mass >= thr.mass
        idx = clusters(k).idx;
        xx = [lags_ms(idx), fliplr(lags_ms(idx))];
        yy = [zeros(size(idx)), fliplr(tvals(idx))];
        patch(xx, yy, col, 'FaceAlpha', 0.2, 'EdgeColor','none');
    end
end
saveas(gcf, fullfile(outdir, ['AOC_alpha_', lower(name), '_tcurve_clusters.png']));
end

function stats = window_stats(S, mask)
% S: subjects × lags; mask: logical over lags
winVals = mean(S(:,mask), 2, 'omitnan');   % mean across lags per subject
z = atanh(winVals);                        % Fisher-z
[~, p, ~, T] = ttest(z);                   % one-sample vs 0
muZ = mean(z); seZ = std(z)/sqrt(numel(z));
ciZ = [muZ - 1.96*seZ, muZ + 1.96*seZ];
stats.r_mean   = tanh(muZ);
stats.r_CI     = tanh(ciZ);
stats.t        = T.tstat;
stats.df       = T.df;
stats.p        = p;
stats.per_pos  = mean(winVals>0);
stats.values_r = winVals;
end

function plot_window_stats(stats, titleStr, col, outdir, fontSize)
vals = stats.values_r(:);
figure('Position',[0 0 700 700]); hold on
% jittered points
x = 0.9 + 0.2*rand(size(vals));
plot(x, vals, 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');
% mean + CI
plot([1 1], stats.r_CI, 'k-', 'LineWidth', 3)
plot(1, stats.r_mean, 's', 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerEdgeColor','k')
yline(0,'k:')
xlim([0.5 1.5])
ylabel('Mean r (Fisher-backtransformed)')
title(sprintf('%s\nr=%.3f, 95%%CI=[%.3f %.3f], t(%d)=%.2f, p=%.4f', ...
    titleStr, stats.r_mean, stats.r_CI(1), stats.r_CI(2), stats.df, stats.t, stats.p), ...
    'FontSize', fontSize-6)
set(gca,'XTick',[], 'FontSize', fontSize-4)
saveas(gcf, fullfile(outdir, ['AOC_alpha_window_', regexprep(lower(titleStr),'[ ()+]','_'), '.png']));
end
