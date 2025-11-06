%% AOC Sternberg — α-envelope ↔ oculomotor couplings (left eye)
% Figures:
%   1) .../AOC_sternberg_alpha_velocity_CCFs_raw.png
%   2) .../AOC_sternberg_alpha_scanpath_CCFs_raw.png
%   3) .../AOC_sternberg_alpha_microsaccades_CCFs_raw.png

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');
% subjects  = subjects(1:10);

set(0,'DefaultFigureColor','w')
fontSize  = 26;

fs        = 500;                 % Hz (EEG & ET native rate)
epochWin  = [-.5 2.5];           % padding for filtering
anaWin    = [0 2];               % retention window
maxLagS   = 0.5;                 % seconds (+/-) for CCF

% 50 ms binning for SPL & MS (20 Hz)
binDur    = 0.05;                % s
fsb       = 1/binDur;            % 20 Hz
maxLag_500= round(maxLagS * fs);   % samples @ 500 Hz
maxLag_20 = round(maxLagS * fsb);  % samples @ 20 Hz

velZthr   = 4;                   % |z|>4 → interpolate
polyOrd   = 3;                   % SG polynomial order

outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions';
if ~exist(outdir,'dir'); mkdir(outdir); end

%% Occipital ROI from labels containing 'O' or 'I'
elecPath = fullfile(path, subjects{1}, 'eeg');
cd(elecPath);
load('power_stern_late_trials.mat'); % to get labels
allLabs = powload2_late.label;
roi_labels = {};
for i = 1:numel(allLabs)
    L = allLabs{i};
    if contains(L, {'O'}) || contains(L, {'I'})
        roi_labels{end+1} = L; %#ok<SAGROW>
    end
end

%% Holders
CCF_speed_allSubs   = {};   % per-subject [nTrials × nLags(500 Hz)]
CCF_spl_allSubs     = {};   % per-subject [nTrials × nLags(20 Hz)]
CCF_ms_allSubs      = {};   % per-subject [nTrials × nLags(20 Hz)]

lags_500   = -maxLag_500:maxLag_500;
lags_20    = -maxLag_20:maxLag_20;
lags_ms_500= (lags_500/fs)*1000;
lags_ms_20 = (lags_20/fsb)*1000;

%% Loop subjects
for s = 1:numel(subjects)
    clc; fprintf('Subject %s (%d/%d)\n', subjects{s}, s, numel(subjects));

    dpeeg  = fullfile(path, subjects{s}, 'eeg');
    dpgaze = fullfile(path, subjects{s}, 'gaze');

    %% EEG
    cd(dpeeg)
    load dataEEG_TFR_sternberg

    % Subject IAF
    IAF_file = fullfile(dpeeg, 'IAF_sternberg_subject.mat');
    if exist(IAF_file, 'file')
        tmp = load(IAF_file, 'IAF_subj');
        IAF_subj = tmp.IAF_subj;
    else
        IAF_subj = NaN;
    end
    if isfinite(IAF_subj)
        bandAlpha = [IAF_subj-4, IAF_subj+2];
    else
        bandAlpha = [8 14];
    end

    % Virtual occipital alpha envelope (trial-aligned, [0 2] s)
    cfg = []; cfg.channel = roi_labels;
    roi = ft_selectdata(cfg, dataTFR);

    cfg = []; cfg.avgoverchan = 'yes';
    roiVS = ft_selectdata(cfg, roi);

    cfg = []; cfg.bpfilter='yes'; cfg.bpfreq=bandAlpha;
    cfg.demean='yes'; cfg.hilbert='complex';
    roiH = ft_preprocessing(cfg, roiVS);

    roiH = ft_selectdata(struct('latency',epochWin), roiH);
    roiH = ft_selectdata(struct('latency',anaWin),   roiH);

    nT_env = numel(roiH.trial);
    aenv_500 = cell(1,nT_env);          % 500 Hz envelope
    aenv_20  = cell(1,nT_env);          % 20 Hz binned means
    for k = 1:nT_env
        ae = abs(roiH.trial{k}(1,:));   % amplitude
        aenv_500{k} = ae;

        % Bin α envelope to 50 ms means for SPL/MS analyses
        aenv_20{k} = bin_series(ae, fs, binDur, 'mean'); % 1×Nbins
    end

    %% Eye-tracking (left eye)
    cd(dpgaze)
    load dataET_sternberg dataETlong

    % Select left-eye channels, pad, crop to [0 2] s
    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y','L-AREA'};
    et  = ft_selectdata(cfg, dataETlong);
    et  = ft_selectdata(struct('latency',epochWin), et);
    et  = ft_selectdata(struct('latency',anaWin),   et);

    nT = min([numel(aenv_500), numel(et.trial), size(dataETlong.trialinfo,1)]);
    aenv_500 = aenv_500(1:nT);
    aenv_20  = aenv_20(1:nT);

    %% Per-trial gaze transforms
    ccf_speed = nan(nT, numel(lags_500));
    ccf_spl   = nan(nT, numel(lags_20));
    ccf_ms    = nan(nT, numel(lags_20));

    for tr = 1:nT
        X = double(et.trial{tr}(1,:));      % px
        Y = double(et.trial{tr}(2,:));      % px
        A = double(et.trial{tr}(3,:));      % area

        % Invert Y to Cartesian up (differences invariant, but we do it for consistency)
        Y = -Y;

        % Blink removal by AREA + linear interpolation with edge padding
        blink = ~isfinite(A) | (A<=0);
        if any(blink)
            X(blink) = NaN; Y(blink) = NaN;
            % edge padding
            if isnan(X(1))
                i1 = find(isfinite(X),1,'first');
                if ~isempty(i1), X(1:i1-1) = X(i1); end
            end
            if isnan(X(end))
                i2 = find(isfinite(X),1,'last');
                if ~isempty(i2), X(i2+1:end) = X(i2); end
            end
            if isnan(Y(1))
                i1 = find(isfinite(Y),1,'first');
                if ~isempty(i1), Y(1:i1-1) = Y(i1); end
            end
            if isnan(Y(end))
                i2 = find(isfinite(Y),1,'last');
                if ~isempty(i2), Y(i2+1:end) = Y(i2); end
            end
            X = fillmissing(X,'linear');
            Y = fillmissing(Y,'linear');
        end

        % --- Fig 1: speed @ 500 Hz (Savitzky–Golay derivative with fallback)
        [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd);
        % Robust outlier handling on velocities
        zvx = (vx - nanmean(vx)) / nanstd(vx);
        zvy = (vy - nanmean(vy)) / nanstd(vy);
        bad = abs(zvx)>velZthr | abs(zvy)>velZthr;
        if any(bad)
            vx(bad) = NaN; vy(bad) = NaN;
            if isnan(vx(1))
                i1 = find(isfinite(vx),1,'first');
                if ~isempty(i1), vx(1:i1-1) = vx(i1); end
            end
            if isnan(vx(end))
                i2 = find(isfinite(vx),1,'last');
                if ~isempty(i2), vx(i2+1:end) = vx(i2); end
            end
            if isnan(vy(1))
                i1 = find(isfinite(vy),1,'first');
                if ~isempty(i1), vy(1:i1-1) = vy(i1); end
            end
            if isnan(vy(end))
                i2 = find(isfinite(vy),1,'last');
                if ~isempty(i2), vy(i2+1:end) = vy(i2); end
            end
            vx = fillmissing(vx,'linear');
            vy = fillmissing(vy,'linear');
        end
        spd = hypot(vx, vy);               % 1×T

        ae = aenv_500{tr};
        % Guard against residual NaNs by local interpolation
        if any(isnan(ae)), ae = fillmissing(ae,'linear','EndValues','nearest'); end
        if any(isnan(spd)), spd= fillmissing(spd,'linear','EndValues','nearest'); end

        aZ  = zscore(ae);
        sZ  = zscore(spd);
        ccf_speed(tr,:) = xcorr(aZ, sZ, maxLag_500, 'coeff');

        % --- Fig 2: scan-path length in 50 ms bins (sum of steps per bin)
        step = sqrt(diff(X).^2 + diff(Y).^2);     % per-sample step (px)
        step = [step, step(end)];                 % pad to length T
        spl_bins = bin_series(step, fs, binDur, 'sum'); % 1×Nbins

        % Match α to 50 ms means
        ae20 = aenv_20{tr};
        if numel(ae20) ~= numel(spl_bins)
            nb = min(numel(ae20), numel(spl_bins));
            ae20 = ae20(1:nb);
            spl_bins = spl_bins(1:nb);
        end
        aZb = zscore(ae20);
        splZ= zscore(spl_bins);
        ccf_spl(tr,:) = xcorr(aZb, splZ, maxLag_20, 'coeff');

        % --- Fig 3: microsaccade rate in 50 ms bins (events/s)
        % Use Engbert detector (left eye, px)
        ms = ms_detect_engbert(X, Y, fs, 6, 6e-3, 0.02); % λ=6, minDur=6 ms, minISI=20 ms
        ms_on = zeros(size(X));
        if ~isempty(ms.onsets)
            on = ms.onsets;
            on = on(on>=1 & on<=numel(ms_on));
            ms_on(on) = 1;
        end
        ms_per_s = bin_series(ms_on, fs, binDur, 'sum') / binDur; % events/s per 50 ms bin

        nb = min(numel(ae20), numel(ms_per_s));
        aZb = zscore(ae20(1:nb));
        msZ = zscore(ms_per_s(1:nb));
        ccf_ms(tr,:) = xcorr(aZb, msZ, maxLag_20, 'coeff');
    end

    CCF_speed_allSubs{end+1} = ccf_speed; %#ok<SAGROW>
    CCF_spl_allSubs{end+1}   = ccf_spl;   %#ok<SAGROW>
    CCF_ms_allSubs{end+1}    = ccf_ms;    %#ok<SAGROW>
end

%% Helper to compute group mean ± SEM
compute_group = @(M) deal( ...
    nanmean(M,1), ...
    nanstd(M,[],1) ./ max(1, sqrt(sum(isfinite(M),1))) );

%% --- Fig 1: α × speed (500 Hz)
CCF_speed = cat(1, CCF_speed_allSubs{:});
CCF_speed = CCF_speed(any(isfinite(CCF_speed),2), :);
[mu, sem] = compute_group(CCF_speed);

close all
figure('Position',[0 0 1200 1000]);
fill([lags_ms_500, fliplr(lags_ms_500)], [mu-sem, fliplr(mu+sem)], 0.85*[1 1 1], ...
    'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_500, mu, 'LineWidth', 3, 'Color', colors(3,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)')
ylabel('Correlation (r)')
xticks(-500:50:500)
ylim(1.25*[-max(abs(mu)) max(abs(mu))])
title('CCF: Alpha Envelope \times Gaze Speed','FontSize',fontSize)
set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_sternberg_alpha_velocity_CCFs_raw.png'));

%% --- Fig 2: α × scan-path length (50 ms bins)
CCF_spl = cat(1, CCF_spl_allSubs{:});
CCF_spl = CCF_spl(any(isfinite(CCF_spl),2), :);
[mu2, sem2] = compute_group(CCF_spl);

figure('Position',[0 0 1200 1000]);
fill([lags_ms_20, fliplr(lags_ms_20)], [mu2-sem2, fliplr(mu2+sem2)], 0.85*[1 1 1], ...
    'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_20, mu2, 'LineWidth', 3, 'Color', colors(2,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)')
ylabel('Correlation (r)')
xticks(-500:50:500)
ylim(1.25*[-max(abs(mu2)) max(abs(mu2))])
title('CCF: Alpha Envelope \times Scan-Path Length (50 ms bins)','FontSize',fontSize)
set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_sternberg_alpha_scanpath_CCFs_raw.png'));

%% --- Fig 3: α × microsaccade rate (50 ms bins)
CCF_ms = cat(1, CCF_ms_allSubs{:});
CCF_ms = CCF_ms(any(isfinite(CCF_ms),2), :);
[mu3, sem3] = compute_group(CCF_ms);

figure('Position',[0 0 1200 1000]);
fill([lags_ms_20, fliplr(lags_ms_20)], [mu3-sem3, fliplr(mu3+sem3)], 0.85*[1 1 1], ...
    'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms_20, mu3, 'LineWidth', 3, 'Color', colors(1,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows \alpha)')
ylabel('Correlation (r)')
xticks(-500:50:500)
ylim(1.25*[-max(abs(mu3)) max(abs(mu3))])
title('CCF: Alpha Envelope \times Microsaccade Rate (50 ms bins)','FontSize',fontSize)
set(gca,'FontSize',fontSize)
saveas(gcf, fullfile(outdir, 'AOC_sternberg_alpha_microsaccades_CCFs_raw.png'));


%% FUNCTIONS
%%


function [vx, vy] = compute_velocity_sg(X, Y, fs, polyOrd)
% Savitzky–Golay smoothing + 1st derivative (fallback: gradient)
% X,Y: 1×T; returns velocities in px/s

Ts = 1/fs;
L  = numel(X);

% choose a reasonable odd frame length
framelen = min(21, L);
if mod(framelen,2)==0, framelen = framelen - 1; end
minLegal = polyOrd + 3;
if mod(minLegal,2)==0, minLegal = minLegal + 1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L, framelen = L - mod(L,2) + 1; end

useFallback = framelen < 5;

if ~useFallback
    Xs = sgolayfilt(X, polyOrd, framelen);
    Ys = sgolayfilt(Y, polyOrd, framelen);
    [~, G] = sgolay(polyOrd, framelen);
    d1 = (factorial(1)/(Ts^1)) * G(:,2)';  % 1st derivative kernel
    vx = conv(Xs, d1, 'same');
    vy = conv(Ys, d1, 'same');
else
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
end
end

function out = bin_series(x, fs, binDur, how)
% Bin a 1×T series 'x' into contiguous bins of length binDur (s)
% how: 'mean' or 'sum'
if iscolumn(x), x = x.'; end
binSamp = round(binDur * fs);
nb = floor(numel(x)/binSamp);
out = nan(1, nb);
for b = 1:nb
    idx1 = (b-1)*binSamp + 1;
    idx2 = b*binSamp;
    seg  = x(idx1:idx2);
    switch lower(how)
        case 'mean'
            out(b) = mean(seg,'omitnan');
        case 'sum'
            out(b) = nansum(seg);
        otherwise
            error('how must be ''mean'' or ''sum''.');
    end
end
end

function ms = ms_detect_engbert(X, Y, fs, lambda, minDurSec, minISISec)
% Minimal Engbert-style microsaccade detector on single-eye data in px.
% Returns onsets, offsets (sample indices).
% Velocity thresholds computed with median-based SD estimator.
% Parameters:
%   lambda     - threshold multiplier (typ. 5–6)
%   minDurSec  - minimum duration (s), e.g. 0.006
%   minISISec  - merge events closer than this ISI (s)

if nargin < 4 || isempty(lambda),    lambda = 6;      end
if nargin < 5 || isempty(minDurSec), minDurSec = 0.006; end
if nargin < 6 || isempty(minISISec), minISISec = 0.02;  end

% Derivatives (px/s)
[vx, vy] = compute_velocity_sg(X, Y, fs, 3);

% Median-based SD (robust)
lambda_x = lambda * 1.4826 * median(abs(vx - median(vx,'omitnan')),'omitnan');
lambda_y = lambda * 1.4826 * median(abs(vy - median(vy,'omitnan')),'omitnan');

% Elliptical threshold
velOK = isfinite(vx) & isfinite(vy);
ex = zeros(size(vx));
ey = zeros(size(vy));
if lambda_x > 0, ex = vx ./ lambda_x; end
if lambda_y > 0, ey = vy ./ lambda_y; end
rad2 = ex.^2 + ey.^2;

cand = velOK & (rad2 > 1);
cand = cand(:);

% Find runs > minDur
minDur = max(1, round(minDurSec * fs));
minISI = max(1, round(minISISec * fs));

d = diff([0; cand; 0]);
on = find(d == 1);
off= find(d == -1) - 1;

dur = off - on + 1;
keep = dur >= minDur;
on  = on(keep);
off = off(keep);

% Merge events closer than minISI
if ~isempty(on)
    merged_on  = on(1);
    merged_off = off(1);
    out_on  = [];
    out_off = [];
    for i = 2:numel(on)
        if on(i) - merged_off <= minISI
            merged_off = off(i);
        else
            out_on  = [out_on;  merged_on]; %#ok<AGROW>
            out_off = [out_off; merged_off]; %#ok<AGROW>
            merged_on  = on(i);
            merged_off = off(i);
        end
    end
    out_on  = [out_on;  merged_on];
    out_off = [out_off; merged_off];
    on  = out_on;
    off = out_off;
end

ms.onsets  = on;
ms.offsets = off;
end
