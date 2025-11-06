%% AOC Sternberg — α-envelope ↔ gaze-speed CCF
% Output: .../AOC_sternberg_alpha_velocity_CCFs_raw.png

%% Setup
startup
[subjects, path, colors, ~] = setup('AOC');
subjects  = subjects(1:10)
set(0,'DefaultFigureColor','w')
fontSize  = 26;

fs        = 500;               % Hz
epochWin  = [-.5 2.5];         % padding for filtering
anaWin    = [0 2];             % retention window
maxLagS   = 0.5;               % seconds (+/-)
maxLag    = round(maxLagS*fs); % samples
velZthr   = 4;                 % |z|>4 → interpolate

outdir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/interactions';
if ~exist(outdir,'dir'); mkdir(outdir); end

%% Electrodes
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
all_lags = -maxLag:maxLag;
lags_ms  = (all_lags/fs)*1000;
CCF_S_allSubs = {};             % per-subject [nTrials × nLags]

%% Loop subjects → α-envelope, gaze speed, trial-wise CCF
for s = 1:numel(subjects)
    clc; fprintf('Subject %s (%d/%d)\n', subjects{s}, s, numel(subjects));

    dpeeg  = fullfile(path, subjects{s}, 'eeg');
    dpgaze = fullfile(path, subjects{s}, 'gaze');

    % EEG
    cd(dpeeg)
    load dataEEG_TFR_sternberg

    % Subject IAF
    IAF_file = fullfile(dpeeg, 'IAF_sternberg_subject.mat');
    if exist(IAF_file, 'file')
        tmp = load(IAF_file, 'IAF_subj'); IAF_subj = tmp.IAF_subj;
    else
        IAF_subj = NaN;
    end
    if isfinite(IAF_subj), bandAlpha = [IAF_subj-4, IAF_subj+2]; else, bandAlpha = [8 14]; end

    % Eye-tracking
    cd(dpgaze)
    load dataET_sternberg dataETlong
    trialinfo = dataETlong.trialinfo;

    % Virtual occipital alpha envelope
    cfg = []; cfg.channel = roi_labels;          roi = ft_selectdata(cfg, dataTFR);
    cfg = []; cfg.avgoverchan = 'yes';           roiVS = ft_selectdata(cfg, roi);
    cfg = []; cfg.bpfilter='yes'; cfg.bpfreq=bandAlpha; cfg.demean='yes'; cfg.hilbert='complex';
    roiH = ft_preprocessing(cfg, roiVS);
    roiH = ft_selectdata(struct('latency',epochWin), roiH);
    roiH = ft_selectdata(struct('latency',anaWin),   roiH);

    nT_env = numel(roiH.trial);
    env    = cell(1,nT_env);
    for k = 1:nT_env
        env{k} = abs(roiH.trial{k}(1,:));
    end

    % Right-eye gaze speed (derivative of X,Y)
    cfg = []; cfg.channel = {'L-GAZE-X','L-GAZE-Y','L-AREA'};
    et  = ft_selectdata(cfg, dataETlong);
    et  = ft_selectdata(struct('latency',epochWin), et);
    et  = ft_selectdata(struct('latency',anaWin),   et);

    nT_et   = numel(et.trial);
    vspd    = cell(1,nT_et);
    valid   = true(nT_et,1);

    polyOrd = 3; Ts = 1/fs;
    for k = 1:nT_et
        X = double(et.trial{k}(1,:));
        Y = double(et.trial{k}(2,:));
        A = double(et.trial{k}(3,:));

        % Fill blinks/invalid area by linear interpolation with edge padding
        blink = ~isfinite(A) | (A<=0);
        if any(blink)
            X(blink)=NaN; Y(blink)=NaN;
            if isnan(X(1)),  i1=find(isfinite(X),1,'first'); if ~isempty(i1), X(1:i1-1)=X(i1); end, end
            if isnan(X(end)), i2=find(isfinite(X),1,'last');  if ~isempty(i2), X(i2+1:end)=X(i2); end, end
            if isnan(Y(1)),  i1=find(isfinite(Y),1,'first'); if ~isempty(i1), Y(1:i1-1)=Y(i1); end, end
            if isnan(Y(end)), i2=find(isfinite(Y),1,'last');  if ~isempty(i2), Y(i2+1:end)=Y(i2); end, end
            X = fillmissing(X,'linear');  Y = fillmissing(Y,'linear');
        end

        % Savitzky–Golay smoothing + derivative (fallback to gradient)
        L = numel(X);
        framelen = min(21, L); if mod(framelen,2)==0, framelen=framelen-1; end
        minLegal = polyOrd+3;  if mod(minLegal,2)==0, minLegal=minLegal+1; end
        if framelen < minLegal, framelen = minLegal; end
        if framelen > L, framelen = L - mod(L,2) + 1; end
        useFallback = framelen < 5;

        if ~useFallback
            Xs = sgolayfilt(X, polyOrd, framelen);
            Ys = sgolayfilt(Y, polyOrd, framelen);
            [~, G] = sgolay(polyOrd, framelen);
            d1 = (factorial(1)/(Ts^1)) * G(:,2)';
            vx = conv(Xs, d1, 'same');
            vy = conv(Ys, d1, 'same');
        else
            vx = gradient(X)*fs;
            vy = gradient(Y)*fs;
        end

        % Robust outlier handling on velocities
        zvx = (vx - nanmean(vx)) / nanstd(vx);
        zvy = (vy - nanmean(vy)) / nanstd(vy);
        bad = abs(zvx)>velZthr | abs(zvy)>velZthr;
        if any(bad)
            vx(bad)=NaN; vy(bad)=NaN;
            if isnan(vx(1)),  i1=find(isfinite(vx),1,'first'); if ~isempty(i1), vx(1:i1-1)=vx(i1); end, end
            if isnan(vx(end)), i2=find(isfinite(vx),1,'last');  if ~isempty(i2), vx(i2+1:end)=vx(i2); end, end
            if isnan(vy(1)),  i1=find(isfinite(vy),1,'first'); if ~isempty(i1), vy(1:i1-1)=vy(i1); end, end
            if isnan(vy(end)), i2=find(isfinite(vy),1,'last');  if ~isempty(i2), vy(i2+1:end)=vy(i2); end, end
            vx = fillmissing(vx,'linear'); vy = fillmissing(vy,'linear');
        end

        vspd{k} = hypot(vx, vy);
        if any(~isfinite(vx)) || any(~isfinite(vy)), valid(k)=false; end
    end

    % Align counts
    nT = min([numel(env), numel(vspd), size(trialinfo,1)]);
    env   = env(1:nT);
    vspd  = vspd(1:nT);
    valid = valid(1:nT);

    % Trial-wise CCFs (α-envelope ↔ speed)
    ccf_spd = nan(nT, numel(all_lags));
    for tr = 1:nT
        if ~valid(tr), continue; end
        a  = zscore(double(env{tr}));
        sp = zscore(double(vspd{tr}));

        if any(isnan(a))
            if isnan(a(1)),  i1=find(isfinite(a),1,'first');  if ~isempty(i1), a(1:i1-1)=a(i1); end, end
            if isnan(a(end)), i2=find(isfinite(a),1,'last');   if ~isempty(i2), a(i2+1:end)=a(i2); end, end
            a = fillmissing(a,'linear');
        end
        if any(isnan(sp))
            if isnan(sp(1)),  i1=find(isfinite(sp),1,'first'); if ~isempty(i1), sp(1:i1-1)=sp(i1); end, end
            if isnan(sp(end)), i2=find(isfinite(sp),1,'last');  if ~isempty(i2), sp(i2+1:end)=sp(i2); end, end
            sp = fillmissing(sp,'linear');
        end

        tmp = xcorr(a, sp, maxLag, 'coeff');
        ccf_spd(tr,:) = tmp(:).';
    end

    CCF_S_allSubs{end+1} = ccf_spd; %#ok<SAGROW>
end

%% Group-level CCF (raw) and figure
CCF_S = cat(1, CCF_S_allSubs{:});
CCF_S = CCF_S(any(isfinite(CCF_S),2),:);

close all
figure('Position',[0 0 1200 1000]);
mu      = nanmean(CCF_S,1);
nPerLag = sum(isfinite(CCF_S),1);
sem     = nanstd(CCF_S,[],1) ./ max(1, sqrt(nPerLag));

fill([lags_ms, fliplr(lags_ms)], [mu-sem, fliplr(mu+sem)], 0.85*[1 1 1], ...
    'FaceAlpha',0.35, 'EdgeColor','none'); hold on
plot(lags_ms, mu, 'LineWidth', 3, 'Color', colors(3,:))
yline(0,'k-'); xline(0,'k--','LineWidth',1.5)
xlabel('Lag [ms] (positive: gaze follows α)')
ylabel('Correlation (r)')
xticks(-500:50:500)
ylim(1.25*[-max(abs(mu)) max(abs(mu))])
title('CCF: Alpha Envelope × Gaze Speed','FontSize',fontSize)
set(gca,'FontSize',fontSize)

saveas(gcf, fullfile(outdir, 'AOC_sternberg_alpha_velocity_CCFs_raw.png'));