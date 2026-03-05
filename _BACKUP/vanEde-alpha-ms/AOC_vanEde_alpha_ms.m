%% AOC_vanEde_alpha_ms.m
% =========================================================================
% Replication of figures from:
%   Liu, Nobre & van Ede (2023). "Microsaccades transiently lateralise
%   EEG alpha activity." Progress in Neurobiology, 224, 102433.
%
% Applied to AOC dataset: Sternberg and N-back tasks.
% Creates Figures 1-4 for both tasks (8 figures total).
%
% Figures:
%   1: Gaze velocity, TFR lateralisation, raw + baselined alpha TC, topographies
%   2: Start vs Return MS: gaze position, TFR lat., raw + baselined alpha TC, topographies
%   3: Ipsilateral vs Contralateral power (start MS only)
%   4: Inter-trial phase coherence / ITPC (start MS only)
%
% Alpha metrics:
%   - Raw lateralisation: ((contra - ipsi) / (contra + ipsi)) * 100
%   - Baselined lateralisation: raw minus pre-MS baseline mean [-500, -150 ms]
%
% Requirements:
%   FieldTrip toolbox, AOC data on network drive
%
% Output:
%   Figures saved to figDir (see below)
% =========================================================================

%% ========================================================================
%  1. SETUP
% ========================================================================
clc; close all;

% --- Initialise FieldTrip and get subjects/paths ---
[subjects, dataPath, ~, headmodel] = setup('AOC');
try
    layANThead = headmodel.layANThead;
catch
    if ispc
        tmp = load('W:\Students\Arne\toolboxes\headmodel\layANThead.mat');
    else
        tmp = load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/layANThead.mat');
    end
    layANThead = tmp.layANThead;
end

% --- Sampling rate ---
fs = 500;

% --- Microsaccade detection parameters (Engbert & Kliegl 2003) ---
velThreshold = 6;        % velocity threshold (multiples of median-based SD)
minMS_dur_s  = 0.006;    % minimum microsaccade duration (s)
minMS_ISI_s  = 0.020;    % minimum inter-saccade interval (s)
maxMS_mag_px = 35;       % max horizontal displacement in pixels (~1 visual deg)
minMS_hor_px = 1.7;      % min horizontal displacement (pixels, ~3.4 arcmin)

% --- Retention period for MS detection ---
msWin = [0.5, 2.0];     % detect MS in this window (s relative to trial onset)

% --- EEG epoch around microsaccade onset ---
epochLim  = [-1.2, 1.2]; % s
epochSamp = round(epochLim * fs);       % [-600, 600]
nEpochSamp = epochSamp(2) - epochSamp(1) + 1;
epochTime  = linspace(epochLim(1), epochLim(2), nEpochSamp);

% --- STFT parameters (matching van Ede 2023) ---
stft_winLen  = 0.300;    % Hanning window length (s)
stft_step    = 0.020;    % time step (s)
stft_foi     = 1:1:50;   % frequencies of interest (Hz)
stft_toi     = -1.0:stft_step:1.0; % TFR time vector (s)

stft_winSamp = round(stft_winLen * fs);  % 150 samples
hann_win     = hanning(stft_winSamp)';   % [1 x winSamp]
nFreq        = numel(stft_foi);
nTime        = numel(stft_toi);

% Zero-pad FFT for 1 Hz frequency resolution
nfft = max(stft_winSamp, fs); % 500 points -> 1 Hz bins
freq_bins = stft_foi + 1;     % MATLAB 1-indexed: 1 Hz -> bin 2, etc.

% --- Alpha band ---
alphaBand = [8 12];
alphaIdx  = stft_foi >= alphaBand(1) & stft_foi <= alphaBand(2);

% --- Key electrodes (following van Ede 2023) ---
elec_L = 'PO7';
elec_R = 'PO8';

% --- Baseline window for Figure 3 ---
blWin = [-0.5, -0.15];  % s relative to MS onset
blIdx = stft_toi >= blWin(1) & stft_toi <= blWin(2);

% --- Task definitions ---
taskName = {'Sternberg', 'Nback'};
eegFile  = {'dataEEG_TFR_sternberg', 'dataEEG_TFR_nback'};
etFile   = {'dataET_sternberg', 'dataET_nback'};

% --- Gaze velocity epoch ---
velEpochLim  = [-0.5, 0.5]; % s
velEpochSamp = round(velEpochLim * fs);
nVelSamp     = velEpochSamp(2) - velEpochSamp(1) + 1;
velTime      = linspace(velEpochLim(1), velEpochLim(2), nVelSamp) * 1000; % ms

% --- Gaze position epoch ---
posEpochLim  = [-0.5, 0.5]; % s
posEpochSamp = round(posEpochLim * fs);
nPosSamp     = posEpochSamp(2) - posEpochSamp(1) + 1;
posTime      = linspace(posEpochLim(1), posEpochLim(2), nPosSamp) * 1000; % ms

% --- Output directory ---
if ispc
    figDir = 'W:\Students\Arne\AOC\figures\tests\vanEde-alpha-ms';
else
    figDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/vanEde-alpha-ms';
end
if ~exist(figDir, 'dir'), mkdir(figDir); end

% --- Figure settings ---
set(0, 'DefaultFigureColor', 'w');
fntSz = 14;
lnWd  = 2.5;

% Diverging blue-white-red colormap
nCmap = 256;
halfC = nCmap / 2;
cmap_blue = [linspace(0.15, 1, halfC)', linspace(0.25, 1, halfC)', ones(halfC, 1)];
cmap_red  = [ones(halfC, 1), linspace(1, 0.25, halfC)', linspace(1, 0.15, halfC)'];
cmap_div  = [cmap_blue; cmap_red];

% Sequential colormaps for power/ITPC
cmap_hot = hot(256);

% TFR time axis in ms for plotting
stft_toi_ms = stft_toi * 1000;

fprintf('\n=== AOC van Ede Alpha-MS Replication ===\n');
fprintf('Subjects: %d\n', numel(subjects));
fprintf('Output: %s\n\n', figDir);

%% ========================================================================
%  2. PROCESSING LOOP (per task)
% ========================================================================
for taskIdx = 1:2
    fprintf('\n%s\n  TASK: %s\n%s\n', repmat('=',1,60), upper(taskName{taskIdx}), repmat('=',1,60));

    nSubj = numel(subjects);

    % --- Per-subject storage ---
    % Power per direction [nChan x nFreq x nTime]
    sub_pow_L  = cell(nSubj, 1); sub_pow_R  = cell(nSubj, 1);
    sub_pow_SL = cell(nSubj, 1); sub_pow_SR = cell(nSubj, 1);
    sub_pow_RL = cell(nSubj, 1); sub_pow_RR = cell(nSubj, 1);
    % ITPC per direction [nChan x nFreq x nTime]
    sub_itpc_L  = cell(nSubj, 1); sub_itpc_R  = cell(nSubj, 1);
    sub_itpc_SL = cell(nSubj, 1); sub_itpc_SR = cell(nSubj, 1);
    sub_itpc_RL = cell(nSubj, 1); sub_itpc_RR = cell(nSubj, 1);
    % Gaze epochs (velocity and position, subject-level means)
    sub_vel      = cell(nSubj, 1); % [1 x nVelSamp]
    sub_pos_SL   = cell(nSubj, 1); sub_pos_SR = cell(nSubj, 1);
    sub_pos_RL   = cell(nSubj, 1); sub_pos_RR = cell(nSubj, 1);
    % Counts and validity
    sub_nMS    = zeros(nSubj, 6); % [nL, nR, nSL, nSR, nRL, nRR]
    sub_valid  = false(nSubj, 1);
    chanLabels = {};

    % ====================================================================
    %  SUBJECT LOOP
    % ====================================================================
    for s = 1:nSubj
        fprintf('  Subject %s (%d/%d): ', subjects{s}, s, nSubj);

        try
            % --- Load EEG (FieldTrip raw) ---
            eegPath = fullfile(dataPath, subjects{s}, 'eeg');
            cd(eegPath); %#ok<*MCCD>
            load(eegFile{taskIdx}, 'dataTFR');

            % --- Load ET (FieldTrip-like) ---
            etPath = fullfile(dataPath, subjects{s}, 'gaze');
            cd(etPath);
            S = load(etFile{taskIdx});
            if isfield(S, 'dataETlong')
                dataET = S.dataETlong;
            elseif isfield(S, 'dataet')
                dataET = S.dataet;
            else
                fn = fieldnames(S); dataET = S.(fn{1});
            end
            clear S;
        catch ME
            fprintf('LOAD ERROR: %s\n', ME.message);
            continue;
        end

        % --- Verify key electrodes exist ---
        idxPO7 = find(strcmpi(dataTFR.label, elec_L), 1);
        idxPO8 = find(strcmpi(dataTFR.label, elec_R), 1);
        if isempty(idxPO7) || isempty(idxPO8)
            fprintf('PO7/PO8 not found, skipping.\n');
            clear dataTFR dataET; continue;
        end
        nChan = numel(dataTFR.label);
        if isempty(chanLabels), chanLabels = dataTFR.label; end

        % --- Match trial count ---
        nTrials = min(numel(dataTFR.trial), numel(dataET.trial));

        % ================================================================
        %  MICROSACCADE DETECTION (across all trials)
        % ================================================================
        ms_all = struct('trial',{}, 'eeg_samp',{}, 'ms_time',{}, ...
                        'dir',{}, 'mstype',{}, 'mag',{});

        for tr = 1:nTrials
            t_eeg  = dataTFR.time{tr};
            t_gaze = dataET.time{tr};

            % Gaze samples in retention window
            gIdx = t_gaze >= msWin(1) & t_gaze <= msWin(2);
            if sum(gIdx) < 100, continue; end

            X_ret = double(dataET.trial{tr}(1, gIdx));
            Y_ret = double(dataET.trial{tr}(2, gIdx));

            % Blink removal
            if size(dataET.trial{tr}, 1) >= 3
                A_ret = double(dataET.trial{tr}(3, gIdx));
                blink = ~isfinite(A_ret) | A_ret <= 0;
            else
                blink = ~isfinite(X_ret) | ~isfinite(Y_ret);
            end
            X_ret(blink) = NaN; Y_ret(blink) = NaN;
            if all(isnan(X_ret)), continue; end
            X_ret = fillmissing(X_ret, 'linear', 'EndValues', 'nearest');
            Y_ret = fillmissing(Y_ret, 'linear', 'EndValues', 'nearest');

            % Savitzky-Golay velocity
            [vx, vy] = compute_velocity_sg_local(X_ret, Y_ret, fs, 3);

            % Engbert & Kliegl threshold
            sx = 1.4826 * median(abs(vx - median(vx,'omitnan')),'omitnan');
            sy = 1.4826 * median(abs(vy - median(vy,'omitnan')),'omitnan');
            tx = velThreshold * max(sx, eps);
            ty = velThreshold * max(sy, eps);
            rad2 = (vx./tx).^2 + (vy./ty).^2;
            above = isfinite(rad2) & rad2 > 1;

            % Find onset/offset
            d = diff([0, above, 0]);
            on = find(d == 1); off = find(d == -1) - 1;

            % Duration filter
            minDurSamp = max(1, round(minMS_dur_s * fs));
            keep = (off - on + 1) >= minDurSamp;
            on = on(keep); off = off(keep);

            % Merge close events
            [on, off] = merge_close_local(on, off, max(1, round(minMS_ISI_s * fs)));

            % Retention gaze time vector and fixation center
            t_ret = t_gaze(gIdx);
            centerX = median(X_ret, 'omitnan');

            % Process each detected saccade
            for k = 1:numel(on)
                onset_r = on(k); % sample in retention window

                % Pre-saccade position (-50 to 0 ms)
                ps = max(1, onset_r - round(0.050*fs));
                preX = mean(X_ret(ps:onset_r), 'omitnan');

                % Post-saccade position (50-100 ms)
                ps1 = min(numel(X_ret), onset_r + round(0.050*fs));
                ps2 = min(numel(X_ret), onset_r + round(0.100*fs));
                if ps2 > numel(X_ret), continue; end
                postX = mean(X_ret(ps1:ps2), 'omitnan');

                % Magnitude
                mag = abs(postX - preX);
                if mag > maxMS_mag_px || mag < minMS_hor_px, continue; end

                % Direction (1=left, 2=right)
                msDir = 2; if postX < preX, msDir = 1; end

                % Start(1) vs Return(2)
                preDist  = abs(preX - centerX);
                postDist = abs(postX - centerX);
                msType = 2; if postDist > preDist, msType = 1; end

                % MS onset in trial time
                ms_time = t_ret(onset_r);

                % Find corresponding EEG sample
                [~, eeg_samp] = min(abs(t_eeg - ms_time));
                eeg_lo = eeg_samp + epochSamp(1);
                eeg_hi = eeg_samp + epochSamp(2);
                if eeg_lo < 1 || eeg_hi > size(dataTFR.trial{tr}, 2), continue; end

                ms_all(end+1) = struct('trial', tr, 'eeg_samp', eeg_samp, ...
                    'ms_time', ms_time, 'dir', msDir, 'mstype', msType, 'mag', mag); %#ok<AGROW>
            end
        end

        nMS     = numel(ms_all);
        dirs    = [ms_all.dir];
        types   = [ms_all.mstype];
        leftIdx = dirs == 1;  rightIdx = dirs == 2;
        startIdx = types == 1; retIdx   = types == 2;
        nL = sum(leftIdx); nR = sum(rightIdx);
        nSL = sum(startIdx & leftIdx); nSR = sum(startIdx & rightIdx);
        nRL = sum(retIdx & leftIdx);   nRR = sum(retIdx & rightIdx);

        fprintf('%d MS (L:%d R:%d | Start-L:%d Start-R:%d Ret-L:%d Ret-R:%d)\n', ...
            nMS, nL, nR, nSL, nSR, nRL, nRR);

        if nL < 5 || nR < 5
            fprintf('    Too few directional MS, skipping.\n');
            clear dataTFR dataET; continue;
        end

        sub_nMS(s,:) = [nL, nR, nSL, nSR, nRL, nRR];
        sub_valid(s) = true;

        % ================================================================
        %  EXTRACT EEG & GAZE EPOCHS
        % ================================================================
        eeg_epochs = zeros(nMS, nChan, nEpochSamp);
        vel_epochs = nan(nMS, nVelSamp);
        pos_epochs = nan(nMS, nPosSamp);

        for mi = 1:nMS
            tr     = ms_all(mi).trial;
            center = ms_all(mi).eeg_samp;

            % EEG epoch [nChan x nEpochSamp]
            eeg_epochs(mi, :, :) = dataTFR.trial{tr}(:, center+epochSamp(1):center+epochSamp(2));

            % Gaze epochs
            t_gaze = dataET.time{tr};
            [~, gc] = min(abs(t_gaze - ms_all(mi).ms_time));

            % Velocity epoch
            vs = gc + velEpochSamp(1); ve = gc + velEpochSamp(2);
            if vs >= 1 && ve <= numel(t_gaze)
                Xf = double(dataET.trial{tr}(1,:));
                Xf = fillmissing(Xf,'linear','EndValues','nearest');
                vel = abs([0, diff(Xf)] * fs);
                vel = smoothdata(vel, 'gaussian', max(3, round(0.007*fs)));
                vel_epochs(mi,:) = vel(vs:ve);
            end

            % Position epoch (horizontal, centred on fixation)
            ps = gc + posEpochSamp(1); pe = gc + posEpochSamp(2);
            if ps >= 1 && pe <= numel(t_gaze)
                Xf = double(dataET.trial{tr}(1,:));
                Xf = fillmissing(Xf,'linear','EndValues','nearest');
                fixCenter = median(Xf, 'omitnan');
                pos_epochs(mi,:) = Xf(ps:pe) - fixCenter;
            end
        end

        % Store subject-mean gaze velocity (all MS)
        valid_vel = ~any(isnan(vel_epochs), 2);
        if any(valid_vel)
            sub_vel{s} = mean(vel_epochs(valid_vel,:), 1, 'omitnan');
        end

        % Store subject-mean gaze position (per direction x type)
        if nSL > 0
            sub_pos_SL{s} = mean(pos_epochs(startIdx & leftIdx, :), 1, 'omitnan');
        end
        if nSR > 0
            sub_pos_SR{s} = mean(pos_epochs(startIdx & rightIdx, :), 1, 'omitnan');
        end
        if nRL > 0
            sub_pos_RL{s} = mean(pos_epochs(retIdx & leftIdx, :), 1, 'omitnan');
        end
        if nRR > 0
            sub_pos_RR{s} = mean(pos_epochs(retIdx & rightIdx, :), 1, 'omitnan');
        end

        % ================================================================
        %  STFT — VECTORIZED COMPUTATION
        % ================================================================
        halfWin = floor(stft_winSamp / 2);

        % Accumulators: power sum and phase-vector sum per direction
        powS_L = zeros(nChan,nFreq,nTime); phS_L = complex(zeros(nChan,nFreq,nTime));
        powS_R = zeros(nChan,nFreq,nTime); phS_R = complex(zeros(nChan,nFreq,nTime));
        powS_SL = zeros(nChan,nFreq,nTime); phS_SL = complex(zeros(nChan,nFreq,nTime));
        powS_SR = zeros(nChan,nFreq,nTime); phS_SR = complex(zeros(nChan,nFreq,nTime));
        powS_RL = zeros(nChan,nFreq,nTime); phS_RL = complex(zeros(nChan,nFreq,nTime));
        powS_RR = zeros(nChan,nFreq,nTime); phS_RR = complex(zeros(nChan,nFreq,nTime));

        for ti = 1:nTime
            [~, cSamp] = min(abs(epochTime - stft_toi(ti)));
            wLo = cSamp - halfWin;
            wHi = wLo + stft_winSamp - 1;
            if wLo < 1 || wHi > nEpochSamp, continue; end

            % Extract & window [nMS x nChan x winSamp]
            seg = eeg_epochs(:, :, wLo:wHi);
            seg = seg .* reshape(hann_win, [1, 1, stft_winSamp]);

            % FFT with zero-padding [nMS x nChan x nfft]
            F_all = fft(seg, nfft, 3);

            % Select frequencies of interest [nMS x nChan x nFreq]
            F_foi = F_all(:, :, freq_bins);
            P = abs(F_foi).^2;
            F_norm = F_foi ./ (abs(F_foi) + eps);

            % Helper to safely squeeze [1 x nChan x nFreq] -> [nChan x nFreq]
            sq = @(X) reshape(X, [nChan, nFreq]);

            % --- Left MS ---
            if nL > 0
                powS_L(:,:,ti) = sq(mean(P(leftIdx,:,:),1));
                phS_L(:,:,ti)  = sq(mean(F_norm(leftIdx,:,:),1));
            end
            % --- Right MS ---
            if nR > 0
                powS_R(:,:,ti) = sq(mean(P(rightIdx,:,:),1));
                phS_R(:,:,ti)  = sq(mean(F_norm(rightIdx,:,:),1));
            end
            % --- Start Left ---
            if nSL > 0
                powS_SL(:,:,ti) = sq(mean(P(startIdx&leftIdx,:,:),1));
                phS_SL(:,:,ti)  = sq(mean(F_norm(startIdx&leftIdx,:,:),1));
            end
            % --- Start Right ---
            if nSR > 0
                powS_SR(:,:,ti) = sq(mean(P(startIdx&rightIdx,:,:),1));
                phS_SR(:,:,ti)  = sq(mean(F_norm(startIdx&rightIdx,:,:),1));
            end
            % --- Return Left ---
            if nRL > 0
                powS_RL(:,:,ti) = sq(mean(P(retIdx&leftIdx,:,:),1));
                phS_RL(:,:,ti)  = sq(mean(F_norm(retIdx&leftIdx,:,:),1));
            end
            % --- Return Right ---
            if nRR > 0
                powS_RR(:,:,ti) = sq(mean(P(retIdx&rightIdx,:,:),1));
                phS_RR(:,:,ti)  = sq(mean(F_norm(retIdx&rightIdx,:,:),1));
            end
        end

        % Store per-subject results
        sub_pow_L{s}  = powS_L;  sub_pow_R{s}  = powS_R;
        sub_pow_SL{s} = powS_SL; sub_pow_SR{s} = powS_SR;
        sub_pow_RL{s} = powS_RL; sub_pow_RR{s} = powS_RR;
        sub_itpc_L{s}  = abs(phS_L);  sub_itpc_R{s}  = abs(phS_R);
        sub_itpc_SL{s} = abs(phS_SL); sub_itpc_SR{s} = abs(phS_SR);
        sub_itpc_RL{s} = abs(phS_RL); sub_itpc_RR{s} = abs(phS_RR);

        clear eeg_epochs dataTFR dataET F_all seg P F_foi F_norm;
    end % subject loop

    % ====================================================================
    %  3. GRAND AVERAGES
    % ====================================================================
    validIdx = find(sub_valid);
    nValid   = numel(validIdx);
    fprintf('\n  Valid subjects: %d / %d\n', nValid, nSubj);
    if nValid < 3
        fprintf('  Too few subjects for grand average, skipping task.\n');
        continue;
    end

    iPO7 = find(strcmpi(chanLabels, elec_L), 1);
    iPO8 = find(strcmpi(chanLabels, elec_R), 1);

    % --- Helper: stack and average ---
    stack3  = @(C) cat(4, C{:});          % [nChan x nFreq x nTime x nSubj]
    stack1  = @(C) cat(1, C{:});          % [nSubj x nSamp]

    % ---- All MS: lateralisation at PO7/PO8 ----
    lat_all_subj = nan(nFreq, nTime, nValid);
    contra_all_subj = nan(nFreq, nTime, nValid);
    ipsi_all_subj   = nan(nFreq, nTime, nValid);
    lat_topo_all_subj = nan(numel(chanLabels), nFreq, nTime, nValid);
    itpc_lat_all_subj = nan(nFreq, nTime, nValid);
    itpc_contra_all_subj = nan(nFreq, nTime, nValid);
    itpc_ipsi_all_subj   = nan(nFreq, nTime, nValid);
    itpc_topo_L_all_subj = nan(numel(chanLabels), nFreq, nTime, nValid);
    itpc_topo_R_all_subj = nan(numel(chanLabels), nFreq, nTime, nValid);

    % ---- Start MS ----
    lat_start_subj = nan(nFreq, nTime, nValid);
    contra_start_subj = nan(nFreq, nTime, nValid);
    ipsi_start_subj   = nan(nFreq, nTime, nValid);
    lat_topo_start_subj = nan(numel(chanLabels), nFreq, nTime, nValid);
    itpc_contra_start_subj = nan(nFreq, nTime, nValid);
    itpc_ipsi_start_subj   = nan(nFreq, nTime, nValid);
    itpc_topo_L_start_subj = nan(numel(chanLabels), nFreq, nTime, nValid);
    itpc_topo_R_start_subj = nan(numel(chanLabels), nFreq, nTime, nValid);

    % ---- Return MS ----
    lat_ret_subj = nan(nFreq, nTime, nValid);
    lat_topo_ret_subj = nan(numel(chanLabels), nFreq, nTime, nValid);

    % ---- Gaze ----
    vel_subj = nan(nValid, nVelSamp);
    pos_SL_subj = nan(nValid, nPosSamp); pos_SR_subj = nan(nValid, nPosSamp);
    pos_RL_subj = nan(nValid, nPosSamp); pos_RR_subj = nan(nValid, nPosSamp);

    for vi = 1:nValid
        s = validIdx(vi);
        pL = sub_pow_L{s}; pR = sub_pow_R{s};

        % Contra/ipsi at PO7/PO8 (all MS)
        contra_po = (reshape(pL(iPO8,:,:),[nFreq,nTime]) + reshape(pR(iPO7,:,:),[nFreq,nTime])) / 2;
        ipsi_po   = (reshape(pL(iPO7,:,:),[nFreq,nTime]) + reshape(pR(iPO8,:,:),[nFreq,nTime])) / 2;
        lat = ((contra_po - ipsi_po) ./ (contra_po + ipsi_po + eps)) * 100;
        lat_all_subj(:,:,vi)     = lat;
        contra_all_subj(:,:,vi)  = contra_po;
        ipsi_all_subj(:,:,vi)    = ipsi_po;

        % Topo lateralisation (left vs right MS)
        lat_topo_all_subj(:,:,:,vi) = ((pL - pR) ./ (pL + pR + eps)) * 100;

        % ITPC contra/ipsi (all MS)
        iL = sub_itpc_L{s}; iR = sub_itpc_R{s};
        itpc_contra = (reshape(iL(iPO8,:,:),[nFreq,nTime]) + reshape(iR(iPO7,:,:),[nFreq,nTime])) / 2;
        itpc_ipsi   = (reshape(iL(iPO7,:,:),[nFreq,nTime]) + reshape(iR(iPO8,:,:),[nFreq,nTime])) / 2;
        itpc_contra_all_subj(:,:,vi) = itpc_contra;
        itpc_ipsi_all_subj(:,:,vi)   = itpc_ipsi;
        itpc_lat_all_subj(:,:,vi)    = itpc_ipsi - itpc_contra;
        itpc_topo_L_all_subj(:,:,:,vi) = iL;
        itpc_topo_R_all_subj(:,:,:,vi) = iR;

        % ---- Start MS ----
        if sub_nMS(s,3) >= 3 && sub_nMS(s,4) >= 3
            pSL = sub_pow_SL{s}; pSR = sub_pow_SR{s};
            c_s = (reshape(pSL(iPO8,:,:),[nFreq,nTime]) + reshape(pSR(iPO7,:,:),[nFreq,nTime])) / 2;
            i_s = (reshape(pSL(iPO7,:,:),[nFreq,nTime]) + reshape(pSR(iPO8,:,:),[nFreq,nTime])) / 2;
            lat_start_subj(:,:,vi)     = ((c_s - i_s)./(c_s + i_s + eps))*100;
            contra_start_subj(:,:,vi)  = c_s;
            ipsi_start_subj(:,:,vi)    = i_s;
            lat_topo_start_subj(:,:,:,vi) = ((pSL - pSR)./(pSL + pSR + eps))*100;

            iSL = sub_itpc_SL{s}; iSR = sub_itpc_SR{s};
            itpc_contra_start_subj(:,:,vi) = (reshape(iSL(iPO8,:,:),[nFreq,nTime])+reshape(iSR(iPO7,:,:),[nFreq,nTime]))/2;
            itpc_ipsi_start_subj(:,:,vi)   = (reshape(iSL(iPO7,:,:),[nFreq,nTime])+reshape(iSR(iPO8,:,:),[nFreq,nTime]))/2;
            itpc_topo_L_start_subj(:,:,:,vi) = iSL;
            itpc_topo_R_start_subj(:,:,:,vi) = iSR;
        end

        % ---- Return MS ----
        if sub_nMS(s,5) >= 3 && sub_nMS(s,6) >= 3
            pRL = sub_pow_RL{s}; pRR = sub_pow_RR{s};
            c_r = (reshape(pRL(iPO8,:,:),[nFreq,nTime]) + reshape(pRR(iPO7,:,:),[nFreq,nTime])) / 2;
            i_r = (reshape(pRL(iPO7,:,:),[nFreq,nTime]) + reshape(pRR(iPO8,:,:),[nFreq,nTime])) / 2;
            lat_ret_subj(:,:,vi) = ((c_r - i_r)./(c_r + i_r + eps))*100;
            lat_topo_ret_subj(:,:,:,vi) = ((pRL - pRR)./(pRL + pRR + eps))*100;
        end

        % ---- Gaze ----
        if ~isempty(sub_vel{s}), vel_subj(vi,:) = sub_vel{s}; end
        if ~isempty(sub_pos_SL{s}), pos_SL_subj(vi,:) = sub_pos_SL{s}; end
        if ~isempty(sub_pos_SR{s}), pos_SR_subj(vi,:) = sub_pos_SR{s}; end
        if ~isempty(sub_pos_RL{s}), pos_RL_subj(vi,:) = sub_pos_RL{s}; end
        if ~isempty(sub_pos_RR{s}), pos_RR_subj(vi,:) = sub_pos_RR{s}; end
    end

    % --- Compute grand averages and SEM ---
    ga = @(X, dim) mean(X, dim, 'omitnan');
    gaSEM = @(X, dim) std(X, [], dim, 'omitnan') ./ sqrt(sum(~isnan(X), dim));

    % All MS
    ga_lat_all      = ga(lat_all_subj, 3);      sem_lat_all = gaSEM(lat_all_subj, 3);
    ga_contra_all   = ga(contra_all_subj, 3);
    ga_ipsi_all     = ga(ipsi_all_subj, 3);
    ga_lat_topo_all = ga(lat_topo_all_subj, 4);

    % Alpha TC (all MS)
    ga_lat_alpha_all    = ga(squeeze(mean(lat_all_subj(alphaIdx,:,:),1)), 2);
    sem_lat_alpha_all   = gaSEM(squeeze(mean(lat_all_subj(alphaIdx,:,:),1)), 2);

    % Start MS
    ga_lat_start      = ga(lat_start_subj, 3);
    ga_contra_start   = ga(contra_start_subj, 3);
    ga_ipsi_start     = ga(ipsi_start_subj, 3);
    ga_lat_topo_start = ga(lat_topo_start_subj, 4);
    ga_lat_alpha_start = ga(squeeze(mean(lat_start_subj(alphaIdx,:,:),1)), 2);
    sem_lat_alpha_start = gaSEM(squeeze(mean(lat_start_subj(alphaIdx,:,:),1)), 2);

    % Return MS
    ga_lat_ret      = ga(lat_ret_subj, 3);
    ga_lat_topo_ret = ga(lat_topo_ret_subj, 4);
    ga_lat_alpha_ret = ga(squeeze(mean(lat_ret_subj(alphaIdx,:,:),1)), 2);
    sem_lat_alpha_ret = gaSEM(squeeze(mean(lat_ret_subj(alphaIdx,:,:),1)), 2);

    % --- Baselined alpha lateralisation (subtract pre-MS baseline mean) ---
    % All MS
    alpha_lat_all_perSubj = squeeze(mean(lat_all_subj(alphaIdx,:,:), 1)); % [nTime x nValid]
    bl_alpha_all = mean(alpha_lat_all_perSubj(blIdx, :), 1);             % [1 x nValid]
    alpha_lat_all_bl_perSubj = alpha_lat_all_perSubj - bl_alpha_all;
    ga_lat_alpha_all_bl  = mean(alpha_lat_all_bl_perSubj, 2, 'omitnan');
    sem_lat_alpha_all_bl = std(alpha_lat_all_bl_perSubj, [], 2, 'omitnan') ./ ...
        sqrt(sum(~isnan(alpha_lat_all_bl_perSubj), 2));

    % Start MS
    alpha_lat_start_perSubj = squeeze(mean(lat_start_subj(alphaIdx,:,:), 1));
    bl_alpha_start = mean(alpha_lat_start_perSubj(blIdx, :), 1);
    alpha_lat_start_bl_perSubj = alpha_lat_start_perSubj - bl_alpha_start;
    ga_lat_alpha_start_bl  = mean(alpha_lat_start_bl_perSubj, 2, 'omitnan');
    sem_lat_alpha_start_bl = std(alpha_lat_start_bl_perSubj, [], 2, 'omitnan') ./ ...
        sqrt(sum(~isnan(alpha_lat_start_bl_perSubj), 2));

    % Return MS
    alpha_lat_ret_perSubj = squeeze(mean(lat_ret_subj(alphaIdx,:,:), 1));
    bl_alpha_ret = mean(alpha_lat_ret_perSubj(blIdx, :), 1);
    alpha_lat_ret_bl_perSubj = alpha_lat_ret_perSubj - bl_alpha_ret;
    ga_lat_alpha_ret_bl  = mean(alpha_lat_ret_bl_perSubj, 2, 'omitnan');
    sem_lat_alpha_ret_bl = std(alpha_lat_ret_bl_perSubj, [], 2, 'omitnan') ./ ...
        sqrt(sum(~isnan(alpha_lat_ret_bl_perSubj), 2));

    % ITPC (all MS)
    ga_itpc_contra_all = ga(itpc_contra_all_subj, 3);
    ga_itpc_ipsi_all   = ga(itpc_ipsi_all_subj, 3);
    ga_itpc_lat_all    = ga(itpc_lat_all_subj, 3);
    ga_itpc_topo_L_all = ga(itpc_topo_L_all_subj, 4);
    ga_itpc_topo_R_all = ga(itpc_topo_R_all_subj, 4);

    % ITPC (start MS)
    ga_itpc_contra_start = ga(itpc_contra_start_subj, 3);
    ga_itpc_ipsi_start   = ga(itpc_ipsi_start_subj, 3);
    ga_itpc_topo_L_start = ga(itpc_topo_L_start_subj, 4);
    ga_itpc_topo_R_start = ga(itpc_topo_R_start_subj, 4);

    % Figure 3: log-power baseline correction (start MS)
    ga_contra_start_log = log(ga_contra_start + eps);
    ga_ipsi_start_log   = log(ga_ipsi_start + eps);
    bl_contra = mean(ga_contra_start_log(:, blIdx), 2);
    bl_ipsi   = mean(ga_ipsi_start_log(:, blIdx), 2);
    ga_contra_start_bl = ga_contra_start_log - bl_contra;
    ga_ipsi_start_bl   = ga_ipsi_start_log   - bl_ipsi;

    % Alpha TC for Fig 3c
    ga_contra_alpha_bl = mean(ga_contra_start_bl(alphaIdx,:), 1);
    ga_ipsi_alpha_bl   = mean(ga_ipsi_start_bl(alphaIdx,:), 1);
    % Per-subject for SEM
    contra_alpha_bl_subj = nan(nValid, nTime);
    ipsi_alpha_bl_subj   = nan(nValid, nTime);
    for vi = 1:nValid
        cs = contra_start_subj(:,:,vi); is_v = ipsi_start_subj(:,:,vi);
        if all(isnan(cs(:))), continue; end
        cl = log(cs + eps); il = log(is_v + eps);
        blc = mean(cl(:,blIdx),2); bli = mean(il(:,blIdx),2);
        contra_alpha_bl_subj(vi,:) = mean(cl(alphaIdx,:) - blc(alphaIdx), 1);
        ipsi_alpha_bl_subj(vi,:)   = mean(il(alphaIdx,:) - bli(alphaIdx), 1);
    end
    sem_contra_alpha_bl = gaSEM(contra_alpha_bl_subj, 1);
    sem_ipsi_alpha_bl   = gaSEM(ipsi_alpha_bl_subj, 1);

    % Gaze velocity
    ga_vel     = ga(vel_subj, 1);
    sem_vel    = gaSEM(vel_subj, 1);

    % Gaze position
    ga_pos_SL = ga(pos_SL_subj, 1); ga_pos_SR = ga(pos_SR_subj, 1);
    ga_pos_RL = ga(pos_RL_subj, 1); ga_pos_RR = ga(pos_RR_subj, 1);

    % ====================================================================
    %  4. FIGURE 1 — All microsaccades
    % ====================================================================
    fprintf('  Plotting Figure 1 ...\n');
    fig1 = figure('Position', [50 50 1400 1000], 'Color', 'w');
    tsk = taskName{taskIdx};

    % --- Panel a: Gaze velocity ---
    ax1a = axes('Position', [0.07 0.76 0.38 0.20]);
    ci95 = 1.96 * sem_vel;
    fill([velTime, fliplr(velTime)], [ga_vel+ci95, fliplr(ga_vel-ci95)], ...
        [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(velTime, ga_vel, 'k-', 'LineWidth', lnWd);
    xline(0, 'k--', 'LineWidth', 1);
    xlabel('Time (ms)'); ylabel('Velocity (px/s)');
    title('a) Gaze velocity'); set(gca, 'FontSize', fntSz);

    % --- Panel b: TFR lateralisation ---
    ax1b = axes('Position', [0.07 0.42 0.55 0.28]);
    imagesc(stft_toi_ms, stft_foi, ga_lat_all);
    set(gca, 'YDir', 'normal'); colormap(ax1b, cmap_div);
    caxis([-max(abs(ga_lat_all(:))) max(abs(ga_lat_all(:)))]);
    cb = colorbar; cb.Label.String = 'Lateralisation (%)';
    xline(0, 'k--', 'LineWidth', 1.5); hold on;
    % Alpha band markers
    yline(alphaBand(1), 'k:', 'LineWidth', 0.8);
    yline(alphaBand(2), 'k:', 'LineWidth', 0.8);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    title('b) Spectral lateralisation (contra - ipsi)');
    set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel c: Alpha time course (raw) ---
    ax1c = axes('Position', [0.07 0.08 0.27 0.26]);
    ci95_a = 1.96 * sem_lat_alpha_all;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_all'+ci95_a', fliplr(ga_lat_alpha_all'-ci95_a')], ...
        [0.5 0.5 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_all, 'b-', 'LineWidth', lnWd);
    xline(0, 'k--', 'LineWidth', 1.5);
    yline(0, 'k:', 'LineWidth', 1);
    xlabel('Time (ms)'); ylabel('8-12 Hz Lat. (%)');
    title('c) Alpha lateralisation (raw)');
    set(gca, 'FontSize', fntSz);

    % --- Panel e: Alpha time course (baselined) ---
    ax1e = axes('Position', [0.37 0.08 0.27 0.26]);
    ci95_bl = 1.96 * sem_lat_alpha_all_bl;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_all_bl'+ci95_bl', fliplr(ga_lat_alpha_all_bl'-ci95_bl')], ...
        [0.5 0.8 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_all_bl, '-', 'Color', [0.1 0.6 0.2], 'LineWidth', lnWd);
    xline(0, 'k--', 'LineWidth', 1.5);
    yline(0, 'k:', 'LineWidth', 1);
    xlabel('Time (ms)'); ylabel('\Delta 8-12 Hz Lat. (%)');
    title('e) Alpha lateralisation (baselined)');
    set(gca, 'FontSize', fntSz);

    % --- Panel d: Topographical maps ---
    topo_times_ms = [-200, -100, 0, 100, 200]; % ms
    nTopos = numel(topo_times_ms);
    freq_topo = make_ft_freq_local(chanLabels, stft_foi, stft_toi, ga_lat_topo_all);

    for ti = 1:nTopos
        ax_d = axes('Position', [0.66 + (ti-1)*0.065, 0.42, 0.06, 0.28]);
        cfg = [];
        cfg.figure    = ax_d;
        cfg.layout    = layANThead;
        cfg.xlim      = topo_times_ms(ti)/1000 * [1 1] + [-0.01 0.01];
        cfg.ylim      = alphaBand;
        cfg.zlim      = 'maxabs';
        cfg.colorbar  = 'no';
        cfg.comment   = 'no';
        cfg.style     = 'straight';
        cfg.shading   = 'interp';
        cfg.colormap  = cmap_div;
        cfg.marker    = 'off';
        cfg.gridscale = 100;
        ft_topoplotTFR(cfg, freq_topo);
        title(sprintf('%d ms', topo_times_ms(ti)), 'FontSize', fntSz-2);
    end

    sgtitle(sprintf('%s — Figure 1: Microsaccades lateralise alpha (N=%d)', tsk, nValid), ...
        'FontSize', fntSz+2, 'FontWeight', 'bold');
    saveas(fig1, fullfile(figDir, sprintf('Fig1_allMS_%s.png', lower(tsk))));
    close(fig1);

    % ====================================================================
    %  5. FIGURE 2 — Start vs Return microsaccades
    % ====================================================================
    fprintf('  Plotting Figure 2 ...\n');
    fig2 = figure('Position', [50 50 1600 1400], 'Color', 'w');

    colL = [0.2 0.4 0.8]; colR = [0.8 0.2 0.2]; % left=blue, right=red

    % --- Panel a: Gaze position ---
    % Start
    ax2a1 = subplot(4, 2, 1);
    plot(posTime, ga_pos_SL, '-', 'Color', colL, 'LineWidth', lnWd); hold on;
    plot(posTime, ga_pos_SR, '-', 'Color', colR, 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('Gaze pos. (px)');
    title('a) Start microsaccades — Gaze position');
    legend('Left MS', 'Right MS', 'Location', 'best');
    set(gca, 'FontSize', fntSz);

    % Return
    ax2a2 = subplot(4, 2, 2);
    plot(posTime, ga_pos_RL, '-', 'Color', colL, 'LineWidth', lnWd); hold on;
    plot(posTime, ga_pos_RR, '-', 'Color', colR, 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('Gaze pos. (px)');
    title('a) Return microsaccades — Gaze position');
    legend('Left MS', 'Right MS', 'Location', 'best');
    set(gca, 'FontSize', fntSz);

    % --- Panel b: TFR lateralisation ---
    ax2b1 = subplot(4, 2, 3);
    imagesc(stft_toi_ms, stft_foi, ga_lat_start);
    set(gca, 'YDir', 'normal'); colormap(ax2b1, cmap_div);
    caxis([-max(abs(ga_lat_start(:))) max(abs(ga_lat_start(:)))]);
    colorbar; xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('b) Start MS — Lateralisation'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    ax2b2 = subplot(4, 2, 4);
    imagesc(stft_toi_ms, stft_foi, ga_lat_ret);
    set(gca, 'YDir', 'normal'); colormap(ax2b2, cmap_div);
    caxis([-max(abs(ga_lat_ret(:))) max(abs(ga_lat_ret(:)))]);
    colorbar; xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('b) Return MS — Lateralisation'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel c: Alpha lateralisation time course (raw) ---
    ax2c1 = subplot(4, 2, 5);
    ci_s = 1.96 * sem_lat_alpha_start;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_start'+ci_s', fliplr(ga_lat_alpha_start'-ci_s')], ...
        [0.5 0.5 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_start, 'b-', 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('8-12 Hz Lat. (%)');
    title('c) Start MS — Alpha lateralisation (raw)');
    set(gca, 'FontSize', fntSz);

    ax2c2 = subplot(4, 2, 6);
    ci_r = 1.96 * sem_lat_alpha_ret;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_ret'+ci_r', fliplr(ga_lat_alpha_ret'-ci_r')], ...
        [0.5 0.5 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_ret, 'b-', 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('8-12 Hz Lat. (%)');
    title('c) Return MS — Alpha lateralisation (raw)');
    set(gca, 'FontSize', fntSz);

    % --- Panel d: Baselined alpha lateralisation ---
    ax2d1 = subplot(4, 2, 7);
    ci_s_bl = 1.96 * sem_lat_alpha_start_bl;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_start_bl'+ci_s_bl', fliplr(ga_lat_alpha_start_bl'-ci_s_bl')], ...
        [0.5 0.8 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_start_bl, '-', 'Color', [0.1 0.6 0.2], 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('\Delta 8-12 Hz Lat. (%)');
    title('d) Start MS — Alpha lateralisation (baselined)');
    set(gca, 'FontSize', fntSz);

    ax2d2 = subplot(4, 2, 8);
    ci_r_bl = 1.96 * sem_lat_alpha_ret_bl;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_lat_alpha_ret_bl'+ci_r_bl', fliplr(ga_lat_alpha_ret_bl'-ci_r_bl')], ...
        [0.5 0.8 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(stft_toi_ms, ga_lat_alpha_ret_bl, '-', 'Color', [0.1 0.6 0.2], 'LineWidth', lnWd);
    xline(0, 'k--'); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('\Delta 8-12 Hz Lat. (%)');
    title('d) Return MS — Alpha lateralisation (baselined)');
    set(gca, 'FontSize', fntSz);

    sgtitle(sprintf('%s — Figure 2: Start vs Return microsaccades (N=%d)', tsk, nValid), ...
        'FontSize', fntSz+2, 'FontWeight', 'bold');
    saveas(fig2, fullfile(figDir, sprintf('Fig2_startReturn_%s.png', lower(tsk))));
    close(fig2);

    % --- Figure 2d: Topographical maps (start / return) ---
    fig2d = figure('Position', [50 50 1200 350], 'Color', 'w');
    topo_t_ms = [0, 100, 200]; % ms
    freq_topo_start = make_ft_freq_local(chanLabels, stft_foi, stft_toi, ga_lat_topo_start);
    freq_topo_ret   = make_ft_freq_local(chanLabels, stft_foi, stft_toi, ga_lat_topo_ret);
    for ti = 1:numel(topo_t_ms)
        % Start
        ax_s = subplot(2, numel(topo_t_ms), ti);
        cfg = [];
        cfg.figure = ax_s;
        cfg.layout = layANThead; cfg.xlim = topo_t_ms(ti)/1000*[1 1]+[-0.01 0.01];
        cfg.ylim = alphaBand; cfg.zlim = 'maxabs'; cfg.colorbar = 'no';
        cfg.comment = 'no'; cfg.style = 'straight'; cfg.shading = 'interp';
        cfg.colormap = cmap_div; cfg.marker = 'off'; cfg.gridscale = 100;
        ft_topoplotTFR(cfg, freq_topo_start);
        title(sprintf('Start %dms', topo_t_ms(ti)), 'FontSize', fntSz-3);
        % Return
        ax_r = subplot(2, numel(topo_t_ms), ti + numel(topo_t_ms));
        cfg.figure = ax_r;
        ft_topoplotTFR(cfg, freq_topo_ret);
        title(sprintf('Return %dms', topo_t_ms(ti)), 'FontSize', fntSz-3);
    end
    sgtitle(sprintf('%s — Fig 2d: Topographies (8-12 Hz, start vs return)', tsk), ...
        'FontSize', fntSz, 'FontWeight', 'bold');
    saveas(fig2d, fullfile(figDir, sprintf('Fig2d_topos_%s.png', lower(tsk))));
    close(fig2d);

    % ====================================================================
    %  6. FIGURE 3 — Ipsilateral vs Contralateral power (start MS)
    % ====================================================================
    fprintf('  Plotting Figure 3 ...\n');
    fig3 = figure('Position', [50 50 1600 1000], 'Color', 'w');

    % --- Panel a: Raw power (contra / ipsi) ---
    ax3a1 = subplot(2, 3, 1);
    imagesc(stft_toi_ms, stft_foi, ga_contra_start);
    set(gca, 'YDir', 'normal'); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('a) Contralateral — Raw power'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    ax3a2 = subplot(2, 3, 4);
    imagesc(stft_toi_ms, stft_foi, ga_ipsi_start);
    set(gca, 'YDir', 'normal'); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('a) Ipsilateral — Raw power'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel b: Baseline-corrected log power ---
    ax3b1 = subplot(2, 3, 2);
    imagesc(stft_toi_ms, stft_foi, ga_contra_start_bl);
    set(gca, 'YDir', 'normal'); colormap(ax3b1, cmap_div);
    cl = max(abs(ga_contra_start_bl(:)));
    caxis([-cl cl]); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    yline(alphaBand(1), 'k:'); yline(alphaBand(2), 'k:');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('b) Contralateral — log power (BL corr)'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    ax3b2 = subplot(2, 3, 5);
    imagesc(stft_toi_ms, stft_foi, ga_ipsi_start_bl);
    set(gca, 'YDir', 'normal'); colormap(ax3b2, cmap_div);
    cl2 = max(abs(ga_ipsi_start_bl(:)));
    caxis([-cl2 cl2]); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    yline(alphaBand(1), 'k:'); yline(alphaBand(2), 'k:');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('b) Ipsilateral — log power (BL corr)'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel c: Alpha time courses ---
    ax3c = subplot(2, 3, [3 6]);
    ci_c = 1.96 * sem_contra_alpha_bl;
    ci_i = 1.96 * sem_ipsi_alpha_bl;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_contra_alpha_bl+ci_c, fliplr(ga_contra_alpha_bl-ci_c)], ...
        colL, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    fill([stft_toi_ms, fliplr(stft_toi_ms)], ...
        [ga_ipsi_alpha_bl+ci_i, fliplr(ga_ipsi_alpha_bl-ci_i)], ...
        colR, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h1 = plot(stft_toi_ms, ga_contra_alpha_bl, '-', 'Color', colL, 'LineWidth', lnWd);
    h2 = plot(stft_toi_ms, ga_ipsi_alpha_bl, '-', 'Color', colR, 'LineWidth', lnWd);
    xline(0, 'k--', 'LineWidth', 1.5); yline(0, 'k:');
    xlabel('Time (ms)'); ylabel('8-12 Hz log power change');
    title('c) Alpha power time course (start MS)');
    legend([h1 h2], {'Contralateral', 'Ipsilateral'}, 'Location', 'best');
    set(gca, 'FontSize', fntSz);

    sgtitle(sprintf('%s — Figure 3: Ipsi/Contra power — start MS (N=%d)', tsk, nValid), ...
        'FontSize', fntSz+2, 'FontWeight', 'bold');
    saveas(fig3, fullfile(figDir, sprintf('Fig3_ipsiContra_%s.png', lower(tsk))));
    close(fig3);

    % ====================================================================
    %  7. FIGURE 4 — ITPC (start MS)
    % ====================================================================
    fprintf('  Plotting Figure 4 ...\n');
    fig4 = figure('Position', [50 50 1600 1000], 'Color', 'w');

    % --- Panel a: ITPC contralateral ---
    ax4a1 = subplot(2, 3, 1);
    imagesc(stft_toi_ms, stft_foi, ga_itpc_contra_start);
    set(gca, 'YDir', 'normal'); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('a) Contralateral ITPC'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel a: ITPC ipsilateral ---
    ax4a2 = subplot(2, 3, 4);
    imagesc(stft_toi_ms, stft_foi, ga_itpc_ipsi_start);
    set(gca, 'YDir', 'normal'); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('a) Ipsilateral ITPC'); set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel b: ITPC difference (ipsi - contra) ---
    ax4b = subplot(2, 3, [2 5]);
    itpc_diff_start = ga_itpc_ipsi_start - ga_itpc_contra_start;
    imagesc(stft_toi_ms, stft_foi, itpc_diff_start);
    set(gca, 'YDir', 'normal'); colormap(ax4b, cmap_div);
    cl_itpc = max(abs(itpc_diff_start(:)));
    caxis([-cl_itpc cl_itpc]); colorbar;
    xline(0, 'k--', 'LineWidth', 1.5);
    yline(alphaBand(1), 'k:'); yline(alphaBand(2), 'k:');
    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title('b) ITPC difference (ipsi - contra)');
    set(gca, 'FontSize', fntSz); ylim([1 30]);

    % --- Panel c/d: ITPC topographies ---
    % 8-12 Hz, 0-250 ms
    itpc_topoWin = [0 0.25]; % s
    itpc_freqWin = alphaBand;

    % Left MS ITPC topo
    ax4c1 = subplot(2, 3, 3);
    freq_itpc_L = make_ft_freq_local(chanLabels, stft_foi, stft_toi, ga_itpc_topo_L_start);
    cfg = [];
    cfg.figure = ax4c1;
    cfg.layout = layANThead; cfg.xlim = itpc_topoWin; cfg.ylim = itpc_freqWin;
    cfg.zlim = 'maxabs'; cfg.colorbar = 'yes'; cfg.comment = 'no';
    cfg.style = 'straight'; cfg.shading = 'interp'; cfg.marker = 'off';
    cfg.gridscale = 100;
    ft_topoplotTFR(cfg, freq_itpc_L);
    title('c) Left MS — ITPC (8-12 Hz, 0-250ms)', 'FontSize', fntSz-2);

    % Right MS ITPC topo
    ax4c2 = subplot(2, 3, 6);
    freq_itpc_R = make_ft_freq_local(chanLabels, stft_foi, stft_toi, ga_itpc_topo_R_start);
    cfg.figure = ax4c2;
    ft_topoplotTFR(cfg, freq_itpc_R);
    title('d) Right MS — ITPC (8-12 Hz, 0-250ms)', 'FontSize', fntSz-2);

    sgtitle(sprintf('%s — Figure 4: ITPC — start MS (N=%d)', tsk, nValid), ...
        'FontSize', fntSz+2, 'FontWeight', 'bold');
    saveas(fig4, fullfile(figDir, sprintf('Fig4_ITPC_%s.png', lower(tsk))));
    close(fig4);

    % --- Figure 4d: ITPC difference topography ---
    fig4d = figure('Position', [50 50 500 400], 'Color', 'w');
    ax4d = axes;
    itpc_topo_diff = ga_itpc_topo_L_start - ga_itpc_topo_R_start;
    freq_itpc_diff = make_ft_freq_local(chanLabels, stft_foi, stft_toi, itpc_topo_diff);
    cfg = [];
    cfg.figure = ax4d;
    cfg.layout = layANThead; cfg.xlim = itpc_topoWin; cfg.ylim = itpc_freqWin;
    cfg.zlim = 'maxabs'; cfg.colorbar = 'yes'; cfg.comment = 'no';
    cfg.style = 'straight'; cfg.shading = 'interp'; cfg.marker = 'off';
    cfg.gridscale = 100; cfg.colormap = cmap_div;
    ft_topoplotTFR(cfg, freq_itpc_diff);
    title(sprintf('%s — ITPC difference (left-right MS, 8-12 Hz, 0-250ms)', tsk), ...
        'FontSize', fntSz);
    saveas(fig4d, fullfile(figDir, sprintf('Fig4d_ITPC_topo_diff_%s.png', lower(tsk))));
    close(fig4d);

    fprintf('  Figures saved for %s.\n', tsk);

end % task loop

fprintf('\n=== All done. Figures saved to: ===\n  %s\n', figDir);

%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [vx, vy] = compute_velocity_sg_local(X, Y, fs, polyOrd)
%COMPUTE_VELOCITY_SG_LOCAL  Savitzky-Golay velocity estimation.
Ts = 1 / fs;
L  = numel(X);
framelen = min(21, L);
if mod(framelen, 2) == 0, framelen = framelen - 1; end
minLegal = polyOrd + 3;
if mod(minLegal, 2) == 0, minLegal = minLegal + 1; end
if framelen < minLegal, framelen = minLegal; end
if framelen > L
    framelen = L;
    if mod(framelen, 2) == 0, framelen = framelen - 1; end
end
if framelen < 5
    vx = gradient(X) * fs;
    vy = gradient(Y) * fs;
    return;
end
Xs = sgolayfilt(double(X), polyOrd, framelen);
Ys = sgolayfilt(double(Y), polyOrd, framelen);
[~, G] = sgolay(polyOrd, framelen);
d1 = (factorial(1) / (Ts^1)) * G(:, 2)';
vx = conv(Xs, d1, 'same');
vy = conv(Ys, d1, 'same');
end

function [on, off] = merge_close_local(on, off, minISI)
%MERGE_CLOSE_LOCAL  Merge saccade events closer than minISI samples.
if numel(on) <= 1, return; end
O = on(1); F = off(1);
onM = []; offM = [];
for i = 2:numel(on)
    if on(i) - F <= minISI
        F = off(i);
    else
        onM = [onM, O]; offM = [offM, F]; %#ok<AGROW>
        O = on(i); F = off(i);
    end
end
onM = [onM, O]; offM = [offM, F];
on = onM; off = offM;
end

function freq = make_ft_freq_local(labels, foi, toi, powspctrm)
%MAKE_FT_FREQ_LOCAL  Create a FieldTrip-compatible freq structure.
%   powspctrm: [nChan x nFreq x nTime]
freq = [];
freq.label      = labels(:);
freq.dimord     = 'chan_freq_time';
freq.freq       = foi;
freq.time       = toi;
freq.powspctrm  = powspctrm;
end
