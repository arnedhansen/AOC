%% AOC EEG Feature Extraction N-back
%
% Extracted features:
%   Power Spectrum
%   IAF, Power at IAF, and Lateralization Index
%   POWER TRIAL-BY-TRIAL
%   TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM

%% POWER Spectrum
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(dataEEG.trialinfo(:, 1) == 21);
        ind2 = find(dataEEG.trialinfo(:, 1) == 22);
        ind3 = find(dataEEG.trialinfo(:, 1) == 23);

        % Select data
        cfg = [];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg, dataEEG);

        % Frequency analysis settings
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'hanning';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'no';% do not keep single trials in output
        cfg.pad = 5;

        % Frequency analysis settings
        cfg.trials = ind1;
        powload1 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind3;
        powload3 = ft_freqanalysis(cfg,dat);

        % Save raw power spectra
        cd(datapath)
        save power_nback powload1 powload2 powload3

    catch ME
        ME.message
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define channels
datapath = strcat(path, subjects{1}, filesep, 'eeg');
cd(datapath);
load('power_nback.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2.label)
    label = powload2.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

% Left and right channels
left_channels = {};
right_channels = {};
for i = 1:length(channels)
    try
        ch = channels{i};
        % Find the first numeric part in the channel name
        numStr = regexp(ch, '\d+', 'match');
        % Convert the first numerical token to a number
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch;
        else
            right_channels{end+1} = ch;
        end
    catch ME
        ME.message
        disp(['Midline channel: ', ch])
    end
end

% Load data and calculate power, IAF and lateralization index
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_nback.mat');

        % Channels selection based on CBPT
        channelIdx = find(ismember(powload1.label, channels));

        % Extract power spectra for selected channels
        powspctrm1 = mean(powload1.powspctrm(channelIdx, :), 1);
        powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
        powspctrm3 = mean(powload3.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));

        % Calculate IAF for 1-back
        alphaPower1 = powspctrm1(alphaIndices);
        [pks1,locs] = findpeaks(alphaPower1);
        [~, ind] = max(pks1);
        IAF1 = powload1.freq(alphaIndices(locs(ind)));
        IAF_range1 = find(powload1.freq > (IAF1-4) & powload1.freq < (IAF1+2));

        % Calculate IAF for 2-back
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        [~, ind] = max(pks2);
        IAF2 = powload2.freq(alphaIndices(locs(ind)));
        IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));

        % Calculate IAF for 3-back
        alphaPower3 = powspctrm3(alphaIndices);
        [pks3,locs] = findpeaks(alphaPower3);
        [~, ind] = max(pks3);
        IAF3 = powload3.freq(alphaIndices(locs(ind)));
        IAF_range3 = find(powload3.freq > (IAF3-4) & powload3.freq < (IAF3+2));

        % Store the power values at the calculated IAFs
        powerIAF1 = mean(powspctrm1(IAF_range1));
        powerIAF2 = mean(powspctrm2(IAF_range2));
        powerIAF3 = mean(powspctrm3(IAF_range3));

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF1 == alphaRange(1) || IAF1 == alphaRange(2)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF3 == alphaRange(1) || IAF3 == alphaRange(2)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF1 > max(pks1)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF3 > max(pks3)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Compute lateralization index
        % as done in Stroganova et al., 2007
        powloads = {powload1, powload2, powload3};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx1 = LatIdx(1);
        LatIdx2 = LatIdx(2);
        LatIdx3 = LatIdx(3);

        % Create a structure array for this subject
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'AlphaPower', num2cell([powerIAF1; powerIAF2; powerIAF3]), 'IAF', num2cell([IAF1; IAF2; IAF3]), 'Lateralization', num2cell([LatIdx1; LatIdx2; LatIdx3]));

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        save lateralization_nback LatIdx1 LatIdx2 LatIdx3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), ' ...
            '3-back: %f Hz (Power: %f) |Lateralization: %f %f %f \n'], subjects{subj}, IAF1, ...
            powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3, LatIdx1, LatIdx2, LatIdx3);
    catch ME
        ME.message
        error(['ERROR calculating alpha power and IAF for Subject ' num2str(subjects{subj}) '!'])
    end
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\eeg_matrix_nback eeg_data_nback
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_nback eeg_data_nback
end

%% POWER TRIAL-BY-TRIAL
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

for subj = 1:length(subjects)
    clc
    disp(['Processing Subject ', subjects{subj}])
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        %% Frequency analysis
        % Identify indices of trials belonging to conditions
        ind1 = find(dataEEG.trialinfo == 21);
        ind1 = find(dataEEG.trialinfo == 22);
        ind3 = find(dataEEG.trialinfo == 23);

        % Select data
        cfg = [];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg,data);

        % Frequency analysis settings
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'yes';
        cfg.pad = 5;

        % Frequency analysis settings
        cfg.trials = ind1;
        powload1_trials = ft_freqanalysis(cfg,dat);
        powload1_trials.trialinfo = ones(1,length(powload1_trials.powspctrm(:, 1, 1)));
        cfg.trials = ind1;
        powload2_trials = ft_freqanalysis(cfg,dat);
        powload2_trials.trialinfo = ones(1,length(powload2_trials.powspctrm(:, 1, 1)))*2;
        cfg.trials = ind3;
        powload3_trials = ft_freqanalysis(cfg,dat);
        powload3_trials.trialinfo = ones(1,length(powload1_trials.powspctrm(:, 1, 1)))*3;

        %% Save
        cd(datapath)
        save power_nback_trials powload1_trials powload2_trials powload3_trials
    catch ME
        ME.message
        error(['ERROR extracting trial-by-trial power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM (NBACK) | AOC
% FOOOF output here is: (model fit - aperiodic) in FOOOF/log space (NOT peaks-only)

clear; close all; clc

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

%% Loop subjects
for subj = 1:length(subjects)

    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    close all
    clc

    disp(['Processing TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM (NBACK) for Subject AOC ', num2str(subjects{subj})])
    disp('Loading EEG TFR data (NBACK)')
    load dataEEG_TFR_nback

    %% Identify indices of trials belonging to conditions (trialinfo codes: 21, 22, 23)
    ind1 = find(dataTFR.trialinfo(:, 1) == 21);
    ind2 = find(dataTFR.trialinfo(:, 1) == 22);
    ind3 = find(dataTFR.trialinfo(:, 1) == 23);

    %% Trial-averaged TFR (for visualisation / later use)
    cfg            = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = 3:1:30;                            % N-back frequency range
    cfg.t_ftimwin  = ones(length(cfg.foi), 1) .* 0.5;   % 0.5 s windows
    cfg.toi        = -1:0.05:2.25;                      % N-back time axis
    cfg.keeptrials = 'no';
    cfg.pad        = 'nextpow2';

    cfg.trials = ind1; tfr1 = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind3; tfr3 = ft_freqanalysis(cfg, dataTFR);

    %% Baselined TFR (raw)
    cfg              = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'db';
    tfr1_bl          = ft_freqbaseline(cfg, tfr1);
    tfr2_bl          = ft_freqbaseline(cfg, tfr2);
    tfr3_bl          = ft_freqbaseline(cfg, tfr3);

    disp(upper('Raw TFR + baseline done...'))

    %%%%%%%%%%%%%%%%%%%%
    %%%%%%  FOOOF %%%%%%
    %%%%%%%%%%%%%%%%%%%%

    %% Sliding-window FOOOF over trial-averaged spectra (mtmfft)
    winLen  = 0.5;     % 500 ms
    stepLen = 0.05;    % 50 ms

    % Use data time axis (robust), then ensure full windows fit
    tMin = dataTFR.time{1}(1);
    tMax = dataTFR.time{1}(end);

    toi_start   = tMin:stepLen:(tMax - winLen);
    toi_centres = toi_start + winLen/2;
    nTimePnts   = numel(toi_centres);

    tfr_fooof = cell(1, 3);

    % FOOOF config
    cfg_fooof            = [];
    cfg_fooof.method     = 'mtmfft';
    cfg_fooof.taper      = 'hanning';
    cfg_fooof.foilim     = [3 30];
    cfg_fooof.pad        = 5;
    cfg_fooof.output     = 'fooof';
    cfg_fooof.keeptrials = 'no';   % average across trials before FOOOF

    %% Conditions
    for tfr_conds = 1:3

        if tfr_conds == 1
            trlIdx   = ind1;
            condName = '1back';
        elseif tfr_conds == 2
            trlIdx   = ind2;
            condName = '2back';
        else
            trlIdx   = ind3;
            condName = '3back';
        end

        disp(' ')
        disp(['Running sliding-window FOOOF (NBACK) for condition ', condName])

        % One test window to get sizes / master freq grid
        cfg_sel0         = [];
        cfg_sel0.latency = [toi_start(1) (toi_start(1) + winLen)];
        cfg_sel0.trials  = trlIdx;
        datTFR_win0      = ft_selectdata(cfg_sel0, dataTFR);

        fooof_test  = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win0);

        nChan       = numel(fooof_test.label);
        freq_master = fooof_test.freq(:);
        nFreq       = numel(freq_master);

        fooof_powspctrm = nan(nChan, nFreq, nTimePnts);   % (model fit - aperiodic) in FOOOF/log space
        fooof_powspec   = nan(nChan, nFreq, nTimePnts);   % input power spectrum in FOOOF/log space
        fooof_aperiodic = nan(nChan, 4,    nTimePnts);    % [offset exponent error r^2]

        for timePnt = 1:nTimePnts

            % Select sliding window + trials
            cfg_sel         = [];
            cfg_sel.latency = [toi_start(timePnt) (toi_start(timePnt) + winLen)];
            cfg_sel.trials  = trlIdx;
            datTFR_win      = ft_selectdata(cfg_sel, dataTFR);

            % Tracker
            clc
            disp(['Running FOOOF for Subject ', num2str(subjects{subj})])
            disp(['Subject:    ', num2str(subj), ' / ', num2str(length(subjects))])
            disp(['Condition:  ', num2str(tfr_conds), ' / 3 (', condName, ')'])
            disp(['Time Point: ', num2str(timePnt), ' / ', num2str(nTimePnts)])

            % Run FOOOF on averaged spectrum
            out = evalc('fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win);');

            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};
            else
                repdata = fooof_out.fooofparams;
            end

            freq_now = fooof_out.freq(:);

            % Map onto master grid (robust to tiny numeric shifts)
            [tf, loc] = ismembertol(freq_now, freq_master, 1e-10);

            local_ps     = nan(nChan, nFreq);
            local_pkcomp = nan(nChan, nFreq);  % (model fit - aperiodic) on master grid

            local_err    = nan(nChan, 1);
            local_rsq    = nan(nChan, 1);
            local_offset = nan(nChan, 1);
            local_expo   = nan(nChan, 1);

            for ch = 1:nChan

                % input spectrum (FOOOF/log space)
                ps_tmp = repdata(ch).power_spectrum(:);
                local_ps(ch, loc(tf)) = ps_tmp(tf).';

                % aperiodic fit (FOOOF/log space)
                ap = repdata(ch).aperiodic_params(:);

                local_err(ch) = repdata(ch).error;
                local_rsq(ch) = repdata(ch).r_squared;

                if numel(ap) == 2
                    offset = ap(1);
                    expo   = ap(2);
                    ap_fit = offset - expo .* log10(freq_now);
                else
                    offset = ap(1);
                    knee   = ap(2);
                    expo   = ap(3);
                    ap_fit = offset - log10(knee + freq_now.^expo);
                end

                local_offset(ch) = offset;
                local_expo(ch)   = expo;

                % full model fit (aperiodic + peaks) in FOOOF/log space
                if isfield(repdata(ch), 'fooofed_spectrum') && ~isempty(repdata(ch).fooofed_spectrum)
                    model_fit = repdata(ch).fooofed_spectrum(:);
                else
                    % fallback: reconstruct via Gaussian peaks + aperiodic fit
                    pk = repdata(ch).peak_params;
                    gauss_sum = zeros(numel(freq_now), 1);
                    if ~isempty(pk)
                        for p = 1:size(pk, 1)
                            cf  = pk(p, 1);
                            amp = pk(p, 2);
                            bw  = pk(p, 3);
                            gauss_sum = gauss_sum + amp .* exp(-(freq_now - cf).^2 ./ (2*bw.^2));
                        end
                    end
                    model_fit = ap_fit(:) + gauss_sum(:);
                end

                % model fit - aperiodic
                pkcomp = model_fit(:) - ap_fit(:);
                local_pkcomp(ch, loc(tf)) = pkcomp(tf).';
            end

            % Save
            fooof_powspctrm(:, :, timePnt) = local_pkcomp;
            fooof_powspec(:,   :, timePnt) = local_ps;

            fooof_aperiodic(:, 1, timePnt) = local_offset;
            fooof_aperiodic(:, 2, timePnt) = local_expo;
            fooof_aperiodic(:, 3, timePnt) = local_err;
            fooof_aperiodic(:, 4, timePnt) = local_rsq;
        end

        % Construct condition-specific FieldTrip freq-like struct
        tfr_ff                = [];
        tfr_ff.label          = fooof_test.label;
        tfr_ff.freq           = fooof_test.freq;
        tfr_ff.time           = toi_centres;
        tfr_ff.powspctrm      = fooof_powspctrm;   % (model fit - aperiodic) in FOOOF/log space
        tfr_ff.power_spectrum = fooof_powspec;     % input spectrum in FOOOF/log space
        tfr_ff.fooofparams    = fooof_aperiodic;   % chan x 4 x time
        tfr_ff.dimord         = 'chan_freq_time';

        tfr_fooof{tfr_conds} = tfr_ff;
    end

    % Unpack per condition
    tfr1_fooof = tfr_fooof{1};
    tfr2_fooof = tfr_fooof{2};
    tfr3_fooof = tfr_fooof{3};

    disp(upper('FOOOF done on trial-averaged spectra (sliding windows): model - aperiodic ...'))

    %% Sanity check (3-panel): input, input-aperiodic, model-aperiodic (stored)
    time_point = 0.5;
    [~, tim]   = min(abs(tfr1_fooof.time - time_point));

    latWin = [toi_start(tim) (toi_start(tim) + winLen)];

    cond_titles = {'1-back','2-back','3-back'};
    trl_all     = {ind1, ind2, ind3};

    figure('Position', [0 0 1500 420], 'Color', 'w');

    for c = 1:3

        trlIdx = trl_all{c};

        cfg_sel         = [];
        cfg_sel.latency = latWin;
        cfg_sel.trials  = trlIdx;
        datTFR_win_sc   = ft_selectdata(cfg_sel, dataTFR);

        fooof_sc = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win_sc);

        if iscell(fooof_sc.fooofparams)
            repdata_sc = fooof_sc.fooofparams{1};
        else
            repdata_sc = fooof_sc.fooofparams;
        end

        freq_sc = fooof_sc.freq(:);
        nChan   = numel(repdata_sc);

        ps_all     = nan(numel(freq_sc), nChan);
        ap_all     = nan(numel(freq_sc), nChan);
        model_all  = nan(numel(freq_sc), nChan);
        pkcomp_all = nan(numel(freq_sc), nChan);

        for ch = 1:nChan

            ps_tmp = repdata_sc(ch).power_spectrum(:);
            ps_all(:, ch) = ps_tmp;

            ap = repdata_sc(ch).aperiodic_params(:);
            if numel(ap) == 2
                offset = ap(1);
                expo   = ap(2);
                ap_fit = offset - expo .* log10(freq_sc);
            else
                offset = ap(1);
                knee   = ap(2);
                expo   = ap(3);
                ap_fit = offset - log10(knee + freq_sc.^expo);
            end
            ap_all(:, ch) = ap_fit(:);

            if isfield(repdata_sc(ch), 'fooofed_spectrum') && ~isempty(repdata_sc(ch).fooofed_spectrum)
                model_fit = repdata_sc(ch).fooofed_spectrum(:);
            else
                pk = repdata_sc(ch).peak_params;
                gauss_sum = zeros(numel(freq_sc), 1);
                if ~isempty(pk)
                    for p = 1:size(pk, 1)
                        cf  = pk(p, 1);
                        amp = pk(p, 2);
                        bw  = pk(p, 3);
                        gauss_sum = gauss_sum + amp .* exp(-(freq_sc - cf).^2 ./ (2*bw.^2));
                    end
                end
                model_fit = ap_fit(:) + gauss_sum(:);
            end

            model_all(:, ch)  = model_fit(:);
            pkcomp_all(:, ch) = model_fit(:) - ap_fit(:);
        end

        ps_in     = mean(ps_all,     2, 'omitnan');
        ap_fit    = mean(ap_all,     2, 'omitnan');
        pk_comp   = mean(pkcomp_all, 2, 'omitnan');
        ps_corr   = ps_in - ap_fit;

        if c == 1
            subplot(1,3,1); hold on
            plot(freq_sc, ps_in, 'LineWidth', 2)
            title('Original Powspctrm')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

            subplot(1,3,2); hold on
            plot(freq_sc, ps_corr, 'LineWidth', 2)
            title('Original - Aperiodic')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

            subplot(1,3,3); hold on
            plot(freq_sc, pk_comp, 'LineWidth', 2)
            title('Model fit - Aperiodic (stored)')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])
        end
    end

    sgtitle(sprintf('FOOOF sanity check (AOC NBACK): Subject %s | Window [%.2f %.2f] s | t = %.2f s', ...
        subjects{subj}, latWin(1), latWin(2), tfr1_fooof.time(tim)), 'FontSize', 18)

    if ispc
        savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
    else
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
    end
    if ~exist(savePathControls, 'dir')
        mkdir(savePathControls)
    end

    saveName = sprintf('AOC_controls_FOOOF_modelMinusAperiodic_nback_subj%s.png', subjects{subj});
    saveas(gcf, fullfile(savePathControls, saveName));

    %% Baselined FOOOF TFR
    cfg              = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';  % FOOOF already in log space
    tfr1_fooof_bl    = ft_freqbaseline(cfg, tfr1_fooof);
    tfr2_fooof_bl    = ft_freqbaseline(cfg, tfr2_fooof);
    tfr3_fooof_bl    = ft_freqbaseline(cfg, tfr3_fooof);

    disp(upper('FOOOF baseline done...'))

    %% Save TFR data
    cd(datapath)
    save tfr_nback ...
        tfr1 tfr2 tfr3 ...
        tfr1_fooof tfr2_fooof tfr3_fooof ...
        tfr1_bl tfr2_bl tfr3_bl ...
        tfr1_fooof_bl tfr2_fooof_bl tfr3_fooof_bl

    %% Convert TFR data to POWSPCTRM (channels x frequency)
    analysisPeriodFull  = [0 2];
    analysisPeriodEarly = [0 1];
    analysisPeriodLate  = [1 2];
    freq_range          = [3 30];

    pow1_fooof          = select_data(analysisPeriodFull,  freq_range, tfr1_fooof);
    pow1_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr1_fooof_bl);
    pow1_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr1_fooof_bl);
    pow1_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr1_fooof_bl);

    pow2_fooof          = select_data(analysisPeriodFull,  freq_range, tfr2_fooof);
    pow2_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr2_fooof_bl);

    pow3_fooof          = select_data(analysisPeriodFull,  freq_range, tfr3_fooof);
    pow3_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr3_fooof_bl);
    pow3_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr3_fooof_bl);
    pow3_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr3_fooof_bl);

    pow1_fooof          = remove_time_dimension(pow1_fooof);
    pow1_fooof_bl       = remove_time_dimension(pow1_fooof_bl);
    pow1_fooof_bl_early = remove_time_dimension(pow1_fooof_bl_early);
    pow1_fooof_bl_late  = remove_time_dimension(pow1_fooof_bl_late);

    pow2_fooof          = remove_time_dimension(pow2_fooof);
    pow2_fooof_bl       = remove_time_dimension(pow2_fooof_bl);
    pow2_fooof_bl_early = remove_time_dimension(pow2_fooof_bl_early);
    pow2_fooof_bl_late  = remove_time_dimension(pow2_fooof_bl_late);

    pow3_fooof          = remove_time_dimension(pow3_fooof);
    pow3_fooof_bl       = remove_time_dimension(pow3_fooof_bl);
    pow3_fooof_bl_early = remove_time_dimension(pow3_fooof_bl_early);
    pow3_fooof_bl_late  = remove_time_dimension(pow3_fooof_bl_late);

    save power_nback_fooof ...
        pow1_fooof pow2_fooof pow3_fooof ...
        pow1_fooof_bl pow2_fooof_bl pow3_fooof_bl ...
        pow1_fooof_bl_early pow2_fooof_bl_early pow3_fooof_bl_early ...
        pow1_fooof_bl_late pow2_fooof_bl_late pow3_fooof_bl_late

    clc
    fprintf('Subject AOC %s (%.3d/%.3d) DONE (sliding-window FOOOF: model - aperiodic) \n', ...
        num2str(subjects{subj}), subj, length(subjects))
end
