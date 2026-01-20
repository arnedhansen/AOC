%% AOC EEG Feature Extraction Sternberg
%
% Extracted features:
%   Power Spectrum (Retention)
%   Power Spectrum (Baseline)
%   IAF, Power at IAF, and Lateralization Index
%   TFR (Raw, FOOOF and Baselined)

%% POWSPCTRM (Retention)
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');
for subj = 1:length(subjects)
    clc
    disp(['Processing Retention POWSPCTRM for Subject AOC ', num2str(subjects{subj})])

    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_sternberg

        % Identify indices of trials belonging to conditions
        ind2 = find(dataEEG.trialinfo(:, 1) == 22); % WM load 2
        ind4 = find(dataEEG.trialinfo(:, 1) == 24); % WM load 4
        ind6 = find(dataEEG.trialinfo(:, 1) == 26); % WM load 6

        % Select data
        cfg = [];
        cfg.latency = [1 2];               % Segmentation for retention interval 1000ms - 2000ms
        dat = ft_selectdata(cfg, dataEEG); % Select data

        % Frequency analysis settings
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'no';% do not keep single trials in output
        cfg.pad = 5;

        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg, dat);
        cfg.trials = ind4;
        powload4 = ft_freqanalysis(cfg, dat);
        cfg.trials = ind6;
        powload6 = ft_freqanalysis(cfg, dat);

        % Save raw power spectra
        cd(datapath)
        save('power_stern_raw.mat', 'powload2', 'powload4', 'powload6')

    catch ME
        disp(ME.message)
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end

%% POWSPCTRM (Baseline & Long Retention Interval)
% Setup
% startup
% [subjects, path, ~ , ~] = setup('AOC');
%
% for subj = 1:length(subjects)
%     clc
%     disp(['Processing Baseline POWSPCTRM for Subject AOC ', num2str(subjects{subj})])
%
%     try
%         % Load data
%         datapath = strcat(path, subjects{subj}, filesep, 'eeg');
%         cd(datapath)
%         close all
%         load dataEEG_TFR_sternberg
%
%         % Identify indices of trials belonging to conditions
%         ind2 = find(dataTFR.trialinfo(:, 1) == 22); % WM load 2
%         ind4 = find(dataTFR.trialinfo(:, 1) == 24); % WM load 4
%         ind6 = find(dataTFR.trialinfo(:, 1) == 26); % WM load 6
%
%         % Frequency analysis
%         % Select data
%         cfg = [];                      % Empty configuration
%         cfg.latency = [-.5 0];     % Segmentation for retention interval
%         datalong = ft_selectdata(cfg, dataTFR);
%
%         % Analysis settings
%         cfg = [];                      % Empty configuration
%         cfg.output = 'pow';            % Estimate power only
%         cfg.method = 'mtmfft';         % Multi-taper FFT method
%         cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
%         cfg.tapsmofrq = 1;             % Smoothening frequency around foi
%         cfg.foilim = [3 30];           % Frequencies of interest
%         cfg.keeptrials = 'no';         % Discard trial information
%         cfg.pad = 5;                   % Add zero-padding
%
%         % Conduct frequency analysis for each condition separately
%         cfg.trials = ind2;
%         powload2_baseline_period = ft_freqanalysis(cfg, datalong);
%         cfg.trials = ind4;
%         powload4_baseline_period = ft_freqanalysis(cfg, datalong);
%         cfg.trials = ind6;
%         powload6_baseline_period = ft_freqanalysis(cfg, datalong);
%
%         % Save baselined power spectra
%         cd(datapath)
%         save power_stern_baseline_period powload2_baseline_period powload4_baseline_period powload6_baseline_period
%
%         % Frequency analysis for retention interval = 200 ms - 2000ms after stimulus presentation
%         disp(['Processing Long Retention Interval (200ms - 2000ms) POWSPCTRM for Subject AOC ', num2str(subjects{subj})])
%
%         % Select data
%         cfg = [];                              % Empty configuration
%         cfg.latency = [0.2 2];                 % Segmentation for retention interval 200ms - 2000ms
%         datLong = ft_selectdata(cfg, dataTFR); % Select data
%
%         % Analysis settings
%         cfg = [];                      % Empty configuration
%         cfg.output = 'pow';            % Estimate power only
%         cfg.method = 'mtmfft';         % Multi-taper FFT method
%         cfg.taper = 'dpss';            % Multiple tapers (discrete prolate spheroidal sequences)
%         cfg.tapsmofrq = 1;             % Smoothening frequency around foi
%         cfg.foilim = [3 30];           % Frequencies of interest
%         cfg.keeptrials = 'no';         % Discard trial information
%         cfg.pad = 5;                   % Add zero-padding
%
%         % Conduct frequency analysis for each condition separately
%         cfg.trials = ind2;
%         powload2long = ft_freqanalysis(cfg, datLong);
%         cfg.trials = ind4;
%         powload4long = ft_freqanalysis(cfg, datLong);
%         cfg.trials = ind6;
%         powload6long = ft_freqanalysis(cfg, datLong);
%
%         % Save data
%         cd(datapath)
%         save power_stern_long powload2long powload4long powload6long
%
%     catch ME
%         ME.message
%         error(['ERROR extracting baslined power for Subject ' num2str(subjects{subj}) '!'])
%     end
% end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern_raw.mat');
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

% Load data and calculate alpha power, IAF and lateralization index
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {});

for subj = 1:length(subjects)
    clc
    disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])

    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_stern_raw.mat');

        % Channel selection
        channelIdx = find(ismember(powload2.label, channels));

        % Extract power spectra for selected channels
        powspctrm2 = mean(powload2.powspctrm(channelIdx, :), 1);
        powspctrm4 = mean(powload4.powspctrm(channelIdx, :), 1);
        powspctrm6 = mean(powload6.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(powload2.freq >= alphaRange(1) & powload2.freq <= alphaRange(2));

        % Calculate IAF for WM load 2
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = powload2.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(powload2.freq > (IAF2-4) & powload2.freq < (IAF2+2));
            powerIAF2 = mean(powspctrm2(IAF_range2));
        end

        % Calculate IAF for WM load 4
        alphaPower4 = powspctrm4(alphaIndices);
        [pks4,locs] = findpeaks(alphaPower4);
        if isempty(pks4)
            IAF4 = NaN;
            IAF_range4 = NaN;
            powerIAF4 = NaN;
        else
            [~, ind] = max(pks4);
            IAF4 = powload4.freq(alphaIndices(locs(ind)));
            IAF_range4 = find(powload4.freq > (IAF4-4) & powload4.freq < (IAF4+2));
            powerIAF4 = mean(powspctrm4(IAF_range4));
        end

        % Calculate IAF for WM load 6
        alphaPower6 = powspctrm6(alphaIndices);
        [pks6,locs] = findpeaks(alphaPower6);
        if isempty(pks6)
            IAF6 = NaN;
            IAF_range6 = NaN;
            powerIAF6 = NaN;
        else
            [~, ind] = max(pks6);
            IAF6 = powload6.freq(alphaIndices(locs(ind)));
            IAF_range6 = find(powload6.freq > (IAF6-4) & powload6.freq < (IAF6+2));
            powerIAF6 = mean(powspctrm6(IAF_range6));
        end

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF4 == alphaRange(1) || IAF4 == alphaRange(2)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if IAF6 == alphaRange(1) || IAF6 == alphaRange(2)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF4 > max(pks4)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if powerIAF6 > max(pks6)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Compute lateralization index as done in Stroganova et al., 2007
        powloads = {powload2, powload4, powload6};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx2 = LatIdx(1);
        LatIdx4 = LatIdx(2);
        LatIdx6 = LatIdx(3);

        % Create a structure array for this subject
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'AlphaPower', num2cell([powerIAF2; powerIAF4; powerIAF6]), 'IAF', num2cell([IAF2; IAF4; IAF6]), ...
            'Lateralization', num2cell([LatIdx2; LatIdx4; LatIdx6]));

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_sternberg_subj subj_data_eeg
        save alpha_power_sternberg powerIAF2 powerIAF4 powerIAF6
        save IAF_sternberg IAF2 IAF4 IAF6
        save lateralization_sternberg LatIdx2 LatIdx4 LatIdx6
        eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: WM2: %f Hz (Power: %f), WM4: %f Hz (Power: %f), ' ...
            'WM6: %f Hz (Power: %f) | Lateralization: %f %f %f \n'], subjects{subj}, ...
            IAF2, powerIAF2, IAF4, powerIAF4, IAF6, powerIAF6, LatIdx2, LatIdx4, LatIdx6);
    catch ME
        ME.message
        error(['ERROR calculating alpha power and IAF for Subject ' num2str(subjects{subj}) '!'])
    end
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg eeg_data_sternberg
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/eeg_matrix_sternberg eeg_data_sternberg
end

%% TFR (Raw, FOOOF and Baselined) | Sternberg (AOC)
clear; close all; clc

startup
[subjects, path, ~, ~] = setup('AOC');

for subj = 1:length(subjects)

    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    close all
    clc

    disp(['Processing TFR (Raw, FOOOF and Baselined) for Subject AOC ', num2str(subjects{subj})])
    disp('Loading EEG TFR data')
    load dataEEG_TFR_sternberg

    %% Identify indices of trials belonging to conditions
    ind2 = find(dataTFR.trialinfo(:, 1) == 22);
    ind4 = find(dataTFR.trialinfo(:, 1) == 24);
    ind6 = find(dataTFR.trialinfo(:, 1) == 26);

    %% Trial-averaged TFR (for visualisation / later use)
    cfg            = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = 2:1:40;
    cfg.t_ftimwin  = ones(length(cfg.foi), 1) .* 0.5;
    cfg.toi        = -1.5:0.05:3;
    cfg.keeptrials = 'no';
    cfg.pad        = 'nextpow2';

    cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind4; tfr4 = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind6; tfr6 = ft_freqanalysis(cfg, dataTFR);

    %% Baselined TFR (raw)
    cfg              = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'db';
    tfr2_bl          = ft_freqbaseline(cfg, tfr2);
    tfr4_bl          = ft_freqbaseline(cfg, tfr4);
    tfr6_bl          = ft_freqbaseline(cfg, tfr6);

    disp(upper('Raw TFR + baseline done...'))

    %% Sliding-window FOOOF over trial-averaged spectra (mtmfft)
    winLen  = 0.5;     % 500 ms
    stepLen = 0.05;    % 50 ms

    % Use data time axis (more robust than hardcoding)
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
    cfg_fooof.foilim     = [2 40];
    cfg_fooof.pad        = 5;
    cfg_fooof.output     = 'fooof';
    cfg_fooof.keeptrials = 'no';   % average across trials before FOOOF

    for tfr_conds = 1:3

        if tfr_conds == 1
            trlIdx   = ind2;
            condName = 'set2';
        elseif tfr_conds == 2
            trlIdx   = ind4;
            condName = 'set4';
        else
            trlIdx   = ind6;
            condName = 'set6';
        end

        disp(' ')
        disp(['Running sliding-window FOOOF for condition ', condName])

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
            disp(['Condition:  ', num2str(tfr_conds), ' / 3'])
            disp(['Time Point: ', num2str(timePnt), ' / ', num2str(nTimePnts)])

            % Run FOOOF on averaged spectrum
            out = evalc('fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win);');

            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};
            else
                repdata = fooof_out.fooofparams;
            end

            freq_now = fooof_out.freq(:);

            % Map onto master grid (robust to any tiny numeric shifts)
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

    tfr2_fooof = tfr_fooof{1};
    tfr4_fooof = tfr_fooof{2};
    tfr6_fooof = tfr_fooof{3};

    disp(upper('FOOOF done on trial-averaged spectra (sliding windows)...'))

    %% Sanity check (3-panel): input, input-aperiodic, model-aperiodic (stored)
    time_point = 0.5;
    [~, tim]   = min(abs(tfr2_fooof.time - time_point));

    latWin = [toi_start(tim) (toi_start(tim) + winLen)];

    cond_titles = {'Set size 2','Set size 4','Set size 6'};
    trl_all     = {ind2, ind4, ind6};

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

    sgtitle(sprintf('FOOOF sanity check (AOC): Subject %s | Window [%.2f %.2f] s | t = %.2f s', ...
        subjects{subj}, latWin(1), latWin(2), tfr2_fooof.time(tim)), 'FontSize', 18)

    if ispc
        savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
    else
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
    end
    if ~exist(savePathControls, 'dir')
        mkdir(savePathControls)
    end

    saveName = sprintf('AOC_controls_FOOOF_modelMinusAperiodic_stern_subj%s.png', subjects{subj});
    saveas(gcf, fullfile(savePathControls, saveName));

    %% Baselined FOOOF TFR
    cfg              = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr2_fooof_bl    = ft_freqbaseline(cfg, tfr2_fooof);
    tfr4_fooof_bl    = ft_freqbaseline(cfg, tfr4_fooof);
    tfr6_fooof_bl    = ft_freqbaseline(cfg, tfr6_fooof);

    disp(upper('FOOOF baseline done...'))

    %% Save TFR data
    cd(datapath)
    save tfr_stern ...
        tfr2 tfr4 tfr6 ...
        tfr2_fooof tfr4_fooof tfr6_fooof ...
        tfr2_bl tfr4_bl tfr6_bl ...
        tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl

    %% Convert TFR data to POWSPCTRM (channels x frequency)
    analysisPeriodFull  = [0 2];
    analysisPeriodEarly = [0 1];
    analysisPeriodLate  = [1 2];
    freq_range          = [2 40];

    pow2_fooof          = select_data(analysisPeriodFull,  freq_range, tfr2_fooof);
    pow2_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr2_fooof_bl);

    pow4_fooof          = select_data(analysisPeriodFull,  freq_range, tfr4_fooof);
    pow4_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr4_fooof_bl);
    pow4_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr4_fooof_bl);
    pow4_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr4_fooof_bl);

    pow6_fooof          = select_data(analysisPeriodFull,  freq_range, tfr6_fooof);
    pow6_fooof_bl       = select_data(analysisPeriodFull,  freq_range, tfr6_fooof_bl);
    pow6_fooof_bl_early = select_data(analysisPeriodEarly, freq_range, tfr6_fooof_bl);
    pow6_fooof_bl_late  = select_data(analysisPeriodLate,  freq_range, tfr6_fooof_bl);

    pow2_fooof          = remove_time_dimension(pow2_fooof);
    pow2_fooof_bl       = remove_time_dimension(pow2_fooof_bl);
    pow2_fooof_bl_early = remove_time_dimension(pow2_fooof_bl_early);
    pow2_fooof_bl_late  = remove_time_dimension(pow2_fooof_bl_late);

    pow4_fooof          = remove_time_dimension(pow4_fooof);
    pow4_fooof_bl       = remove_time_dimension(pow4_fooof_bl);
    pow4_fooof_bl_early = remove_time_dimension(pow4_fooof_bl_early);
    pow4_fooof_bl_late  = remove_time_dimension(pow4_fooof_bl_late);

    pow6_fooof          = remove_time_dimension(pow6_fooof);
    pow6_fooof_bl       = remove_time_dimension(pow6_fooof_bl);
    pow6_fooof_bl_early = remove_time_dimension(pow6_fooof_bl_early);
    pow6_fooof_bl_late  = remove_time_dimension(pow6_fooof_bl_late);

    save power_stern_fooof ...
        pow2_fooof pow4_fooof pow6_fooof ...
        pow2_fooof_bl pow4_fooof_bl pow6_fooof_bl ...
        pow2_fooof_bl_early pow4_fooof_bl_early pow6_fooof_bl_early ...
        pow2_fooof_bl_late pow4_fooof_bl_late pow6_fooof_bl_late

    clc
    fprintf('Subject AOC %s (%.3d/%.3d) DONE (sliding-window FOOOF: model - aperiodic) \n', ...
        num2str(subjects{subj}), subj, length(subjects))
end
