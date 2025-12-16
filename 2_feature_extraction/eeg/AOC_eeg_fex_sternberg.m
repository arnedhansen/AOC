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

%% TFR (Raw, FOOOF and Baselined) POWSPCTRM
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Define parfor output function
D = parallel.pool.DataQueue;
afterEach(D, @printProgress);
function printProgress(s)
    fprintf('Subj %d | Cond %d | Time %d/%d finished\n', ...
        s.subj, s.cond, s.time, s.nTimePnts);
end

% Read data, segment and convert to FieldTrip data structure
for subj = 1 %%%%% : length(subjects)
    D = parallel.pool.DataQueue;
    afterEach(D, @(x) fprintf('Timepoint %d finished\n', x));

    % Check existing data
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    clc
    disp(['Processing TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM for Subject AOC ', num2str(subjects{subj})])
    cd(datapath)
    close all
    disp('Loading EEG TFR data')
    load dataEEG_TFR_sternberg

    % Identify indices of trials belonging to conditions
    ind2 = find(dataTFR.trialinfo(:, 1) == 22);
    ind4 = find(dataTFR.trialinfo(:, 1) == 24);
    ind6 = find(dataTFR.trialinfo(:, 1) == 26);

    % Time frequency analysis (trial-averaged TFR for visualisation / later use)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:1:40;                         % analysis 2 to 40 Hz in steps of 1 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1.5:0.05:3;
    cfg.keeptrials   = 'no';
    cfg.pad          = 'nextpow2';

    cfg.trials = ind2;
    tfr2       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind4;
    tfr4       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind6;
    tfr6       = ft_freqanalysis(cfg, dataTFR);

    % Baselined TFR
    % Raw powspctrm baselined
    cfg                      = [];
    cfg.baseline             = [-.5 -.25];
    cfg.baselinetype         = 'db';
    tfr2_bl                  = ft_freqbaseline(cfg, tfr2);
    tfr4_bl                  = ft_freqbaseline(cfg, tfr4);
    tfr6_bl                  = ft_freqbaseline(cfg, tfr6);

    %%%%%%  FOOOF %%%%%%
    % Sliding-window FOOOF over trial-averaged spectra
    startWin_FOOOF = [-1.5 -1];    % 500 ms window
    steps_FOOOF    = 0.05;         % 50 ms step
    toi_FOOOF      = -1.5:0.05:3;  % full time range for 50 ms steps
    nTimePnts      = round((abs(toi_FOOOF(1)) + toi_FOOOF(end)) / steps_FOOOF) + 1;

    toi_centres    = startWin_FOOOF(1):steps_FOOOF:startWin_FOOOF(1) + steps_FOOOF*(nTimePnts-1);
    toi_centres    = toi_centres + 0.25;  % shift by half window

    % Container for each condition
    tfr_fooof = cell(1, 3);

    % Conditions
    for tfr_conds = 1 : 3

        if tfr_conds == 1
            trlIdx = ind2;
        elseif tfr_conds == 2
            trlIdx = ind4;
        elseif tfr_conds == 3
            trlIdx = ind6;
        end

        disp(' ')
        disp(['Running FOOOF on trial-averaged spectra for condition ' num2str(tfr_conds)])

        % Prepare FOOOF config
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = [2 40];
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';        % average across trials before FOOOF

        % One test window to get sizes
        cfg_sel0         = [];
        cfg_sel0.latency = startWin_FOOOF;      % first 500 ms window
        cfg_sel0.trials  = trlIdx;
        datTFR_win0      = ft_selectdata(cfg_sel0, dataTFR);
        fooof_test       = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win0);
        nChan            = numel(fooof_test.label);
        nFreq            = numel(fooof_test.freq);
        fooof_powspctrm = nan(nChan, nFreq, nTimePnts);
        fooof_powspec   = nan(nChan, nFreq, nTimePnts);
        fooof_aperiodic = nan(nChan, 4, nTimePnts);  % [offset slope error r^2]

        % parfor over timepoints
        parfor timePnt = 1 : nTimePnts

            % Select data window and trials for this condition
            cfg_sel         = [];
            cfg_sel.latency = startWin_FOOOF + steps_FOOOF * (timePnt-1);
            cfg_sel.trials  = trlIdx;
            datTFR_win      = ft_selectdata(cfg_sel, dataTFR);

            % Run FOOOF on the averaged spectrum
            fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win);

            % Store FOOOFed power (chan x freq)
            local_pow = fooof_out.powspctrm;      % chan x freq

            % Extract FOOOF parameters per channel
            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};   % single averaged "trial"
            else
                repdata = fooof_out.fooofparams;
            end

            tmpaperdiodic = {repdata.aperiodic_params};
            tmperror      = {repdata.error};
            tmpr_sq       = {repdata.r_squared};
            tmp_pwr_spec  = {repdata.power_spectrum};

            local_aper = nan(nChan, 2);  % [offset slope]
            local_err  = nan(nChan, 1);
            local_rsq  = nan(nChan, 1);
            local_ps   = nan(nChan, nFreq);

            for electrode = 1 : nChan
                aper_params              = tmpaperdiodic{electrode};   % [offset slope]
                local_aper(electrode, :) = aper_params(:).';

                local_err(electrode, 1)  = tmperror{electrode};
                local_rsq(electrode, 1)  = tmpr_sq{electrode};
                local_ps(electrode, :)   = tmp_pwr_spec{electrode};
            end

            % Write into sliced arrays (parfor-friendly)
            fooof_powspctrm(:, :, timePnt) = local_pow;
            fooof_powspec(:, :, timePnt)   = local_ps;

            % Pack aperiodic parameters into one [chan x 4] matrix
            local_aper_all        = nan(nChan, 4);
            local_aper_all(:, 1)  = local_aper(:, 1);   % intercept
            local_aper_all(:, 2)  = local_aper(:, 2);   % slope
            local_aper_all(:, 3)  = local_err;          % error
            local_aper_all(:, 4)  = local_rsq;          % r squared

            % Single consistent sliced write for parfor
            fooof_aperiodic(:, :, timePnt) = local_aper_all;

            % Output
            s           = struct();
            s.subj      = subj;
            s.cond      = tfr_conds;
            s.time      = timePnt;
            s.nTimePnts = nTimePnts;
            send(D, s);
        end

        % Construct condition-specific FOOOF struct
        tfr_ff                    = [];
        tfr_ff.label              = fooof_test.label;
        tfr_ff.freq               = fooof_test.freq;
        tfr_ff.time               = toi_centres(1:nTimePnts);
        tfr_ff.powspctrm          = fooof_powspctrm;   % FOOOFed TFR: chan x freq x time
        tfr_ff.power_spectrum     = fooof_powspec;     % full model spectrum: chan x freq x time
        tfr_ff.fooofparams        = fooof_aperiodic;   % chan x 4 x time
        tfr_ff.dimord             = 'chan_freq_time';

        tfr_fooof{tfr_conds}      = tfr_ff;
    end

    % Assign condition-specific outputs
    tfr2_fooof = tfr_fooof{1};   % averages (aperiodic + spectrum) cond 1
    tfr4_fooof = tfr_fooof{2};   % averages cond 2
    tfr6_fooof = tfr_fooof{3};   % averages cond 3
    disp(upper('FOOOF done on trial-averaged spectra...'))

    %% Sanity Check: averaged FOOOF output, all channels, all three conditions
    time_point = 0.5;
    [~, tim] = min(abs(tfr2_fooof.time - time_point));

    tfr_all     = {tfr2_fooof, tfr4_fooof, tfr6_fooof};
    cond_titles = {'Set size 2','Set size 4','Set size 6'};

    figure('Position', [0 0 1512 500], 'Color', 'w');

    for c = 1:3
        tfr_cond = tfr_all{c};
        freq     = tfr_cond.freq;

        % pick a channel to avoid mixing spaces across channels at first
        ch = 1;

        % recover the repdata for THIS condition/timepoint:
        % (best is to store repdata per timepoint, but for a quick sanity plot
        %  you can rerun one window or just plot what you already stored if you saved repdata)
        %
        % Here: you need repdata from the same time window that produced tim.
        % So easiest is: rerun ONE window here, for ONE condition, ONE subject.
        %
        % For now assume you still have repdata in workspace from the last fooof_out call.
        ps_in = repdata(ch).power_spectrum(:);

        ap = repdata(ch).aperiodic_params(:);
        if numel(ap) == 2
            offset = ap(1);
            expo   = ap(2);
            ap_fit = offset - expo .* log10(freq(:));
        elseif numel(ap) == 3
            offset = ap(1);
            knee   = ap(2);
            expo   = ap(3);
            ap_fit = offset - log10(knee + freq(:).^expo);
        else
            ap_fit = nan(numel(freq), 1);
        end

        pk = repdata(ch).peak_params;
        gauss_sum = zeros(numel(freq), 1);

        if ~isempty(pk)
            for p = 1:size(pk, 1)
                cf  = pk(p, 1);
                amp = pk(p, 2);
                bw  = pk(p, 3);

                % Option A: treat bw as Gaussian std (common)
                gauss = amp .* exp(-(freq(:) - cf).^2 ./ (2*bw.^2));

                % Option B (fallback): treat bw as FWHM -> convert to std
                % std = bw / (2*sqrt(2*log(2)));
                % gauss = amp .* exp(-(freq(:) - cf).^2 ./ (2*std.^2));

                gauss_sum = gauss_sum + gauss;
            end
        end

        model_fit = ap_fit + gauss_sum;

        subplot(1,3,c); hold on
        plot(freq, ps_in,     'LineWidth', 3)
        plot(freq, model_fit, 'LineWidth', 3)
        plot(freq, ap_fit,    '--', 'LineWidth', 3)

        xlabel('Frequency (Hz)')
        if c == 1
            ylabel('Power (FOOOF space)')
            legend({'Input spectrum','Model fit','Aperiodic fit'}, 'Location', 'best')
        end
        title(sprintf('%s | t = %.2f s', cond_titles{c}, tfr_cond.time(tim)))
        set(gca, 'FontSize', 15)
    end

    sgtitle(sprintf('FOOOF sanity check: Subject %s', ...
        subjects{subj}), 'FontSize', 20)

    % Save figure using same path logic
    if ispc
        savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
    else
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
    end
    saveName = sprintf('AOC_controls_FOOOF_powspctrm_subj%s.png', subjects{subj});
    saveas(gcf, fullfile(savePathControls, saveName));

    %% FOOOFed powspctrm baselined
    cfg                      = [];
    cfg.baseline             = [-.5 -.25];
    cfg.baselinetype         = 'absolute';   % FOOOF already sets log scale, so no 'dB' here
    tfr2_fooof_bl            = ft_freqbaseline(cfg, tfr2_fooof);
    tfr4_fooof_bl            = ft_freqbaseline(cfg, tfr4_fooof);
    tfr6_fooof_bl            = ft_freqbaseline(cfg, tfr6_fooof);

    % Save data
    cd(datapath)
    save tfr_stern ...
        tfr2 tfr4 tfr6 ...
        tfr2_fooof tfr4_fooof tfr6_fooof ...
        tfr2_bl tfr4_bl tfr6_bl ...
        tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl

    %% Convert TFR data to POWSCPTRM (channels x frequency)
    analysisPeriodFull  = [0 2];
    analysisPeriodEarly = [0 1];
    analysisPeriodLate  = [1 2];
    freq_range          = [2 40];

    % Select data
    pow2_fooof                                = select_data(analysisPeriodFull, freq_range, tfr2_fooof);
    pow2_fooof_bl                             = select_data(analysisPeriodFull, freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_early                       = select_data(analysisPeriodEarly, freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_late                        = select_data(analysisPeriodLate, freq_range, tfr2_fooof_bl);

    pow4_fooof                                = select_data(analysisPeriodFull, freq_range, tfr4_fooof);
    pow4_fooof_bl                             = select_data(analysisPeriodFull, freq_range, tfr4_fooof_bl);
    pow4_fooof_bl_early                       = select_data(analysisPeriodEarly, freq_range, tfr4_fooof_bl);
    pow4_fooof_bl_late                        = select_data(analysisPeriodLate, freq_range, tfr4_fooof_bl);

    pow6_fooof                                = select_data(analysisPeriodFull, freq_range, tfr6_fooof);
    pow6_fooof_bl                             = select_data(analysisPeriodFull, freq_range, tfr6_fooof_bl);
    pow6_fooof_bl_early                       = select_data(analysisPeriodEarly, freq_range, tfr6_fooof_bl);
    pow6_fooof_bl_late                        = select_data(analysisPeriodLate, freq_range, tfr6_fooof_bl);

    % Remove time dimension for POWSCPTRM (channels x frequency)
    pow2_fooof                                = remove_time_dimension(pow2_fooof);
    pow2_fooof_bl                             = remove_time_dimension(pow2_fooof_bl);
    pow2_fooof_bl_early                       = remove_time_dimension(pow2_fooof_bl_early);
    pow2_fooof_bl_late                        = remove_time_dimension(pow2_fooof_bl_late);

    pow4_fooof                                = remove_time_dimension(pow4_fooof);
    pow4_fooof_bl                             = remove_time_dimension(pow4_fooof_bl);
    pow4_fooof_bl_early                       = remove_time_dimension(pow4_fooof_bl_early);
    pow4_fooof_bl_late                        = remove_time_dimension(pow4_fooof_bl_late);

    pow6_fooof                                = remove_time_dimension(pow6_fooof);
    pow6_fooof_bl                             = remove_time_dimension(pow6_fooof_bl);
    pow6_fooof_bl_early                       = remove_time_dimension(pow6_fooof_bl_early);
    pow6_fooof_bl_late                        = remove_time_dimension(pow6_fooof_bl_late);

    save power_stern_fooof ...
        pow2_fooof pow4_fooof pow6_fooof ...
        pow2_fooof_bl pow4_fooof_bl pow6_fooof_bl ...
        pow2_fooof_bl_early pow4_fooof_bl_early pow6_fooof_bl_early ...
        pow2_fooof_bl_late pow4_fooof_bl_late pow6_fooof_bl_late

    clc
    disp('TFR and FOOOFed POWSPCTRM COMPUTED...');
end
