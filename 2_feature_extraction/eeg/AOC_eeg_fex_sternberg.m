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

% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)

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

    % Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:1:40;                         % analysis 2 to 40 Hz in steps of 1 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1.5:0.05:3;
    cfg.keeptrials   = 'no';

    cfg.trials = ind2;
    tfr2       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind4;
    tfr4       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind6;
    tfr6       = ft_freqanalysis(cfg, dataTFR);

    % Baselined TFR
    % Raw powspctrm baselined
    cfg                              = [];
    cfg.baseline                     = [-.5 -.25];
    cfg.baselinetype                 = 'db';
    tfr2_bl                          = ft_freqbaseline(cfg, tfr2);
    tfr4_bl                          = ft_freqbaseline(cfg, tfr4);
    tfr6_bl                          = ft_freqbaseline(cfg, tfr6);

    %%%%%%%%%%%%%%%%%%%%
    %%%%%%  FOOOF %%%%%%
    %%%%%%%%%%%%%%%%%%%%
    startWin_FOOOF = [-1.5 -1]; % window of 500ms
    steps_FOOOF = 0.05; % 50 ms
    toi_FOOOF = -1.5:0.05:3; % -1.5 to 3 s with in 50s steps
    nTimePnts = round((abs(toi_FOOOF(1))+toi_FOOOF(end))/steps_FOOOF)+1;

    % Parallel processing of conditions
    for tfr_conds = 1 : 3 %parfor
        tmpfreq        = [];
        tmpfooofparams = [];
        allff          = [];
        if tfr_conds == 1
            trlIdx = ind2;
            cfg.trials = trlIdx;
        elseif tfr_conds == 2
            trlIdx = ind4;
            cfg.trials = trlIdx;
        elseif tfr_conds == 3
            trlIdx = ind6;
            cfg.trials = trlIdx;
        end

        % Loop over timepoints
        for timePnt = 1 : nTimePnts

            % Select data
            cfg         = [];
            cfg.latency = startWin_FOOOF + 0.05 * (timePnt-1);
            datTFR      = ft_selectdata(cfg, dataTFR);

            % Loop over trials
            for trl = 1 : numel(trlIdx)
                clc
                disp(['Subject    ' num2str(subj)])
                disp(['Condition  ' num2str(tfr_conds)])
                disp(['Timepoint  ' num2str(timePnt)])
                disp(['Trial      ' num2str(trl)])

                % FOOOF FieldTrip configs
                cfg            = [];
                cfg.method     = 'mtmfft';
                cfg.taper      = 'hanning';
                cfg.foilim     = [2 40];
                cfg.pad        = 5;
                cfg.output     = 'fooof';
                cfg.trials     = trlIdx(trl);

                % Run FOOOF (settings Marius)
                tmpfreq{trl} = ft_freqanalysis_Arne_FOOOF(cfg, datTFR); 
                tmpfooofparams{trl, 1}  =  tmpfreq{trl}.fooofparams; % save fooof params

                % Check fit
                %mean(cell2mat({tmpfreq{trl}.fooofparams.r_squared}))
            end
            % Compute avg over trials
            cfg                = [];
            cfg.keepindividual = 'yes';
            ff                 = ft_freqgrandaverage(cfg, tmpfreq{:});
            ff.fooofparams     = tmpfooofparams;
            ff.cfg             = [];
            allff{timePnt}     = ff;
        end

        % Concatenate data over time points
        tfr_ff_trl    = tmpfreq{1};
        fooofed_power = [];
        ff_foof       = {};
        for timePnt = 1 : length(allff)
            fooofed_power(:, :, :, timePnt) = allff{timePnt}.powspctrm;
            ff_foof{timePnt} = allff{timePnt}.fooofparams;
        end

        % Extract aperiodic signal
        tmp_aperiodic = [];
        power_spectrum = [];
        for timePnt = 1 : length(ff_foof)
            tdata = ff_foof{timePnt};
            for trl = 1 : length(tdata)
                repdata       = tdata{trl};
                tmpaperdiodic = {repdata.aperiodic_params};
                tmperror      = {repdata.error};
                tmpr_sq       = {repdata.r_squared};
                tmp_pwr_spec  = {repdata.power_spectrum};
                elec_data     = [];
                datafit       = [];
                pwr_spec      = [];
                for electrode = 1 : size(tfr_ff_trl.label, 2)
                    elec_data(1, electrode, :) = tmpaperdiodic{electrode};
                    datafit(1, electrode, 1)   = tmperror{electrode};
                    datafit(1, electrode, 2)   = tmpr_sq{electrode};
                    pwr_spec(electrode, :)     = tmp_pwr_spec{electrode};
                end
                tmp_aperiodic(trl, :, 1, timePnt)  = elec_data(1, :, 1); % intercept
                tmp_aperiodic(trl, :, 2, timePnt)  = elec_data(1, :, 2); % slope
                tmp_aperiodic(trl, :, 3, timePnt)  = datafit(1, :, 1);   % error
                tmp_aperiodic(trl, :, 4, timePnt)  = datafit(1, :, 2);   % r squared
                power_spectrum(trl, :, :, timePnt) = pwr_spec;
            end
        end
        % Save FOOOF output: trl x chan x freq x time
        tfr_ff_trl.dimord         = 'rpt_chan_freq_time';
        tfr_ff_trl.powspctrm      = fooofed_power;
        tfr_ff_trl.power_spectrum = power_spectrum;
        tfr_ff_trl.trialinfo      = trlIdx;
        tfr_ff_trl.fooofparams    = tmp_aperiodic;
        tfr_ff_trl.time           = toi_FOOOF;

        % Average over trials: chan x freq x time
        cfg                 = [];
        cfg.parameter       = 'powspctrm';
        tfr_ff_avg          = ft_freqdescriptives(cfg, tfr_ff_trl);

        % Average aperiodic parameters across trials as well
        aperiodic_avg       = squeeze(mean(tfr_ff_trl.fooofparams, 1, 'omitnan'));   % chan x 4 x time
        power_spectrum_avg  = squeeze(mean(tfr_ff_trl.power_spectrum, 1, 'omitnan')); % chan x freq x time
        tfr_ff                    = [];
        tfr_ff.avg_fooofparams    = aperiodic_avg;
        tfr_ff.avg_power_spectrum = power_spectrum_avg;
        tfr_ff.avg_powspctrm      = tfr_ff_avg.powspctrm;

        % Assign to sliced outputs for parfor
        tfr_fooof_trl{tfr_conds} = tfr_ff_trl;
        tfr_fooof_avg{tfr_conds} = tfr_ff;
        aperiodic_avg{tfr_conds} = aperiodic_mean;
        powspec_avg{tfr_conds}   = powspec_mean;

        % Save
        saveName = sprintf('AOC_controls_FOOOF_powspctrm_subj%s_cond1_ch%s_t%d.png', ...
            subjects{subj}, tfr_cond.label{chan}, tim);
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
        saveas(gcf, fullfile(savePathControls, saveName));
    end
    % after parfor: unpack
    tfr2_fooof = tfr_fooof_trl{1};
    tfr4_fooof = tfr_fooof_trl{2};
    tfr6_fooof = tfr_fooof_trl{3};
    disp(upper('FOOOF done...'))

    %% Sanity Check
    % example sanity check for condition 1 (set size 2) after parfor
    tfr_cond   = tfr2_fooof;         % rpt_chan_freq_time
    chan_label = 'Oz';               % whatever channel you want
    chan       = find(strcmp(tfr_cond.label, chan_label));

    [~, tim]   = min(abs(tfr_cond.time - 0.5));
    freq       = tfr_cond.freq;

    % average across trials for that channel & time
    raw_spec   = squeeze(mean(tfr_cond.powspctrm(:, chan, :, tim), 1));           % freq
    model_spec = squeeze(mean(tfr_cond.power_spectrum(:, chan, :, tim), 1));      % freq

    offset     = squeeze(mean(tfr_cond.fooofparams(:, chan, 1, tim), 1));         % intercept
    slope      = squeeze(mean(tfr_cond.fooofparams(:, chan, 2, tim), 1));         % slope

    aperiodic_fit = offset - slope .* log10(freq);

    figure('Position', [0 0 1512/2 982], 'Color', 'w');
    plot(freq, log10(raw_spec), 'LineWidth', 3)
    hold on
    plot(freq, model_spec, 'LineWidth', 3)
    plot(freq, aperiodic_fit, 'LineWidth', 3, 'LineStyle', '--');
    ylabel('Power (log_{10})')
    xlabel('Frequency (Hz)')
    set(gca, 'FontSize', 15)
    legend({'Raw Power', 'Final Fit', 'Aperiodic Fit'}, 'Location', 'best')
    title(sprintf('Powspctrm: Subject %s | Cond 1 (set size 2) | t = %.2f s', ...
        subjects{subj}, tfr_cond.time(tim)), 'FontSize', 20)

    %% FOOOFed powspctrm baselined
    cfg                              = [];
    cfg.baseline                     = [-.5 -.25];
    cfg.baselinetype                 = 'absolute';   % FOOOF already sets log scale, so no 'dB' here
    tfr2_fooof_bl                    = ft_freqbaseline(cfg, tfr2_fooof);
    tfr4_fooof_bl                    = ft_freqbaseline(cfg, tfr4_fooof);
    tfr6_fooof_bl                    = ft_freqbaseline(cfg, tfr6_fooof);

    % Save data
    cd(datapath)
    save tfr_stern ...
        tfr2 tfr4 tfr6 ...
        tfr2_fooof tfr4_fooof tfr6_fooof ...
        tfr2_bl tfr4_bl tfr6_bl ...
        tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl

    %% Convert TFR data to POWSCPTRM (channels x frequency)
    analysisPeriodFull = [0 2];
    analysisPeriodEarly = [0 1];
    analysisPeriodLate = [1 2];
    freq_range = [2 40];

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