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
        ind1 = find(data.trialinfo == 21);
        ind1 = find(data.trialinfo == 22);
        ind3 = find(data.trialinfo == 23);

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

%% TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM
% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

% Read data, segment and convert to FieldTrip data structure
for subj = 1 : length(subjects)

    % Check existing data
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    clc
    disp(['Processing TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM (NBACK) for Subject AOC ', num2str(subjects{subj})])
    cd(datapath)
    close all
    disp('Loading EEG TFR data (NBACK)')
    load dataEEG_TFR_nback

    % Identify indices of trials belonging to conditions
    % trialinfo codes: 21, 22, 23 (as in your old script)
    ind1 = find(dataTFR.trialinfo(:, 1) == 21);
    ind2 = find(dataTFR.trialinfo(:, 1) == 22);
    ind3 = find(dataTFR.trialinfo(:, 1) == 23);

    % Time-frequency analysis (keep nback timing, adapt style)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 3:1:30;                         % keep nback frequency range
    cfg.t_ftimwin    = ones(length(cfg.foi), 1).*0.5;  % 0.5 s windows
    cfg.toi          = -1:0.05:2.25;                   % keep nback time axis
    cfg.keeptrials   = 'no';
    cfg.pad          = 'nextpow2';

    cfg.trials = ind1;
    tfr1       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind2;
    tfr2       = ft_freqanalysis(cfg, dataTFR);
    cfg.trials = ind3;
    tfr3       = ft_freqanalysis(cfg, dataTFR);

    % Baselined TFR (raw powspctrm)
    cfg                      = [];
    cfg.baseline             = [-.5 -.25];             % keep nback baseline
    cfg.baselinetype         = 'db';
    tfr1_bl                  = ft_freqbaseline(cfg, tfr1);
    tfr2_bl                  = ft_freqbaseline(cfg, tfr2);
    tfr3_bl                  = ft_freqbaseline(cfg, tfr3);

    %%%%%%%%%%%%%%%%%%%%
    %%%%%%  FOOOF %%%%%%
    %%%%%%%%%%%%%%%%%%%%

    % Match nback epoch: -1 to 2.25 s, 500 ms windows, 50 ms steps
    startWin_FOOOF = [-1 -0.5];    % 500 ms window
    steps_FOOOF    = 0.05;         % 50 ms
    toi_FOOOF      = -1:steps_FOOOF:2.25;
    nTimePnts      = round((abs(toi_FOOOF(1)) + toi_FOOOF(end)) / steps_FOOOF) + 1;

    toi_centres    = startWin_FOOOF(1):steps_FOOOF:startWin_FOOOF(1) + steps_FOOOF*(nTimePnts-1);
    toi_centres    = toi_centres + 0.25;  % shift by half window

    % Prepare containers per condition
    tfr_fooof_trl = cell(1, 3);
    tfr_fooof_avg = cell(1, 3);

    % Conditions
    for tfr_conds = 1 : 3
        cfg   = [];
        allff = {};

        if tfr_conds == 1
            trlIdx = ind1;
        elseif tfr_conds == 2
            trlIdx = ind2;
        elseif tfr_conds == 3
            trlIdx = ind3;
        end

        % Loop over timepoints
        for timePnt = 1 : nTimePnts
            nTrl           = numel(trlIdx);
            tmpfreq        = cell(nTrl, 1);
            tmpfooofparams = cell(nTrl, 1);

            % Select data segment for FOOOF window
            cfg         = [];
            cfg.latency = startWin_FOOOF + steps_FOOOF * (timePnt-1);
            datTFR      = ft_selectdata(cfg, dataTFR);

            % Loop over trials (parallel)
            parfor trl = 1 : numel(trlIdx)
                disp(['Subject    ' num2str(subj)])
                disp(['Condition  ' num2str(tfr_conds)])
                disp(['Timepoint  ' num2str(timePnt)])
                disp(['Trial      ' num2str(trl)])

                % FOOOF FieldTrip configs
                cfg            = [];
                cfg.method     = 'mtmfft';
                cfg.taper      = 'hanning';
                cfg.foilim     = [3 30];   % match nback frequency range
                cfg.pad        = 5;
                cfg.output     = 'fooof';
                cfg.trials     = trlIdx(trl);

                % Run FOOOF via your wrapper
                tmpfreq{trl}        = ft_freqanalysis_Arne_FOOOF(cfg, datTFR);
                tmpfooofparams{trl} = tmpfreq{trl}.fooofparams; % save fooofparams
            end

            % Compute avg over trials (keepindividual = yes)
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
            ff_foof{timePnt}                = allff{timePnt}.fooofparams;
        end

        % Extract aperiodic signal and full power spectrum
        tmp_aperiodic  = [];
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
                tmp_aperiodic(trl, :, 4, timePnt)  = datafit(1, :, 2);   % r-squared
                power_spectrum(trl, :, :, timePnt) = pwr_spec;
            end
        end

        % Save FOOOF output: trl x chan x freq x time
        tfr_ff_trl.dimord         = 'rpt_chan_freq_time';
        tfr_ff_trl.powspctrm      = fooofed_power;
        tfr_ff_trl.power_spectrum = power_spectrum;
        tfr_ff_trl.trialinfo      = trlIdx;
        tfr_ff_trl.fooofparams    = tmp_aperiodic;
        tfr_ff_trl.time           = toi_centres;

        % Average over trials: chan x freq x time
        tmp_for_avg = tfr_ff_trl;
        tmp_for_avg = rmfield(tmp_for_avg, {'fooofparams', 'power_spectrum'});

        cfg           = [];
        cfg.parameter = 'powspctrm';
        tfr_ff_avg    = ft_freqdescriptives(cfg, tmp_for_avg);

        % Average aperiodic parameters and power_spectrum across trials
        aperiodic_mean = squeeze(mean(tfr_ff_trl.fooofparams, 1, 'omitnan'));    % chan x 4 x time
        powspec_mean   = squeeze(mean(tfr_ff_trl.power_spectrum, 1, 'omitnan')); % chan x freq x time

        % Struct with only averages
        tfr_ff                    = [];
        tfr_ff.fooofparams        = aperiodic_mean;
        tfr_ff.power_spectrum     = powspec_mean;
        tfr_ff.powspctrm          = tfr_ff_avg.powspctrm;
        tfr_ff.label              = tfr_ff_trl.label;
        tfr_ff.freq               = tfr_ff_trl.freq;
        tfr_ff.time               = tfr_ff_trl.time;
        tfr_ff.dimord             = 'chan_freq_time';

        % Assign to sliced outputs
        tfr_fooof_trl{tfr_conds} = tfr_ff_trl;   % full rpt_chan_freq_time struct
        tfr_fooof_avg{tfr_conds} = tfr_ff;       % averages only
    end

    % Unpack per condition
    tfr1_fooof_trl   = tfr_fooof_trl{1};
    tfr2_fooof_trl   = tfr_fooof_trl{2};
    tfr3_fooof_trl   = tfr_fooof_trl{3};

    tfr1_fooof       = tfr_fooof_avg{1};
    tfr2_fooof       = tfr_fooof_avg{2};
    tfr3_fooof       = tfr_fooof_avg{3};
    disp(upper('FOOOF done...'))

    %% Sanity Check: averaged FOOOF output, all channels, all three conditions
    time_point = 0.5; % time (s) to inspect

    % Collect condition data (already averaged over trials)
    tfr_all     = {tfr1_fooof, tfr2_fooof, tfr3_fooof};
    cond_titles = {'Cond 1', 'Cond 2', 'Cond 3'};  % adapt labels if you like (0-back, 1-back, 2-back)

    % Use time axis from first condition
    [~, tim] = min(abs(tfr1_fooof.time - time_point));
    freq     = tfr1_fooof.freq;

    figure('Position', [0 0 1512 500], 'Color', 'w');

    for c = 1 : 3
        tfr_cond = tfr_all{c};

        % powspctrm: chan x freq x time
        raw_spec = squeeze( ...
                        mean( ...
                            tfr_cond.powspctrm(:, :, tim), ...
                        1, 'omitnan') ...
                    );                     % -> freq

        % power_spectrum: chan x freq x time (already log10 power from FOOOF)
        model_spec = squeeze( ...
                        mean( ...
                            tfr_cond.power_spectrum(:, :, tim), ...
                        1, 'omitnan') ...
                      );                   % -> freq

        % fooofparams: chan x 4 x time, dim 2: 1 = intercept, 2 = slope
        offset_vec = squeeze(tfr_cond.fooofparams(:, 1, tim)); % chan
        slope_vec  = squeeze(tfr_cond.fooofparams(:, 2, tim)); % chan
        offset     = mean(offset_vec, 'omitnan');
        slope      = mean(slope_vec,  'omitnan');

        % log-transform raw spectrum so everything lives in log10(power)
        raw_log       = log10(raw_spec);
        model_log     = model_spec;                 % already log10(power)
        aperiodic_fit = offset - slope .* log10(freq);

        subplot(1, 3, c);
        plot(freq, raw_log, 'LineWidth', 3)
        hold on
        plot(freq, model_log, 'LineWidth', 3)
        plot(freq, aperiodic_fit, 'LineWidth', 3, 'LineStyle', '--')

        xlabel('Frequency (Hz)')
        if c == 1
            ylabel('Power (log_{10})')
            legend({'Raw Power', 'Final Fit', 'Aperiodic Fit'}, 'Location', 'best')
        end
        set(gca, 'FontSize', 15)
        title(sprintf('%s | t = %.2f s', ...
              cond_titles{c}, tfr_cond.time(tim)), 'FontSize', 16)
    end

    sgtitle(sprintf('FOOOF sanity check (NBACK): Subject %s', ...
            subjects{subj}), 'FontSize', 20)

    % Save figure (adapt path if you want a separate NBACK folder)
    if ispc
        savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
    else
        savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
    end
    saveName = sprintf('AOC_controls_FOOOF_powspctrm_nback_subj%s.png', subjects{subj});
    saveas(gcf, fullfile(savePathControls, saveName));

    %% FOOOFed powspctrm baselined
    cfg                      = [];
    cfg.baseline             = [-.5 -.25];      % same as raw baseline
    cfg.baselinetype         = 'absolute';      % FOOOF already in log space
    tfr1_fooof_bl            = ft_freqbaseline(cfg, tfr1_fooof);
    tfr2_fooof_bl            = ft_freqbaseline(cfg, tfr2_fooof);
    tfr3_fooof_bl            = ft_freqbaseline(cfg, tfr3_fooof);

    % Save data
    cd(datapath)
    save tfr_nback ...
        tfr1 tfr2 tfr3 ...
        tfr1_fooof tfr2_fooof tfr3_fooof ...
        tfr1_fooof_trl tfr2_fooof_trl tfr3_fooof_trl ...
        tfr1_bl tfr2_bl tfr3_bl ...
        tfr1_fooof_bl tfr2_fooof_bl tfr3_fooof_bl

    %% Convert TFR data to POWSPCTRM (channels x frequency)
    analysisPeriodFull  = [0 2];
    analysisPeriodEarly = [0 1];
    analysisPeriodLate  = [1 2];
    freq_range          = [2 40];   % keep as in your previous scripts

    % Select data
    pow1_fooof              = select_data(analysisPeriodFull,  freq_range, tfr1_fooof);
    pow1_fooof_bl           = select_data(analysisPeriodFull,  freq_range, tfr1_fooof_bl);
    pow1_fooof_bl_early     = select_data(analysisPeriodEarly, freq_range, tfr1_fooof_bl);
    pow1_fooof_bl_late      = select_data(analysisPeriodLate,  freq_range, tfr1_fooof_bl);

    pow2_fooof              = select_data(analysisPeriodFull,  freq_range, tfr2_fooof);
    pow2_fooof_bl           = select_data(analysisPeriodFull,  freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_early     = select_data(analysisPeriodEarly, freq_range, tfr2_fooof_bl);
    pow2_fooof_bl_late      = select_data(analysisPeriodLate,  freq_range, tfr2_fooof_bl);

    pow3_fooof              = select_data(analysisPeriodFull,  freq_range, tfr3_fooof);
    pow3_fooof_bl           = select_data(analysisPeriodFull,  freq_range, tfr3_fooof_bl);
    pow3_fooof_bl_early     = select_data(analysisPeriodEarly, freq_range, tfr3_fooof_bl);
    pow3_fooof_bl_late      = select_data(analysisPeriodLate,  freq_range, tfr3_fooof_bl);

    % Remove time dimension for POWSPCTRM (channels x frequency)
    pow1_fooof              = remove_time_dimension(pow1_fooof);
    pow1_fooof_bl           = remove_time_dimension(pow1_fooof_bl);
    pow1_fooof_bl_early     = remove_time_dimension(pow1_fooof_bl_early);
    pow1_fooof_bl_late      = remove_time_dimension(pow1_fooof_bl_late);

    pow2_fooof              = remove_time_dimension(pow2_fooof);
    pow2_fooof_bl           = remove_time_dimension(pow2_fooof_bl);
    pow2_fooof_bl_early     = remove_time_dimension(pow2_fooof_bl_early);
    pow2_fooof_bl_late      = remove_time_dimension(pow2_fooof_bl_late);

    pow3_fooof              = remove_time_dimension(pow3_fooof);
    pow3_fooof_bl           = remove_time_dimension(pow3_fooof_bl);
    pow3_fooof_bl_early     = remove_time_dimension(pow3_fooof_bl_early);
    pow3_fooof_bl_late      = remove_time_dimension(pow3_fooof_bl_late);

    save power_nback_fooof ...
        pow1_fooof pow2_fooof pow3_fooof ...
        pow1_fooof_bl pow2_fooof_bl pow3_fooof_bl ...
        pow1_fooof_bl_early pow2_fooof_bl_early pow3_fooof_bl_early ...
        pow1_fooof_bl_late pow2_fooof_bl_late pow3_fooof_bl_late

    clc
    disp('TFR and FOOOFed POWSPCTRM (NBACK) COMPUTED...');
end