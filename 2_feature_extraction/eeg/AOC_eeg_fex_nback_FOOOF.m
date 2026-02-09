%% TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM (NBACK) | AOC
% FOOOF output here is: (model fit - aperiodic) in FOOOF/log space (NOT peaks-only)

clear; close all; clc

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

% Setup logging
logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
scriptName = 'AOC_eeg_fex_nback_FOOOF';

%% Loop subjects
for subj = 1:length(subjects)
    try
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

        % Exclude bad FOOOF fits (r_squared < 0.90)
        r_sq_thresh = 0.90;
        bad_fits = squeeze(fooof_aperiodic(:, 4, :)) < r_sq_thresh;
        n_bad = sum(bad_fits(:));
        if n_bad > 0
            fprintf('  Excluding %d / %d channel x time FOOOF fits (r^2 < %.2f)\n', ...
                n_bad, numel(bad_fits), r_sq_thresh);
            for tp = 1:nTimePnts
                fooof_powspctrm(bad_fits(:, tp), :, tp) = NaN;
                fooof_powspec(bad_fits(:, tp), :, tp) = NaN;
            end
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

    figure('Position', [0 0 1500 420], 'Color', 'w');

    cfg_sel         = [];
    cfg_sel.latency = latWin;
    cfg_sel.trials  = ind1;
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

    ps_in      = mean(ps_all,     2, 'omitnan');
    ap_mean    = mean(ap_all,     2, 'omitnan');
    model_mean = mean(model_all,  2, 'omitnan');
    pk_comp    = mean(pkcomp_all, 2, 'omitnan');
    ps_corr    = ps_in - ap_mean;

    subplot(1,3,1)
    plot(freq_sc, ps_in, 'k', 'LineWidth', 2); hold on
    plot(freq_sc, ap_mean, 'b--', 'LineWidth', 2)
    plot(freq_sc, model_mean, 'r', 'LineWidth', 2); hold off
    legend('Original', 'Aperiodic', 'Model fit', 'Location', 'best')
    title('Original Powspctrm')
    xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

    subplot(1,3,2)
    plot(freq_sc, ps_corr, 'LineWidth', 2)
    title('Original - Aperiodic')
    xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

    subplot(1,3,3)
    plot(freq_sc, pk_comp, 'LineWidth', 2)
    title('Model fit - Aperiodic (stored)')
    xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
    set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

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

    saveName = sprintf('AOC_controls_FOOOF_nback_subj%s.png', subjects{subj});
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
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end
