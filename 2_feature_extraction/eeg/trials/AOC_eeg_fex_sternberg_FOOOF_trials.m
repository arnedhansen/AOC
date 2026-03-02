%% TFR (Raw, FOOOF and Baselined) | Sternberg (AOC) | Trial-level
startup
[subjects, path, ~, ~] = setup('AOC');

% Setup logging
logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
scriptName = 'AOC_eeg_fex_sternberg_FOOOF_trials';

for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        clc

        disp(['Processing TRIAL-LEVEL TFR (Raw, FOOOF and Baselined) for Subject AOC ', num2str(subjects{subj})])
        disp('Loading EEG TFR data')
        load dataEEG_TFR_sternberg

        %% Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind4 = find(dataTFR.trialinfo(:, 1) == 24);
        ind6 = find(dataTFR.trialinfo(:, 1) == 26);

        %% Trial-level TFR
        cfg            = [];
        cfg.output     = 'pow';
        cfg.method     = 'mtmconvol';
        cfg.taper      = 'hanning';
        cfg.foi        = 3:1:30;
        cfg.t_ftimwin  = ones(length(cfg.foi), 1) .* 0.5;
        cfg.toi        = -1.5:0.05:3;
        cfg.keeptrials = 'yes';
        cfg.pad        = 'nextpow2';

        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR); % rpt_chan_freq_time
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

        %% Sliding-window FOOOF over single-trial spectra (mtmfft)
        winLen  = 0.5;   % 500 ms
        stepLen = 0.05;  % 50 ms

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
        cfg_fooof.keeptrials = 'no';  % each loop passes one trial window

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

            nTrials = numel(trlIdx);

            disp(' ')
            disp(['Running trial-level sliding-window FOOOF for condition ', condName])

            % One test window / trial to get sizes and master frequency grid
            cfg_sel0         = [];
            cfg_sel0.latency = [toi_start(1) (toi_start(1) + winLen)];
            cfg_sel0.trials  = trlIdx(1);
            datTFR_win0      = ft_selectdata(cfg_sel0, dataTFR);

            fooof_test  = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win0);
            nChan       = numel(fooof_test.label);
            freq_master = fooof_test.freq(:);
            nFreq       = numel(freq_master);

            fooof_powspctrm = nan(nTrials, nChan, nFreq, nTimePnts); % (model fit - aperiodic), log space
            fooof_powspec   = nan(nTrials, nChan, nFreq, nTimePnts); % input power spectrum, log space
            fooof_aperiodic = nan(nTrials, nChan, 4, nTimePnts);     % [offset exponent error r^2]

            for tr = 1:nTrials
                for timePnt = 1:nTimePnts

                    % Select one trial + one sliding window
                    cfg_sel         = [];
                    cfg_sel.latency = [toi_start(timePnt) (toi_start(timePnt) + winLen)];
                    cfg_sel.trials  = trlIdx(tr);
                    datTFR_win      = ft_selectdata(cfg_sel, dataTFR);

                    % Tracker
                    clc
                    disp(['Running TRIAL-LEVEL FOOOF for Subject ', num2str(subjects{subj})])
                    disp(['Subject:    ', num2str(subj), ' / ', num2str(length(subjects))])
                    disp(['Condition:  ', num2str(tfr_conds), ' / 3'])
                    disp(['Trial:      ', num2str(tr), ' / ', num2str(nTrials)])
                    disp(['Time Point: ', num2str(timePnt), ' / ', num2str(nTimePnts)])

                    % Run FOOOF on current trial window
                    evalc('fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win);');

                    if iscell(fooof_out.fooofparams)
                        repdata = fooof_out.fooofparams{1};
                    else
                        repdata = fooof_out.fooofparams;
                    end

                    freq_now = fooof_out.freq(:);
                    [tf, loc] = ismembertol(freq_now, freq_master, 1e-10);

                    local_ps     = nan(nChan, nFreq);
                    local_pkcomp = nan(nChan, nFreq);

                    local_err    = nan(nChan, 1);
                    local_rsq    = nan(nChan, 1);
                    local_offset = nan(nChan, 1);
                    local_expo   = nan(nChan, 1);

                    for ch = 1:nChan
                        % Input spectrum (FOOOF/log space)
                        ps_tmp = repdata(ch).power_spectrum(:);
                        local_ps(ch, loc(tf)) = ps_tmp(tf).';

                        % Aperiodic fit (FOOOF/log space)
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

                        if isfield(repdata(ch), 'fooofed_spectrum') && ~isempty(repdata(ch).fooofed_spectrum)
                            model_fit = repdata(ch).fooofed_spectrum(:);
                        else
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

                    fooof_powspctrm(tr, :, :, timePnt) = local_pkcomp;
                    fooof_powspec(tr,   :, :, timePnt) = local_ps;

                    fooof_aperiodic(tr, :, 1, timePnt) = local_offset;
                    fooof_aperiodic(tr, :, 2, timePnt) = local_expo;
                    fooof_aperiodic(tr, :, 3, timePnt) = local_err;
                    fooof_aperiodic(tr, :, 4, timePnt) = local_rsq;
                end
            end

            % Exclude bad FOOOF fits (r_squared < 0.90)
            r_sq_thresh = 0.90;
            bad_fits = squeeze(fooof_aperiodic(:, :, 4, :)) < r_sq_thresh; % rpt x chan x time
            n_bad = sum(bad_fits(:));
            if n_bad > 0
                fprintf('  Excluding %d / %d trial x channel x time FOOOF fits (r^2 < %.2f)\n', ...
                    n_bad, numel(bad_fits), r_sq_thresh);
                for tp = 1:nTimePnts
                    bad_tp = bad_fits(:, :, tp);
                    bad_mask = repmat(bad_tp, [1 1 nFreq]);

                    block_pk = squeeze(fooof_powspctrm(:, :, :, tp));
                    block_ps = squeeze(fooof_powspec(:, :, :, tp));

                    block_pk(bad_mask) = NaN;
                    block_ps(bad_mask) = NaN;

                    fooof_powspctrm(:, :, :, tp) = block_pk;
                    fooof_powspec(:, :, :, tp) = block_ps;
                end
            end

            % Construct condition-specific FieldTrip freq-like struct
            tfr_ff                = [];
            tfr_ff.label          = fooof_test.label;
            tfr_ff.freq           = fooof_test.freq;
            tfr_ff.time           = toi_centres;
            tfr_ff.powspctrm      = fooof_powspctrm;
            tfr_ff.power_spectrum = fooof_powspec;
            tfr_ff.fooofparams    = fooof_aperiodic;
            tfr_ff.trialinfo      = dataTFR.trialinfo(trlIdx, :);
            tfr_ff.dimord         = 'rpt_chan_freq_time';

            tfr_fooof{tfr_conds} = tfr_ff;
        end

        tfr2_fooof = tfr_fooof{1};
        tfr4_fooof = tfr_fooof{2};
        tfr6_fooof = tfr_fooof{3};

        disp(upper('FOOOF done on single-trial spectra (sliding windows)...'))

        %% Sanity check (single trial, 3-panel)
        time_point = 0.5;
        trial_idx = 1;
        [~, tim] = min(abs(tfr2_fooof.time - time_point));
        latWin = [toi_start(tim) (toi_start(tim) + winLen)];

        figure('Position', [0 0 1512 982]);

        cfg_sel         = [];
        cfg_sel.latency = latWin;
        cfg_sel.trials  = ind2(trial_idx);
        datTFR_win_sc   = ft_selectdata(cfg_sel, dataTFR);

        fooof_sc = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win_sc);

        if iscell(fooof_sc.fooofparams)
            repdata_sc = fooof_sc.fooofparams{1};
        else
            repdata_sc = fooof_sc.fooofparams;
        end

        freq_sc = fooof_sc.freq(:);
        nChan = numel(repdata_sc);

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

        subplot(3,1,1)
        plot(freq_sc, ps_in, 'k', 'LineWidth', 2); hold on
        plot(freq_sc, ap_mean, 'b--', 'LineWidth', 2)
        plot(freq_sc, model_mean, 'r', 'LineWidth', 2); hold off
        legend('Original', 'Aperiodic', 'Model fit', 'Location', 'best')
        title('Original Powspctrm')
        xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
        set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

        subplot(3,1,2)
        plot(freq_sc, ps_corr, 'LineWidth', 2)
        title('Original - Aperiodic')
        xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
        set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

        subplot(3,1,3)
        plot(freq_sc, pk_comp, 'LineWidth', 2)
        title('Model fit - Aperiodic (stored)')
        xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
        set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

        sgtitle(sprintf(['FOOOF sanity check TRIAL-level (AOC STERNBERG): Subject %s | Trial %d | ', ...
            'Window [%.2f %.2f] s | t = %.2f s'], ...
            subjects{subj}, trial_idx, latWin(1), latWin(2), tfr2_fooof.time(tim)), 'FontSize', 18)

        if ispc
            savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
        else
            savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
        end
        if ~exist(savePathControls, 'dir')
            mkdir(savePathControls)
        end

        saveName = sprintf('AOC_controls_FOOOF_stern_trials_subj%s.png', subjects{subj});
        saveas(gcf, fullfile(savePathControls, saveName));

        %% Baselined FOOOF TFR
        cfg              = [];
        cfg.baseline     = [-.5 -.25];
        cfg.baselinetype = 'absolute';
        tfr2_fooof_bl    = ft_freqbaseline(cfg, tfr2_fooof);
        tfr4_fooof_bl    = ft_freqbaseline(cfg, tfr4_fooof);
        tfr6_fooof_bl    = ft_freqbaseline(cfg, tfr6_fooof);

        disp(upper('FOOOF baseline done...'))

        %% Save trial-level TFR data
        cd(datapath)
        save tfr_stern_trials_fooof ...
            tfr2 tfr4 tfr6 ...
            tfr2_fooof tfr4_fooof tfr6_fooof ...
            tfr2_bl tfr4_bl tfr6_bl ...
            tfr2_fooof_bl tfr4_fooof_bl tfr6_fooof_bl

        %% Convert trial-level TFR data to trial-level POWSPCTRM (rpt x channels x frequency)
        analysisPeriodFull  = [0 2];
        analysisPeriodEarly = [0 1];
        analysisPeriodLate  = [1 2];
        freq_range          = [3 30];

        pow2_fooof          = triallevel_select_and_collapse(tfr2_fooof,    analysisPeriodFull,  freq_range);
        pow2_fooof_bl       = triallevel_select_and_collapse(tfr2_fooof_bl, analysisPeriodFull,  freq_range);
        pow2_fooof_bl_early = triallevel_select_and_collapse(tfr2_fooof_bl, analysisPeriodEarly, freq_range);
        pow2_fooof_bl_late  = triallevel_select_and_collapse(tfr2_fooof_bl, analysisPeriodLate,  freq_range);

        pow4_fooof          = triallevel_select_and_collapse(tfr4_fooof,    analysisPeriodFull,  freq_range);
        pow4_fooof_bl       = triallevel_select_and_collapse(tfr4_fooof_bl, analysisPeriodFull,  freq_range);
        pow4_fooof_bl_early = triallevel_select_and_collapse(tfr4_fooof_bl, analysisPeriodEarly, freq_range);
        pow4_fooof_bl_late  = triallevel_select_and_collapse(tfr4_fooof_bl, analysisPeriodLate,  freq_range);

        pow6_fooof          = triallevel_select_and_collapse(tfr6_fooof,    analysisPeriodFull,  freq_range);
        pow6_fooof_bl       = triallevel_select_and_collapse(tfr6_fooof_bl, analysisPeriodFull,  freq_range);
        pow6_fooof_bl_early = triallevel_select_and_collapse(tfr6_fooof_bl, analysisPeriodEarly, freq_range);
        pow6_fooof_bl_late  = triallevel_select_and_collapse(tfr6_fooof_bl, analysisPeriodLate,  freq_range);

        save power_stern_fooof_trials ...
            pow2_fooof pow4_fooof pow6_fooof ...
            pow2_fooof_bl pow4_fooof_bl pow6_fooof_bl ...
            pow2_fooof_bl_early pow4_fooof_bl_early pow6_fooof_bl_early ...
            pow2_fooof_bl_late pow4_fooof_bl_late pow6_fooof_bl_late

        clc
        fprintf(['Subject AOC %s (%.3d/%.3d) DONE (trial-level sliding-window FOOOF: ', ...
            'model - aperiodic)\n'], num2str(subjects{subj}), subj, length(subjects))
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end

function out = triallevel_select_and_collapse(in, latency_window, freq_window)
cfg = [];
cfg.latency = latency_window;
cfg.frequency = freq_window;
out = ft_selectdata(cfg, in);

if ndims(out.powspctrm) == 4
    out.powspctrm = mean(out.powspctrm, 4, 'omitnan');
    if isfield(out, 'time')
        out = rmfield(out, 'time');
    end
    out.dimord = 'rpt_chan_freq';
end
end
