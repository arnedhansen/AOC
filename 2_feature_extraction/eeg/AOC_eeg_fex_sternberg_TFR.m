%% AOC EEG FOOOF Extraction — Sternberg (TFR Branch)
% Computes TFR-based FOOOF outputs and builds FOOOF-only EEG products:
% `eeg_data_sternberg_FOOOF` (MAT + CSV).
% IAF is loaded from non-FOOOF CSV (`AOC_eeg_matrix_sternberg.csv`).
% If IAF is missing/invalid for a subject-condition, FOOOF alpha is computed
% with fallback band [8 14] Hz (aligned with normal non-FOOOF scripts).
clear; close all; clc

startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;
featPath = paths.features;

% Setup logging
logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
scriptName = 'AOC_eeg_fex_sternberg_FOOOF';
eeg_data_sternberg_FOOOF = struct('ID', {}, 'Condition', {}, ...
    'AlphaPower_FOOOF', {}, 'AlphaPower_FOOOF_bl', {}, ...
    'AlphaPower_FOOOF_bl_early', {}, 'AlphaPower_FOOOF_bl_late', {});

iaf_csv = fullfile(featPath, 'AOC_eeg_matrix_sternberg.csv');
if ~isfile(iaf_csv)
    error('Missing non-FOOOF EEG CSV for IAF lookup: %s', iaf_csv);
end
iaf_table = readtable(iaf_csv);

for subj = 1:length(subjects)
    try
        datapath = fullfile(path, subjects{subj}, 'eeg');
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
        cfg.foi        = 3:1:30;
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
        cfg_fooof.foilim     = [3 30];   % match methods & N-back for comparable alpha
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

        tfr2_fooof = tfr_fooof{1};
        tfr4_fooof = tfr_fooof{2};
        tfr6_fooof = tfr_fooof{3};

        disp(upper('FOOOF done on trial-averaged spectra (sliding windows)...'))

        %% Sanity check (3-panel): input, input-aperiodic, model-aperiodic (stored)
        time_point = 0.5;
        [~, tim]   = min(abs(tfr2_fooof.time - time_point));

        latWin = [toi_start(tim) (toi_start(tim) + winLen)];

        figure('Position', [0 0 1500 420], 'Color', 'w');

        cfg_sel         = [];
        cfg_sel.latency = latWin;
        cfg_sel.trials  = ind2;
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

        sgtitle(sprintf('FOOOF sanity check (AOC STERNBERG): Subject %s | Window [%.2f %.2f] s | t = %.2f s', ...
            subjects{subj}, latWin(1), latWin(2), tfr2_fooof.time(tim)), 'FontSize', 18)

        if ispc
            savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
        else
            savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
        end
        if ~exist(savePathControls, 'dir')
            mkdir(savePathControls)
        end

        saveName = sprintf('AOC_controls_FOOOF_stern_subj%s.png', subjects{subj});
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
        freq_range          = [3 30];   % match FOOOF foilim

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

        save power_stern_fooof_TFR ...
            pow2_fooof pow4_fooof pow6_fooof ...
            pow2_fooof_bl pow4_fooof_bl pow6_fooof_bl ...
            pow2_fooof_bl_early pow4_fooof_bl_early pow6_fooof_bl_early ...
            pow2_fooof_bl_late pow4_fooof_bl_late pow6_fooof_bl_late

        % Build FOOOF-only subject EEG rows using IAF from non-FOOOF CSV.
        occ_channels_fooof = occ_channels_from_labels(pow2_fooof.label);
        occ_idx = find(ismember(pow2_fooof.label, occ_channels_fooof));
        condVals = [2; 4; 6];
        subID = str2double(subjects{subj});
        if isnan(subID)
            subID = str2num(subjects{subj}); %#ok<ST2NM>
        end
        AlphaPower_FOOOF = nan(3, 1);
        AlphaPower_FOOOF_bl = nan(3, 1);
        AlphaPower_FOOOF_bl_early = nan(3, 1);
        AlphaPower_FOOOF_bl_late = nan(3, 1);
        condPows = {pow2_fooof, pow4_fooof, pow6_fooof};
        condPows_bl = {pow2_fooof_bl, pow4_fooof_bl, pow6_fooof_bl};
        condPows_bl_early = {pow2_fooof_bl_early, pow4_fooof_bl_early, pow6_fooof_bl_early};
        condPows_bl_late = {pow2_fooof_bl_late, pow4_fooof_bl_late, pow6_fooof_bl_late};
        nIAF_fallback = 0;
        for c = 1:3
            IAF_now = get_iaf_from_table(iaf_table, subID, condVals(c));
            band = [IAF_now - 4, IAF_now + 2];
            if ~isfinite(IAF_now) || any(~isfinite(band)) || band(1) >= band(2)
                band = [8 14];
                nIAF_fallback = nIAF_fallback + 1;
            end
            AlphaPower_FOOOF(c) = robust_roi_pow(condPows{c}, occ_idx, band);
            AlphaPower_FOOOF_bl(c) = robust_roi_pow(condPows_bl{c}, occ_idx, band);
            AlphaPower_FOOOF_bl_early(c) = robust_roi_pow(condPows_bl_early{c}, occ_idx, band);
            AlphaPower_FOOOF_bl_late(c) = robust_roi_pow(condPows_bl_late{c}, occ_idx, band);
        end
        if nIAF_fallback > 0
            fprintf('Subject %s: IAF fallback [8 14] used in %d/3 conditions.\n', subjects{subj}, nIAF_fallback);
        end
        subj_data_fooof = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell(condVals), ...
            'AlphaPower_FOOOF', num2cell(AlphaPower_FOOOF), ...
            'AlphaPower_FOOOF_bl', num2cell(AlphaPower_FOOOF_bl), ...
            'AlphaPower_FOOOF_bl_early', num2cell(AlphaPower_FOOOF_bl_early), ...
            'AlphaPower_FOOOF_bl_late', num2cell(AlphaPower_FOOOF_bl_late));
        save eeg_matrix_sternberg_FOOOF_subj subj_data_fooof
        eeg_data_sternberg_FOOOF = [eeg_data_sternberg_FOOOF; subj_data_fooof];

        %% Additional sanity check: WM-load differences (subject-level alpha, FOOOF-bl)
        % Uses split-style cloud + box + jitter aesthetics.
        occ_channels = occ_channels_from_labels(pow2_fooof_bl.label);
        alphaRange = [8 14];
        a2 = extract_subject_alpha_samples(pow2_fooof_bl, occ_channels, alphaRange);
        a4 = extract_subject_alpha_samples(pow4_fooof_bl, occ_channels, alphaRange);
        a6 = extract_subject_alpha_samples(pow6_fooof_bl, occ_channels, alphaRange);

        figure('Position', [0 0 1512 982]);
        hold on
        cols = [0.60 0.78 0.88; 0.63 0.82 0.61; 0.93 0.70 0.78];
        draw_one_cloud_sanity(a2, 1, cols(1, :), 0.30, 120, 0.75);
        draw_one_cloud_sanity(a4, 2, cols(2, :), 0.30, 120, 0.75);
        draw_one_cloud_sanity(a6, 3, cols(3, :), 0.30, 120, 0.75);
        yline(0, '--', 'Color', [0.6 0.6 0.6]);
        set(gca, 'XTick', 1:3, 'XTickLabel', {'WM load 2', 'WM load 4', 'WM load 6'}, 'FontSize', 20);
        ylabel('AlphaPower\_FOOOF\_bl');
        title(sprintf('Sanity check WM-load differences: Subject %s', subjects{subj}), 'Interpreter', 'none');
        box off

        saveName = sprintf('AOC_controls_FOOOF_stern_subj%s_alphaLoadSanity.png', subjects{subj});
        saveas(gcf, fullfile(savePathControls, saveName));
        close(gcf);

        clc
        fprintf('Subject AOC %s (%.3d/%.3d) DONE (sliding-window FOOOF: model - aperiodic) \n', ...
            num2str(subjects{subj}), subj, length(subjects))
    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end

if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_sternberg_FOOOF eeg_data_sternberg_FOOOF
    writetable(struct2table(eeg_data_sternberg_FOOOF), 'W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_sternberg_FOOOF.csv')
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_sternberg_FOOOF eeg_data_sternberg_FOOOF
    writetable(struct2table(eeg_data_sternberg_FOOOF), '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_sternberg_FOOOF.csv')
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end+1} = lab; %#ok<AGROW>
    end
end
if isempty(ch)
    ch = labels;
end
end

function alpha_vals = extract_subject_alpha_samples(S, channels, band)
ch_idx = find(ismember(S.label, channels));
f_idx = S.freq >= band(1) & S.freq <= band(2);
if isempty(ch_idx) || ~any(f_idx)
    alpha_vals = NaN;
    return
end
X = S.powspctrm(ch_idx, f_idx);
alpha_vals = X(:);
alpha_vals = alpha_vals(isfinite(alpha_vals));
if isempty(alpha_vals), alpha_vals = NaN; end
end

function iaf_val = get_iaf_from_table(T, id_val, cond_val)
row = T.ID == id_val & T.Condition == cond_val;
if ~any(row) || ~ismember('IAF', T.Properties.VariableNames)
    iaf_val = NaN;
    return
end
vals = T.IAF(row);
iaf_val = vals(1);
end

function v = robust_roi_pow(S, channelIdx, band)
fmask = S.freq >= band(1) & S.freq <= band(2);
if ~any(fmask)
    v = NaN;
    return
end
x = S.powspctrm(channelIdx, fmask);
x = x(:);
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = mean(x, 'omitnan');
end
end

function draw_one_cloud_sanity(yvals, xpos, col, box_w, dot_size, dot_alpha)
y = yvals(isfinite(yvals));
if isempty(y)
    return
end
[f, xi] = ksdensity(y, 'NumPoints', 120);
if max(f) > 0
    f = f / max(f) * 0.30;
else
    f = zeros(size(f));
end
x_den = xpos - 0.08;
x_box = xpos + 0.03;
fill([x_den - f, fliplr(repmat(x_den, 1, numel(f)))], [xi, fliplr(xi)], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1);
q1 = prctile(y, 25);
q3 = prctile(y, 75);
med = median(y);
p5 = prctile(y, 5);
p95 = prctile(y, 95);
plot([x_box x_box], [p5 q1], '-k', 'LineWidth', 1.2);
plot([x_box x_box], [q3 p95], '-k', 'LineWidth', 1.2);
rectangle('Position', [x_box-box_w/2, q1, box_w, q3-q1], ...
    'FaceColor', [col 0.08], 'EdgeColor', 'k', 'LineWidth', 1.2);
plot(x_box + [-box_w/2 box_w/2], [med med], '-k', 'LineWidth', 2);
jit = box_w * (rand(numel(y), 1) - 0.5);
scatter(x_box + jit, y, dot_size, col, 'filled', ...
    'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
end
