%% FOOOF Feature Extraction — Resting State (AOC)
% Sliding-window FOOOF on continuous resting-state EEG.
% Extracts oscillatory power (model fit - aperiodic) for the full recording.
%
% FOOOF output: (model fit - aperiodic) in FOOOF/log space
%   Positive values = oscillatory peak above aperiodic background
%   Zero / negative  = no oscillatory peak at that frequency
%
% Input:
%   dataEEG_resting.mat (from AOC_preprocessing_resting.m)
%
% Output per subject:
%   fooof_resting.mat — FieldTrip-like struct (chan x freq x time)
%     .powspctrm      = oscillatory component (model fit - aperiodic)
%     .power_spectrum  = input spectrum (FOOOF/log space)
%     .fooofparams     = [offset exponent error r²] per chan per window
%     .time            = window centre times (s)
%     .freq            = frequency axis (Hz)
%     .label           = channel labels
%     .dimord          = 'chan_freq_time'
%   segInfo            = fixcross/blank sample ranges (passed through)

clear; close all; clc

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

% Setup logging
logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
scriptName = 'AOC_eeg_fex_resting_FOOOF';

% FOOOF sliding-window parameters (match interaction script)
winLen  = 2;       % 2 s windows
stepLen = 1;       % 1 s step

% FOOOF config
cfg_fooof            = [];
cfg_fooof.method     = 'mtmfft';
cfg_fooof.taper      = 'hanning';
cfg_fooof.foilim     = [3 30];    % match Sternberg/N-back FOOOF
cfg_fooof.pad        = 5;         % pad to 5 s → 0.2 Hz resolution
cfg_fooof.output     = 'fooof';
cfg_fooof.keeptrials = 'no';

% Quality threshold
r_sq_thresh = 0.90;

%% Loop subjects
for subj = 1:length(subjects)
    try
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        clc

        disp(['FOOOF Resting State — Subject AOC ', subjects{subj}, ...
              ' (', num2str(subj), '/', num2str(length(subjects)), ')'])

        % Check if already done
        if exist('fooof_resting.mat', 'file')
            disp(['Subject ', subjects{subj}, ' already done. SKIPPING...'])
            continue
        end

        % Load preprocessed resting EEG
        if ~exist('dataEEG_resting.mat', 'file')
            fprintf('  No resting data for subject %s. SKIPPING...\n', subjects{subj})
            continue
        end
        load('dataEEG_resting.mat', 'dataEEG', 'segInfo')

        %% Define time axis for sliding windows
        tMin = dataEEG.time{1}(1);
        tMax = dataEEG.time{1}(end);

        toi_start   = tMin:stepLen:(tMax - winLen);
        toi_centres = toi_start + winLen / 2;
        nTimePnts   = numel(toi_centres);

        fprintf('  Recording: %.1f s, %d FOOOF windows\n', tMax - tMin, nTimePnts)

        %% Test window to get frequency grid and channel count
        cfg_sel0         = [];
        cfg_sel0.latency = [toi_start(1), toi_start(1) + winLen];
        datWin0          = ft_selectdata(cfg_sel0, dataEEG);

        fooof_test  = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datWin0);

        nChan       = numel(fooof_test.label);
        freq_master = fooof_test.freq(:);
        nFreq       = numel(freq_master);

        %% Pre-allocate
        fooof_powspctrm = nan(nChan, nFreq, nTimePnts);   % oscillatory (model - aperiodic)
        fooof_powspec   = nan(nChan, nFreq, nTimePnts);   % input spectrum
        fooof_aperiodic = nan(nChan, 4,    nTimePnts);    % [offset exponent error r²]

        %% Sliding-window FOOOF
        for tp = 1:nTimePnts

            % Select window
            cfg_sel         = [];
            cfg_sel.latency = [toi_start(tp), toi_start(tp) + winLen];
            datWin          = ft_selectdata(cfg_sel, dataEEG);

            % Progress tracker
            if mod(tp, 25) == 0 || tp == 1
                clc
                fprintf('FOOOF Resting — Subject %s (%d/%d)\n', subjects{subj}, subj, length(subjects))
                fprintf('Window %d / %d  (%.0f%%)\n', tp, nTimePnts, 100 * tp / nTimePnts)
            end

            % Run FOOOF (suppress console output)
            try
                out = evalc('fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datWin);'); %#ok<ASGLU>
            catch
                continue;  % skip window on FOOOF error
            end

            if iscell(fooof_out.fooofparams)
                repdata = fooof_out.fooofparams{1};
            else
                repdata = fooof_out.fooofparams;
            end

            freq_now = fooof_out.freq(:);

            % Map onto master frequency grid
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

                % Aperiodic parameters
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

                % Full model fit (aperiodic + peaks) in FOOOF/log space
                if isfield(repdata(ch), 'fooofed_spectrum') && ~isempty(repdata(ch).fooofed_spectrum)
                    model_fit = repdata(ch).fooofed_spectrum(:);
                else
                    % Fallback: reconstruct via Gaussian peaks + aperiodic
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

                % Oscillatory component = model fit - aperiodic
                pkcomp = model_fit(:) - ap_fit(:);
                local_pkcomp(ch, loc(tf)) = pkcomp(tf).';
            end

            % Store
            fooof_powspctrm(:, :, tp) = local_pkcomp;
            fooof_powspec(:,   :, tp) = local_ps;

            fooof_aperiodic(:, 1, tp) = local_offset;
            fooof_aperiodic(:, 2, tp) = local_expo;
            fooof_aperiodic(:, 3, tp) = local_err;
            fooof_aperiodic(:, 4, tp) = local_rsq;
        end

        %% Exclude bad FOOOF fits (r_squared < threshold)
        bad_fits = squeeze(fooof_aperiodic(:, 4, :)) < r_sq_thresh;
        n_bad = sum(bad_fits(:));
        if n_bad > 0
            fprintf('  Excluding %d / %d channel x time FOOOF fits (r^2 < %.2f)\n', ...
                n_bad, numel(bad_fits), r_sq_thresh);
            for tp = 1:nTimePnts
                fooof_powspctrm(bad_fits(:, tp), :, tp) = NaN;
                fooof_powspec(bad_fits(:, tp), :, tp)   = NaN;
            end
        end

        %% Build output struct
        fooof_resting                = [];
        fooof_resting.label          = fooof_test.label;
        fooof_resting.freq           = fooof_test.freq;
        fooof_resting.time           = toi_centres;
        fooof_resting.powspctrm      = fooof_powspctrm;     % oscillatory (model - aperiodic)
        fooof_resting.power_spectrum = fooof_powspec;        % input spectrum
        fooof_resting.fooofparams    = fooof_aperiodic;      % chan x 4 x time
        fooof_resting.dimord         = 'chan_freq_time';

        %% Sanity check figure (mid-recording window)
        midTP = round(nTimePnts / 2);
        latWin = [toi_start(midTP), toi_start(midTP) + winLen];

        figure('Position', [0 0 1500 420], 'Color', 'w');

        cfg_sel_sc         = [];
        cfg_sel_sc.latency = latWin;
        datWin_sc          = ft_selectdata(cfg_sel_sc, dataEEG);

        try
            fooof_sc = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datWin_sc);

            if iscell(fooof_sc.fooofparams)
                repdata_sc = fooof_sc.fooofparams{1};
            else
                repdata_sc = fooof_sc.fooofparams;
            end

            freq_sc    = fooof_sc.freq(:);
            nChan_sc   = numel(repdata_sc);

            ps_all     = nan(numel(freq_sc), nChan_sc);
            ap_all     = nan(numel(freq_sc), nChan_sc);
            model_all  = nan(numel(freq_sc), nChan_sc);
            pkcomp_all = nan(numel(freq_sc), nChan_sc);

            for ch = 1:nChan_sc
                ps_all(:, ch) = repdata_sc(ch).power_spectrum(:);

                ap = repdata_sc(ch).aperiodic_params(:);
                if numel(ap) == 2
                    ap_fit = ap(1) - ap(2) .* log10(freq_sc);
                else
                    ap_fit = ap(1) - log10(ap(2) + freq_sc.^ap(3));
                end
                ap_all(:, ch) = ap_fit(:);

                if isfield(repdata_sc(ch), 'fooofed_spectrum') && ~isempty(repdata_sc(ch).fooofed_spectrum)
                    mf = repdata_sc(ch).fooofed_spectrum(:);
                else
                    pk = repdata_sc(ch).peak_params;
                    gs = zeros(numel(freq_sc), 1);
                    if ~isempty(pk)
                        for p = 1:size(pk, 1)
                            gs = gs + pk(p,2) .* exp(-(freq_sc - pk(p,1)).^2 ./ (2*pk(p,3).^2));
                        end
                    end
                    mf = ap_fit(:) + gs(:);
                end
                model_all(:, ch)  = mf(:);
                pkcomp_all(:, ch) = mf(:) - ap_fit(:);
            end

            subplot(1,3,1)
            plot(freq_sc, mean(ps_all,2,'omitnan'), 'k', 'LineWidth', 2); hold on
            plot(freq_sc, mean(ap_all,2,'omitnan'), 'b--', 'LineWidth', 2)
            plot(freq_sc, mean(model_all,2,'omitnan'), 'r', 'LineWidth', 2); hold off
            legend('Original', 'Aperiodic', 'Model fit', 'Location', 'best')
            title('Original Powspctrm')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

            subplot(1,3,2)
            plot(freq_sc, mean(ps_all,2,'omitnan') - mean(ap_all,2,'omitnan'), 'LineWidth', 2)
            title('Original - Aperiodic')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

            subplot(1,3,3)
            plot(freq_sc, mean(pkcomp_all,2,'omitnan'), 'LineWidth', 2)
            title('Model fit - Aperiodic (stored)')
            xlabel('Frequency (Hz)'); ylabel('Power (FOOOF/log space)')
            set(gca, 'FontSize', 13); xlim([min(freq_sc) max(freq_sc)])

            sgtitle(sprintf('FOOOF sanity check (AOC RESTING): Subject %s | Window [%.1f %.1f] s', ...
                subjects{subj}, latWin(1), latWin(2)), 'FontSize', 18)
        catch
            text(0.5, 0.5, 'Sanity check failed', 'HorizontalAlignment', 'center');
        end

        if ispc
            savePathControls = 'W:\Students\Arne\AOC\data\controls\FOOOF\';
        else
            savePathControls = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/FOOOF/';
        end
        if ~exist(savePathControls, 'dir'), mkdir(savePathControls); end
        saveName = sprintf('AOC_controls_FOOOF_resting_subj%s.png', subjects{subj});
        saveas(gcf, fullfile(savePathControls, saveName));

        %% Save
        cd(datapath)
        save('fooof_resting', 'fooof_resting', 'segInfo', '-v7.3')

        clear dataEEG fooof_resting
        clc
        fprintf('Subject AOC %s (%d/%d) DONE (sliding-window FOOOF resting state)\n', ...
            subjects{subj}, subj, length(subjects))

    catch ME
        log_error(scriptName, subjects{subj}, subj, length(subjects), ME, logDir);
        fprintf('Continuing to next subject...\n');
    end
end
disp('RESTING STATE FOOOF FEATURE EXTRACTION COMPLETE.')
