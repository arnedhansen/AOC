%% AOC EEG FOOOF Feature Extraction - Sternberg
% Computes subject-level FOOOF alpha features and IAF_specParam for sternberg.
% Saves `AOC_eeg_matrix_sternberg_FOOOF` (MAT + CSV) and per-subject `IAF_specParam_sternberg.mat`.
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;
featPath = paths.features;

if exist('ft_freqanalysis_Arne_FOOOF', 'file') ~= 2
    error(['ft_freqanalysis_Arne_FOOOF not found on MATLAB path. ', ...
        'Add it first, for example from /Users/Arne/Documents/GitHub/functions.']);
end

iaf_csv = fullfile(featPath, 'AOC_eeg_matrix_sternberg.csv');
if ~isfile(iaf_csv)
    error('Missing non-FOOOF EEG CSV for IAF lookup: %s', iaf_csv);
end
iaf_table = readtable(iaf_csv);

alphaRange = [8 14];
winIAF = [1 2];
eeg_data_sternberg_FOOOF = struct('ID', {}, 'Condition', {}, ...
    'AlphaPower_FOOOF', {}, 'AlphaPower_FOOOF_bl', {}, ...
    'AlphaPower_FOOOF_bl_early', {}, 'AlphaPower_FOOOF_bl_late', {}, ...
    'IAF_specParam', {});

for subj = 1:length(subjects)
    try
        clc
        fprintf('[EEG FOOOF - STERNBERG] Subject %d / %d\n', subj, length(subjects))
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        load dataEEG_TFR_sternberg

        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind4 = find(dataTFR.trialinfo(:, 1) == 24);
        ind6 = find(dataTFR.trialinfo(:, 1) == 26);
        condTrials = {ind2, ind4, ind6};
        condVals = [2; 4; 6];
        condNames = {'set2', 'set4', 'set6'};

        occ_channels = occ_channels_from_labels(dataTFR.label);

        winLen = 0.5;
        stepLen = 0.05;
        tMin = dataTFR.time{1}(1);
        tMax = dataTFR.time{1}(end);
        toi_start = tMin:stepLen:(tMax - winLen);
        toi_centres = toi_start + winLen / 2;
        nTimePnts = numel(toi_centres);

        cfg_fooof = [];
        cfg_fooof.method = 'mtmfft';
        cfg_fooof.taper = 'hanning';
        cfg_fooof.foilim = [2 40];
        cfg_fooof.pad = 5;
        cfg_fooof.output = 'fooof';
        cfg_fooof.keeptrials = 'no';

        tfr_fooof = cell(1, 3);
        tfr_fooof_bl = cell(1, 3);
        for c = 1:3
            trlIdx = condTrials{c};
            if isempty(trlIdx)
                continue
            end
            clc
            fprintf('[EEG FOOOF - STERNBERG] Subject %d / %d, condition %s\n', subj, length(subjects), condNames{c})

            cfg_sel0 = [];
            cfg_sel0.latency = [toi_start(1) (toi_start(1) + winLen)];
            cfg_sel0.trials = trlIdx;
            datTFR_win0 = ft_selectdata(cfg_sel0, dataTFR);
            fooof_test = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win0);
            nChan = numel(fooof_test.label);
            freq_master = fooof_test.freq(:);
            nFreq = numel(freq_master);

            fooof_powspctrm = nan(nChan, nFreq, nTimePnts);
            fooof_aperiodic = nan(nChan, 4, nTimePnts);

            for timePnt = 1:nTimePnts
                cfg_sel = [];
                cfg_sel.latency = [toi_start(timePnt) (toi_start(timePnt) + winLen)];
                cfg_sel.trials = trlIdx;
                datTFR_win = ft_selectdata(cfg_sel, dataTFR);
                fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, datTFR_win);

                if iscell(fooof_out.fooofparams)
                    repdata = fooof_out.fooofparams{1};
                else
                    repdata = fooof_out.fooofparams;
                end

                freq_now = fooof_out.freq(:);
                [tf, loc] = ismembertol(freq_now, freq_master, 1e-10);
                local_pkcomp = nan(nChan, nFreq);
                local_err = nan(nChan, 1);
                local_rsq = nan(nChan, 1);
                local_offset = nan(nChan, 1);
                local_expo = nan(nChan, 1);

                for ch = 1:nChan
                    ap = repdata(ch).aperiodic_params(:);
                    local_err(ch) = repdata(ch).error;
                    local_rsq(ch) = repdata(ch).r_squared;
                    if numel(ap) == 2
                        offset = ap(1);
                        expo = ap(2);
                        ap_fit = offset - expo .* log10(freq_now);
                    else
                        offset = ap(1);
                        knee = ap(2);
                        expo = ap(3);
                        ap_fit = offset - log10(knee + freq_now .^ expo);
                    end
                    local_offset(ch) = offset;
                    local_expo(ch) = expo;

                    if isfield(repdata(ch), 'fooofed_spectrum') && ~isempty(repdata(ch).fooofed_spectrum)
                        model_fit = repdata(ch).fooofed_spectrum(:);
                    else
                        pk = repdata(ch).peak_params;
                        gauss_sum = zeros(numel(freq_now), 1);
                        if ~isempty(pk)
                            for p = 1:size(pk, 1)
                                cf = pk(p, 1);
                                amp = pk(p, 2);
                                bw = pk(p, 3);
                                gauss_sum = gauss_sum + amp .* exp(-(freq_now - cf) .^ 2 ./ (2 * bw .^ 2));
                            end
                        end
                        model_fit = ap_fit(:) + gauss_sum(:);
                    end
                    pkcomp = model_fit(:) - ap_fit(:);
                    local_pkcomp(ch, loc(tf)) = pkcomp(tf).';
                end

                fooof_powspctrm(:, :, timePnt) = local_pkcomp;
                fooof_aperiodic(:, 1, timePnt) = local_offset;
                fooof_aperiodic(:, 2, timePnt) = local_expo;
                fooof_aperiodic(:, 3, timePnt) = local_err;
                fooof_aperiodic(:, 4, timePnt) = local_rsq;
            end

            bad_fits = squeeze(fooof_aperiodic(:, 4, :)) < 0.90;
            for tp = 1:nTimePnts
                fooof_powspctrm(bad_fits(:, tp), :, tp) = NaN;
            end

            tfr_ff = [];
            tfr_ff.label = fooof_test.label;
            tfr_ff.freq = fooof_test.freq;
            tfr_ff.time = toi_centres;
            tfr_ff.powspctrm = fooof_powspctrm;
            tfr_ff.fooofparams = fooof_aperiodic;
            tfr_ff.dimord = 'chan_freq_time';
            tfr_fooof{c} = tfr_ff;

            cfgb = [];
            cfgb.baseline = [-1.5 -0.5];
            cfgb.baselinetype = 'absolute';
            tfr_fooof_bl{c} = ft_freqbaseline(cfgb, tfr_ff);
        end

        if any(cellfun(@isempty, tfr_fooof))
            error('Missing FOOOF TFR output for subject %s', subjects{subj});
        end

        window.full = [0 2];
        window.early = [0 1];
        window.late = [1 2];
        freq_range = [2 40];

        pow2_fooof_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof{1}));
        pow4_fooof_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof{2}));
        pow6_fooof_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof{3}));
        pow2_fooof_bl_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof_bl{1}));
        pow4_fooof_bl_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof_bl{2}));
        pow6_fooof_bl_full = remove_time_dimension(select_data(window.full, freq_range, tfr_fooof_bl{3}));
        pow2_fooof_bl_early = remove_time_dimension(select_data(window.early, freq_range, tfr_fooof_bl{1}));
        pow4_fooof_bl_early = remove_time_dimension(select_data(window.early, freq_range, tfr_fooof_bl{2}));
        pow6_fooof_bl_early = remove_time_dimension(select_data(window.early, freq_range, tfr_fooof_bl{3}));
        pow2_fooof_bl_late = remove_time_dimension(select_data(window.late, freq_range, tfr_fooof_bl{1}));
        pow4_fooof_bl_late = remove_time_dimension(select_data(window.late, freq_range, tfr_fooof_bl{2}));
        pow6_fooof_bl_late = remove_time_dimension(select_data(window.late, freq_range, tfr_fooof_bl{3}));

        occ_idx = find(ismember(pow2_fooof_full.label, occ_channels));
        if isempty(occ_idx)
            occ_idx = 1:numel(pow2_fooof_full.label);
        end

        IAF_specParam2 = iaf_specparam_fooof(dataTFR, ind2, winIAF, occ_channels, alphaRange);
        IAF_specParam4 = iaf_specparam_fooof(dataTFR, ind4, winIAF, occ_channels, alphaRange);
        IAF_specParam6 = iaf_specparam_fooof(dataTFR, ind6, winIAF, occ_channels, alphaRange);
        IAF_specParam = [IAF_specParam2; IAF_specParam4; IAF_specParam6];

        subID = str2double(subjects{subj});
        if isnan(subID)
            subID = str2num(subjects{subj}); %#ok<ST2NM>
        end

        AlphaPower_FOOOF = nan(3, 1);
        AlphaPower_FOOOF_bl = nan(3, 1);
        AlphaPower_FOOOF_bl_early = nan(3, 1);
        AlphaPower_FOOOF_bl_late = nan(3, 1);
        condPows_full = {pow2_fooof_full, pow4_fooof_full, pow6_fooof_full};
        condPows_bl_full = {pow2_fooof_bl_full, pow4_fooof_bl_full, pow6_fooof_bl_full};
        condPows_bl_early = {pow2_fooof_bl_early, pow4_fooof_bl_early, pow6_fooof_bl_early};
        condPows_bl_late = {pow2_fooof_bl_late, pow4_fooof_bl_late, pow6_fooof_bl_late};
        for c = 1:3
            IAF_now = get_iaf_from_table(iaf_table, subID, condVals(c));
            band = [IAF_now - 4, IAF_now + 2];
            if ~isfinite(IAF_now) || any(~isfinite(band)) || band(1) >= band(2)
                band = [8 14];
            end
            AlphaPower_FOOOF(c) = robust_roi_pow(condPows_full{c}, occ_idx, band);
            AlphaPower_FOOOF_bl(c) = robust_roi_pow(condPows_bl_full{c}, occ_idx, band);
            AlphaPower_FOOOF_bl_early(c) = robust_roi_pow(condPows_bl_early{c}, occ_idx, band);
            AlphaPower_FOOOF_bl_late(c) = robust_roi_pow(condPows_bl_late{c}, occ_idx, band);
        end

        subj_data_fooof = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell(condVals), ...
            'AlphaPower_FOOOF', num2cell(AlphaPower_FOOOF), ...
            'AlphaPower_FOOOF_bl', num2cell(AlphaPower_FOOOF_bl), ...
            'AlphaPower_FOOOF_bl_early', num2cell(AlphaPower_FOOOF_bl_early), ...
            'AlphaPower_FOOOF_bl_late', num2cell(AlphaPower_FOOOF_bl_late), ...
            'IAF_specParam', num2cell(IAF_specParam));
        fooof_specparam = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell(condVals), ...
            'IAF_specParam', num2cell(IAF_specParam));

        save(fullfile(datapath, 'IAF_specParam_sternberg.mat'), 'fooof_specparam')
        save(fullfile(datapath, 'power_sternberg_fooof_TFR.mat'), ...
            'pow2_fooof_full', 'pow4_fooof_full', 'pow6_fooof_full', ...
            'pow2_fooof_bl_full', 'pow4_fooof_bl_full', 'pow6_fooof_bl_full', ...
            'pow2_fooof_bl_early', 'pow4_fooof_bl_early', 'pow6_fooof_bl_early', ...
            'pow2_fooof_bl_late', 'pow4_fooof_bl_late', 'pow6_fooof_bl_late')
        eeg_data_sternberg_FOOOF = [eeg_data_sternberg_FOOOF; subj_data_fooof];
    catch ME
        fprintf('STERNBERG FOOOF failed for subject %s: %s\n', subjects{subj}, ME.message);
    end
end

save_or_merge_sternberg_fooof_matrix(featPath, eeg_data_sternberg_FOOOF);

function cf = iaf_specparam_fooof(dataTFR, trialinds, winSec, chLabs, alphaRange)
cf = NaN;
if isempty(trialinds) || isempty(chLabs)
    return
end
try
    cfg_fooof = [];
    cfg_fooof.method = 'mtmfft';
    cfg_fooof.taper = 'hanning';
    cfg_fooof.foilim = [2 40];
    cfg_fooof.pad = 5;
    cfg_fooof.output = 'fooof';
    cfg_fooof.keeptrials = 'no';

    cfg_sel = [];
    cfg_sel.latency = winSec;
    cfg_sel.trials = trialinds(:)';
    cfg_sel.channel = chLabs(:);
    dat = ft_selectdata(cfg_sel, dataTFR);
    fooof_out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat);

    if iscell(fooof_out.fooofparams)
        rep = fooof_out.fooofparams{1};
    else
        rep = fooof_out.fooofparams;
    end

    cfs = [];
    for k = 1:numel(rep)
        if ~isfield(rep(k), 'peak_params') || isempty(rep(k).peak_params)
            continue
        end
        pk = rep(k).peak_params;
        inb = pk(:, 1) >= alphaRange(1) & pk(:, 1) <= alphaRange(2);
        if ~any(inb)
            continue
        end
        pk = pk(inb, :);
        [~, im] = max(pk(:, 2));
        cfs(end + 1) = pk(im, 1); %#ok<AGROW>
    end
    if ~isempty(cfs)
        cf = median(cfs, 'omitnan');
    end
catch
    cf = NaN;
end
end

function save_or_merge_sternberg_fooof_matrix(featPath, fooofData)
fooofMat = fullfile(featPath, 'AOC_eeg_matrix_sternberg_FOOOF.mat');
fooofCsv = fullfile(featPath, 'AOC_eeg_matrix_sternberg_FOOOF.csv');
if isempty(fooofData)
    warning('No sternberg FOOOF data to save.');
    return
end
Tout = dedupe_by_key(struct2table(fooofData), {'ID', 'Condition'});
eeg_data_sternberg_FOOOF = table2struct(Tout); %#ok<NASGU>
save(fooofMat, 'eeg_data_sternberg_FOOOF');
writetable(Tout, fooofCsv);
end

function T = dedupe_by_key(T, keyVars)
[~, ia] = unique(T(:, keyVars), 'rows', 'stable');
T = T(sort(ia), :);
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end + 1} = lab;
    end
end
if isempty(ch)
    ch = labels;
end
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
