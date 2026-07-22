%% AOC EEG Feature Extraction N Back Trial Level ERSD
% Compute trial level ERSD features for N back only.
% ERSD band is condition-wise IAF (IAF-4, IAF+2); fallback [8 14] when IAF undefined.
% Output: `AOC_eeg_matrix_nback_trials.mat` with `eeg_data_nback_trials`.
startup
[subjects, paths, ~, ~] = setup('AOC');
featPath = paths.features;
alphaRange = [8 14];
winIAF = [0 2];

eeg_data_nback_trials = struct( ...
    'ID', {}, ...
    'Trial', {}, ...
    'Condition', {}, ...
    'ERSD_early', {}, ...
    'ERSD_late', {}, ...
    'ERSD_full', {});

for subj = 1:length(subjects)
    try
        clc
        fprintf('[EEG TRIAL ERSD NBACK] Subject %d/%d (%s)\n', subj, length(subjects), subjects{subj});
        datapath = fullfile(featPath, subjects{subj}, 'eeg');
        cd(datapath)
        load dataEEG_TFR_nback

        condCodesRaw = dataTFR.trialinfo(:, 1);
        condSetRaw = unique(condCodesRaw(:))';
        if ~all(ismember(condSetRaw, [21 22 23]))
            error('Unexpected N back condition code(s): %s', mat2str(condSetRaw));
        end

        ind1 = find(condCodesRaw == 21);
        ind2 = find(condCodesRaw == 22);
        ind3 = find(condCodesRaw == 23);
        condIndices = {ind1, ind2, ind3};
        condOutVals = [1, 2, 3];

        occ_ch = occ_channels_from_labels(dataTFR.label);
        [IAF1, ~] = iaf_from_retention_mtmfft(dataTFR, ind1, winIAF, occ_ch, alphaRange);
        [IAF2, ~] = iaf_from_retention_mtmfft(dataTFR, ind2, winIAF, occ_ch, alphaRange);
        [IAF3, ~] = iaf_from_retention_mtmfft(dataTFR, ind3, winIAF, occ_ch, alphaRange);
        iafVals = [IAF1, IAF2, IAF3];

        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.output = 'pow';
        cfg.taper = 'hanning';
        cfg.foi = 2:2:40;
        cfg.t_ftimwin = ones(size(cfg.foi)) * 0.5;
        cfg.toi = -1.5:0.05:2.25;
        cfg.pad = 'nextpow2';
        cfg.keeptrials = 'yes';

        cfgb = [];
        cfgb.baseline = [-1.5 -0.5];
        cfgb.baselinetype = 'db';

        subID = str2double(subjects{subj});
        subj_data_eeg_trials = repmat(eeg_data_nback_trials, 0, 1);

        for c = 1:3
            trlIdx = condIndices{c};
            if isempty(trlIdx)
                continue
            end

            cfg.trials = trlIdx;
            tf = ft_freqanalysis(cfg, dataTFR);
            tf_bl = ft_freqbaseline(cfgb, tf);
            [ersdEarly, ersdLate, ersdFull] = compute_trial_ersd(tf_bl, iafVals(c), alphaRange);
            trialIDs = tf_bl.trialinfo(:, 2);
            nTrials = numel(trialIDs);

            condField = num2cell(repmat(condOutVals(c), nTrials, 1));
            subField = num2cell(repmat(subID, nTrials, 1));
            trialField = num2cell(trialIDs);

            condStruct = struct( ...
                'ID', subField, ...
                'Trial', trialField, ...
                'Condition', condField, ...
                'ERSD_early', num2cell(ersdEarly), ...
                'ERSD_late', num2cell(ersdLate), ...
                'ERSD_full', num2cell(ersdFull));
            subj_data_eeg_trials = [subj_data_eeg_trials; condStruct]; %#ok<AGROW>
        end

        save(fullfile(datapath, 'eeg_matrix_nback_subj_trials.mat'), 'subj_data_eeg_trials');
        eeg_data_nback_trials = [eeg_data_nback_trials; subj_data_eeg_trials]; %#ok<AGROW>
    catch ME
        fprintf('[EEG TRIAL ERSD NBACK] Failed for subject %s: %s\n', subjects{subj}, ME.message);
    end
end

save(fullfile(featPath, 'AOC_eeg_matrix_nback_trials.mat'), 'eeg_data_nback_trials');
writetable(struct2table(eeg_data_nback_trials), fullfile(featPath, 'AOC_eeg_matrix_nback_trials.csv'));

function [ersdEarly, ersdLate, ersdFull] = compute_trial_ersd(tf_bl, iafVal, alphaRange)
chUse = occ_channels_from_labels(tf_bl.label);
chIdx = find(ismember(tf_bl.label, chUse));
if isempty(chIdx)
    chIdx = 1:numel(tf_bl.label);
end

nTrials = size(tf_bl.powspctrm, 1);
ersdEarly = nan(nTrials, 1);
ersdLate = nan(nTrials, 1);
ersdFull = nan(nTrials, 1);

band = ersd_alpha_band(iafVal, alphaRange);
fMask = tf_bl.freq >= band(1) & tf_bl.freq <= band(2);
tMaskEarly = tf_bl.time >= 0 & tf_bl.time <= 1;
tMaskLate = tf_bl.time >= 1 & tf_bl.time <= 2;
tMaskFull = tf_bl.time >= 0 & tf_bl.time <= 2;

ersdEarly = reduce_trial_ersd(tf_bl.powspctrm, chIdx, fMask, tMaskEarly);
ersdLate = reduce_trial_ersd(tf_bl.powspctrm, chIdx, fMask, tMaskLate);
ersdFull = reduce_trial_ersd(tf_bl.powspctrm, chIdx, fMask, tMaskFull);
end

function out = reduce_trial_ersd(pow4d, chIdx, fMask, tMask)
nTrials = size(pow4d, 1);
out = nan(nTrials, 1);
if ~any(fMask) || ~any(tMask)
    return
end
for tr = 1:nTrials
    x = squeeze(pow4d(tr, chIdx, fMask, tMask));
    x = x(:);
    out(tr) = mean(x, 'omitnan');
end
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end + 1} = lab; %#ok<AGROW>
    end
end
if isempty(ch)
    ch = labels;
end
end

function [IAF, powerIAF] = iaf_from_retention_mtmfft(dataTFR, trialinds, winSec, chLabs, alphaRange)
IAF = NaN;
powerIAF = NaN;
if isempty(trialinds) || isempty(chLabs)
    return
end
try
    cfg_sel = [];
    cfg_sel.latency = winSec;
    cfg_sel.trials = trialinds(:);
    cfg_sel.channel = chLabs(:);
    dw = ft_selectdata(cfg_sel, dataTFR);
    if isempty(dw.trial)
        return
    end
    cfgf = [];
    cfgf.method = 'mtmfft';
    cfgf.output = 'pow';
    cfgf.taper = 'dpss';
    cfgf.tapsmofrq = 2;
    cfgf.foilim = [6 18];
    cfgf.pad = 'nextpow2';
    cfgf.keeptrials = 'no';
    fr = ft_freqanalysis(cfgf, dw);
    ps = mean(fr.powspctrm, 1);
    [IAF, powerIAF] = iaf_peak_rules(fr.freq(:), ps(:), alphaRange);
catch
    IAF = NaN;
    powerIAF = NaN;
end
end

function [IAF, powerIAF] = iaf_peak_rules(freq, spec, alphaRange)
freq = freq(:);
spec = spec(:);
IAF = NaN;
powerIAF = NaN;
amask = freq >= alphaRange(1) & freq <= alphaRange(2);
alphaFreqs = freq(amask);
alphaSpec = spec(amask);
if numel(alphaSpec) < 3
    return
end
[pks, locs] = findpeaks(alphaSpec);
if isempty(pks)
    return
end
[~, ind] = max(pks);
IAF = alphaFreqs(locs(ind));
bandIdx = freq > (IAF - 4) & freq < (IAF + 2);
if any(bandIdx)
    powerIAF = mean(spec(bandIdx));
end
if locs(ind) == 1 || locs(ind) == numel(alphaFreqs)
    powerIAF = NaN;
end
end

function band = ersd_alpha_band(IAF, alphaRange)
% ERSD band: (IAF-4, IAF+2) when IAF is valid; otherwise fixed alphaRange ([8 14]).
if ~isfinite(IAF)
    band = alphaRange;
else
    band = [IAF - 4, IAF + 2];
end
if any(~isfinite(band)) || band(1) >= band(2)
    band = alphaRange;
end
end
