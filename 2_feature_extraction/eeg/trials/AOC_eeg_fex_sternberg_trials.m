%% AOC EEG Feature Extraction Sternberg Trial Level ERSD
% Compute trial level ERSD features for Sternberg only.
% Output: `AOC_eeg_matrix_sternberg_trials.mat` with `eeg_data_sternberg_trials`.
startup
[subjects, paths, ~, ~] = setup('AOC');
featPath = paths.features;

eeg_data_sternberg_trials = struct( ...
    'ID', {}, ...
    'Trial', {}, ...
    'Condition', {}, ...
    'ERSD_early', {}, ...
    'ERSD_late', {}, ...
    'ERSD_full', {});

for subj = 1:length(subjects)
    try
        clc
        fprintf('[EEG TRIAL ERSD STERNBERG] Subject %d/%d (%s)\n', subj, length(subjects), subjects{subj});
        datapath = fullfile(featPath, subjects{subj}, 'eeg');
        cd(datapath)
        load dataEEG_TFR_sternberg

        condCodesRaw = dataTFR.trialinfo(:, 1);
        condSetRaw = unique(condCodesRaw(:))';
        if ~all(ismember(condSetRaw, [22 24 26]))
            error('Unexpected Sternberg condition code(s): %s', mat2str(condSetRaw));
        end

        ind2 = find(condCodesRaw == 22);
        ind4 = find(condCodesRaw == 24);
        ind6 = find(condCodesRaw == 26);
        condIndices = {ind2, ind4, ind6};
        condOutVals = [2, 4, 6];

        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.output = 'pow';
        cfg.taper = 'hanning';
        cfg.foi = 2:2:40;
        cfg.t_ftimwin = ones(size(cfg.foi)) * 0.5;
        cfg.toi = -1.5:0.05:3;
        cfg.pad = 'nextpow2';
        cfg.keeptrials = 'yes';

        cfgb = [];
        cfgb.baseline = [-1.5 -0.5];
        cfgb.baselinetype = 'db';

        subID = str2double(subjects{subj});
        subj_data_eeg_trials = repmat(eeg_data_sternberg_trials, 0, 1);

        for c = 1:3
            trlIdx = condIndices{c};
            if isempty(trlIdx)
                continue
            end

            cfg.trials = trlIdx;
            tf = ft_freqanalysis(cfg, dataTFR);
            tf_bl = ft_freqbaseline(cfgb, tf);
            [ersdEarly, ersdLate, ersdFull] = compute_trial_ersd(tf_bl);
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

        save(fullfile(datapath, 'eeg_matrix_sternberg_subj_trials.mat'), 'subj_data_eeg_trials');
        eeg_data_sternberg_trials = [eeg_data_sternberg_trials; subj_data_eeg_trials]; %#ok<AGROW>
    catch ME
        fprintf('[EEG TRIAL ERSD STERNBERG] Failed for subject %s: %s\n', subjects{subj}, ME.message);
    end
end

save(fullfile(featPath, 'AOC_eeg_matrix_sternberg_trials.mat'), 'eeg_data_sternberg_trials');
writetable(struct2table(eeg_data_sternberg_trials), fullfile(featPath, 'AOC_eeg_matrix_sternberg_trials.csv'));

function [ersdEarly, ersdLate, ersdFull] = compute_trial_ersd(tf_bl)
chUse = occ_channels_from_labels(tf_bl.label);
chIdx = find(ismember(tf_bl.label, chUse));
if isempty(chIdx)
    chIdx = 1:numel(tf_bl.label);
end

fMask = tf_bl.freq >= 8 & tf_bl.freq <= 14;
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
