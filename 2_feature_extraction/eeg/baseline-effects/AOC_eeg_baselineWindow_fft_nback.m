%% AOC EEG Baseline Window FFT Feature Extraction — N-Back
% Baseline-window multitaper FFT [-1.5 -0.5] s per n-back level.
% Occipital ROI: channel labels containing O or I.
% Each subject×condition ROI spectrum is sum-normalized over foilim so that
% absolute scale is removed and AlphaPower_baselineWindow is a relative
% measure comparable across subjects.
% IAF is peaked in [8 14] Hz; alpha power is the mean of the normalized
% spectrum in (IAF-4, IAF+2).
%
% Key outputs:
%   AOC_eeg_matrix_nback_baselineWindow.mat / .csv
%     ID, Condition, IAF_baselineWindow, AlphaPower_baselineWindow

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
featPath = paths.features;

winBase = [-1.5 -0.5];
alphaRange = [8 14];
foilim = [6 18];

eeg_data_nback_baselineWindow = struct( ...
    'ID', {}, 'Condition', {}, ...
    'IAF_baselineWindow', {}, 'AlphaPower_baselineWindow', {});

%% Per-subject extraction
for subj = 1:length(subjects)
    try
        clc; fprintf('[EEG BASELINEWINDOW FFT - NBACK] Subject %d / %d (%s)\n', ...
            subj, length(subjects), subjects{subj});
        datapath = fullfile(featPath, subjects{subj}, 'eeg');
        cd(datapath)
        load dataEEG_TFR_nback

        ind1 = find(dataTFR.trialinfo(:, 1) == 21);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind3 = find(dataTFR.trialinfo(:, 1) == 23);
        trialPack = {ind1, ind2, ind3};
        condVals = [1; 2; 3];

        occ_channels = occ_channels_from_labels(dataTFR.label);
        if isempty(occ_channels)
            error('No O/I occipital channels found for subject %s', subjects{subj});
        end

        IAF_b = nan(3, 1);
        Pow_b = nan(3, 1);
        for c = 1:3
            [IAF_b(c), Pow_b(c)] = iaf_alpha_from_baselineWindow_mtmfft( ...
                dataTFR, trialPack{c}, winBase, occ_channels, alphaRange, foilim);
        end

        subID = str2double(subjects{subj});
        subj_data = struct( ...
            'ID', num2cell([subID; subID; subID]), ...
            'Condition', num2cell(condVals), ...
            'IAF_baselineWindow', num2cell(IAF_b), ...
            'AlphaPower_baselineWindow', num2cell(Pow_b));

        save(fullfile(datapath, 'eeg_matrix_nback_baselineWindow_subj.mat'), 'subj_data');
        eeg_data_nback_baselineWindow = [eeg_data_nback_baselineWindow; subj_data];

        fprintf('  IAF_blWin 1/2/3-back: %.2f / %.2f / %.2f | Alpha_blWin: %.4f / %.4f / %.4f\n', ...
            IAF_b(1), IAF_b(2), IAF_b(3), Pow_b(1), Pow_b(2), Pow_b(3));
    catch ME
        fprintf('Continuing to next subject (%s): %s\n', subjects{subj}, ME.message);
    end
end

%% Save group matrices
outMat = fullfile(featPath, 'AOC_eeg_matrix_nback_baselineWindow.mat');
outCsv = fullfile(featPath, 'AOC_eeg_matrix_nback_baselineWindow.csv');
save(outMat, 'eeg_data_nback_baselineWindow');
writetable(struct2table(eeg_data_nback_baselineWindow), outCsv);
fprintf('Saved %s\n', outMat);
fprintf('Saved %s\n', outCsv);

%% Local helpers
function [IAF, alphaPow] = iaf_alpha_from_baselineWindow_mtmfft(dataTFR, trialinds, winSec, chLabs, alphaRange, foilim)
IAF = NaN;
alphaPow = NaN;
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
    cfgf.foilim = foilim;
    cfgf.pad = 'nextpow2';
    cfgf.keeptrials = 'no';
    fr = ft_freqanalysis(cfgf, dw);

    % ROI-mean spectrum, then sum-normalize so total power over foilim = 1.
    % This removes absolute scale differences between subjects.
    ps = mean(fr.powspctrm, 1);
    ps = ps(:);
    denom = sum(ps(isfinite(ps)));
    if ~(isfinite(denom) && denom > 0)
        return
    end
    ps_norm = ps ./ denom;

    [IAF, alphaPow] = iaf_peak_rules(fr.freq(:), ps_norm, alphaRange);
catch
    IAF = NaN;
    alphaPow = NaN;
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

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end + 1} = lab; %#ok<AGROW>
    end
end
end
