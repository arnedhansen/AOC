%% AOC EEG FOOOF Feature Extraction - Sternberg
% Computes subject-level IAF_specParam (FOOOF alpha centre frequency) for sternberg.
% Output per subject: IAF_specParam_sternberg.mat with variable `fooof_specparam`.
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;

if exist('ft_freqanalysis_Arne_FOOOF', 'file') ~= 2
    error(['ft_freqanalysis_Arne_FOOOF not found on MATLAB path. ', ...
        'Add it first, for example from /Users/Arne/Documents/GitHub/functions.']);
end

alphaRange = [8 14];
winIAF = [1 2];

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

        occ_channels = {};
        for i = 1:length(dataTFR.label)
            label = dataTFR.label{i};
            if contains(label, {'O'}) || contains(label, {'I'})
                occ_channels{end + 1} = label; %#ok<AGROW>
            end
        end
        if isempty(occ_channels)
            error('No occipital channels found for subject %s.', subjects{subj});
        end

        IAF_specParam2 = iaf_specparam_fooof(dataTFR, ind2, winIAF, occ_channels, alphaRange);
        IAF_specParam4 = iaf_specparam_fooof(dataTFR, ind4, winIAF, occ_channels, alphaRange);
        IAF_specParam6 = iaf_specparam_fooof(dataTFR, ind6, winIAF, occ_channels, alphaRange);

        subID = str2double(subjects{subj});
        fooof_specparam = struct( ...
            'ID', num2cell([subID; subID; subID]), ...
            'Condition', num2cell([2; 4; 6]), ...
            'IAF_specParam', num2cell([IAF_specParam2; IAF_specParam4; IAF_specParam6]));

        save(fullfile(datapath, 'IAF_specParam_sternberg.mat'), 'fooof_specparam')
    catch ME
        fprintf('STERNBERG FOOOF failed for subject %s: %s\n', subjects{subj}, ME.message);
    end
end

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
