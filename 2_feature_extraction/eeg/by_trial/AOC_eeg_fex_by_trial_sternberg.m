%% AOC EEG Feature Extraction Sternberg TRIAL-BY-TRIAL
%
% Extracted features:
%   Power Spectrum (Retention)
%   IAF, Power at IAF, and Lateralization Index

%% AOC Trial-by-Trial Alpha Power + Lateralisation (Retention Interval)
startup
[subjects, path, ~ , ~] = setup('AOC');
eeg_data_sternberg_trials = struct('ID', {}, 'Trial', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralisation', {});
alphaRange = [8 14];  % Hz

for subj = 1:length(subjects)
    clc
    disp(['Processing trial-by-trial features for AOC ', num2str(subjects{subj})])
    try
        % Load data
        datapath = fullfile(path, subjects{subj}, 'eeg');
        cd(datapath)
        load dataEEG_sternberg

        % Select retention interval: 1–2s after stimulus
        cfg = [];
        cfg.latency = [1 2];
        dat = ft_selectdata(cfg, dataEEG);

        % Frequency analysis: trial-by-trial
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper = 'dpss';
        cfg.tapsmofrq = 1;
        cfg.foilim = [3 30];
        cfg.keeptrials = 'yes';
        cfg.pad = 5;
        PowAll = ft_freqanalysis(cfg, dat);

        freqs = PowAll.freq;
        alphaIdx = find(freqs >= alphaRange(1) & freqs <= alphaRange(2));
        sid = str2double(subjects{subj});
        nTrials = size(PowAll.powspctrm,1);

        % Channel selection: occipital
        occIdx = find(contains(PowAll.label, 'O') & ~contains(PowAll.label, 'EOG') | contains(PowAll.label, {'I'}));

        % Separate into left/right occipital for lateralisation
        left_channels = {};
        right_channels = {};
        for i = 1:length(PowAll.label)
            ch = PowAll.label{i};
            if contains(ch, 'O') || contains(ch, 'I')
                numStr = regexp(ch, '\d+', 'match');
                if ~isempty(numStr)
                    numVal = str2double(numStr{1});
                    if mod(numVal,2)==1  % odd = left
                        left_channels{end+1} = ch;
                    else                 % even = right
                        right_channels{end+1} = ch;
                    end
                end
            end
        end
        left_idx = find(ismember(PowAll.label, left_channels));
        right_idx = find(ismember(PowAll.label, right_channels));

        % Loop through trials
        for t = 1:nTrials
            cond = PowAll.trialinfo(t);  % 22, 24, 26

            % IAF estimation
            spec = squeeze(mean(PowAll.powspctrm(t, occIdx, alphaIdx), 2));
            [pks, locs] = findpeaks(spec);

            if isempty(pks)
                IAF = NaN;
                powerIAF = NaN;
            else
                [~, maxIdx] = max(pks);
                IAF = freqs(alphaIdx(locs(maxIdx)));
                IAF_range = find(freqs > (IAF - 4) & freqs < (IAF + 2));
                powerIAF = mean(squeeze(mean(PowAll.powspctrm(t, occIdx, IAF_range), 2)));

                % Validity checks
                if IAF == alphaRange(1) || IAF == alphaRange(2) || powerIAF > max(pks)
                    IAF = NaN;
                    powerIAF = NaN;
                end
            end

            % Lateralisation index (Stroganova)
            left_alpha = mean(mean(PowAll.powspctrm(t, left_idx, alphaIdx), 2));
            right_alpha = mean(mean(PowAll.powspctrm(t, right_idx, alphaIdx), 2));
            if (right_alpha + left_alpha) == 0
                latIdx = NaN;
            else
                latIdx = (right_alpha - left_alpha) / (right_alpha + left_alpha);
            end

            % Add to struct
            eeg_data_sternberg_trials(end+1) = struct( ...
                'ID', sid, ...
                'Trial', t, ...
                'Condition', cond - 20, ...
                'AlphaPower', powerIAF, ...
                'IAF', IAF, ...
                'Lateralisation', latIdx );
        end

        % Save per subject
        save(fullfile(datapath, 'eeg_matrix_sternberg_subj_trials.mat'), 'eeg_data_sternberg_trials');
        save(fullfile(datapath, 'power_stern_trials.mat'), 'PowAll');

        fprintf('Subject %s done\n', subjects{subj});

    catch ME
        warning(['⚠️ ERROR with subject ', subjects{subj}, ': ', ME.message])
    end
end

% Save full group-level table
if ispc
    save('W:\Students\Arne\AOC\data\features\eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
else
    save('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat', 'eeg_data_sternberg_trials');
end

fprintf('Trial-by-trial alpha power, IAF and lateralisation COMPLETED. \n');
