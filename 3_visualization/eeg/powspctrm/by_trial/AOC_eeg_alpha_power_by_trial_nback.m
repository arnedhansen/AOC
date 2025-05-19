%% AOC Alpha Power Nback (trial-by-trial)

%% Setup
startup
[subjects, path, ~ , ~] = setup('AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
load('power_nback_trials.mat');
occ_channels = {};
for i = 1:length(powload1_trials.label)
    label = powload1_trials.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data and calculate alpha power and IAF
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load('power_nback_trials.mat');
    channelIdx = find(ismember(powload1_trials.label, channels));

    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload1_trials.freq >= alphaRange(1) & powload1_trials.freq <= alphaRange(2));

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = length(powload1_trials.trialinfo)+length(powload2_trials.trialinfo)+length(powload3_trials.trialinfo);
    condition = [];
    alpha = [];
    IAF = [];

    %% Extract power spectra
    for trl = 1:num_trials
        if trl <= length(powload1_trials.trialinfo)
            %% Extract power spectra for selected channels and calculate IAF for 1-back
            powspctrm1 = powload1_trials.powspctrm(trl, channelIdx, :);
            alphaPower1 = powspctrm1(alphaIndices);
            [~, maxIndex1] = max(alphaPower1);
            IAFsub = powload1_trials.freq(alphaIndices(maxIndex1));
            cond = 1;
            alphaPower = alphaPower1(maxIndex1);
        elseif trl <= length(powload2_trials.trialinfo)+length(powload1_trials.trialinfo)
            %% Extract power spectra for selected channels and calculate IAF for 2-back
            powspctrm2 = powload2_trials.powspctrm(trl-length(powload1_trials.trialinfo), channelIdx, :);
            alphaPower2 = powspctrm2(alphaIndices);
            [~, maxIndex2] = max(alphaPower2);
            IAFsub = powload2_trials.freq(alphaIndices(maxIndex2));
            cond = 2;
            alphaPower = alphaPower2(maxIndex2);
        else 
            %% Extract power spectra for selected channels and calculate IAF for 3-back
            powspctrm3 = powload3_trials.powspctrm(trl-(length(powload2_trials.trialinfo)+length(powload1_trials.trialinfo)), channelIdx, :);
            alphaPower3 = powspctrm3(alphaIndices);
            [~, maxIndex3] = max(alphaPower3);
            IAFsub = powload3_trials.freq(alphaIndices(maxIndex3));
            cond = 3;
            alphaPower = alphaPower3(maxIndex3);
        end

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; trl];
        condition = [condition; cond];
        alpha = [alpha; alphaPower];
        IAF = [IAF; IAFsub];

        %% Create a structure array for this subject
        subj_data_eeg_trial = struct('ID', num2cell(subject_id), 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), 'AlphaPower', num2cell(alpha), 'IAF', num2cell(IAF));

        %% Calculate subject-specific AlphaPower by condition
        l1 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 1);
        l1alphapow = mean([l1.AlphaPower], 'omitnan');
        l2 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 2);
        l2alphapow = mean([l2.AlphaPower], 'omitnan');
        l3 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 3);
        l3alphapow = mean([l3.AlphaPower], 'omitnan');
    end

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_nback_subj_trial subj_data_eeg_trial
    save alpha_power_nback_trial l1alphapow l2alphapow l3alphapow
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

end