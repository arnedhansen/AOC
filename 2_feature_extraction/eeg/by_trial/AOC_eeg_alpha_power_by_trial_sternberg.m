%% AOC Alpha Power Sternberg (trial-by-trial)

%% Setup
clear
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/features/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, '/eeg');
cd(datapath);
load('power_stern_trials.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload2_trials.label)
    label = powload2_trials.label{i};
    if contains(label, {'O'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data and calculate alpha power and IAF
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load('power_stern_trials.mat');
    channelIdx = find(ismember(powload2_trials.label, channels));
    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload2_trials.freq >= alphaRange(1) & powload2_trials.freq <= alphaRange(2));

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = length(powload2_trials.trialinfo)+length(powload4_trials.trialinfo)+length(powload6_trials.trialinfo);
    condition = [];
    alpha = [];
    IAF = [];

    %% Extract power spectra
    for trl = 1:num_trials
        if trl <= length(powload2_trials.trialinfo)
            %% Extract power spectra for selected channels and calculate IAF for WM load 2
            powspctrm2 = powload2_trials.powspctrm(trl, channelIdx, :);
            alphaPower2 = powspctrm2(alphaIndices);
            [~, maxIndex2] = max(alphaPower2);
            IAFsub = powload2_trials.freq(alphaIndices(maxIndex2));
            cond = 2;
            alphaPower = alphaPower2(maxIndex2);
        elseif trl <= length(powload4_trials.trialinfo)+length(powload2_trials.trialinfo)
            %% Extract power spectra for selected channels and calculate IAF for WM load 4
            powspctrm4 = powload4_trials.powspctrm(trl-length(powload2_trials.trialinfo), channelIdx, :);
            alphaPower4 = powspctrm4(alphaIndices);
            [~, maxIndex4] = max(alphaPower4);
            IAFsub = powload4_trials.freq(alphaIndices(maxIndex4));
            cond = 4;
            alphaPower = alphaPower4(maxIndex4);
        else 
            %% Extract power spectra for selected channels and calculate IAF for WM load 6
            powspctrm6 = powload6_trials.powspctrm(trl-(length(powload4_trials.trialinfo)+length(powload2_trials.trialinfo)), channelIdx, :);
            alphaPower6 = powspctrm6(alphaIndices);
            [~, maxIndex6] = max(alphaPower6);
            IAFsub = powload6_trials.freq(alphaIndices(maxIndex6));
            cond = 6;
            alphaPower = alphaPower6(maxIndex6);
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
        l2 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 2);
        l2alphapow = mean([l2.AlphaPower], 'omitnan');
        l4 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 4);
        l4alphapow = mean([l4.AlphaPower], 'omitnan');
        l6 = subj_data_eeg_trial([subj_data_eeg_trial.Condition] == 6);
        l6alphapow = mean([l6.AlphaPower], 'omitnan');
    end

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/eeg/');
    mkdir(savepath)
    cd(savepath)
    save eeg_matrix_sternberg_subj_trial subj_data_eeg_trial
    save alpha_power_sternberg_trial l2alphapow l4alphapow l6alphapow
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

end