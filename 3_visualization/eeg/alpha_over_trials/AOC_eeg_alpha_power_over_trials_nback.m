%% AOC Alpha Power change over trials N-back

%% Setup
startup
clear
addEEGLab
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\merged\';
else
    path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');

%% Define channels
subj = 1;
datapath = strcat(path, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_nback_trials.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(powload1_trials.label)
    label = powload1_trials.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

%% Load data
for subj = 1:length(subjects)
    try
        % Load data
        datapath = strcat(path, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(data.trialinfo == 21);
        ind2 = find(data.trialinfo == 22);
        ind3 = find(data.trialinfo == 23);

        % Select data
        cfg = [];
        cfg.latency = [0 2]; % Segment from 0 to 2 [seconds]
        dat = ft_selectdata(cfg,data);

        % Frequency analysis settings
        cfg = [];% empty config
        cfg.output = 'pow';% estimates power only
        cfg.method = 'mtmfft';% multi taper fft method
        cfg.taper = 'dpss';% multiple tapers
        cfg.tapsmofrq = 1;% smoothening frequency around foi
        cfg.foilim = [3 30];% frequencies of interest (foi)
        cfg.keeptrials = 'no';% do not keep single trials in output
        cfg.pad = 10;

        % Frequency analysis settings
        cfg.trials = ind1;
        powload1 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind2;
        powload2 = ft_freqanalysis(cfg,dat);
        cfg.trials = ind3;
        powload3 = ft_freqanalysis(cfg,dat);

        % Save raw power spectra
        cd(datapath)
        save power_nback powload1 powload2 powload3

    catch ME
        ME.message
        error(['ERROR extracting power for Subject ' num2str(subjects{subj}) '!'])
    end
end


%% Load data
% Load powspctrm data
for subj = 10%1:3%%%%%%%%%%length(subjects)
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load power_nback_trials
    
    % Compute power over time from ft_freqanalysis structure
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.channel = channels;
    cfg.frequency = [8 14];
    % powload1_trials in form: trials x electrodes x frequency bins
    alphapowl1{subj} = ft_freqdescriptives(cfg, powload1_trials);
    alphapowl2{subj} = ft_freqdescriptives(cfg, powload2_trials);
    alphapowl3{subj} = ft_freqdescriptives(cfg, powload3_trials);
    disp(['Subject ', num2str(subjects{subj}), ' loaded.'])
end

%% Plot INDIVIDUAL alpha power over time




%%%%%%%%%% PUT PLOTS AFTER ANOTHER LIKE PSF in Condition colors from
%%%%%%%%%% colors(1, :)





% Compute alpha powers over trials
alphapow1 = mean(alphapowl1{1, subj}.powspctrm, 3);
alphapow1 = mean(alphapow1, 2);
alphapow2 = mean(alphapowl2{1, subj}.powspctrm, 3);
alphapow2 = mean(alphapow2, 2);
alphapow3 = mean(alphapowl3{1, subj}.powspctrm, 3);
alphapow3 = mean(alphapow3, 2);

% Plot
figure;
sgtitle('Alpha Power over ')
set(gcf, 'Position', [0 0 600 1000], 'Color', 'W');
subplot(3, 1, 1)
bar(alphapow1, 'b')
legend('1-back');
xlabel('Trials');
ylabel('Alpha Power (8-14 Hz)');
title('1-back Grand Average Alpha Power over Trials');

subplot(3, 1, 2)
bar(alphapow2, 'g')
legend('2-back');
xlabel('Trials');
ylabel('Alpha Power (8-14 Hz)');
title('2-back Grand Average Alpha Power over Trials');

subplot(3, 1, 3)
bar(alphapow3, 'r')
legend('3-back');
xlabel('Trials');
ylabel('Alpha Power (8-14 Hz)');
title('3-back Grand Average Alpha Power over Trials');

