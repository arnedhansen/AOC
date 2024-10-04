%% AOC Behavioral Feature Extraction Sternberg
% 
% Extracted features:
%   Accuracy
%   Reaction Time

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
behav_data_sternberg = struct('ID', {}, 'Condition', {}, 'Accuracy', {}, 'ReactionTime', {});

%% Read data
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    
    % Initialize subject-specific arrays
    subject_id = [];
    trial_num = [];
    condition = [];
    accuracy = [];
    reaction_time = [];
    stimuli = [];
    probe = [];
    match = [];
    
    %% Read blocks
    trial_counter = 1;
    for block = 1:6
        load(strcat(subjects{subj}, '_AOC_Sternberg_block', num2str(block), '_task.mat'))
        num_trials = length(saves.data.correct);
        
        % Append data for this block
        subject_id = [subject_id; repmat({saves.subject.ID}, num_trials, 1)];
        trial_num = [trial_num; (trial_counter:(trial_counter + num_trials - 1))'];
        condition = [condition; saves.data.trialSetSize'];
        accuracy = [accuracy; saves.data.correct'];
        reaction_time = [reaction_time; saves.data.reactionTime'];
        stimuli = [stimuli; saves.data.stimuli'];
        probe = [probe; saves.data.probe'];
        match = [match; saves.data.match'];
        trial_counter = trial_counter + num_trials;
    end
    % Convert ASCII numbers to letters
    probe = char(probe);
    % Set RT > 2 to NaN
    reaction_time(reaction_time > 2) = NaN;
    
    %% Create a trial-by-trial structure array for this subject
    subj_data_behav_trial = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', num2cell(condition), ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time), 'Stimuli', num2cell(stimuli), ...
        'Probe', num2cell(probe), 'Match', num2cell(match));

    %% Calculate subject-specific Acc and RT by condition
    l2 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 2);
    l2acc = sum([l2.Accuracy])/length(l2)*100;
    l2rt = mean([l2.ReactionTime], 'omitnan');
    l4 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 4);
    l4acc = sum([l4.Accuracy])/length(l4)*100;
    l4rt = mean([l4.ReactionTime], 'omitnan');
    l6 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 6);
    l6acc = sum([l6.Accuracy])/length(l6)*100;
    l6rt = mean([l6.ReactionTime], 'omitnan');

    %% Create across condition structure
    subj_data_behav = struct('ID', num2cell([subject_id{1}; subject_id{1}; subject_id{1}]), 'Condition', num2cell([2; 4; 6]), 'Accuracy', num2cell([l2acc; l4acc; l6acc]), 'ReactionTime', num2cell([l2rt; l4rt; l6rt]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/behavioral/');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_sternberg_subj_trial subj_data_behav_trial
    save behavioral_matrix_sternberg_subj subj_data_behav
    save acc_sternberg l2acc l4acc l6acc
    save rt_sternberg l2rt l4rt l6rt
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])
    
    % Append to the final structure array
    behav_data_sternberg = [behav_data_sternberg; subj_data_behav];
end
save /Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_sternberg behav_data_sternberg