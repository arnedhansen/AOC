%% AOC Behavioral Feature Extraction N-back
% 
% Extracted features:
%   Accuracy
%   Reaction Time

%% Setup
clear
clc
close all
path = '/Volumes/methlab_data/OCC/AOC/data/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
behav_data_nback = struct('ID', {}, 'Condition', {}, 'Accuracy', {}, 'ReactionTime', {});

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
        fileNback = dir(strcat(subjects{subj}, '_AOC_Nback_block', num2str(block), '_*_task.mat'));
        load(fileNback.name)
        num_trials = length(saves.data.correct);

        % Extract condition from file name
        tokens = regexp(fileNback.name, '_\d+back_', 'match');
        token = tokens{1}; % There should be only one match
        cond = str2double(token(2)); % Convert the number to double
        
        % Append data for this block
        subject_id = [subject_id; repmat({saves.subject.ID}, num_trials, 1)];
        trial_num = [trial_num; (trial_counter:(trial_counter + num_trials - 1))'];
        condition = [condition; repmat({cond}, num_trials, 1)];
        accuracy = [accuracy; saves.data.correct'];
        reaction_time = [reaction_time; saves.data.reactionTime'];
        stimuli = [stimuli; saves.data.stimulus'];
        match = [match; saves.data.match'];
        trial_counter = trial_counter + num_trials;
    end
    % Convert ASCII numbers to letters
    stimuli = char(stimuli);
    % Set RT > 2 to NaN
    reaction_time(reaction_time > 2) = NaN;
    
    %% Create a structure array for this subject
    subj_data_behav_trial = struct('ID', subject_id, 'Trial', num2cell(trial_num), 'Condition', condition, ...
        'Accuracy', num2cell(accuracy), 'ReactionTime', num2cell(reaction_time), 'Stimuli', num2cell(stimuli), 'Match', num2cell(match));

    %% Calculate subject-specific Acc and RT by condition
    l1 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 1);
    l1acc = sum([l1.Accuracy], 'omitnan')/length(l1)*100;
    l1rt = mean([l1.ReactionTime], 'omitnan');
    l2 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 2);
    l2acc = sum([l2.Accuracy], 'omitnan')/length(l2)*100;
    l2rt = mean([l2.ReactionTime], 'omitnan');
    l3 = subj_data_behav_trial([subj_data_behav_trial.Condition] == 3);
    l3acc = sum([l3.Accuracy], 'omitnan')/length(l3)*100;
    l3rt = mean([l3.ReactionTime], 'omitnan');

    %% Create across condition structure
    subj_data_behav = struct('ID', num2cell([subject_id{1}; subject_id{1}; subject_id{1}]), 'Condition', num2cell([1; 2; 3]), ...
        'Accuracy', num2cell([l1acc; l2acc; l3acc]), 'ReactionTime', num2cell([l1rt; l2rt; l3rt]));

    %% Save
    savepath = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/',subjects{subj}, '/behavioral/');
    mkdir(savepath)
    cd(savepath)
    save behavioral_matrix_nback_subj_trial subj_data_behav_trial
    save behavioral_matrix_nback_subj subj_data_behav
    save acc_nback l1acc l2acc l3acc
    save rt_nback l1rt l2rt l3rt
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])
    
    % Append to the final structure array
    behav_data_nback = [behav_data_nback; subj_data_behav];
end
save /Volumes/methlab/Students/Arne/AOC/data/features/behavioral_matrix_nback behav_data_nback