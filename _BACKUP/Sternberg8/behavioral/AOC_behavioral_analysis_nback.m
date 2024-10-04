%% AOC Behavioral Measures N-back

%% Setup
clear;
addpath('/Users/Arne/Documents/matlabtools/eeglab2024.0');
eeglab;
clc;
close all;
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Initialize results storage
results = struct('subject', {}, 'condition', {}, 'trial', {}, 'correct', {}, 'reaction_time', {});

%% Read data, segment and convert to FieldTrip data struct
for subj=1:length(subjects)
    datapath = strcat(path, subjects{subj});
    addpath(datapath);

    %% Read blocks
    conditions = {'1back', '2back', '3back'};
    event_codes = {'21', '22', '23'};
    
    for cond=1:length(conditions)
        filename = strcat(subjects{subj}, '_EEG_ET_Nback_', conditions{cond}, '_merged.mat');
        load(filename);
        EEG = pop_epoch(EEG, {event_codes{cond}}, [-1.5 2.5]);
        
        %% Process events for accuracy and reaction time
        for trial=1:length(EEG.epoch)
            events = EEG.epoch(trial).eventtype;
            latencies = EEG.epoch(trial).eventlatency;
            
            correct = false;
            reaction_time = NaN;
            
            % Check for matching trial
            match_event = find(strcmp(events, '4'));
            match_response = find(strcmp(events, '87'));
            nonmatch_event = find(strcmp(events, '5'));
            nonmatch_no_response = find(strcmp(events, '88'));
            wrong_button = find(strcmp(events, '89'));
            stimulus_present = find(strcmp(events, event_codes{cond}));
            
            % Determine correctness and reaction time
            if ~isempty(match_event) && ~isempty(match_response)
                correct = true;
                reaction_time = latencies{match_response} - latencies{stimulus_present};
            elseif ~isempty(nonmatch_event) && ~isempty(nonmatch_no_response)
                correct = true;
                reaction_time = NaN; % No reaction time for correct non-responses
            elseif ~isempty(match_event) && isempty(match_response)
                correct = false;
                reaction_time = NaN;
            elseif ~isempty(nonmatch_event) && isempty(nonmatch_no_response)
                correct = false;
                reaction_time = NaN;
            elseif ~isempty(wrong_button)
                correct = false;
                reaction_time = NaN;
            end
            
            % Store results
            results(end+1) = struct('subject', subjects{subj}, 'condition', conditions{cond}, 'trial', trial, 'correct', correct, 'reaction_time', reaction_time);
        end
    end
end

%% Calculate accuracy for each condition and save
accuracy_results = struct('subject', {}, 'condition', {}, 'accuracy', {}, 'mean_reaction_time', {});
for subj=1:length(subjects)
    for cond=1:length(conditions)
        subject_trials = results(strcmp({results.subject}, subjects{subj}) & strcmp({results.condition}, conditions{cond}));
        correct_trials = [subject_trials.correct];
        reaction_times = [subject_trials.reaction_time];
        
        accuracy = sum(correct_trials) / length(correct_trials) * 100;
        mean_reaction_time = mean(reaction_times(correct_trials & ~isnan(reaction_times))); % Calculate mean reaction time for correct trials only
        
        accuracy_results(end+1) = struct('subject', subjects{subj}, 'condition', conditions{cond}, 'accuracy', accuracy, 'mean_reaction_time', mean_reaction_time); %#ok<AGROW>
    end
end

%% Save
save('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral/accuracy.mat', 'accuracy_results');
save('/Volumes/methlab/Students/Arne/AOC/data/features/behavioral/reaction_times.mat', 'results');
