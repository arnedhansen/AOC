%% AOC Alpha Power Nback (trial-by-trial)

%% Setup
clear
[subjects, path, ~ , ~] = setup('AOC');

%% Load data
nback_trial_matrix = [];
for subj = 1:3%%%%%:length(subjects)
    datapath = strcat(path, subjects{subj}, '/eeg');
    cd(datapath);
    load eeg_matrix_nback_subj_trial subj_data_eeg_trial
    nback_trial_matrix = [nback_trial_matrix; subj_data_eeg_trial];
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' loaded.'])
end

%% Add condition Order info
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
for subj = 1:length(subjects)
    subjID = subjects{subj};
    datapath = strcat(path, num2str(subjID));
    cd(datapath)
    files = dir(['/Volumes/methlab_data/OCC/AOC/data/', num2str(subjID), filesep, num2str(subjID), '_AOC_Nback_block*_task.mat']);
    for i = 1:numel(files)
        fname = files(i).name;
        tok = regexp( fname, ...
            '_NBack_block(\d+)_(\d)back_task\.mat$', ...
            'tokens', 'once' );
        blockIdx = str2double( tok{1} );
        backLev  = str2double( tok{2} );
        condOrder(blockIdx) = backLev;
    end
    disp(condOrder)

    % Save in data
    nback_trial_matrix.CondOrder([nback_trial_matrix.ID] == str2num(subjID)) = [{condOrder}; {condOrder}; {condOrder}];
end
disp('Finding Condition Orders DONE')

%% Visualize averaged trial-by-trial alpha power for each condition

%% Visualize for each indiividual participant the trial-by-trial alpha power over blocks