%% AOC Preprocessing for N-back task

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

%% Read data, segment and convert to FieldTrip data structure
tic;
for subj = 1:length(subjects)
    clearvars -except subjects subj path
    clear data dataTFR data1 data2 data3
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    fprintf('Processing Subject %s\n', subjects{subj})

    % Only process new data
    if ispc == 1
        newDataFolder = dir(['W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\eeg\dataEEG_nback.mat']);
    else
        newDataFolder = dir(['/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/dataEEG_nback.mat']);
    end
    %if isempty(newDataFolder)
        clear alleeg
        %% Read blocks
        for block = 1:6
            try % Do not load emtpy blocks
                load(strcat(subjects{subj}, '_EEG_ET_Nback_block',num2str(block),'_merged.mat'))
                alleeg{block} = EEG;
                clear EEG
                fprintf('Block %.1d loaded \n', block)
            catch ME
                disp(['ERROR loading Block ' num2str(block) '!'])
            end
        end

        % Skip subject if there is no N-back data
        if exist('alleeg', 'var') == 0
            fprintf('No N-back data... SKIPPING processing of Subject %s\n....', subjects{subj})
            continue;
        end

        %% Segment data into epochs -1.5s before and 2.5s after stim onset and
        %  convert to Fieldtrip data structure AND extract gaze metrics from raw EEG data
        epoch_window = [-1.5 2.5];
        analysis_window = [0 2]; % Analysis window for N-back eye metric extraction
        for block = 1:6
            % 21 = PRESENTATION1 (Trigger for letter presentation (1-back))
            % 22 = PRESENTATION2 (Trigger for letter presentation (2-back))
            % 23 = PRESENTATION3 (Trigger for letter presentation (3-back))

            %% Segment 1-back data
            try
                EEG1back = pop_epoch(alleeg{block}, {'21'}, epoch_window);

                % Exclude matching trials to avoid including motor responses
                matching_trials = find(strcmp({EEG1back.event.type}, '4'));
                exclude_epochs = unique([EEG1back.event(matching_trials).epoch]);
                EEG1back = pop_select(EEG1back, 'notrial', exclude_epochs);
                data1{block} = eeglab2fieldtrip(EEG1back, 'raw');

                % 1-back gaze metrics extraction
                gaze_metrics_1back = pop_epoch(alleeg{block}, {'21'}, analysis_window);
                trl_1back(block) = gaze_metrics_1back.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_1back.event(strcmp({gaze_metrics_1back.event.type}, 'L_blink') | strcmp({gaze_metrics_1back.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_1back.event(strcmp({gaze_metrics_1back.event.type}, 'L_saccade') | strcmp({gaze_metrics_1back.event.type}, 'R_saccade'));

                % Exclude saccades around blinks
                valid_saccades = 0;
                for s = 1:length(saccade_events)
                    saccade_time = saccade_events(s).latency;
                    % Check if this saccade is within 100 ms of any blink
                    near_blink = any(abs(saccade_time - blink_times) <= 50); % 50 samples = 100 ms
                    if ~near_blink
                        valid_saccades = valid_saccades + 1;
                    end
                end

                % Count 1-back gaze metrics
                sacc_1back(block) = valid_saccades;
                fix_1back(block) = sum(ismember({gaze_metrics_1back.event.type}, {'L_fixation', 'R_fixation'}));
                blink_1back(block) = numel(blink_times);

                fprintf('1-back data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing 1-back for block %d: %s\n', block, ME.message);
                data1{block} = struct;
            end

            %% Segment 2-back data
            try
                EEG2back = pop_epoch(alleeg{block}, {'22'}, epoch_window);

                % Exclude matching trials to avoid including motor responses
                matching_trials = find(strcmp({EEG2back.event.type}, '4'));
                exclude_epochs = unique([EEG2back.event(matching_trials).epoch]);
                EEG2back = pop_select(EEG2back, 'notrial', exclude_epochs);
                data2{block} = eeglab2fieldtrip(EEG2back, 'raw');

                % 2-back gaze metrics extraction
                gaze_metrics_2back = pop_epoch(alleeg{block}, {'22'}, analysis_window);
                trl_2back(block) = gaze_metrics_2back.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_2back.event(strcmp({gaze_metrics_2back.event.type}, 'L_blink') | strcmp({gaze_metrics_2back.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_2back.event(strcmp({gaze_metrics_2back.event.type}, 'L_saccade') | strcmp({gaze_metrics_2back.event.type}, 'R_saccade'));

                % Exclude saccades around blinks
                valid_saccades = 0;
                for s = 1:length(saccade_events)
                    saccade_time = saccade_events(s).latency;
                    % Check if this saccade is within 100 ms of any blink
                    near_blink = any(abs(saccade_time - blink_times) <= 50); % 50 samples = 100 ms
                    if ~near_blink
                        valid_saccades = valid_saccades + 1;
                    end
                end

                % Count 2-back gaze metrics
                sacc_2back(block) = valid_saccades;
                fix_2back(block) = sum(ismember({gaze_metrics_2back.event.type}, {'L_fixation', 'R_fixation'}));
                blink_2back(block) = numel(blink_times);

                fprintf('2-back data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing 2-back for block %d: %s\n', block, ME.message);
                data2{block} = struct;
            end

            %% Segment 3-back data
            try
                EEG3back = pop_epoch(alleeg{block}, {'23'}, epoch_window);
                % Exclude matching trials to avoid including motor responses
                matching_trials = find(strcmp({EEG3back.event.type}, '4'));
                exclude_epochs = unique([EEG3back.event(matching_trials).epoch]);
                EEG3back = pop_select(EEG3back, 'notrial', exclude_epochs);
                data3{block} = eeglab2fieldtrip(EEG3back, 'raw');

                % 3-back gaze metrics extraction
                gaze_metrics_3back = pop_epoch(alleeg{block}, {'23'}, analysis_window);
                trl_3back(block) = gaze_metrics_3back.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_3back.event(strcmp({gaze_metrics_3back.event.type}, 'L_blink') | strcmp({gaze_metrics_3back.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_3back.event(strcmp({gaze_metrics_3back.event.type}, 'L_saccade') | strcmp({gaze_metrics_3back.event.type}, 'R_saccade'));

                % Exclude saccades around blinks
                valid_saccades = 0;
                for s = 1:length(saccade_events)
                    saccade_time = saccade_events(s).latency;
                    % Check if this saccade is within 100 ms of any blink
                    near_blink = any(abs(saccade_time - blink_times) <= 50); % 50 samples = 100 ms
                    if ~near_blink
                        valid_saccades = valid_saccades + 1;
                    end
                end

                % Count 3-back gaze metrics
                sacc_3back(block) = valid_saccades;
                fix_3back(block) = sum(ismember({gaze_metrics_3back.event.type}, {'L_fixation', 'R_fixation'}));
                blink_3back(block) = numel(blink_times);

                fprintf('3-back data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing 3-back for block %d: %s\n', block, ME.message);
                data3{block} = struct;
            end
        end

        %% Remove empty blocks
        data1 = data1(~cellfun(@(x) isempty(fieldnames(x)), data1));
        data2 = data2(~cellfun(@(x) isempty(fieldnames(x)), data2));
        data3 = data3(~cellfun(@(x) isempty(fieldnames(x)), data3));

        %% Equalize labels
        update_labels(data1);
        update_labels(data2);
        update_labels(data3);

        % Fix labels for subjects with only left eye recording (AOC419)
        if strcmp(subjects{subj}, '419')
            data2{1}.label = data2{2}.label; % Reduce to 132 electrodes
            for i = 1:length(data2{1}.trial) % Reduce to 132 electrodes
                data2{1}.trial{i} = data2{1}.trial{i}(1:132, :);
            end
        end

        %% Append data for conditions
        cfg = [];
        cfg.keepsampleinfo = 'no';
        data1 = ft_appenddata(cfg,data1{:});
        data2 = ft_appenddata(cfg,data2{:});
        data3 = ft_appenddata(cfg,data3{:});

        %% Add trialinfo
        data1.trialinfo = zeros(numel(data1.trial), 1) + 21;
        disp(['Adding trialinfo 1-back in Block ' num2str(block) '!']);
        data2.trialinfo = zeros(numel(data2.trial), 1) + 22;
        disp(['Adding trialinfo 2-back in Block ' num2str(block) '!']);
        data3.trialinfo = zeros(numel(data3.trial), 1) + 23;
        disp(['Adding trialinfo 3-back in Block ' num2str(block) '!']);

        %% Append all data into single data file with appropriate trialinfo
        cfg = [];
        cfg.keepsampleinfo = 'no';
        data = ft_appenddata(cfg, data1, data2, data3);
        trialinfo = [data1.trialinfo; data2.trialinfo; data3.trialinfo];
        data.trialinfo = trialinfo;

        %% Pre-stim fixation check
        preStimWindow = [-0.5 0];
        fixThresh = 0.8; % 80% of trials should be within fixation box
        distOK = 45;     % 1 degree (dva) from the center
        [trialsToKeep, excludedTrialIdx, distL, invalidTrials] = fixCheck(data, preStimWindow, fixThresh, distOK);

        % Save excluded trials info
        preStimFixInfo.subject = subjects{subj};
        preStimFixInfo.excludedTrials = find(~trialsToKeep);
        preStimFixInfo.invalidTrials = invalidTrials;
        preStimFixInfo.totalTrials = numel(trialsToKeep);
        preStimFixInfo.keptTrials = find(trialsToKeep);
        preStimFixInfo.distL = distL;
        if ispc == 1
            savepathControlsFix = (['W:\Students\Arne\AOC\data\controls\preStimFixation\', subjects{subj}]);
            mkdir(savepathControlsFix)
            save([savepathControlsFix, filesep, 'AOC_preStimFixation_', subjects{subj}, '_nback'], "preStimFixInfo");
        else
            savepathControlsFix = ['/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/', subjects{subj}];
            mkdir(savepathControlsFix)
            save([savepathControlsFix, filesep, 'AOC_preStimFixation_', subjects{subj}, '_nback'], "preStimFixInfo");
        end

        % Continue analyses with correct fix trials
        data = ft_selectdata(struct('trials', trialsToKeep), data);
        trialinfo = data.trialinfo; 

        %% Get EyeTracking data
        cfg = [];
        cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA', 'R-GAZE-X'  'R-GAZE-Y' 'R-AREA'};
        dataet = ft_selectdata(cfg, data);
        dataETlong = dataet;
        dataETlong.trialinfo = trialinfo;

        %% Get EEG data (excl. ET and EOG data)
        cfg = [];
        cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA', '-R-GAZE-X'  '-R-GAZE-Y' '-R-AREA'};
        data = ft_selectdata(cfg, data);

        %% Resegment data to avoid filter ringing
        cfg = [];
        dataTFR = ft_selectdata(cfg, data); % TRF data long
        dataTFR.trialinfo = trialinfo;
        cfg = [];
        cfg.latency = [0 2]; % Time window for N-back task
        data = ft_selectdata(cfg, data);
        dataet = ft_selectdata(cfg, dataet);
        dataet.trialinfo = trialinfo;

        %% Re-reference data to average or common reference
        cfg = [];
        cfg.reref   = 'yes';
        cfg.refchannel = 'all';
        data = ft_preprocessing(cfg,data);
        data.trialinfo = trialinfo;

        %% Compute gaze metric data
        % 1-back gaze metrics average across trials
        saccades_1back = sum(sacc_1back(:)) / sum(trl_1back(:));
        fixations_1back = sum(fix_1back(:)) / sum(trl_1back(:));
        blinks_1back = sum(blink_1back(:)) / sum(trl_1back(:));

        % 2-back gaze metrics average across trials
        saccades_2back = sum(sacc_2back(:)) / sum(trl_2back(:));
        fixations_2back = sum(fix_2back(:)) / sum(trl_2back(:));
        blinks_2back = sum(blink_2back(:)) / sum(trl_2back(:));

        % 3-back gaze metrics average across trials
        saccades_3back = sum(sacc_3back(:)) / sum(trl_3back(:));
        fixations_3back = sum(fix_3back(:)) / sum(trl_3back(:));
        blinks_3back = sum(blink_3back(:)) / sum(trl_3back(:));

        %% Save data
        if ispc == 1
            savepathEEG = strcat('W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\eeg\');
        else
            savepathEEG = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepathEEG)
        cd(savepathEEG)
        save dataEEG_nback data
        save dataEEG_TFR_nback dataTFR
        if ispc == 1
            savepathET = strcat('W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\gaze\');
        else
            savepathET = strcat('/Volumes/methlab/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
        end
        mkdir(savepathET)
        cd(savepathET)
        save dataET_nback dataet dataETlong
        save gaze_metrics_nback saccades_1back fixations_1back blinks_1back ...
            saccades_2back fixations_2back blinks_2back ...
            saccades_3back fixations_3back blinks_3back
        clc
        if subj == length(subjects)
            disp(['Subject AOC ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. N-BACK PREPROCESSING FINALIZED.'])
        else
            disp(['Subject AOC ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
        end
    %else
    %    disp(['Subject ', num2str(subjects{subj}), ' already done. SKIPPING...'])
    %end
end
toc;
%finishedScriptMail;