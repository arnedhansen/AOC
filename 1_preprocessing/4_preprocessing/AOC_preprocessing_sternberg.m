%% AOC Preprocessing — Sternberg EEG to FieldTrip
% Reads merged EEG–ET blocks, segments into trials (load 2/4/6), converts to FieldTrip, applies baseline and filtering. Produces dataEEG_sternberg.mat per subject.
%
% Key outputs:
%   dataEEG_sternberg.mat (FieldTrip) with trialinfo, WM load 2/4/6
%   Preprocessed, epoched, and filtered Sternberg EEG

%% Setup
startup
clear
[~, ~, ~ , ~] = setup('AOC');
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\merged\';
else
    path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/merged/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');
runMode = askRunMode();

% Setup logging
logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
scriptName = 'AOC_preprocessing_sternberg';

%% Read data, segment and convert to FieldTrip data structure
tic;
for subj = 1:length(subjects)
    try
    datapath = strcat(path,subjects{subj});
    cd(datapath)

    % Only process new data
    if ispc == 1
        newDataFolder = dir(['W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\eeg\dataEEG_sternberg.mat']);
    else
        newDataFolder = dir(['/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/dataEEG_sternberg.mat']);
    end

    if strcmp(runMode,'all') || isempty(newDataFolder)
        clc
        disp(['Preprocessing Subject AOC ', num2str(subjects{subj})])
        clear alleeg data2 data4 data6
        %% Read blocks
        %  and assign globalID values
        trialCounter = 0;
        for block = 1:6
            try % Do not load emtpy blocks
                load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
                alleeg{block} = EEG;

                % Ensure field exists for all events
                if ~isfield(EEG.event, 'globalTrial')
                    [EEG.event.globalTrial] = deal(0);
                end

                % Identify stimulus events
                stimCodes = {'22','24','26'};
                stimEventIdx = find(ismember({EEG.event.type}, stimCodes));
                [~, sortOrder] = sort([EEG.event(stimEventIdx).latency]);
                stimEventIdx = stimEventIdx(sortOrder);

                % Assign chronological trial IDs
                for t = 1:length(stimEventIdx)
                    EEG.event(stimEventIdx(t)).globalTrial = t + trialCounter;
                end

                % Update counter and save back
                trialCounter = trialCounter + length(stimEventIdx);
                alleeg{block} = EEG;

                clear EEG
                fprintf('Block %.1d loaded \n', block)
            catch ME
                fprintf('[WARNING] Subject %s (iteration %d/%d): Error loading Block %d: %s\n', ...
                    subjects{subj}, subj, length(subjects), block, ME.message);
            end
        end

        % Skip subject if there is no data
        if exist('alleeg', 'var') == 0
            fprintf('No Sternberg data... SKIPPING processing of Subject %s\n....', subjects{subj})
            continue;
        end

        %% Segment data into epochs -2s before and 3.5s after stim onset and
        %  convert to Fieldtrip data structure AND extract gaze metrics from raw EEG data
        epoch_window = [-2 3.5];
        analysis_window = [1 2]; % Analysis window for Sternberg eye metric extraction
        for block = 1:6
            % PRESENTATION2 = 22 (Trigger for stimuli presentation (setSize = 2))
            % PRESENTATION4 = 24 (Trigger for stimuli presentation (setSize = 4))
            % PRESENTATION6 = 26 (Trigger for stimuli presentation (setSize = 6))
            % 52 = RETENTION2 (Trigger for retention (setSize = 2))
            % 54 = RETENTION4 (Trigger for retention (setSize = 4))
            % 56 = RETENTION6 (Trigger for retention (setSize = 6))

            %% Segment WM load 2 data
            try
                EEGload2 = pop_epoch(alleeg{block}, {'22'}, epoch_window);
                data2{block} = eeglab2fieldtrip(EEGload2, 'raw');

                % Add global trial ID info
                global_ids = nan(numel(EEGload2.epoch),1);
                for e = 1:numel(EEGload2.epoch)
                    evTypes = {EEGload2.event(EEGload2.epoch(e).event).type};
                    matchIdx = find(ismember(evTypes, {'22'}), 1, 'first');
                    if ~isempty(matchIdx)
                        evIndex = EEGload2.epoch(e).event(matchIdx);
                        global_ids(e) = EEGload2.event(evIndex).globalTrial;
                    else
                        global_ids(e) = NaN;
                    end
                end
                data2{block}.trialinfo = [repmat(22, numel(global_ids), 1), global_ids(:)];

                % WM load 2 gaze metrics extraction
                gaze_metrics_l2 = pop_epoch(alleeg{block}, {'22'}, analysis_window);
                trl_l2(block) = gaze_metrics_l2.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_l2.event(strcmp({gaze_metrics_l2.event.type}, 'L_blink') | strcmp({gaze_metrics_l2.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_l2.event(strcmp({gaze_metrics_l2.event.type}, 'L_saccade') | strcmp({gaze_metrics_l2.event.type}, 'R_saccade'));

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

                % Count WM load 2 gaze metrics
                sacc_l2(block) = valid_saccades;
                fix_l2(block) = sum(ismember({gaze_metrics_l2.event.type}, {'L_fixation', 'R_fixation'}));
                blink_l2(block) = numel(blink_times);

                fprintf('WM load 2 data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing WM load 2 for block %d: %s\n', block, ME.message);
                data2{block} = struct;
            end

            %% Segment WM load 4 data
            try
                EEGload4 = pop_epoch(alleeg{block}, {'24'}, epoch_window);
                data4{block} = eeglab2fieldtrip(EEGload4, 'raw');

                % Add global trial ID info
                global_ids = nan(numel(EEGload4.epoch),1);
                for e = 1:numel(EEGload4.epoch)
                    evTypes = {EEGload4.event(EEGload4.epoch(e).event).type};
                    matchIdx = find(ismember(evTypes, {'24'}), 1, 'first');
                    if ~isempty(matchIdx)
                        evIndex = EEGload4.epoch(e).event(matchIdx);
                        global_ids(e) = EEGload4.event(evIndex).globalTrial;
                    else
                        global_ids(e) = NaN;
                    end
                end
                data4{block}.trialinfo = [repmat(24, numel(global_ids), 1), global_ids(:)];

                % WM load 4 gaze metrics extraction
                gaze_metrics_l4 = pop_epoch(alleeg{block}, {'24'}, analysis_window);
                trl_l4(block) = gaze_metrics_l4.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_l4.event(strcmp({gaze_metrics_l4.event.type}, 'L_blink') | strcmp({gaze_metrics_l4.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_l4.event(strcmp({gaze_metrics_l4.event.type}, 'L_saccade') | strcmp({gaze_metrics_l4.event.type}, 'R_saccade'));

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

                % Count WM load 4 gaze metrics
                sacc_l4(block) = valid_saccades;
                fix_l4(block) = sum(ismember({gaze_metrics_l4.event.type}, {'L_fixation', 'R_fixation'}));
                blink_l4(block) = numel(blink_times);

                fprintf('WM load 4 data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing WM load 4 for block %d: %s\n', block, ME.message);
                data4{block} = struct;
            end

            %% Segment WM load 6 data
            try
                EEGload6 = pop_epoch(alleeg{block}, {'26'}, epoch_window);
                data6{block} = eeglab2fieldtrip(EEGload6, 'raw');

                % Add global trial ID info
                global_ids = nan(numel(EEGload6.epoch),1);
                for e = 1:numel(EEGload6.epoch)
                    evTypes = {EEGload6.event(EEGload6.epoch(e).event).type};
                    matchIdx = find(ismember(evTypes, {'26'}), 1, 'first');
                    if ~isempty(matchIdx)
                        evIndex = EEGload6.epoch(e).event(matchIdx);
                        global_ids(e) = EEGload6.event(evIndex).globalTrial;
                    else
                        global_ids(e) = NaN;
                    end
                end
                data6{block}.trialinfo = [repmat(26, numel(global_ids), 1), global_ids(:)];

                % WM load 6 gaze metrics extraction
                gaze_metrics_l6 = pop_epoch(alleeg{block}, {'26'}, analysis_window);
                trl_l6(block) = gaze_metrics_l6.trials;
                % Extract blink timepoints
                blink_times = [gaze_metrics_l6.event(strcmp({gaze_metrics_l6.event.type}, 'L_blink') | strcmp({gaze_metrics_l6.event.type}, 'R_blink')).latency];
                % Extract saccades timepoints
                saccade_events = gaze_metrics_l6.event(strcmp({gaze_metrics_l6.event.type}, 'L_saccade') | strcmp({gaze_metrics_l6.event.type}, 'R_saccade'));

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

                % Count WM load 6 gaze metrics
                sacc_l6(block) = valid_saccades;
                fix_l6(block) = sum(ismember({gaze_metrics_l6.event.type}, {'L_fixation', 'R_fixation'}));
                blink_l6(block) = numel(blink_times);

                fprintf('WM load 6 data processed for block %d.\n', block);
            catch ME
                fprintf('ERROR processing WM load 6 for block %d: %s\n', block, ME.message);
                data6{block} = struct;
            end
        end

        %% Remove empty blocks
        data2 = data2(~cellfun(@(x) isempty(fieldnames(x)), data2));
        data4 = data4(~cellfun(@(x) isempty(fieldnames(x)), data4));
        data6 = data6(~cellfun(@(x) isempty(fieldnames(x)), data6));

        %% Append data for conditions
        cfg = [];
        cfg.keepsampleinfo = 'yes';
        data2 = ft_appenddata(cfg,data2{:});
        data4 = ft_appenddata(cfg,data4{:});
        data6 = ft_appenddata(cfg,data6{:});

        %% Append all data into single data file with appropriate trialinfo
        cfg = [];
        cfg.keepsampleinfo = 'yes';
        data = ft_appenddata(cfg,data2, data4, data6);
        trialinfo = [data2.trialinfo; data4.trialinfo; data6.trialinfo];
        data.trialinfo = trialinfo;

        %% Sort appended data by chronological order based on global trial ID
        % sort by chronological order (globalTrial), push NaNs to the end
        [~, order] = sort(data.trialinfo(:,2), 'ascend', 'MissingPlacement', 'last');

        % keep only valid (non-NaN, >0) global IDs in this reordering pass
        valid = ~isnan(data.trialinfo(order,2)) & data.trialinfo(order,2) > 0;
        order = order(valid);

        % reorder the FieldTrip data *manually* (ft_selectdata won't reorder)
        data = reorder_trials_ft(data, order);

        % always refresh trialinfo from the reordered data
        trialinfo = data.trialinfo;

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
            save([savepathControlsFix, filesep, 'AOC_preStimFixation_', subjects{subj}, '_sternberg'], "preStimFixInfo");
        else
            savepathControlsFix = ['/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/preStimFixation/', subjects{subj}];
            mkdir(savepathControlsFix)
            save([savepathControlsFix, filesep, 'AOC_preStimFixation_', subjects{subj}, '_sternberg'], "preStimFixInfo");
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
        cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL' '-VEOGU' '-VEOGL' ...
            '-L-GAZE-X' '-L-GAZE-Y' '-L-AREA' ...
            '-R-GAZE-X' '-R-GAZE-Y' '-R-AREA'};
        dataEEG = ft_selectdata(cfg, data);

        %% Resegment data to avoid filter ringing
        % TFR data
        cfg = [];
        dataTFR = ft_selectdata(cfg, dataEEG);

        % EEG & ET data for Sternberg retention interval
        cfg = [];
        cfg.latency = [1 2]; % Time window for Sternberg task
        dataEEG = ft_selectdata(cfg, dataEEG); % EEG data
        dataet = ft_selectdata(cfg,dataet); % ET data
        dataet.trialinfo = trialinfo;

        %% Re-reference data to average or common reference
        cfg = [];
        cfg.reref   = 'yes';
        cfg.refchannel = 'all';
        dataTFR = ft_preprocessing(cfg,dataTFR);
        dataTFR.trialinfo = trialinfo;
        dataEEG = ft_preprocessing(cfg,dataEEG);
        dataEEG.trialinfo = trialinfo;

        %% Compute gaze metric data
        % WM load 2 gaze metrics average across trials
        saccades_l2 = sum(sacc_l2(:)) / sum(trl_l2(:));
        fixations_l2 = sum(fix_l2(:)) / sum(trl_l2(:));
        blinks_l2 = sum(blink_l2(:)) / sum(trl_l2(:));

        % WM load 4 gaze metrics average across trials
        saccades_l4 = sum(sacc_l4(:)) / sum(trl_l4(:));
        fixations_l4 = sum(fix_l4(:)) / sum(trl_l4(:));
        blinks_l4 = sum(blink_l4(:)) / sum(trl_l4(:));

        % WM load 6 gaze metrics average across trials
        saccades_l6 = sum(sacc_l6(:)) / sum(trl_l6(:));
        fixations_l6 = sum(fix_l6(:)) / sum(trl_l6(:));
        blinks_l6 = sum(blink_l6(:)) / sum(trl_l6(:));

        %% Save data
        clc
        if ispc == 1
            savepathEEG = strcat('W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\eeg\');
        else
            savepathEEG = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepathEEG)
        cd(savepathEEG)
        disp('SAVING dataEEG_sternberg...')
        save dataEEG_sternberg dataEEG
        disp('SAVING dataEEG_TFR_sternberg...')
        save dataEEG_TFR_sternberg dataTFR
        if ispc == 1
            savepathET = strcat('W:\Students\Arne\AOC\data\features\' , subjects{subj}, '\gaze\');
        else
            savepathET = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
        end
        mkdir(savepathET)
        cd(savepathET)
        disp('SAVING dataET_sternberg...')
        save dataET_sternberg dataet dataETlong dataETearly dataETlate
        disp('SAVING gaze_metrics_sternberg...')
        save gaze_metrics_sternberg ...
            saccades_l2 fixations_l2 blinks_l2 ...
            saccades_l4 fixations_l4 blinks_l4 ...
            saccades_l6 fixations_l6 blinks_l6
        clc
        if subj == length(subjects)
            disp(['Subject AOC ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. STERNBERG PREPROCESSING FINALIZED.'])
        else
            disp(['Subject AOC ' num2str(subjects{subj})  ' (' num2str(subj) '/' num2str(length(subjects)) ') done. Loading next subject...'])
        end
    else
        disp(['Subject ', num2str(subjects{subj}), ' already done. SKIPPING...'])
    end
end
toc;

%% CHECK GLOBAL TRIAL ID
% ids_all   = trialinfo(:,2);
% conds_all = trialinfo(:,1);
% ids_nonzero = ids_all(ids_all>0);
% 
% % 6. Visual diagnostic plots
% figure('Name','Trialinfo diagnostics','Color','w','Position',[200 200 1200 400]);
% subplot(1,3,1)
% scatter(1:numel(ids_all), ids_all, 10, conds_all, 'filled');
% xlabel('Trial index'); ylabel('Global ID');
% title('Chronological distribution (coloured by condition)');
% box on
% 
% subplot(1,3,2)
% histogram(ids_nonzero, 'BinWidth',5, 'FaceColor',[0.3 0.3 0.8]);
% xlabel('Global ID'); ylabel('Count');
% title('Global ID histogram'); box on
% 
% subplot(1,3,3)
% plot(ids_nonzero, conds_all(ids_all>0), '.', 'MarkerSize',10);
% xlabel('Global ID'); ylabel('Condition code');
% title('Condition vs. Global ID'); box on