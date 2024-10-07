%% Load, convert, preprocess ET data.
% Uses FieldTrip
function output = ICET_Conversion_v1(pTop,reRun,pConfig);
warning('off')
%version Number
vN = 2;

%Parralelililililise
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     parpool();
% end

%Paths
pDat = fullfile(pTop,'01-data');
pOut = fullfile(pTop,'01-data/tidy_data');

%% Function Code.

%Subj Configs
conFile = dir(fullfile(pConfig,'*subject_configs*')); conFile = conFile.name;
% loads subjConStruct - % subjConStruct.subjConf; -% subjConStruct.subjVersionCtrl;
load(fullfile(pConfig,conFile));

%Conv Checklist
if reRun == 0
    ppIdx = subjConStruct.subjConf.converted ~= 1;
    subj = subjConStruct.subjConf.subjID(ppIdx);
    confK = 0;
else
    subj = subjConStruct.subjConf.subjID;
    confK = 1; %for checking if to update configs.
end

%pull raw ET data from non-processed ET data.
rawFiles = dir(fullfile(pDat,'raw_et','*asc')); rawFiles = {rawFiles.name};
%Pull the unprocessed data from the list of nonconverted subjs

%% Tidy asc Files

for iS = 1:length(subj)

    pSubjRaw = fullfile(pOut,subj{iS},'raw_et');
    mkdir(pSubjRaw)
    pSubjConv = fullfile(pOut,subj{iS},'ET');
    mkdir(pSubjConv)

    %part files
    subjFiles = find(contains(rawFiles, subj{iS}, 'IgnoreCase', true));

    %Move files to subject directories.
    if isempty(subjFiles)
        eK = 1;
        subjFiles = dir(fullfile(pSubjRaw,'*asc'));
        subjFiles = {subjFiles.name};
        compFiles = dir(fullfile(pSubjConv,'*.mat'));
        if isempty(subjFiles)
            disp([subj{iS},' is Empty']); %Need to make a workaround for already converted trials?
            %probably no
            continue
        end
    else
        ek = 0;
    end

    disp(['Tidying ', subj{iS}])
    %parproc batch
    batch{iS}.subj = subj{iS};
    batch{iS}.outPath = fullfile(pSubjConv);
    batch{iS}.rawPath = pSubjRaw;

    for iF = 1:length(subjFiles)

        if eK==1 %just to control if they've already been moved but not converted.
            batch{iS}.cfg.dataset{iF} = fullfile(pSubjRaw,subjFiles{iF});
        else % if not moved.
            batch{iS}.cfg.dataset{iF} = fullfile(pSubjRaw,rawFiles{subjFiles(iF)});  %it'll take from here after it's all moved.
            movefile(fullfile(pDat,'raw_et/',rawFiles{subjFiles(iF)}),pSubjRaw)
        end
    end
end

%% Convert Data.
disp('Converting Data');

% Collate part data for parproc.
parfor iS = 1:length(batch)

% for iS = 1:length(batch)
    % Concatenates all trials into a single file already. Nice!

    % 
    event_eye{iS} = ft_read_event(batch{iS}.cfg.dataset);
% trialExtract(event_eye{iS})
    %Pull the trials out of event eye and feed them into pre-processing to
    %do this automatically - function it below.

    data_eye{iS} = ft_preprocessing(batch{iS}.cfg);
    
end

for iD = 1:length(data_eye)

    et_data = data_eye{iD};
    et_events = event_eye{iD};
    save(fullfile(batch{iD}.outPath,[batch{iD}.subj,'_conv.mat']),'et_data','et_events');
end

%Update config
if length(subj)>0
    for iS = 1:length(subj)
        sIdx = find(contains(subjConStruct.subjConf.subjID,subj{iS}));
        subjConStruct.subjConf.converted(sIdx) = 1;

        fileP = matlab.desktop.editor.getActiveFilename;
        [~ ,scrName] = fileparts(fileP);
        % Save version information for conversion.
        verStr = strcat(scrName,'_',string(datetime('today','Format','yyyy-MM-dd')));

        % Add version control data to subj config.
        subjConStruct.subjVersionCtrl{sIdx,1} = [subjConStruct.subjVersionCtrl{sIdx,1};verStr];
        disp(strcat('Converted ',subj{iS}));
    end

    save(fullfile(pConfig,'subject_configs.mat'),'subjConStruct'); %save the old version.
    output.Subj = subj;
    disp(strcat('Converted ',num2str(length(subj)),' subject data'))
end
warning('on')
end


% 
% function trialOut = trialExtract(event_data)
% 
% msgIdx = find(contains({event_data(:).type},'MSG'));
% 
% tempTrIdx = find(~cellfun(@isempty, regexp({event_data(msgIdx).value}, '^Onset_([1-9]|[1-9][0-9]) $')));
% %index of MSG events - Trial Onsets.
% trialIdx = msgIdx(tempTrIdx);
% 
% %Trials and time onset - Reference for later? (Maybe make all this a
% %table).
% trials = {et_events(trialIdx).value};
% trials(2,:) = {et_events(trialIdx).sample};%timing
% trials(3,:) = num2cell(trialIdx);
% 
% 
% 
% 
% 
% 
% 
% 
% end
% 





