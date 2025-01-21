%% Convert and pre-process eye-tracking data.

function  ICET_PreProcessingv2(pTop,reRun,pConfig)



warning('off')
%version Number
vN = 2;

% Parralelililililise
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(3);
end

%Paths
pDat = fullfile(pTop,'01-data');
pOut = fullfile(pTop,'01-data/tidy_data');

%% Function Code.

% Subj Configs
conFile = dir(fullfile(pConfig,'*subject_configs*')); conFile = conFile.name;
% loads subjConStruct - % subjConStruct.subjConf; -% subjConStruct.subjVersionCtrl;
load(fullfile(pConfig,conFile));

%Conv Checklist
if reRun == 0
    ppIdx = subjConStruct.subjConf.converted ~= 1;
    tmpSubj = subjConStruct.subjConf.subjID(ppIdx);
    confK = 0;
else
    tmpSubj = subjConStruct.subjConf.subjID;
    confK = 1; %for checking if to update configs.
end

% pull raw ET data from non-processed ET data.
rawFiles = dir(fullfile(pDat,'raw_et','*edf')); rawFiles = {rawFiles.name};
%Pull the unprocessed data from the list of nonconverted subjs

%to avoid processing already processed.
subjFiles = find(contains(rawFiles, tmpSubj, 'IgnoreCase', true));
subjFiles = rawFiles(subjFiles);

subjs = cellfun(@(x) x(1:6),subjFiles,'UniformOutput',false);
subjs = unique(subjs);

%% Extract data, collate all of it.
% There are more extractable things from this, already written, just not
% passed out.

for iS = 1:length(subjs)
    datFiles = find(contains(subjFiles, subjs(iS), 'IgnoreCase', true));
    datFiles = subjFiles(datFiles);

    clear subj_ET_data
    for iD = 1:length(datFiles)
        subj_ET_data{iD}  = process_ET_data(fullfile(pDat,'raw_et',datFiles{iD}));
    end

    pSubjOut = fullfile(pOut,subjs{iS})
    mkdir(pSubjOut)

    subjstr = fullfile(pSubjOut,[subjs{iS},'_ET_converted.mat']);
    save(subjstr,'subj_ET_data');

        disp([subjs{iS},' finished']);

end

end



%% The preprocessing function itself.
function [subj_ET_data]  = process_ET_data(dataPath)


etStruct = edfmex(dataPath);

%% Index Trial Onsets

%This is the system time of the first trial, subtract from all values.
% sysTime = etStruct.RECORDINGS(1).time;
et_events = etStruct.FEVENT;
et_data = etStruct.FSAMPLE;

%Save memory.
clear etStruct

% et_events.sttime = cell2mat({et_events.sttime})-sysTime;
msgIdx = find(contains({et_events(:).codestring},'MESSAGEEVENT'));

% I made silly strings and now they need this to be pulled.
%Temp Idx where msgIdx == Onset
tempTrIdx = find(~cellfun(@isempty, regexp({et_events(msgIdx).message}, '^Onset_([1-9]|[1-9][0-9])$')));
%index of MSG events - Trial Onsets.
trialIdx = msgIdx(tempTrIdx);

%Trials and time onset - Reference for later? (Maybe make all this a
%table).
onsets = {et_events(trialIdx).sttime};
onsetTimes = cell2mat({et_events(trialIdx).sttime});%timing
eventIdx = trialIdx; %whre it is in et_events struct.
% Extraction area for saccade rate increase decrease.
winStarts = onsetTimes-2000; %2 seconds before.
winEnds = onsetTimes+4000; % 4 seconds after.

%% Blink Removal and Thresholding
% Take only trials with >75%  with one eye. Exclude blinks +/-100ms

%pull all blinks
tempBlIdx = find(strcmp({et_events.codestring},'ENDBLINK')); %EndBlink sttime = start. Entime = end.

lBlTmp = cell2mat({et_events(tempBlIdx).eye}==0); %Left eye
leftBlIdx = tempBlIdx(lBlTmp);

rBlTmp = cell2mat({et_events(tempBlIdx).eye}==1);
rightBlIdx = tempBlIdx(rBlTmp);

% Construct vector of blink for both eyes. +/- 100ms buffers
%left
lbOn = cell2mat({et_events(leftBlIdx).sttime}); %Raw onsets for blink resps
%Right
rbOn = cell2mat({et_events(rightBlIdx).sttime});


%% Extract Blinks as function of Stimulus.
% Extract all trials with a blink within 500ms of stimuli onset.
% maybe extend to 1000ms due to rising intensity being later?
% Measure blink response time by stimuli.

%% Removable?

% blinkResponseWin = onsetTimes+1500; %gen vector of stim onsets
%
% %Loop this. Makes most sense.
% blinkRespTime = [];
% for iO = 1:length(onsetTimes)
%
%     tV = onsetTimes(iO):onsetTimes(iO)+1500; % Sampling vector for when a blink can occur to be considered a response.3
%
%     rCheck = ismember(tV,rbOn);
%     lCheck = ismember(tV,lbOn);
%
%     rResp = find(rCheck,1,'first'); %find the first value in that period - blink onset in response to trial.
%     lResp = find(lCheck,1,'first');
%
%     if ~isempty(rResp) & ~isempty(lResp)
%         %Time after stimuli onset - Mean of both eyes blink.
%         blinkRespTime(iO) = round(mean([lResp,rResp]));
%     else % If only one eye blinks, it's not considered a blink response.
%         blinkRespTime(iO) = nan;
%     end
% end
%
% % Blink Responses Linked to Trials
% % number 1:18, tie back to stimuli
% trials = {et_events(trialIdx).message}';
% blinkRespTime = blinkRespTime';
%
% % This table contains all trials "onset_N" - which match the randomisation
% % N and indexing. Can be used for blink response analysis.
% blinkTable = table(trials,blinkRespTime);

%% Remove Blinks from time series

% ALSO ADJUST ALL TIMING DATA TO BE RELATIVE TO BEGINNING OF FIRST BLOCK!
% ET DATA HAS WEIRD TIMING SYSTEM AND IT RESETS!

% Note - ET_event.timestamp is system time stamp, et_data.time is
% experiment time adjusted to 0 from trial1, block 1.

%try to trim excess data that is before experiment start.
% This is the start of every block.
startIdx = find(contains({et_events.codestring},'STARTSAMPLES'));
% Mark ends by the using the display_coords messages as indexing.
endIdx = find(contains({et_events.codestring},'ENDSAMPLES'));
endIdx(1) = []; %Just trim round oddities in ET data struct.
endIdx = [endIdx,length(et_events)]-1; %if you don't put in -1, nothing works.

%Times - This lines up exactly with the et_data.time.
startTime = double(cell2mat({et_events(startIdx).sttime}));
endTime = cell2mat({et_events(endIdx).sttime});
samplesLength = endTime-startTime;

% Blink Exclusions +/-
blBuf = 75; %blink buffer

lbOnAdj = lbOn-blBuf; %adjusted blink onsets for excluding time series
lbOffs = cell2mat({et_events(leftBlIdx).entime})+blBuf;

leftBlinkPeriod = arrayfun(@(x, y) x:y, lbOnAdj, lbOffs, 'UniformOutput', false);
leftBlinkPeriod = [leftBlinkPeriod{:}];

rbOnadj = rbOn-blBuf; %ajdusted as above.
rbOffs = cell2mat({et_events(rightBlIdx).entime})+blBuf;

rightBlinkPeriod = arrayfun(@(x, y) x:y, rbOn, rbOffs, 'UniformOutput', false);
rightBlinkPeriod = [rightBlinkPeriod{:}];

% Trim
rightBlinkPeriod(rightBlinkPeriod>endTime | rightBlinkPeriod<startTime) = [];
leftBlinkPeriod(leftBlinkPeriod>endTime | leftBlinkPeriod<startTime) = [];

leftBlinkPeriod(leftBlinkPeriod-startTime==0)=[];
rightBlinkPeriod(leftBlinkPeriod-startTime==0)=[];

% Create logic vector of "blink present"
%can be used as indexing or anythign else really.
blinkMarker = zeros(1,samplesLength);
blinkMarker((leftBlinkPeriod-startTime)) = 1;
blinkMarker((rightBlinkPeriod-startTime)) = 1;

%% Compile Data Table
%pupil x
tmp.time = et_data.time;
tmp.pXL = et_data.px(1,:);
tmp.pXR = et_data.px(2,:);
%pupil y
tmp.pYL = et_data.py(1,:);
tmp.pYR = et_data.py(2,:);
%pupil dialation.
tmp.pDL =  et_data.pa(1,:);
tmp.pDR =  et_data.pa(2,:);
% gaze
tmp.gXL = et_data.gx(1,:);
tmp.gXR = et_data.gx(2,:);
tmp.gYL = et_data.gy(1,:);
tmp.gYR = et_data.gy(2,:);
%pixels per degree
tmp.ppX = et_data.rx(1,:);
tmp.ppY = et_data.ry(1,:);
% eye velocity
tmp.velXL = et_data.rxvel(1,:);
tmp.velXR = et_data.rxvel(2,:);
tmp.velYL = et_data.ryvel(1,:);
tmp.velYR = et_data.ryvel(2,:);

%eye track Data.
tfnames = fieldnames(tmp);
eye_table = table();
for iT = 1:length(tfnames)

    var = double(tmp.(tfnames{iT}));
    eye_table = addvars(eye_table,var');
end

eye_table = renamevars(eye_table,1:width(eye_table),tfnames);
eye_table = addvars(eye_table,blinkMarker');
eye_table = renamevars(eye_table,18,'blink');


%% Microsaccade Detection.
%% From the FT Scripts.
% Engbert R, Kliegl R (2003) Microsaccades uncover the orientation of
% covert attention. Vision Res 43:1035-1045.
fsample = 1000;
kernel = [1 1 0 -1 -1].*(fsample/6); % this is equivalent to Engbert et al (2003) Vis Res, eqn. (1)
velthres = 6;
mindur  =  3; %duration a micsac must last in samples

% The Velocity Data As Extracted by EyeTracker.

% Remove Blinks - Be more subtle in this?
eye_table2 = eye_table(eye_table.blink==0,:);
velDataX = [eye_table2.pXL,eye_table2.pXR];
nsample = length(velDataX);

% Demean the velocities.
% see ft_Detect movement line 150...
velDataX = ft_preproc_polyremoval(velDataX', 0, 1, nsample);

% Detect the microsaccade.
n = size(kernel,2);
pad = ceil(n/2);
%padding is for the sake of the kernal operation I assume.
% Velocity of both eye at once.

%Padding for convolution to work.
dat = ft_preproc_padding(velDataX, 'localmean', pad);

% convolution. See Engbert et al (2003) Vis Res, eqn. (1)
vel = convn(dat,kernel,'same');
% cut the egdes
vel = ft_preproc_padding(vel, 'remove', pad);

%% microsaccade detection
% compute velocity thresholds as in Engbert et al (2003) Vis Res, eqn. (2)
medianstd = sqrt( median(vel.^2,2) - (median(vel,2)).^2 );

% Catch incase eye tracking fails? For some reason it seems to do this for
% some trials maybe?
if medianstd(1) == 0 | medianstd(2) == 0
    subj_ET_data.pupilTab = nan;
    subj_ET_data.saccadeDataTable = nan;
    return
end

% Engbert et al (2003) Vis Res, eqn. (3)
radius = velthres*medianstd;

% compute test criterion: ellipse equation
test = sum((vel./radius(:,ones(1,nsample))).^2,1);
sacsmp = find(test>1); % microsaccade's indexing
j = find(diff(sacsmp)==1);
j1 = [j; j+1];
com = intersect(j,j+1);
cut = ~ismember(j1,com);
sacidx = reshape(j1(cut),2,[]);


microsaccades = [];
for k=1:size(sacidx,2)
    duration = sacidx(1,k):sacidx(2,k);
    if size(duration,2) >= mindur
        % finding peak velocity by Pitagoras
        begtrl = sacsmp(duration(1,1));
        endtrl = sacsmp(duration(1,end));

        [peakvel, smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
        veltrl = sacsmp(duration(1,smptrl)); % peak velocity microsaccade sample -> important for spike conversion

        trlsmp = eye_table.time(eye_table.blink==0);
        begsample = trlsmp( begtrl); % begining microsaccade sample
        endsample = trlsmp( endtrl); % end microsaccade sample
        velsample = trlsmp( veltrl); % velocity peak microsaccade sample
        microsaccades(end+1,:) = [begsample endsample velsample];
    end
end

%adjust msac times
microsaccades = microsaccades-startTime;
%just in case of any 0 values.
microsaccades(microsaccades==0) = 1;

% Mark saccade onsets in vector.
saccMark = zeros(1,samplesLength);
saccMark(microsaccades(:,1)) = 1;
% Save into Eye data table.
eye_table = addvars(eye_table,saccMark','NewVariableNames','mSaccadeOn');
mSacTable = table(microsaccades(:,1),microsaccades(:,2),microsaccades(:,3),'VariableNames',{'Onset','Peak','Offset'});

%% Create a single Table of Event times and labels (i.e. trials and mSaccades)

timeAdjWinStarts = double(winStarts - startTime); % for indexing out which where blinked in.
timeAdjWinEnds = double(winEnds - startTime); % for indexing out which where blinked in.

%loop through the trials, gen indexes for the whole window
% Exclude if >N% of window is covered?

blinkIdx = find(blinkMarker);

%exclusion threshold for blink overlap in this window.
bThresh = 0.25;
adjOnsets = onsetTimes-startTime;
%  Marks a trial vector of exclusions due to blinks
% >10% of trial covered by blink, or >25% of 6 second window covered by
% blink.
excIdx = [];
for iT=1:length(timeAdjWinStarts)
    % First, find out if they blink through most of the actual stim, if
    % >10% exclude.
    stimTIdx = [adjOnsets(iT):adjOnsets(iT)+1000];
    stimTest = sum(ismember(stimTIdx,blinkIdx));
    if stimTest > 100
        % break or record?
        excIdx(iT) = 1;
        continue
    end
    timeIdxs = [timeAdjWinStarts(iT):timeAdjWinEnds(iT)];
    winTest = sum(ismember(timeIdxs,blinkIdx));

    if winTest > length(timeIdxs).*bThresh;
        excIdx(iT) = 1;
    else
        excIdx(iT) = 0;
    end
end

% Take +/-5s either side trials.
% Trim out blink moments.
% Take msacc, pupilom, event data.

extOnsets = adjOnsets(~excIdx)-3000;
extOffsets = adjOnsets(~excIdx)+5000;

trialExtractPeriod = arrayfun(@(x, y) x:y, extOnsets, extOffsets, 'UniformOutput', false);
trialExtractPeriod = [trialExtractPeriod{:}];
trialExtractPeriod = unique(trialExtractPeriod);

% This can lead to overlaps in the period... Is that an issue? It won't be
% overlapping in window, but it just gives more samples for extracting
% saccade rate baselines etc.

%An index for the data file
extractVec = zeros(1,samplesLength);
extractVec(trialExtractPeriod) = 1;
extractVec(~blinkMarker) =1;
extractIdx = find(extractVec);

%Really need to decide whether to standardise things - here it doesn't
%really matter as all msaccades should be included, just in case they fall
%inside an excluded trial due to blinks, uncommon.
%ok yes it'll all be adjusted from block onset.

incSacIdx = ismember(table2array(mSacTable(:,1)),extractIdx);
mSacEvents = table2array(mSacTable(incSacIdx,1));
mSacLabels = repmat({'MSaccade'},length(mSacEvents),1);

trialEvents = adjOnsets(~excIdx)';
trialLabels = {et_events(trialIdx(~excIdx)).message};

combEvents = [mSacEvents;trialEvents];
combLabels = [mSacLabels;trialLabels'];

saccadeDataTable = table(combEvents,combLabels,'VariableNames',{'AdjTime','Event'});

saccadeDataTable = sortrows(saccadeDataTable);

%% Pupilometry Extraction.

pupil_L = eye_table.pDL(extractIdx);
pupil_R = eye_table.pDR(extractIdx);
timePoints = eye_table.time(extractIdx)-startTime;
pupilTab = table(timePoints,pupil_L,pupil_R,'VariableNames',{'adjTime','Left','Right'});

%Extracted Data

subj_ET_data.pupilTab = pupilTab;
subj_ET_data.saccadeDataTable = saccadeDataTable;

end