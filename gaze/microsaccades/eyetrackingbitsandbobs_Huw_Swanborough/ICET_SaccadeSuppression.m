%% For calculating MSaccade suppression



ccc;
rans = load('V:\hswanb\01-InfCry\02-IC_EyeTracking\04-stimuli\randomisations\IC_Randomisations.mat');
ranIdx = rans.randomisation_cells{1}{2};
ranTrials = rans.randomisation_cells{1}{1};

% For assigning randomisation feats later in tables.
[stiRefIdx,b,c] = unique(ranIdx);
stimNames = ranTrials(b);




pDat = 'V:\hswanb\01-InfCry\02-IC_EyeTracking\01-data\tidy_data';

subj = dir(fullfile(pDat,'IC*'));
subj = {subj.name};



studyTracks = zeros(18,5501,length(subj));
studyTally = zeros(18,5501,length(subj));
studyPupils = zeros(18,5501,length(subj));

% Assign to cell arrays for tables
tracks = {};
tally = {}; 
tallyWin = []; %just 0.2-0.7s after stim.
contours = {};
intensity = {};
sourceLang = {};
subjID = {};
stimID = [];

trialTotals = [];

for iS = 1:length(subj)

    disp(subj{iS});
    pSubj = fullfile(pDat,subj{iS});
    datFile = dir(fullfile(pSubj,'IC*converted*'));
    load(fullfile(datFile.folder,datFile.name)); %loads subj_ET_data

    % iterate through the six blocs
    saccadeData = table();
    pupilData = table();
    for iF = 1:length(subj_ET_data)

        tempSacDat = subj_ET_data{iF}.saccadeDataTable;
        pupilData = subj_ET_data{iF}.pupilTab;
        saccadeData = [saccadeData;tempSacDat];

    end
    saccadeData.AdjTime = double(saccadeData.AdjTime);
    pupilData.AdjTime = double(pupilData.adjTime);


    onsets = unique(saccadeData.Event);
    onsets = onsets(2:end); %Get all the trials.


    %Iterate through stimuli, spike train matrices. of -2:4seconds x N trials.
    trialReference = {};
    subjectTimeTracks = [];
    stimuliOnsetIdx = []; %-2 to 4.
    subjectPupilTracks = [];
    ik = 1;
    for iO = 1:length(onsets)

        trialIdx = find(strcmp(saccadeData.Event, onsets{iO}));

        onsetTimes = saccadeData.AdjTime(trialIdx);
        trialTimeTracks =[]; %-2 to 4.
        pupilTimeTracks = [];
        for iT = 1:length(trialIdx)

            tIdx = trialIdx(iT); %

            tempVec = [tIdx-10:tIdx+10];

            % times = saccadeData.AdjTime(tempVec)-saccadeData.AdjTime(tIdx);

            %Incase the indexing goes outside of trial ranges.
            longK = any(tempVec>length(saccadeData.AdjTime));
            shortK = any(tempVec<0);

            kSkip = 0;

            if longK | shortK

                tempVec(tempVec>length(saccadeData.AdjTime)|tempVec<0) = [];
                times = saccadeData.AdjTime(tempVec)-saccadeData.AdjTime(tIdx);
                kSkip = 1; % to force skipping the next if statement due to range issues.
            else
                times = saccadeData.AdjTime(tempVec)-saccadeData.AdjTime(tIdx);

            end

            % times = saccadeData.AdjTime(tempVec)-saccadeData.AdjTime(tIdx);
            if min(times)<-2000 && max(times)>4000 && kSkip == 0;
                % If no values ? -2000 & +4000; then iterate through til it is
                % exceeded.

                k=0;
                itAd = 0;
                % Will iterate through til it samples times outside of the window.
                % Then use that to index sacade events.
                while k == 0;
                    itAd = itAd+1;
                    tempVec = [tIdx-(10+itAd):tIdx+(10+itAd)];
                    times = saccadeData.AdjTime(tempVec)-saccadeData.AdjTime(tIdx);
                    k= min(times)<-2000 && max(times)>4000;
                end
            end

            % Find Samples within the range.
            trialT = saccadeData.AdjTime(tIdx);
            eventTimes = saccadeData.AdjTime(tempVec)-trialT;
            eventIdx = eventTimes > -2000 & eventTimes<4000;
            saccTimes = eventTimes(eventIdx);
            %The indexes for pulling times and events from the list.
            eventIdx = tempVec(eventIdx);
            % Remove trial
            tCheck = find(eventIdx==tIdx);

            eventIdx(eventIdx==tCheck) = [];
            saccTimes(saccTimes==0) = [];
            %Saccade Times
            saccTimes = saccTimes+2000; %Readjust for marking saccade events above!
            saccVec = zeros(1,6000);
            saccVec(saccTimes) = 1;
            trialTimeTracks(iT,:) = saccVec;

            %% Pupil Dilation
            % Find all instances within the window of the trial for pupil
            % % dilation.
            % indexingTimeVec = [saccadeData.AdjTime(tIdx)-2000:saccadeData.AdjTime(tIdx)+4000];
            % indexingTimeVec(1) = []; %remove extar value cos im bad at indexing.
            % ttIdx = find(ismember(pupilData.adjTime, indexingTimeVec));
            %
            % %pull values
            % pupilValues = [pupilData.Left(ttIdx),pupilData.Right(ttIdx)];
            % pupilValues(pupilValues==0)=nan; % ET inputs pupil data as 0 if it loses it.
            % % Adjust times to 0 point.
            % truePupilTimes = (pupilData.adjTime(ttIdx) - (saccadeData.AdjTime(tIdx)-2000));
            % %to not have a 0 indexing point.
            % truePupilTimes(truePupilTimes==0)=[];
            % %Mean pupil Dilation.
            % pupilValues = mean(pupilValues,2,'omitnan');

            %
            % if sum(isnan(pupilValues)) == 6000;
            %     return
            % end
            % %
            % % pupVec = nan(1,6000);
            % % pupVec(truePupilTimes) = pupilValues;
            %
            %
            % timestester{ik} = truePupilTimes;
            %
            % finderror{ik} = [iO,iT];
            %
            % pupilTimeTracks(iT,:) = pupVec;
            %
            %
            %
            % testpupVec(ik) = sum(isnan(pupVec));
            %
            % testval(ik) = sum(isnan(pupilValues));
            % ik = ik+1;
        end

        trialOnset = repmat(onsets(iO),length(trialIdx),1);
        trialReference = [trialReference;trialOnset];
        stimuliOnsetIdx = [stimuliOnsetIdx;repmat(iO,length(trialIdx),1)];
        subjectTimeTracks = [subjectTimeTracks;trialTimeTracks];

        % subjectPupilTracks = [subjectPupilTracks;pupilTimeTracks];
    end

trialTotals = [trialTotals;length(stimuliOnsetIdx)];


    %% Bottom up saccade suppression paper.
    % Sampling frequency (for example, 1000 Hz)
    Fs = 1000;
    % Decay parameter (in s^-1)
    alpha = 50;
    t_max = 0.01;
    % Create the kernel
    t = linspace(0, t_max, t_max * Fs);
    origKernel = alpha^2 * t .* exp(-alpha * t);
    % Normalize the kernel to have a maximum value of 1
    origKernel = origKernel / max(origKernel);

    %Guassian Kernel
    sigma = 25; % Standard deviation in samples
    kSize = 2 * sigma;
    x = linspace(-3*sigma, 3*sigma, kSize);
    gKernel = exp(-x.^2 / (2 * sigma^2));
    gKernel = gKernel / sum(gKernel); % norm

    % Rolling frequency Window - Huw design
    % Ok, I think this is called "event desnsity estimation with
    % convolution"
    tWin = tukeywin(100);


    % Which Kernel To Use?
    kern = origKernel;

    % To plot the tracks in terms of frequency, multiply the value of the
    % avager kernel by hzLim, and then multiply the tracks by this.
    % gives approx frequency.
    hzLim = 1000/length(kern);
    kWeight = mean(kern);
    hzScaler = hzLim/kWeight;


    % % Detrend the Eye tracks
    % subjectTimeTracks = detrend(subjectTimeTracks');
    % subjectTimeTracks = subjectTimeTracks';

    % REmove first 500ms. Artefact issues.
    subjectTimeTracks = subjectTimeTracks(:,500:end);
    %Smoothing whole thing for basweline/
    overallTrack = fconv(sum(stimTrack,1),Fs,convWidth)*10*50/convWidth;
    % Subject Baseline
    baselineMSRate = mean(overallTrack(1:1500));

    %% Stimuli level Analysis - do for each exemplar, then group them after.
    % stimuliTracks = zeros(18,length(overallTrack));
    % stimuliTally = zeros(18,length(overallTrack));
    stimuliTally ={};
    stimuliWinTally = [];
    stimuliTracks = {};

    for iT = 1:max(stimuliOnsetIdx)

        stimIdx = stimuliOnsetIdx==iT;
        stimTrack = subjectTimeTracks(stimIdx,:);

        %Kernal method
        % smoothTrialTrack =  conv(sum(stimTrack,1), kern, 'same');

        convWidth = 50;
        smoothTrialTrack = fconv(sum(stimTrack,1),1000,convWidth)*10*50/convWidth;
        % Correct to baseline - maybe do this better?

        bsLineTrack = smoothTrialTrack-baselineMSRate;
        % mym = mym - nanmean(mym(smoothTrialTrack<0 & smoothTrialTrack>bc));
        % bsLineTrack = detrend(bsLineTrack');
        stimuliTracks(iT) = {bsLineTrack};
        trackMat(iT,:) = bsLineTrack;
        % Tallying
        stimuliTally(iT) = {sum(stimTrack,1)};
        stimuliWinTally(iT) = mean(sum(stimTrack(:,[1700:2200]),2));

    end

%Test Sijia baselinecorr.
bLine = mean(trackMat);

bLine2 = bLine- mean(bLine(10:1500));


    % subjectStimCells = stimNames(stimuliOnsetIdx);
    %Langs
    contCells = cellfun(@(x) x([1,2]),stimNames,'UniformOutput',false);
    intCells = cellfun(@(x) x(4),stimNames,'UniformOutput',false);
    sourceCells = cellfun(@(x) x(6),stimNames,'UniformOutput',false);


    % Assign to cell arrays for tables
    subjID = [subjID; repmat(subj(iS),length(contCells),1)];
    tracks = [tracks;stimuliTracks'];
    tally = [tally;stimuliTally'];
    tallyWin = [tallyWin;stimuliWinTally'];
    contours = [contours;contCells];
    intensity = [intensity;intCells];
    sourceLang = [sourceLang;sourceCells];
    stimID = [stimID,1:18];


    % Assign to Table cos it makes more sense.


    % % Assign pupilometry.
    % % To average on participant level or stimuli group level?
    % % Confusion.....
    % %both?
    % stimuliPupils = nan(18,6000);
    %
    % % normSubjPupils = normalize(subjectPupilTracks,2) ;
    %
    % % De-mean them.
    % normSubjPupils = subjectPupilTracks-mean(subjectPupilTracks,2,'omitnan') ;
    %
    % for iT = 1:max(stimuliOnsetIdx)
    %
    %     % Adjust pupil baseline.
    %
    %
    %
    %     stimIdx = stimuliOnsetIdx==iT;
    %     pupTracks = normSubjPupils(stimIdx,:);
    %
    %     plot(pupTracks(1,:))
    %     hold on
    %     % meanPupils = mean(pupTracks,'omitnan')
    %     % Pupilometry
    %     % stimuliPupils(iT,:) = subjectPupilTracks(stimIdx,:);
    %
    %
    % end
    %


    % %Saccades
    % studyTracks(:,:,iS) = stimuliTracks;
    % studyTally(:,:,iS) = stimuliTally;
    % %Pupil Dilation
    % % studyPupils(:,:,iS) = stimuliPupils;

end
stimID = stimID';
suppression_table = table(subjID,tracks,tally,tallyWin,contours,intensity,sourceLang,stimID);

%% Analysis

%% 50% crossing average
% Take the 50% point between suppression dip and prior peak for the average
% of all trials.

% When other conditions cross this is there suppression rate.

% Do also a return side of it.

%% Inhibition Onset - Language
% Onset is locked to 1500ms.

meanAll = mean(cat(1,suppression_table.tracks{:}));
meanAll = detrend(meanAll);

meanAll([1:100]) = meanAll(100);
meanAll([end-100:end]) = meanAll(end-100);

% to put in density frequency (approx).
% meanAll = meanAll*hzScaler;

[apVals, allPeaks] = findpeaks(meanAll);
[avVals, allValleys] = findpeaks(-meanAll);

% Find peak prior to lowest valley.

%take max valley at suppression
[allSupValue,allSupIdx] = max(avVals);
%accounts for peak/valley orders.
allSuppTime = allValleys(allSupIdx);

if allPeaks(allSupIdx)>allValleys(allSupIdx)
    allPriorPeak = apVals(allSupIdx-1);
    allPostPeak =  apVals(allSupIdx);
else
    allPriorPeak = apVals(allSupIdx);
    allPostPeak =  apVals(allSupIdx+1);
end

%Mean Crossing Points of onset and offset of suppresion
allOnsetXing = allPriorPeak-((allPriorPeak-allSupValue)/2);
allReleaseXing = allPriorPeak-((allPostPeak-allSupValue)/2);

% Mean Release Extent - The release extent seems to vary greatly.
%Time Values of crossing midpoint
[~, allOnsetTime] = min(abs(meanAll(1:allSuppTime)-allOnsetXing));
[~, allOffsetTime] = min(abs(meanAll(allSuppTime:allSuppTime+200)-allReleaseXing));
allOffsetTime = allOffsetTime+allSuppTime;

% Loop through conditions of offsets/onsets dependent on stimuli etc.
conts = unique(suppression_table.contours);
ints = unique(suppression_table.intensity);
nkIdx = nchoosek(1:3,2);

contsXints = combinations(conts,ints);
% contsXints = [conts(nkIdx(:,1)),ints(nkIdx(:,1))];
contsXints = cellfun(@(x,y) {x, y}, contsXints.conts(:), contsXints.ints(:), 'UniformOutput', false);


% BOOTSTRAPPING GOES IN THIS LOOP.
factors = [conts;ints;contsXints];
close all;
plot(meanAll,LineWidth=3)
hold on
for iF = 1:3%length(factors)

    fct = factors{iF};
    %Index
    fCIdx = contains(suppression_table.contours,fct);
    fIIdx = contains(suppression_table.intensity,fct);

    %Combine Indexes
    fIdx = (fCIdx+fIIdx)==max(fCIdx+fIIdx);

    % Mean tracks
    meanTracks = mean(cat(1,suppression_table.tracks{fIdx}));
     meanTracks = detrend(meanTracks);
    meanTracks([1:100]) = meanTracks(100);
    meanTracks([end-100:end]) = meanTracks(end-100);    % meanTracks = detrend(meanTracks);
    % meanTracks = meanTracks.*hzScaler;
    plot(meanTracks)
    %same peak process as above. but altered for when crossing the meean
    %crossing point of all trials.
    [pkValues, peakIdx] = findpeaks(meanTracks);
    [vyValues, valleyIdx] = findpeaks(-meanTracks);

    % Find peak prior to lowest valley.

    %take max valley at suppression
    [suppresionV,suppIdx] = max(vyValues);
    %accounts for peak/valley orders.
    suppresionT = valleyIdx(suppIdx);

    if peakIdx(suppIdx)>valleyIdx(suppIdx)
        prePeak = pkValues(suppIdx-1);
        % postPeak =  pkValues(suppIdx);
    else
        prePeak = pkValues(suppIdx);
        % postPeak =  pkValues(suppIdx+1);
    end


% %Mean Crossing Points of onset and offset of suppresion
% onsetFactMidPoint = prePeak-((prePeak-suppresionV)/2);
% % releaseFactMidPoint = postPeak-((postPeak-suppresionV)/2);
% 
% % Get time that it crosses the allOnset/Offset points.
% [~, factorOnsetXing] = min(abs(meanTracks(1500:suppresionT)-allOnsetXing));
% factorOnsetXing = factorOnsetXing+1500;
% [~, factorOffsetXing] = min(abs(meanTracks(suppresionT:suppresionT+200)-allReleaseXing));
% factorOffsetXing = factorOffsetXing+suppresionT;
% 
% onXTimes(iF) = factorOnsetXing;
% supTimes(iF) = suppresionT;
% releaseXTimes(iF) = factorOffsetXing;
clear factorOnsetXing suppresionT factorOffsetXing

end
legend;

% tempTab = table(onXTimes',supTimes',releaseXTimes');
% supTimeTotal = tempTab.Var3-tempTab.Var1;

return


% BootStrapping Notes

% until we get more subjects, we're going to have to bootstrap by doing
% rand perms of all tracks, and then using that for the SD, and using a one
% directional T-test for the crossing/release points.


% Otherwise, time to first Saccade could be used as a metric
% Or average tallies within a window of X.
% Area under a curve analysis?

% We could remove neutral too.
% Do STDEV across conditions to see which are weirdest?




%% Tally Analysis
% from 1500 to 2500

for iF = 1:length(factors)

    fct = factors{iF};
    %Index
    fCIdx = contains(suppression_table.contours,fct);
    fIIdx = contains(suppression_table.intensity,fct);

    %Combine Indexes
    fIdx = (fCIdx+fIIdx)==max(fCIdx+fIIdx);

    stimTallies = cat(1,suppression_table.tally{fIdx});

    extractSaccades  = stimTallies(:,1700:2000);
    winTally(iF) = mean(sum(extractSaccades,2));



    winSDV(iF) = std(sum(extractSaccades,2));
end

%% Time after stimuli average position of MSaccades.

for iF = 1:length(factors)
    fct = factors{iF};
    %Index
    fCIdx = contains(suppression_table.contours,fct);
    fIIdx = contains(suppression_table.intensity,fct);

    %Combine Indexes
    fIdx = (fCIdx+fIIdx)==max(fCIdx+fIIdx);
    stimTallies = cat(1,suppression_table.tally{fIdx});
    extractSaccades  = stimTallies(:,1550:1800);

    %Convert to Time units
    [r c] = size(extractSaccades);
    esOut = nan(r,c);

    for iR = 1:height(extractSaccades)
        timeLoc = find(extractSaccades(iR,:));
        esOut(iR,timeLoc) = timeLoc;


    end


    winTally(iF) = mean(esOut(:),'omitnan');

    winSDV(iF) = std(esOut(:),'omitnan');
end

% Hacky Raster Plot


rastY = [];
rastX = [];
for iS = 1:18

    tIdx = suppression_table.stimID == iS;

    trialTallies =cat(1,suppression_table.tally{tIdx});

    % rastTrack = [rastTrack;trialTallies];

    % for timeidx

    offset=1/35;


    tempX = [];
    tempY = [];
    for iR = 1:height(trialTallies)

        tx=find(trialTallies(iR,:));

        % tmp = zeros(1,5501);
        % tmp()

        tempX=[tempX,tx];

        tY= repmat(iS-(offset*iR),1,length(tx));
        tempY=[tempY,tY];

    end
rastY = [rastY,tempY];
rastX = [rastX,tempX];

end

scatter(rastX,rastY,".")
hold on

for iS = 1:18
plot(1:5501,repmat(iS,1,5501),LineStyle = '-.',LineWidth=2,Color=[0,0,0,1]);
end
plot(repmat(1500,1,19),0:18,LineStyle ='--',LineWidth=2,Color=[0,0,0,1])
plot(repmat(2500,1,19),0:18,LineStyle ='--',LineWidth=2,Color=[0,0,0,1])

colormap("bone")

a=gca();
a.YTick = [0.5:18.5]
a.YLabel = stimNames;
set(gca,"YTick",[0.5:18.5])
set(gca,"YLabel",stimNames)

ylabel(stimNames)







%% Stimuli ID


close all;
plot(meanAll,LineWidth=3)
hold on
for iF = 1:3%18

    %Combine Indexes
    fIdx = suppression_table.stimID == iF;

    % Mean tracks
    meanTracks = mean(cat(1,suppression_table.tracks{fIdx}));
    meanTracks([1:25]) = mean(meanTracks);
    meanTracks([end-25:end]) = mean(meanTracks);    % meanTracks = detrend(meanTracks);
    % meanTracks = meanTracks.*hzScaler;
    plot(meanTracks)
    %same peak process as above. but altered for when crossing the meean
    %crossing point of all trials.
    [pkValues, peakIdx] = findpeaks(meanTracks);
    [vyValues, valleyIdx] = findpeaks(-meanTracks);

    % Find peak prior to lowest valley.

    %take max valley at suppression
    [suppresionV,suppIdx] = max(vyValues);
    %accounts for peak/valley orders.
    suppresionT = valleyIdx(suppIdx);

    if peakIdx(suppIdx)>valleyIdx(suppIdx)
        prePeak = pkValues(suppIdx-1);
        postPeak =  pkValues(suppIdx);
    else
        prePeak = pkValues(suppIdx);
        postPeak =  pkValues(suppIdx+1);
    end


% %Mean Crossing Points of onset and offset of suppresion
% onsetFactMidPoint = prePeak-((prePeak-suppresionV)/2);
% releaseFactMidPoint = postPeak-((postPeak-suppresionV)/2);
% 
% % Get time that it crosses the allOnset/Offset points.
% [~, factorOnsetXing] = min(abs(meanTracks(1500:suppresionT)-allOnsetXing));
% factorOnsetXing = factorOnsetXing+1500;
% [~, factorOffsetXing] = min(abs(meanTracks(suppresionT:suppresionT+200)-allReleaseXing));
% factorOffsetXing = factorOffsetXing+suppresionT;
% 
% onXTimes(iF) = factorOnsetXing;
% supTimes(iF) = suppresionT;
% releaseXTimes(iF) = factorOffsetXing;
% clear factorOnsetXing suppresionT factorOffsetXing

end

legend;

%% Append MPS to table of each participants tallys. Throw that in R.

% saccadeTotal = cellfun(@sum, suppression_table.tallyWin);
saccadeTotal = suppression_table.tallyWin;


pStim = 'V:\hswanb\01-InfCry\02-IC_EyeTracking\04-stimuli\441_targets';
for iS = 1:length(stimNames)

    fStr = fullfile(pStim,stimNames{iS});
    output(iS) = MPS_analysis_roughout(fStr);
end



%Get demographs
% Load and append demographics data
pDemo = 'V:\hswanb\01-InfCry\02-IC_EyeTracking\01-data\questionnaires';
demFiles = dir(fullfile(pDemo,'*csv'));
dgraphs = fullfile(pDemo,demFiles.name);
dgraphData = readtable(dgraphs);

demTable = readtable('V:\hswanb\01-InfCry\02-IC_EyeTracking\02-analysis\02-behavioural\r_analyis\IC_ET_Behavioural.csv');



resTable = readtable('V:\hswanb\01-InfCry\02-IC_EyeTracking\02-analysis\02-behavioural\r_analyis\IC_ET_Behavioural.csv');


trials = unique(resTable.trials);

aversion = [];
salience = [];

for iS = 1:length(subj)

    sIdx= find(resTable.subjID==iS);
    subjTable = resTable(sIdx,:);

    for iT = 1:length(trials)
        tempTrial = trials{iT};

        trIdx = find(strcmp(subjTable.trials,tempTrial));

        avResp(iT) = mean(subjTable.aversive(trIdx));
        salResp(iT) = mean(subjTable.salient(trIdx));

    end
    

aversion = [aversion;avResp'];
salience = [salience;salResp'];

end



sex = {};
parental = {};
PANAS_Pos = [];
PANAS_Neg = [];
uncondPositive = [];
antAnnoyance = [];
contWilling = [];

for iS = 1:length(subj)

    % sIdx = find(strcmp(subj{iS},demTable.subjID))

    sIdx= find(demTable.subjID==iS);
    sIdx = sIdx(1);
    sSex = repmat(demTable.sex(sIdx),18,1);
    pSat = repmat(demTable.parental(sIdx),18,1);
    panasPos = repmat(demTable.PANAS_Pos(sIdx),18,1);
    panasNeg = repmat(demTable.PANAS_Neg(sIdx),18,1);
    ucPos = repmat(demTable.uncondPositive(sIdx),18,1);
    anAnn = repmat(demTable.antAnnoyance(sIdx),18,1);
    cWill = repmat(demTable.contWilling(sIdx),18,1);


    %demographs
    sex = [sex;sSex];
    parental = [parental;pSat];
    PANAS_Pos = [PANAS_Pos;panasPos];
    PANAS_Neg = [PANAS_Neg;panasNeg];
    uncondPositive = [uncondPositive;ucPos];
    antAnnoyance = [antAnnoyance;anAnn];
    contWilling = [contWilling;cWill];

end


tallyTable = suppression_table;
tallyTable.tracks = [];
tallyTable.tally = [];


stimMPS = repmat(output,1,length(subj));
stimMPS = stimMPS';
tallyTable = addvars(tallyTable,saccadeTotal,stimMPS,sex,parental,PANAS_Pos,PANAS_Neg,uncondPositive,antAnnoyance,contWilling,aversion,salience);

writetable(tallyTable,'V:\hswanb\01-InfCry\02-IC_EyeTracking\02-analysis\01-ET\R_analysis\Tally_Rates.csv');


nVec = strcmp(tallyTable.contours,'Ne');
fVec = strcmp(tallyTable.contours,'Fr');
gVec = strcmp(tallyTable.contours,'Ge');

[r, p] = corr(tallyTable.saccadeTotal(nVec),tallyTable.salience(nVec))
[r, p] = corr(tallyTable.saccadeTotal(fVec),tallyTable.salience(fVec))
[r, p] = corr(tallyTable.saccadeTotal(gVec),tallyTable.salience(gVec))




[r, p] = corr(tallyTable.aversion,tallyTable.salience)

[r, p] = corr(tallyTable.saccadeTotal,tallyTable.salience)

[r, p] = corr(tallyTable.saccadeTotal,tallyTable.aversion)
return





% For bootstrapping, see the method inside r_lmersampler.
% Do the same thing for behavioural and ET - do it for tallies, is it
% possible to do it for the tracks as well? Probably less so.....





% Correlate ratings and tallys per Sec

R  = corr2(output,winTally);





for iF = 1:18

    %Combine Indexes
    fIdx  = suppression_table.stimID == iF;


    stimTallies = cat(1,suppression_table.tally{fIdx});

    extractSaccades  = stimTallies(:,1700:2000);
    winTally(iF) = mean(sum(extractSaccades,2));



    winSDV(iF) = std(sum(extractSaccades,2));
end


newLabels = cellfun(@(x) x([1,2]),stimNames,'UniformOutput',false)

for ii = 1:18

contMPS(ii) = strcat(newLabels(ii),'-',num2str(output(ii)));
end


scatter(output,winTally);

text(output, winTally,contMPS, 'Vert','bottom', 'Horiz','left', 'FontSize',7);


R  = corr2(output,winTally);


[r, p] = corr(output',winTally')

std(tallyTable.saccadeTotal)

std(tallyTable.saccadeTotal(strcmp(tallyTable.contours,'Ge')))
std(tallyTable.saccadeTotal(strcmp(tallyTable.contours,'Fr')))
std(tallyTable.saccadeTotal(strcmp(tallyTable.contours,'Ne')))


%% Quick Descriptives of Surviving Trials


includedTrials = trialTotals/324;

mean(includedTrials);
std(includedTrials);





% 
% 
% 
% 
% 
% 
% 
% 
% 
% % adjust to be placed in frequency domain
% meanTracks = meanTracks*hzScaler;
% 
% meanAll = mean(meanTracks,1);
% freMean = mean(meanTracks(1:6,:));
% gerMean = mean(meanTracks(7:12,:));
% neuMean = mean(meanTracks(13:18,:));
% 
% [fpVals, frPeaks] = findpeaks(freMean);
% [gpVals, gePeaks] = findpeaks(gerMean);
% [npVals, nePeaks] = findpeaks(neuMean);
% 
% [fvVals, frValleys] = findpeaks(-freMean);
% [gvVals, geValleys] = findpeaks(-gerMean);
% [nvVals, neValleys] = findpeaks(-neuMean);
% 
% 
% 
% % Could angle be used as analysis?
% 
% [fSupp2,fsidx2] = max(fvVals)
% 
% 
% % time lock the last peak before suppression.
% [fSupp,allSupIdx] = max(fvVals);
% [gSupp,gsidx] = max(gvVals);
% [nSupp,nsidx] = max(nvVals);
% 
% 
% frOffset = frValleys(allSupIdx) - frPeaks(allSupIdx);
% geOffset = geValleys(gsidx) - gePeaks(gsidx);
% neOffset = neValleys(nsidx) - nePeaks(nsidx);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% % Bootstrapping? 
% % Take N% of trials/participants?
% 
% %% Tally Analysis
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Notes
% % interesting interaction between intensity and contour for speed of
% % suppression and release from suppression at both onset and stimuli
% % offset.
% 
% % Align stimuli with peak amplitude.
% 
% % pStim = 'V:\hswanb\01-InfCry\02-IC_EyeTracking\04-stimuli\441_targets';
% % 
% % files = dir(fullfile(pStim,'*.wav')); files = {files.name};
% % 
% % 
% % close all
% % for iF = 1:length(files)
% %     [wavs fs]= audioread(fullfile(pStim,files{iF}));
% %     env = envelope(wavs,fs*0.1,'rms');
% %     [x mIdx(iF)]= max(env);
% % 
% %     lvl(iF) = rms(wavs);
% % end
% 
% % %sample down to time.
% % peakTimes = round(mIdx/fs*1000)
% % % adjustment = (peakTimes-min(peakTimes))+1;
% % adjustment = peakTimes+1500;
% % adjustment = 1500;
% meanTracks = mean(studyTracks,3);
% meanTally = mean(studyTally,3);
% 
% % adjust to be placed in frequency domain
% meanTracks = meanTracks*hzScaler;
% meanTracks = meanTracks(:,1:end-500);
% holdTracks = meanTracks;
% %Align with peaks
% % 
% % clear meanTracks
% % meanTracks = nan(18,length(holdTracks)+1000);
% 
% 
% % 
% % 
% % 
% % for iR = 1:18
% % 
% % 
% % 
% %     tempVec = holdTracks(iR,adjustment(iR)-500:end);
% % 
% %     meanTracks(iR,:) =  [tempVec,nan(1,length(meanTracks)-length(tempVec))];
% % 
% % end
% 
% meanTracks = meanTracks(:,1000:end)
% 
% meanAll = mean(meanTracks,1);
% 
% fre = mean(meanTracks(1:6,:));
% ger = mean(meanTracks(7:12,:));
% neu = mean(meanTracks(13:18,:));
% 
% 
% 
% 
% figure()
% plot(meanAll,LineWidth=3)
% hold on
% plot(fre,LineWidth=1)
% plot(ger,LineWidth=1)
% plot(neu,LineWidth=1)
% 
% legend all french german neutral
% 
% 
% 
% 
% 
% figure()
% 
% hgh = mean(meanTracks([1,2,7,8,13,14],:));
% med = mean(meanTracks([5,6,11,12,17,18],:));
% low = mean(meanTracks([3,4,9,10,15,16],:));
% 
% plot(meanAll,LineWidth=3)
% hold on
% plot(hgh)
% plot(med)
% plot(low)
% 
% legend all high medium low
% 
% 
% 
% 
% 
% 
% hghFr = mean(meanTracks([1,2],:));% - min(mean(meanTracks([1,2],:)));
% medFr = mean(meanTracks([5,6],:));% - min(mean(meanTracks([5,6],:)));
% lowFr = mean(meanTracks([3,4],:));% - min(mean(meanTracks([3,4],:)));
% 
% 
% hghGe = mean(meanTracks([7,8],:));% - min(mean(meanTracks([7,8],:)))
% medge = mean(meanTracks([11,12],:));% - min(mean(meanTracks([11,12],:)));
% lowge = mean(meanTracks([9,10],:));% - min(mean(meanTracks([9,10],:)));
% 
% 
% figure()
% plot(meanAll,LineWidth=3)
% hold on
% plot(hghFr)
% plot(medFr)
% plot(lowFr)
% legend a hf mf lf
% 
% 
% figure()
% plot(meanAll,LineWidth=3)
% hold on
% 
% plot(hghGe)
% plot(medge)
% plot(lowge)
% 
% legend a hg mg lg
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Peaks and Valleys Analysis
% 
% % Align all conditionsby the last Peak before their min valley.
% 
% % 
% % [fpVals, frPeaks] = findpeaks(fre);
% % [gpVals, gePeaks] = findpeaks(ger);
% % [npVals, nePeaks] = findpeaks(neu);
% % 
% % [fvVals, frValleys] = findpeaks(-fre);
% % [gvVals, geValleys] = findpeaks(-ger);
% % [nvVals, neValleys] = findpeaks(-neu);
% % 
% % % time lock the last peak before suppression.
% % [fSupp,fsidx] = max(fvVals);
% % [gSupp,gsidx] = max(gvVals);
% % [nSupp,nsidx] = max(nvVals);
% % 
% % 
% % frOffset = frValleys(fsidx) - frPeaks(fsidx);
% % geOffset = geValleys(gsidx) - gePeaks(gsidx);
% % neOffset = neValleys(nsidx) - nePeaks(nsidx);
% % 
% % 
% % osFr = fre;
% % % osFr(1:24) = []
% % 
% % osGe = ger;
% % osGe(1:12) = [];
% % 
% % osNe = neu;
% % % osNe(1:26) = [];
% % 
% % 
% % figure()
% % plot(osFr)
% % hold on
% plot(osGe)
% plot(osNe)
% 
% legend french german neutral
% 
% 
% 
% %
% % sumTal = sum(studyTally,3);
% % sumTal = sumTal(:,500:end-2000);
% %
% imagesc(sum(sumTal,3))
% colormap('bone')
%
% a=1


% Derive baseline firing from summary and smoothing of ALL trials
% pre-detection.

%Sum and smooth for each stimuli presentation.

%
% plot(mean(subjectTimeTracks,1))
%
% % Set a baseline, pre-trial microsaccade rate prior to stimuli onset.
% Correct to baseline(?) I guess just subtract each subjects baseline
% average.



%permute through these for stats testing I guess?






%
% end































