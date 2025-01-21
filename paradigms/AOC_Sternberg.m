% #AOC Sternberg Arne
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2023b.
%
% Sternberg Task with 6 blocks of WM load 2, 4, and 6.

%% Initialize EEG and ET
% Calibrate ET (Tobii Pro Fusion)
disp('CALIBRATING ET...');
calibrateET;

if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;

    % Wait ten seconds to initialize EEG
    disp('INITIALIZING EEG... PLEASE WAIT 10 SECONDS')
    for i=1:10
        waitbar(i/10, 'INITIALIZING EEG');
        pause(1);
    end
    wbar = findall(0,'type','figure','tag','TMWWaitbar');
    delete(wbar)
    disp('EEG INITIALIZED!')
end

% Hide cursor on participant screen
HideCursor(whichScreen);

%% Define triggers
MATCH = 4; % trigger for probe stimulus
NO_MATCH = 5; % trigger for probe stimulus
TASK_START = 10; % trigger for ET cutting
FIXATION = 15; % trigger for fixation cross
PRESENTATION2 = 22; % trigger for letter presentation
PRESENTATION4 = 24; % trigger for letter presentation
PRESENTATION6 = 26; % trigger for letter presentation
DIGITOFF = 20; % trigger for change of digit to blank
BLOCK0 = 29; % trigger for start training block
BLOCK1 = 31; % trigger for start of block 1
BLOCK2 = 32; % trigger for start of block 2
BLOCK3 = 33; % trigger for start of block 3
BLOCK4 = 34; % trigger for start of block 4
BLOCK5 = 35; % trigger for start of block 5
BLOCK6 = 36; % trigger for start of block 6
ENDBLOCK0 = 39; % trigger for end training block
ENDBLOCK1 = 41; % trigger for end of block 1
ENDBLOCK2 = 42; % trigger for end of block 2
ENDBLOCK3 = 43; % trigger for end of block 3
ENDBLOCK4 = 44; % trigger for end of block 4
ENDBLOCK5 = 45; % trigger for end of block 5
ENDBLOCK6 = 46; % trigger for end of block 6
RETENTION2 = 52; % trigger for retention (setSize = 2)
RETENTION4 = 54; % trigger for retention (setSize = 4)
RETENTION6 = 56; % trigger for retention (setSize = 6)
RESP_YES = 87; % trigger for response yes (depends on changing key bindings)
RESP_NO = 88; % trigger for response no (depends on changing key bindings)
BAD_RESP = 89; % trigger for wrong keyboard input (any key apart from 'A', 'L' or 'Space')
TASK_END = 90; % trigger for ET cutting

%% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    exp.nTrials = 10;
else
    exp.nTrials = 50;            % 6 blocks x 50 trials = 300 trials
end
exp.setSizes = [2, 4, 6];     % Number of items presented on the screen

% Enable (=1) or disable (=0) screenshots
enableScreenshots = 0;

%% Set up equipment parameters
equipment.viewDist = 680;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
equipment.greyVal = .5;
equipment.blackVal = 0;
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

% Set up stimulus parameters Fixation
stimulus.fixationSize_dva = .5;         % Size of fixation cross in degress of visual angle
stimulus.fixationColor = 1;             % Color of fixation cross (1 = white)
stimulus.fixationLineWidth = 3;         % Line width of fixation cross

% Location
stimulus.regionHeight_dva = 7.3;         % Height of the region
stimulus.regionWidth_dva = 4;            % Width of the region
stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

% Set up color parameters
color.black = 0;
color.white = 1;

%% Set up texts
if TRAINING == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n On each trial, you will be shown a number of letters in a row. \n\n' ...
        'The sides will be filled with ''Xs''. These do not count! \n\n' ...
        'Example:  X N L + T P X \n\n' ...
        '\n\n' ...
        'Please ALWAYS fixate the center of the screen! \n\n' ...
        'Afterwards, you will be presented with one letter. \n\n' ...
        'Your task is to determine if this single letter was included previously. \n\n' ...
        'In this training session you''ll get feedback about the correctness of your responses. \n\n' ...
        '\n\n' ...
        'Press any key to continue.'];
else
    if BLOCK == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['On each trial, you will be shown a number of letters in a row. \n\n' ...
            'The sides will be filled with ''Xs''. These do not count! \n\n' ...
            'Example:  X D S + M T X \n\n' ...
            '\n\n' ...
            'Please ALWAYS fixate the center of the screen! \n\n' ...
            'Afterwards, you will be presented with one letter. \n\n' ...
            'Your task is to determine if this single letter was included previously. \n\n' ...
            'There will be no feedback (e.g., ''Correct!'') anymore. \n\n' ...
            'Press any key to continue.'];
    else
        loadingText = 'Loading actual task...';
        startExperimentText = ['Block ' num2str(BLOCK) ' / 6 \n\n' ...
            'Press any key to continue.'];
    end
end

performanceBonusText = ['In the following task there is a performance bonus! \n\n' ...
    'Try to be as accurate as possible. \n\n \n\n' ...
    'Press any key to continue.'];

startBlockText = 'Press any key to begin the next block.';

%% Set up temporal parameters (in seconds)
timing.letterPresentation = 0.2;            % Duration of letter presentation (200ms)
timing.rest = 2;                            % Duration of blank resting interval (2000ms)
timing.retentionInterval = 2.8;             % Duration of blank retention interval (2800ms)

%% Shuffle rng for random elements
rng('default');
rng('shuffle');                     % Use MATLAB twister for rng

%% Set up Psychtoolbox Pipeline
AssertOpenGL;

% Imaging set up
screenID = whichScreen;
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference', 'SkipSyncTests', 0); % For linux (can be 0)

% Set verbosity to disallow CW output
Screen('Preference','Verbosity', 0);

% Window setup
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
exp.runPriority = MaxPriority(ptbWindow);

% Set font and font size for instructions and stimuli
Screen('TextFont', ptbWindow, 'Liberation Mono');
Screen('TextSize', ptbWindow, 20);

global psych_default_colormode;                     % Sets colormode to be unclamped 0-1 range.
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text
ptb_drawformattedtext_disableClipping = 1;

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',color.black);
Screen('Flip',ptbWindow);

%% Retrieve response keys
KeyCodeA = KbName('A');                     % Retrieve key code for the A button
KeyCodeL = KbName('L'); % KbName('/?');     % Retrieve key code for the L button
spaceKeyCode = KbName('Space');             % Retrieve key code for spacebar

% Assign response keys
if mod(subject.ID,2) == 0 % Use subject ID for assignment to ensure counterbalancing of button presses
    YesIsL = true;        % L is YES, A is NO
    responseInstructionText = ['If you think the letter was included previously, press L. \n\n' ...
        'If you think the letter was not included previously, press A. \n\n'...
        ' \n\n' ...
        ' \n\n' ...
        'Press any key to continue.'];
elseif mod(subject.ID,2) == 1
    YesIsL = false;       % L is NO, A is YES
    responseInstructionText = ['If you think the letter was included previously, press A. \n\n' ...
        'If you think the letter was not included previously, press L. \n\n'...
        ' \n\n' ...
        ' \n\n' ...
        'Press any key to continue.'];
end

%% Calculate equipment parameters
equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
equipment.ppd = equipment.ppm*equipment.mpd;    % Pixels per degree

% Fix coordiantes for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];
fixPos = [screenCentreX, screenCentreY];

%% Define stimulus letter pool
consonants = char(setdiff('A':'Z', 'AEIOUXQ')); % Exclude vowels and 'X' and 'Q'

%% Create data structure for preallocating data
data = struct;
data.stimuli{1, exp.nTrials} = NaN;
data.trialSetSize(1, exp.nTrials) = NaN;
data.probe(1, exp.nTrials) = NaN;
data.match(1, exp.nTrials) = NaN;
data.responses(1, exp.nTrials) = NaN;
data.correct(1, exp.nTrials) = NaN;
data.fixation(1, exp.nTrials) = NaN;
data.reactionTime(1:exp.nTrials) = NaN;
data.trlDuration(1:exp.nTrials) = NaN;
count5trials = NaN;

%% Fixate randomized setSizes for each block
setS2 = ones(1, 12)*exp.setSizes(1);
setS4 = ones(1, 12)*exp.setSizes(2);
setS6 = ones(1, 12)*exp.setSizes(3);
chance = randsample(1:3, 1);
if chance == 1
    extraNums = [exp.setSizes(1) exp.setSizes(1)];
elseif chance == 2
    extraNums = [exp.setSizes(2) exp.setSizes(2)];
elseif chance == 3
    extraNums = [exp.setSizes(3) exp.setSizes(3)];
end
setALL = [setS2, setS4, setS6, extraNums];
for trialSetSizes = 1:exp.nTrials
    data.trialSetSize(trialSetSizes) = randsample(setALL, 1);
end
% Count occurences of set sizes in trials
data.setSizeOccurences(1) = numel(find(data.trialSetSize == exp.setSizes(1)));
data.setSizeOccurences(2) = numel(find(data.trialSetSize == exp.setSizes(2)));
data.setSizeOccurences(3) = numel(find(data.trialSetSize == exp.setSizes(3)));


%% Display info texts
clc
if TRAINING == 0
    % Show performance bonus incentive text
    DrawFormattedText(ptbWindow,performanceBonusText,'center','center',color.black);
    Screen('Flip',ptbWindow);
    disp('Participant is reading the performance bonus text');
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end
end

% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.black);
startExperimentTime = Screen('Flip',ptbWindow);
disp('Participant is reading the instructions');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end
screenshot('AOC_Sternberg_screenshot_startExperimentText.png', ptbWindow, enableScreenshots);

% Show response instruction text
DrawFormattedText(ptbWindow,responseInstructionText,'center','center',color.black);
Screen('Flip',ptbWindow);
disp('Participant is reading the response instructions');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end
screenshot('AOC_Sternberg_screenshot_responseInstructionText.png', ptbWindow, enableScreenshots);
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
Screen('Flip',ptbWindow);

% Send triggers for start of task
if TRAINING == 1
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
else
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
    sendtrigger(TASK_START,port,SITE,stayup);
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = BLOCK1;
elseif BLOCK == 2
    TRIGGER = BLOCK2;
elseif BLOCK == 3
    TRIGGER = BLOCK3;
elseif BLOCK == 4
    TRIGGER = BLOCK4;
elseif BLOCK == 5
    TRIGGER = BLOCK5;
elseif BLOCK == 6
    TRIGGER = BLOCK6;
else
    TRIGGER = BLOCK0;
end

if TRAINING == 1
    disp('Start of Block 0 (Training)');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end

% Experiment prep
HideCursor(whichScreen); % Make sure to hide cursor from participant screen
timing.startTime = datestr(now, 'dd/mm/yy-HH:MM:SS'); % Measure duration

%% Experiment Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trl = 1:exp.nTrials
    tic;
    % Randomize letterSequence
    lettersRand = randperm(length(consonants));
    sequenceIdx = lettersRand(1:data.trialSetSize(trl)); % Pick # of letters (# up to length of trialSetSize)
    clear sequenceLetters
    for numLoc = 1:data.trialSetSize(trl)
        sequenceLetters(numLoc) = consonants(sequenceIdx(numLoc));
    end
    % Save sequence of letters of this trial in data
    data.stimuli{trl} = sequenceLetters;

    %% Central fixation interval (jittered 500 - 1500ms)
    Screen('DrawLines', ptbWindow, fixCoords, stimulus.fixationLineWidth, stimulus.fixationColor, fixPos, 2); % Draw fixation cross
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    if TRAINING == 1
        Eyelink('Message', num2str(FIXATION));
        Eyelink('command', 'record_status_message "FIXATION"');
    else
        Eyelink('Message', num2str(FIXATION));
        Eyelink('command', 'record_status_message "FIXATION"');
        sendtrigger(FIXATION,port,SITE,stayup);
    end
    screenshot('AOC_Sternberg_screenshot_fixcross.png', ptbWindow, enableScreenshots);
    timing.cfi(trl) = (randsample(500:1500, 1))/1000; % Duration of the jittered inter-trial interval
    WaitSecs(timing.cfi(trl));

    %% Check fixation just before stimulus presentation
    noFixation = checkFixation(screenWidth, screenHeight, screenCentreX, screenCentreY);

    %% Presentation of stimuli (200ms)
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 50);
    % Define stimulus
    if data.trialSetSize(trl) == exp.setSizes(1)
        stimulusText = ['X ', 'X ', num2str(sequenceLetters(1)), '   ', ...
            num2str(sequenceLetters(2)), ' X', ' X'];
    elseif data.trialSetSize(trl) == exp.setSizes(2)
        stimulusText = ['X ', num2str(sequenceLetters(1)), ' ', num2str(sequenceLetters(2)), '   ', ...
            num2str(sequenceLetters(3)), ' ', num2str(sequenceLetters(4)), ' X'];
    elseif data.trialSetSize(trl) == exp.setSizes(3)
        stimulusText = [num2str(sequenceLetters(1)), ' ', num2str(sequenceLetters(2)), ' ', ...
            num2str(sequenceLetters(3)), '   ', num2str(sequenceLetters(4)), ' ', ...
            num2str(sequenceLetters(5)), ' ', num2str(sequenceLetters(6))];
    end
    data.stimulusLetters(trl) = {sequenceLetters(1:data.trialSetSize(trl))};
    data.stimulusText(trl) = {stimulusText};

    % Present stimuli
    DrawFormattedText(ptbWindow, stimulusText, 'center', 'center', color.white); % Draw stimuli in white
    Screen('DrawLines', ptbWindow, fixCoords, stimulus.fixationLineWidth, stimulus.fixationColor, fixPos, 2); % Draw fixation cross
    Screen('DrawDots', ptbWindow, backPos, backDiameter, backColor, [], 1);
    Screen('DrawDots', ptbWindow, stimPos, stimDiameter, stimColor, [], 1);
    Screen('Flip', ptbWindow);
    screenshot(sprintf('AOC_Sternberg_screenshot_WMload%d_trl%d.png', data.trialSetSize(trl), trl), ptbWindow, enableScreenshots);

    % Send triggers for Presentation
    if data.trialSetSize(trl) == exp.setSizes(1)
        TRIGGER = PRESENTATION2;
    elseif data.trialSetSize(trl) == exp.setSizes(2)
        TRIGGER = PRESENTATION4;
    elseif data.trialSetSize(trl) == exp.setSizes(3)
        TRIGGER = PRESENTATION6;
    end

    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMULUS"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMULUS"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.letterPresentation);

    % Return size of text to default
    Screen('TextSize', ptbWindow, 20);

    if TRAINING == 1
        Eyelink('Message', num2str(DIGITOFF));
        Eyelink('command', 'record_status_message "DIGITOFF"');
    else
        Eyelink('Message', num2str(DIGITOFF));
        Eyelink('command', 'record_status_message "DIGITOFF"');
        sendtrigger(DIGITOFF,port,SITE,stayup);
    end

    %% Retention interval (2800ms)
    Screen('DrawLines', ptbWindow, fixCoords, stimulus.fixationLineWidth, stimulus.fixationColor, fixPos, 2); % Draw fixation cross
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);

    if data.trialSetSize(trl) == exp.setSizes(1)
        TRIGGER = RETENTION2;
    elseif data.trialSetSize(trl) == exp.setSizes(2)
        TRIGGER = RETENTION4;
    elseif data.trialSetSize(trl) == exp.setSizes(3)
        TRIGGER = RETENTION6;
    end

    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "RETENTION"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "RETENTION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.retentionInterval);

    %% Randomize matching between letters in sequence, present probe stimulus and draw (probeLetter)
    chance = randsample(1:2, 1);
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 50);
    if chance == 1
        % Pick random matching probe stimulus from letterSequence
        probeLetter = randsample(sequenceLetters, 1);
        % Draw probe stimulus
        DrawFormattedText(ptbWindow,[num2str(probeLetter)],'center','center',color.white);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
        Screen('Flip', ptbWindow);
        match = 1;
        TRIGGER = MATCH;
    else
        % Pick random NON-matching probe stimulus from letters
        tmpAlphabet = consonants;
        for removeIdx = 1:data.trialSetSize(trl)
            tmpAlphabet = erase(tmpAlphabet, sequenceLetters(removeIdx));
        end
        sequenceLetters_excluded = tmpAlphabet;
        probeLetter = randsample(sequenceLetters_excluded, 1);
        % Draw probe stimulus
        DrawFormattedText(ptbWindow,[num2str(probeLetter)],'center','center',color.white);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
        Screen('Flip', ptbWindow);
        match = 0;
        TRIGGER = NO_MATCH;
    end
    probePresentationTime = GetSecs;
    screenshot(sprintf('AOC_Sternberg_screenshot_probe_trl%d.png', trl), ptbWindow, enableScreenshots);

    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PROBE STIMULUS"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PROBE STIMULUS"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end

    % Return size of text to 20 pts
    Screen('TextSize', ptbWindow, 20);

    % Save probe letter
    data.probe(trl) = probeLetter;
    % Save match/no match
    data.match(trl) = match;

    %% Get response (max 2000ms)
    getResponse = true;
    badResponseFlag = false;
    maxResponseTime = GetSecs + 2;
    responseTime = NaN;
    while getResponse
        [time,keyCode] = KbWait(-1, 2, maxResponseTime); % Wait 2 seconds for response, continue afterwards if there is no input
        whichKey = find(keyCode);

        % Subject pressed button A or L
        if ~isempty(whichKey)
            responseTime = GetSecs;
            if whichKey == KeyCodeA || whichKey == KeyCodeL
                getResponse = false;
                data.responses(trl) = whichKey;

                % Send triggers
                if whichKey == KeyCodeA && YesIsL == true
                    TRIGGER = RESP_NO;
                elseif whichKey == KeyCodeA && YesIsL == false
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL && YesIsL == true
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL && YesIsL == false
                    TRIGGER = RESP_NO;
                end

                if TRAINING == 1
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "RESPONSE"');
                else
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "RESPONSE"');
                    sendtrigger(TRIGGER,port,SITE,stayup)
                end

            else
                % Subject pressed other button than A or L
                responseTime = probePresentationTime + 4; % +4 to check; will be set to NaN
                TRIGGER = BAD_RESP;
                if TRAINING == 1
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "BAD RESPONSE"');
                else
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "BAD RESPONSE"');
                    sendtrigger(TRIGGER,port,SITE,stayup)
                end
                badResponseFlag = true;
            end

            % No input by participant
        elseif isempty(whichKey)
            data.responses(trl) = 0;
            responseTime = probePresentationTime + 5; % +5 to check; will be set to NaN
        end
        if ~isempty(whichKey)
            if time < maxResponseTime
                WaitSecs(maxResponseTime - time);
            end
        end
        if time > 1
            getResponse = false;
        end
    end

    % Save reaction time for each trial
    data.reactionTime(trl) = responseTime - probePresentationTime;

    %% Check if response was correct
    if YesIsL == 1       % L is YES, A is NO
        if data.responses(trl) == 0
            data.correct(trl) = 0;
        elseif match == 1     % Matched trial
            data.correct(trl) = data.responses(trl) == KeyCodeL;
        elseif match == 0     % Unmatched trial
            data.correct(trl) = data.responses(trl) == KeyCodeA;
        end
    elseif YesIsL == 0   % L is NO, A is YES
        if data.responses(trl) == 0
            data.correct(trl) = 0;
        elseif match == 1     % Matched trial
            data.correct(trl) = data.responses(trl) == KeyCodeA;
        elseif match == 0     % Unmatched trial
            data.correct(trl) = data.responses(trl) == KeyCodeL;
        end
    end

    %% Feedback for training block and CW output
    % CW Feedback
    if data.correct(trl) == 1
        feedbackText = 'Correct!  ';
    elseif data.correct(trl) == 0 && data.responses(trl) == 0
        feedbackText = 'NO RESPONSE';
    elseif data.correct(trl) == 0 && badResponseFlag == false
        feedbackText = 'Incorrect!';
    elseif data.correct(trl) == 0 && badResponseFlag == true
        feedbackText = ['Wrong button! \n\n' ...
            'Use only A or L.'];
    end

    % Give feedback in training block
    if TRAINING == 1
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        screenshot(sprintf('AOC_Sternberg_screenshot_feedback_trl%d.png', trl), ptbWindow, enableScreenshots);
        WaitSecs(2);
        % Give feedback for no response (too slow)
    elseif TRAINING == 0 && data.correct(trl) == 0 && data.responses(trl) == 0
        feedbackText = 'TOO SLOW! ';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
    end

    %% Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 60%
    responsesLastTrials = 0;
    if trl >= 10
        responsesLastTrials = data.correct(trl-9 : trl);
        percentLastTrialsCorrect = sum(responsesLastTrials)*10;
        if percentLastTrialsCorrect < 60 && count5trials <= trl-5
            count5trials = trl;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n ' ...
                '\n\n Please stay focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %. [' num2str(responsesLastTrials) ']']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.white);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            screenshot(sprintf('AOC_Sternberg_screenshot_accuracy_reminder_trl%d.png', trl), ptbWindow, enableScreenshots);
            WaitSecs(3);
        end
    end

    %% Fixation reminder
    % noFixation = 0;
    if noFixation > 0
        Screen('TextSize', ptbWindow, 30);
        fixText = 'ALWAYS LOOK AT THE CENTER OF THE SCREEN!';
        DrawFormattedText(ptbWindow, fixText, 'center', 'center', color.white);
        Screen('DrawDots', ptbWindow, backPos, backDiameter, backColor, [], 1);
        Screen('Flip', ptbWindow);
        disp('FIXATION REMINDER')
        screenshot(sprintf('AOC_Sternberg_screenshot_fixation_reminder_trl%d.png', trl), ptbWindow, enableScreenshots);
        WaitSecs(3);
        data.fixation(trl) = 0;
        Screen('TextSize', ptbWindow, 20);
    else
        data.fixation(trl) = 1;
    end

    %% Trial Info CW output
    overall_accuracy = round((sum(data.correct(1:trl))/trl)*100);
    reactionTime = num2str(round(data.reactionTime(trl), 2), '%.2f');
    if trl < 10
        disp(['Response to Trial  ' num2str(trl) '/' num2str(exp.nTrials) ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (WM load ' num2str(data.trialSetSize(trl)) ' | Acc: ' num2str(overall_accuracy) '% | RT: ' reactionTime 's)']);
    else
        disp(['Response to Trial ' num2str(trl) '/' num2str(exp.nTrials) ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (WM load ' num2str(data.trialSetSize(trl)) ' | Acc: ' num2str(overall_accuracy) '% | RT: ' reactionTime 's)']);
    end

    %% Blank screen for resting interval (2000ms)
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    screenshot('AOC_Sternberg_screenshot_resting_blank.png', ptbWindow, enableScreenshots);
    WaitSecs(timing.rest);

    % Save trial duration in seconds
    data.trlDuration(trl) = toc;
end

%% Send triggers for end of block and output
if BLOCK == 1
    TRIGGER = ENDBLOCK1;
elseif BLOCK == 2
    TRIGGER = ENDBLOCK2;
elseif BLOCK == 3
    TRIGGER = ENDBLOCK3;
elseif BLOCK == 4
    TRIGGER = ENDBLOCK4;
elseif BLOCK == 5
    TRIGGER = ENDBLOCK5;
elseif BLOCK == 6
    TRIGGER = ENDBLOCK6;
else
    TRIGGER = ENDBLOCK0;
end

disp(['End of Block ' num2str(BLOCK)]);

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end

% Send triggers for end of task (ET cutting)
if TRAINING == 1
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
else
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
    sendtrigger(TASK_END,port,SITE,stayup);
end

%% Record block duration
timing.endTime = datestr(now, 'dd/mm/yy-HH:MM:SS');
% Convert to datetime objects
startTime = datetime(timing.startTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
endTime = datetime(timing.endTime, 'InputFormat', 'dd/MM/yy-HH:mm:ss');
% Calculate block duration in seconds
timing.duration = seconds(endTime - startTime);

%% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)
if TRAINING == 1
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_training.mat'];
else
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_task.mat'];
end

% Save data structure
saves = struct;
saves.data = data;
saves.KeyCodeA = KeyCodeA;
saves.KeyCodeL = KeyCodeL;
saves.KeyBindingsYesIsL = YesIsL;
saves.experiment = experiment;
saves.screen.screenWidth = screenWidth;
saves.screen.screenHeight = screenHeight;
saves.screen.screenCentreX = screenCentreX;
saves.screen.screenCentreY = screenCentreY;
saves.startBlockText = startBlockText;
saves.startExperimentText = startExperimentText;
saves.subjectID = subjectID;
saves.subject = subject;
saves.text = text;
saves.timing = timing;
saves.flipInterval = flipInterval;

%% Save triggers
trigger = struct;
trigger.MATCH = MATCH;
trigger.NO_MATCH = NO_MATCH;
trigger.FIXATION = FIXATION;
trigger.FIX_RETENTION = FIX_RETENTION;
trigger.TASK_START = TASK_START;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.PRESENTATION4 = PRESENTATION4;
trigger.PRESENTATION6 = PRESENTATION6;
trigger.DIGITOFF = DIGITOFF;
trigger.BLOCK0 = BLOCK0;
trigger.BLOCK1 = BLOCK1;
trigger.BLOCK2 = BLOCK2;
trigger.BLOCK3 = BLOCK3;
trigger.BLOCK4 = BLOCK4;
trigger.BLOCK5 = BLOCK5;
trigger.BLOCK6 = BLOCK6;
trigger.ENDBLOCK0 = ENDBLOCK0;
trigger.ENDBLOCK1 = ENDBLOCK1;
trigger.ENDBLOCK2 = ENDBLOCK2;
trigger.ENDBLOCK3 = ENDBLOCK3;
trigger.ENDBLOCK4 = ENDBLOCK4;
trigger.ENDBLOCK5 = ENDBLOCK5;
trigger.ENDBLOCK6 = ENDBLOCK6;
trigger.RETENTION2 = RETENTION2;
trigger.RETENTION4 = RETENTION4;
trigger.RETENTION6 = RETENTION6;
trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;
trigger.BAD_RESP = BAD_RESP;
trigger.TASK_END = TASK_END;

% Stop and close EEG and ET recordings
disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger');
closeEEGandET;

try
    PsychPortAudio('Close');
catch
end

%% Compute accuracy and report after each block (no additional cash for training task)
if BLOCK == 0
    totalCorrect = sum(data.correct);
    totalTrials = trl;
    percentTotalCorrect = totalCorrect / totalTrials * 100;

    feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) '%. '];

    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.black);
    disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif BLOCK == 6
    totalCorrect = sum(data.correct);
    totalTrials = trl;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.0125;

    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional CHF ' num2str(amountCHFextra(BLOCK)) '.'];

    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.black);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif BLOCK > 0
    totalCorrect = sum(data.correct);
    totalTrials = trl;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.0125;

    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional CHF ' num2str(amountCHFextra(BLOCK)) '.' ...
        '\n\n Keep it up!'];

    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.black);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
end

% Show break instruction text
if TRAINING == 1
    if percentTotalCorrect >= 60
        breakInstructionText = 'Well done! \n\n Press any key to start the actual task.';
    else
        breakInstructionText = ['Score too low! ' num2str(percentTotalCorrect) ' % correct. ' ...
            '\n\n Press any key to repeat the training task.'];
    end
elseif BLOCK == 6
    breakInstructionText = ['End of the Task! ' ...
        '\n\n Press any key to view your stats.'];
else
    breakInstructionText = ['Break! Rest for a while... ' ...
        '\n\n Press any key to start the mandatory break of at least 10 seconds.'];
end
DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.black);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Wait at least 10 Seconds between Blocks (only after Block 1 has finished, not after Block 6)
if TRAINING == 1 && percentTotalCorrect < 60
    waitTime = 10;
    intervalTime = 1;
    timePassed = 0;
    printTime = 10;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n You can repeat the training task afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.black);
    Screen('Flip',ptbWindow);

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n You can repeat the training task afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.black);
        Screen('Flip',ptbWindow);
    end
elseif BLOCK >= 1 && BLOCK < 6
    waitTime = 10;
    intervalTime = 1;
    timePassed = 0;
    printTime = 10;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.black);
    Screen('Flip',ptbWindow);
    disp('Break started');

    while timePassed < waitTime
        waitbar(timePassed/10, 'INITIALIZING EEG');
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.black);
        Screen('Flip',ptbWindow);
    end
    wbar = findall(0,'type','figure','tag','TMWWaitbar');
    delete(wbar)
    disp('Break done!')
end

%% Quit
Screen('CloseAll');