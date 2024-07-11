% #AOC Nback
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2023a.
%
% N-back task with 6 blocks of 1-, 2- and 3-back conditions in randomized order.

%% Initialize EEG and ET
if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;
end

% Calibrate ET
disp('CALIBRATING ET...');
calibrateET

% Hide cursor
HideCursor(whichScreen);

%% Define triggers
TASK_START = 10; % Trigger for start of task
MATCH = 4; % trigger for matching condition
NO_MATCH = 5; % trigger for non-matching condition
FIXATION = 15; % trigger for fixation cross
PRESENTATION0 = 19; % trigger for letter presentation (training task; 1-back)
PRESENTATION1 = 21; % trigger for letter presentation (1-back)
PRESENTATION2 = 22; % trigger for letter presentation (2-back)
PRESENTATION3 = 23; % trigger for letter presentation (3-back)
STIMOFF = 28; % trigger for change of letter to cfi
COND1 = 31; % Trigger for block with 1-back
COND2 = 32; % Trigger for block with 2-back
COND3 = 33; % Trigger for block with 3-back
BLOCK0 = 60; % trigger for start of training block
BLOCK1 = 61; % trigger for start of block 1
BLOCK2 = 62; % trigger for start of block 2
BLOCK3 = 63; % trigger for start of block 3
BLOCK4 = 64; % trigger for start of block 4
BLOCK5 = 65; % trigger for start of block 5
BLOCK6 = 66; % trigger for start of block 6
ENDBLOCK0 = 70; % trigger for end of training block
ENDBLOCK1 = 71; % trigger for end of block 1
ENDBLOCK2 = 72; % trigger for end of block 2
ENDBLOCK3 = 73; % trigger for end of block 3
ENDBLOCK4 = 74; % trigger for end of block 4
ENDBLOCK5 = 75; % trigger for end of block 5
ENDBLOCK6 = 76; % trigger for end of block 6
RESP_YES = 87; % trigger for response yes (spacebar)
RESP_NO = 88; % trigger for response no (no input)
RESP_WRONG = 89;% trigger for wrong keyboard input response
TASK_END = 90; % Trigger for end of task

%% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    experiment.nTrials = 12;
else
    experiment.nTrials = 100;           % 6 blocks x 100 trials = 600 trials
end

% Set up equipment parameters
equipment.viewDist = 680;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
equipment.greyVal = .5;
equipment.blackVal = 0;
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

% Set up stimulus parameters Fixation
stimulus.fixationOn = 1;                % Toggle fixation on (1) or off (0)
stimulus.fixationSize_dva = .5;         % Size of fixation cross in degress of visual angle
stimulus.fixationColor = 1;             % Color of fixation cross (1 = white)
stimulus.fixationLineWidth = 3;         % Line width of fixation cross

% Location
stimulus.regionHeight_dva = 7.3;         % Height of the region
stimulus.regionWidth_dva = 4;            % Width of the region
stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

% Set up text parameters
text.instructionFont = 'Menlo';         % Font of instruction text
text.instructionPoints = 12;            % Size of instruction text (This if overwritten by )
text.instructionStyle = 0;              % Styling of instruction text (0 = normal)
text.instructionWrap = 80;              % Number of characters at which to wrap instruction text

% Set up color parameters
color.Black = 0;                      
color.White = 1;

%% Define startExperimentText
if TRAINING == 1 && COND == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
        'You will see a series of random letters. \n\n' ...
        'Your task is to press SPACE if you see the same letter twice in a row. \n\n' ...
        'Otherwise, do not press any button. \n\n' ...
        'Please always use your right hand.' ...
        '\n\n Don''t worry, you can do a training sequence in the beginning. \n\n' ...
        '\n\n Press any key to continue.'];
elseif TRAINING == 1 && COND == 2
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
        'You will see a series of random letters. \n\n' ...
        'Your task is to press SPACE if the letter you see \n\n' ...
        'is the same letter as the one two letters before. \n\n' ...
        'Example: K  -  Q  -  K \n\n' ...
        'Otherwise, do not press any button. \n\n' ...
        'Please always use your right hand.' ...
        '\n\n Press any key to continue.'];
elseif TRAINING == 1 && COND == 3
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
        'You will see a series of random letters. \n\n' ...
        'Your task is to press SPACE if the letter you see \n\n' ...
        'is the same letter as the one three letters before. \n\n' ...
        'Example: S - Q - P - S \n\n' ...
        'Otherwise, do not press any button. \n\n' ...
        'Please always use your right hand.' ...
        '\n\n Press any key to continue.'];
else
    if COND == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if you see the same letter twice in a row. \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    elseif COND == 2
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if the letter you see \n\n' ...
            'is the same letter as the one two letters before. \n\n' ...
            'Example: K  -  Q  -  K \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    elseif COND == 3
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if the letter you see \n\n' ...
            'is the same letter as the one three letters before. \n\n' ...
            'Example: S - Q - P - S \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    end
end

performanceBonusText = ['In the following task there is a performance bonus! \n\n' ...
    'Try to be as accurate as possible. \n\n \n\n' ...
    'Press any key to continue.'];

%% Setup
% Set up temporal parameters (all in seconds)
timing.blank = 1;                   % Duration of blank screen

% Shuffle rng for random elements
rng('default');
rng('shuffle');                     % Use MATLAB twister for rng

% Set up Psychtoolbox Pipeline
AssertOpenGL;

% Imaging set up
screenID = whichScreen;
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference', 'SkipSyncTests', 0); % For linux

% Window set-up
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
experiment.runPriority = MaxPriority(ptbWindow);

% Set font size for instructions and stimuli
Screen('TextSize', ptbWindow, 20);

global psych_default_colormode;                     % Sets colormode to be unclamped 0-1 range.
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text
ptb_drawformattedtext_disableClipping = 1;

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',color.Black);
Screen('Flip',ptbWindow);

% Retrieve response key
spaceKeyCode = KbName('Space'); % Retrieve key code for spacebar

% Calculate equipment parameters
equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
equipment.ppd = equipment.ppm*equipment.mpd;    % Pixels per degree

% Fix coordiantes for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];

%% Create data structure for preallocating data
data = struct;
data.letterSequence = NaN; % Stimuli
data.match(1:experiment.nTrials) = NaN; % Matching
data.responses(1:experiment.nTrials) = NaN; % Button press response
data.correct(1:experiment.nTrials) = NaN; % Correct answer
data.stimulus(1:experiment.nTrials) = NaN; % Probe stimulus
data.reactionTime(1:experiment.nTrials) = NaN; % Reaction time
data.fixation(1:experiment.nTrials) = NaN; % Reaction time 
count5trials = NaN; % Feedback variable

%% Define letterSequence depending on block iteration
letterSequence = createLetterSequence(TRAINING, COND);
data.letterSequence = letterSequence; % 1x102 double

%% Display instructive texts
if TRAINING == 0
    % Show performance bonus incentive text
    DrawFormattedText(ptbWindow,performanceBonusText,'center','center',color.Black);
    Screen('Flip',ptbWindow);
    disp('Participant is reading the performance bonus text');
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end
end

% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.Black);
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
Screen('Flip',ptbWindow);
disp('Participant is reading the instructions.');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Send triggers
% Send triggers for task start
if TRAINING == 1
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "TASK_START"');
else
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "TASK_START"');
    sendtrigger(TASK_START,port,SITE,stayup);
end

% Send triggers for condition
if COND == 1
    TRIGGER = COND1;
elseif COND == 2
    TRIGGER = COND2;
elseif COND == 3
    TRIGGER = COND3;
end

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
else
    Eyelink('Message', num2str(TRIGGER));
    sendtrigger(TRIGGER,port,SITE,stayup); % EEG
end

% Send triggers for block 
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
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end

if TRAINING == 1
    disp('Start of Training Block.');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end
HideCursor(whichScreen);
tic;

%% Experiment Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trl = 1:experiment.nTrials
    %% Jittered CFI before presentation of letter (3000ms +/- 1000ms)
    Screen('DrawLines',ptbWindow,fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2); % Draw fixation cross
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    TRIGGER = STIMOFF;
    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMOFF"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMOFF"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    TRIGGER = FIXATION;
    timing.cfi(trl) = (randsample(2000:4000, 1))/1000;    % Randomize the jittered central fixation interval on trial
    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "FIXATION"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "FIXATION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.cfi(trl));                            % Wait duration of the jittered central fixation interval

    %% Check fixation just before stimulus presentation
    noFixation = checkFixation(screenWidth, screenHeight, screenCentreX, screenCentreY);
    
    %% Present stimulus from letterSequence (2000ms)
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 60);
    DrawFormattedText(ptbWindow,[letterSequence(trl)],'center','center',color.White);
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
    Screen('Flip', ptbWindow);
    probePresentationTime = GetSecs;
    % Return size of text to default
    Screen('TextSize', ptbWindow, 20);
    % Send triggers for Presentation
    if TRAINING == 1
        TRIGGER = PRESENTATION0;
    elseif BLOCK == 1
        TRIGGER = PRESENTATION1;
    elseif BLOCK == 2
        TRIGGER = PRESENTATION2;
    elseif BLOCK == 3
        TRIGGER = PRESENTATION3;
    end

    if TRAINING == 1
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
    else
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    data.stimulus(trl) = letterSequence(trl);

    %% Get response
    getResponse = true;
    badResponseFlag = false;
    maxResponseTime = GetSecs + 2;
    responseTime = NaN;
    while getResponse
        [time,keyCode] = KbWait(-1, 2, maxResponseTime); % Wait 2 seconds for response, continue afterwards if there is no input.
        whichKey = find(keyCode);

        if ~isempty(whichKey)
            responseTime = GetSecs;
            if whichKey == spaceKeyCode
                getResponse = false;
                data.responses(trl) = whichKey;
                TRIGGER = RESP_YES;
            else
                TRIGGER = RESP_WRONG;
                data.responses(trl) = whichKey;
                badResponseFlag = true;
            end

        elseif isempty(whichKey)
            data.responses(trl) = 0;
            TRIGGER = RESP_NO;
            responseTime = probePresentationTime + 5;
        end

        % send triggers
        if TRAINING == 1
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "RESPONSE"');
        else
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "RESPONSE"');
            sendtrigger(TRIGGER,port,SITE,stayup)
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

    %% Save match/no match
    if COND == 1 && trl > 1
        if letterSequence(trl-1) == letterSequence(trl)
            trlMatch = 1;
            TRIGGER = MATCH;
            el_msg = 'MATCH';
        else
            trlMatch = 0;
            TRIGGER = NO_MATCH;
            el_msg = 'NO MATCH';
        end
        data.match(trl) = trlMatch;
        % Send triggers: trial matching. If training, send only ET triggers
        if TRAINING == 1
            Eyelink('Message', el_msg);
        else
            Eyelink('Message', el_msg);
            sendtrigger(TRIGGER,port,SITE,stayup); % EEG
        end
    elseif COND == 2 && trl > 2
        if letterSequence(trl-2) == letterSequence(trl)
            trlMatch = 1;
            TRIGGER = MATCH;
            el_msg = 'MATCH';
        else
            trlMatch = 0;
            TRIGGER = NO_MATCH;
            el_msg = 'NO MATCH';
        end
        data.match(trl) = trlMatch;
        % Send triggers: trial matching. If training, send only ET triggers
        if TRAINING == 1
            Eyelink('Message', el_msg);
        else
            Eyelink('Message', el_msg);
            sendtrigger(TRIGGER,port,SITE,stayup); % EEG
        end
    elseif COND == 3 && trl > 3
        if letterSequence(trl-3) == letterSequence(trl)
            trlMatch = 1;
            TRIGGER = MATCH;
            el_msg = 'MATCH';
        else
            trlMatch = 0;
            TRIGGER = NO_MATCH;
            el_msg = 'NO MATCH';
        end
        data.match(trl) = trlMatch;
        % Send triggers: trial matching. If training, send only ET triggers
        if TRAINING == 1
            Eyelink('Message', el_msg);
        else
            Eyelink('Message', el_msg);
            sendtrigger(TRIGGER,port,SITE,stayup); % EEG
        end
    end

    % Check if response was correct
    if COND == 1 && trl > 1
        if trlMatch == 1 && data.responses(trl) == spaceKeyCode  % Correct matched trial
            data.correct(trl) = 1;
        elseif trlMatch == 1 && data.responses(trl) == 0  % Incorrect matched trial
            data.correct(trl) = 0;
        elseif trlMatch == 0 && data.responses(trl) == 0  % Correct unmatched trial
            data.correct(trl) = 1;
        elseif trlMatch == 0 && data.responses(trl) == spaceKeyCode  % Incorrect unmatched trial
            data.correct(trl) = 0;
        elseif data.responses(trl) ~= spaceKeyCode
            data.correct(trl) = 0;
        end
    elseif COND == 2 && trl > 2
        if trlMatch == 1 && data.responses(trl) == spaceKeyCode  % Correct matched trial
            data.correct(trl) = 1;
        elseif trlMatch == 1 && data.responses(trl) == 0  % Incorrect matched trial
            data.correct(trl) = 0;
        elseif trlMatch == 0 && data.responses(trl) == 0  % Correct unmatched trial
            data.correct(trl) = 1;
        elseif trlMatch == 0 && data.responses(trl) == spaceKeyCode  % Incorrect unmatched trial
            data.correct(trl) = 0;
        elseif data.responses(trl) ~= spaceKeyCode
            data.correct(trl) = 0;
        end
    elseif COND == 3 && trl > 3
        if trlMatch == 1 && data.responses(trl) == spaceKeyCode  % Correct matched trial
            data.correct(trl) = 1;
        elseif trlMatch == 1 && data.responses(trl) == 0  % Incorrect matched trial
            data.correct(trl) = 0;
        elseif trlMatch == 0 && data.responses(trl) == 0  % Correct unmatched trial
            data.correct(trl) = 1;
        elseif trlMatch == 0 && data.responses(trl) == spaceKeyCode  % Incorrect unmatched trial
            data.correct(trl) = 0;
        elseif data.responses(trl) ~= spaceKeyCode
            data.correct(trl) = 0;
        end
    end

    % Display (in-)correct response in CW
    if data.correct(trl) == 1 && trl > 1
        feedbackText = 'Correct!';
    elseif data.correct(trl) == 0 && badResponseFlag == false && trl > 1
        feedbackText = 'Incorrect!';
    elseif data.correct(trl) == 0 && badResponseFlag == true && trl > 1
        feedbackText = 'Wrong button! Use only SPACE.';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.Black);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(3);
    elseif trl == 1
        disp('No Response to Trial 1 in N-Back Task');
    elseif COND == 2 && trl == 2
        disp('No Response to Trial 2 in 3-Back Task');
    elseif COND == 3 && trl == 2
        disp('No Response to Trial 2 in 3-Back Task');
    elseif COND == 3 && trl == 3
        disp('No Response to Trial 3 in 3-Back Task');
    end

    %% Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 60%
    responsesLastTrials = 0;
    if COND == 1 && trl > 11
        % Get 10 last trials, but ignore last data point
        responsesLastTrials = data.correct(end-10 : end-1);
        percentLastTrialsCorrect = (sum(responsesLastTrials)/length(responsesLastTrials))*100;
        if percentLastTrialsCorrect < 60 && count5trials <= trl-5
            count5trials = trl;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n Of the last 10 trials ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                '\n\n You can earn more if you perform better.' ...
                '\n\n Please keep focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) '%.']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.Black);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(5);
        end
    elseif COND == 2 && trl > 12 || COND == 3 && trl > 13
        % Get 10 last trials, but ignore first two and last data point
        responsesLastTrials = data.correct(end-9 : end-1);
        percentLastTrialsCorrect = (sum(responsesLastTrials)/length(responsesLastTrials))*100;
        if percentLastTrialsCorrect < 70 && count5trials <= trl-5
            count5trials = trl;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n Of the last 10 trials ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                '\n\n You can earn more if you perform better.' ...
                '\n\n Please keep focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) '%.']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.Black);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(5);
        end
    end

    %% Fixation reminder
    if noFixation > 0
        Screen('TextSize', ptbWindow, 30);
        fixText = 'ALWAYS LOOK AT THE CENTER OF THE SCREEN!';
        DrawFormattedText(ptbWindow, fixText, 'center', 'center', color.White);
        Screen('DrawDots', ptbWindow, backPos, backDiameter, backColor, [], 1);
        Screen('Flip', ptbWindow);
        disp('FIXATION REMINDER')
        WaitSecs(3);
        data.fixation(trl) = 0;
        Screen('TextSize', ptbWindow, 20);
    else
        data.fixation(trl) = 1;
    end

    %% Trialf Info CW output
    overall_accuracy = round((sum(data.correct(1:trl), 'omitnan')/(trl-COND))*100);
    if trl == 1 || COND == 2 && trl == 2 || COND == 3 && trl == 2 || COND == 3 && trl == 3
        feedbackText = ('Correct! (no response)');
    end
    reactionTime = num2str(round(data.reactionTime(trl), 2), '%.2f');
    disp(['Response to Trial ' num2str(trl) ' in Block ' num2str(BLOCK) ' is ' feedbackText ' (' num2str(COND) '-back | Acc: ' num2str(overall_accuracy) '% | RT: ' num2str(reactionTime) 's)']);

end

%% Send triggers to end task
endT = Screen('Flip',ptbWindow);
if TRAINING == 1
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
else
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
    sendtrigger(TASK_END,port,SITE,stayup)
end

% Send triggers for block and output
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

if TRAINING == 1
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    disp('End of Training Block.');
else
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
    disp(['End of Block ' num2str(BLOCK)]);
end

% Record block duration
timing.duration = toc;

%% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)
if TRAINING == 1
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_' num2str(COND) 'back_training.mat'];
else
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_' num2str(COND) 'back_task.mat'];
end

% Save data structure
saves = struct;
saves.data = data;
saves.spaceKeyCode = spaceKeyCode;
saves.data.timing = timing;
saves.experiment = experiment;
saves.screen.screenWidth = screenWidth;
saves.screen.screenHeight = screenHeight;
saves.screen.screenCentreX = screenCentreX;
saves.screen.screenCentreY = screenCentreY;
saves.startExperimentText = startExperimentText;
saves.stimulus = stimulus;
saves.subjectID = subjectID;
saves.subject = subject;
saves.text = text;
saves.flipInterval = flipInterval;

%% Save triggers
trigger = struct;
trigger.TASK_START = TASK_START;
trigger.MATCH = MATCH;
trigger.NO_MATCH = NO_MATCH;
trigger.FIXATION = FIXATION;
trigger.PRESENTATION0 = PRESENTATION0;
trigger.PRESENTATION1 = PRESENTATION1;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.PRESENTATION3 = PRESENTATION3;
trigger.STIMOFF = STIMOFF;
trigger.COND1 = COND1;
trigger.COND2 = COND2;
trigger.COND3 = COND3;
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
trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;
trigger.RESP_WRONG = RESP_WRONG;
trigger.TASK_END = TASK_END;

if BLOCK == 6
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
end

%% Stop and close EEG and ET recordings and SAVE FILE
if TRAINING == 1
    disp('TRAINING FINISHED...');
else
    disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
end
disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger');
closeEEGandET;

try
    PsychPortAudio('Close');
catch
end

%% Compute accuracy and report after each block (no additional cash for training task)
if TRAINING == 1
    % Get sum of correct responses, but ignore first and last data point
    totalCorrect = sum(data.correct(1:end-1), 'omitnan');
    totalTrials = trl-2;
    percentTotalCorrect = round(totalCorrect / totalTrials * 100);
    format bank % Change format for display
    feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) '%. '];
    disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.Black);
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif COND == 1
    % Get sum of correct responses, but ignore first and last data point
    totalCorrect = sum(data.correct(1, 2:end-1));
    totalTrials = trl-2;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    format bank % Change format for display
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.01;
    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) '%. ' ...
        '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
        '\n\n Keep it up!'];

    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.Black);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif COND > 1
    % Get sum of correct responses, but ignore first 2/3/4 and last data point
    if COND == 2
        totalCorrect = sum(data.correct(1, 3:end-1));
        totalTrials = trl-3;
    elseif COND == 3
        totalCorrect = sum(data.correct(1, 4:end-1));
        totalTrials = trl-4;
    end
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    format bank % Change format for display
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.01;
    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) '%. ' ...
        '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
        '\n\n Keep it up!'];
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.Black);
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
elseif BLOCK == 3
    breakInstructionText = ['End of the Task! ' ...
        '\n\n Press any key to view your stats.'];
else
    breakInstructionText = ['Break! Rest for a while... ' ...
        '\n\n Press any key to start the mandatory break of at least 15 seconds.'];
end
DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.Black);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Wait at least 15 Seconds between Blocks (only after Block 1 has finished, not after Block 2)
if TRAINING == 1 && percentTotalCorrect < 60
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n You can repeat the training task afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
    Screen('Flip',ptbWindow);

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n You can repeat the training task afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
        Screen('Flip',ptbWindow);
    end
elseif BLOCK == 1 && TRAINING == 1
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n Block 1 of the N-back task will start afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
    Screen('Flip',ptbWindow);
    disp('Break started');

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n Block 1 of the N-back task will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
        Screen('Flip',ptbWindow);
        disp(printTime);
    end
elseif BLOCK > 1 && TRAINING == 0
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
    Screen('Flip',ptbWindow);
    disp('Break started');

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.Black);
        Screen('Flip',ptbWindow);
        disp(printTime);
    end
end

%% Save total amount earned and display
if BLOCK == 6
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
    format bank % Change format for display
    endTextCash = ['Well done! You have completed the task.' ...
        ' \n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextraTotal) ' CHF in total.' ...
        ' \n\n ' ...
        ' \n\n Block 1: ' num2str(percentTotalCorrect(1)) ' % accuracy earned you ' num2str(amountCHFextra(1)) ' CHF.' ...
        ' \n\n Block 2: ' num2str(percentTotalCorrect(2)) ' % accuracy earned you ' num2str(amountCHFextra(2)) ' CHF.' ...
        ' \n\n Block 3: ' num2str(percentTotalCorrect(3)) ' % accuracy earned you ' num2str(amountCHFextra(3)) ' CHF.' ...
        ' \n\n Block 4: ' num2str(percentTotalCorrect(4)) ' % accuracy earned you ' num2str(amountCHFextra(4)) ' CHF.' ...
        ' \n\n Block 5: ' num2str(percentTotalCorrect(5)) ' % accuracy earned you ' num2str(amountCHFextra(5)) ' CHF.' ...
        ' \n\n Block 6: ' num2str(percentTotalCorrect(6)) ' % accuracy earned you ' num2str(amountCHFextra(6)) ' CHF.' ...
        ' \n\n ' ...
        ' \n\n ' ...
        ' \n\n Press any key to end the task.'];
    DrawFormattedText(ptbWindow,endTextCash,'center','center',color.Black); % Display output for participant
    disp(['End of Block ' num2str(BLOCK) '. Participant ' num2str(subjectID) ' has earned CHF ' num2str(amountCHFextraTotal) ' extra in total.']);
    statsCW = ['Block 1: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(1)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(1)) '%' ...
        ' \n\n Block 2: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(2)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(2)) '%' ...
        ' \n\n Block 3: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(3)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(3)) '%' ...
        ' \n\n Block 4: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(4)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(4)) '%' ...
        ' \n\n Block 5: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(5)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(5)) '%' ...
        ' \n\n Block 6: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(6)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(6)) '%'];
    disp(statsCW)
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end
end

% Quit
Screen('CloseAll');