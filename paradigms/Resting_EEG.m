%% Resting state EEG for the AOC and GCP studies
% 2.5 minutes fixation cross
% 2.5 minutes of blank

%% Resting
TASK = 'Resting';
TRAINING = 0; % no training for Resting

% Run Resting
disp('RESTING EEG...');

%% start recording EEG
disp('STARTING EEG RECORDING...');
initEEG;

% Calibrate ET (Tobii Pro Fusion)
disp('CALIBRATING ET...');
calibrateET;

%% Settings
testmode = 0;
monitorwidth_cm = 53;
dist_cm = 80;

% Text
tSize1 = 15;
tSize2 = 20;
tSize3 = 30;
colorText = 0;
colorBG = [];
colorBrightBG = [255,255,255];
colorInfo = [255,0,0];
colorBrightGray = [];
colorDarkGray = [];

%% Instructions
ins=struct();
ins.misc=struct();
ins.misc.mouse = [...
    'Press any key to start the task'...
    ];
ins.misc.finished = [...
    'Finished!'...
    ];
ins.resting=struct();
ins.resting.inst = [...
    'Experiment: Resting EEG' ...
    '\n\n\n You will see a cross in the middle of the screen. '...
    '\n\n Focus your gaze on this cross. \n'...
    ];
ins.resting.end = [...
    'Thank you! Now we will go on to other tasks. '...
    ];

%% Setting the Trigger codes
START = 10;
FLIPTONOFIXCROSS = 50;
END  = 90;

%% Screen Calculations
[scresw, scresh]=Screen('WindowSize',whichScreen);  % Get screen resolution
center = [scresw scresh]/2;     % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)
fixRect = [center-2 center+2];  % fixation dot
hz=Screen('FrameRate', whichScreen, 1);
cm2px = scresw/monitorwidth_cm;     % multiplication factor to convert cm to pixels
deg2px = dist_cm*cm2px*pi/180;      % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).
load gammafnCRT;   % load the gamma function parameters for this monitor - or some other CRT and hope they're similar! (none of our questions rely on precise quantification of physical contrast)
maxLum = GrayLevel2Lum(255,Cg,gam,b0);
par.BGcolor = Lum2GrayLevel(maxLum/2,Cg,gam,b0);

%% Experiment ptbWindow
clc;
ptbWindow=Screen('OpenWindow', whichScreen, par.BGcolor); % dont need to open a screen again
Screen('TextSize', ptbWindow, tSize2);
DrawFormattedText(ptbWindow, ins.resting.inst, scresw / 3, scresh / 3, colorText);
DrawFormattedText(ptbWindow, ins.misc.mouse,'center', 0.9*scresh, colorText);
Screen('Flip', ptbWindow);
HideCursor(whichScreen);
clc;
disp('RESTING STATE. THE SUBJECT IS READING THE INSTRUCTIONS...');

% Wait for participant to start the resting
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Send start triggers
Eyelink('Message', num2str(START));
Eyelink('command', 'record_status_message "START"');
sendtrigger(START,port,SITE,stayup)

disp('STARTING RESTING STATE');
% Display fixation cross for 2.5 minutes
Screen('DrawLine', ptbWindow,[0 0 0],center(1)-12,center(2), center(1)+12,center(2));
Screen('DrawLine', ptbWindow,[0 0 0],center(1),center(2)-12, center(1),center(2)+12);
Screen('Flip',ptbWindow);
WaitSecs(150);

%% Switch to blank
% Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
Screen('Flip', ptbWindow);
% Send blank triggers
Eyelink('Message', num2str(FLIPTONOFIXCROSS));
Eyelink('command', 'record_status_message "FLIPTONOFIXCROSS"');
sendtrigger(FLIPTONOFIXCROSS,port,SITE,stayup)
disp('2.5 MINUTES DONE. FLIP TO BLANK...');
WaitSecs(150);

%% Send end triggers
Eyelink('Message', num2str(END));
Eyelink('command', 'record_status_message "END"');
sendtrigger(END,port,SITE,stayup)
disp('RESTING EEG FINISHED');

%% Display end texts
Screen('TextSize', ptbWindow, tSize3);
DrawFormattedText(ptbWindow, ins.misc.finished,'center', 0.4*scresh, colorText);
Screen('TextSize', ptbWindow, tSize2);
DrawFormattedText(ptbWindow, ins.resting.end,'center', 0.5*scresh, colorText);
Screen('Flip', ptbWindow);
ShowCursor(whichScreen);
WaitSecs(5);

%% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)
save(fullfile(filePath, [subjectID,'_', TASK, '.mat']),'par','eyeO','eyeC');

%% Close and save EEG and ET
disp('SAVING DATA...');
closeEEGandET;

sca; % If Eyetracker wasn't used, close the Screens now
