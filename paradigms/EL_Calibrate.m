try
    % Setup for Calibration
    % Close any open audio devices to avoid conflicts
    try
        Snd('Close');
        PsychPortAudio('Close');
    catch
        % Ignore errors if there are no open devices
    end
    
    fprintf('Starting Calibration \n');
    
    % Configure Eyelink settings
    Eyelink('Command', 'saccade_velocity_threshold = 35'); % Set saccade velocity threshold
    Eyelink('Command', 'saccade_acceleration_threshold = 9500'); % Set saccade acceleration threshold
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA,PUPIL,HREF'); % Request left and right eye data, gaze position, pupil size, and head-reflected distance (HREF)
    %Eyelink('Command', 'active_eye = LEFT'); % Set the active eye to left
    Eyelink('Command', 'calibration_type = HV9'); % Set calibration type to 9-point
    Eyelink('Command', 'enable_automatic_calibration = YES'); % Enable automatic calibration
    Eyelink('Command', 'automatic_calibration_pacing = 500'); % Set calibration pacing to 500ms
    Eyelink('Command', 'set_idle_mode'); % Set the tracker to idle mode
    
    % Run Calibration
    HideCursor(whichScreen); % Hide the cursor during calibration
    Eyelink('StartSetup', 1); % Start the setup and calibration procedure
    disp('Calibration done');
    
    % Close any open audio devices again to ensure clean state
    try
        Snd('Close');
        PsychPortAudio('Close');
    catch
        % Ignore errors if there are no open devices
    end
    Screen('CloseAll'); % Close all open screens
    
catch ME
    % Display error message if something goes wrong
    disp('Error running the calibration');
    disp(ME.message);
end
