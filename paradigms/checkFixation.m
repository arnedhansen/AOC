function noFixation = checkFixation(screenCentreX, screenCentreY, fixCheckDuration)
    numSamples = (fixCheckDuration*1000)/2; % 1 sample every 2 ms at 500Hz
    numSamples = ceil(numSamples);
    fixThresh = numSamples * 0.80; % 80% threshold
    distOK = 45 + 45; % 1 degree from the center (+ 0.5 deg of ET error)

    % Initialize noFixation counter
    noFixation = 0;

    % Collect gaze data for a specified number of samples
    samples = zeros(numSamples, 2); % Initialize matrix for gaze samples
    i = 1;

    % Collect gaze data for 4 ms
    while GetSecs < startTimeSample + 0.004
        if i <= numSamples
            % Fetch gaze data sample
            evt = Eyelink('NewestFloatSample');
            if evt.time > 0 % Ensure a valid sample is received
                gaze_x = evt.gx(1);
                gaze_y = evt.gy(1);
                samples(i, :) = [gaze_x, gaze_y]; % Store gaze sample
                i = i + 1;
            end
        else
            break; % Stop if numSamples is reached
        end
    end

    % Check fixation
    validSamples = sum(samples(:, 1) > 0 & samples(:, 2) > 0); % Count valid samples
    if validSamples > fixThresh % Ensure enough valid samples before checking fixation
        xFix = sum(samples(:, 1) > screenCentreX - distOK & samples(:, 1) < screenCentreX + distOK);
        yFix = sum(samples(:, 2) > screenCentreY - distOK & samples(:, 2) < screenCentreY + distOK);

        if xFix <= fixThresh || yFix <= fixThresh
            % No fixation detected
            noFixation = 1;
        end
    else
        disp('Not enough valid samples for fixation check.');
    end
end