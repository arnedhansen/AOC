function [trialsToKeep, excludedTrialIdx, distL, invalidTrials] = fixCheck(dataET, window, fixThresh, distOK)
% fixCheck checks if participants were fixating the centre during pre-stimulus period
%
% Inputs:
%   - dataET: FieldTrip data structure with eye-tracking channels
%   - window: [start, end] in seconds (e.g., [-0.5 0])
%   - fixThresh: proportion of samples that must fall within allowed distance
%   - distOK: maximum Euclidean distance from centre (in pixels, same unit as gaze)
%
% Outputs:
%   - trialsToKeep: logical index vector for good trials
%   - excludedTrialIdx: indices of trials excluded due to poor fixation (with valid gaze)
%   - distL: cell array of distance vectors per trial (valid samples only)
%   - invalidTrials: indices of trials with no valid gaze data (e.g. outside screen)

    screenCentreX = 400;  % half of screen width (800)
    screenCentreY = 300;  % half of screen height (600)

    nTrials = numel(dataET.trial);
    trialsToKeep = true(nTrials, 1);  % Initialise all to true
    invalidTrials = [];              % Trials with no valid samples
    distL = cell(nTrials, 1);         % Store distances from centre

    for t = 1:nTrials
        timeVec = dataET.time{t};
        idx = find(timeVec >= window(1) & timeVec <= window(2));

        if isempty(idx)
            trialsToKeep(t) = false;
            invalidTrials(end+1) = t;
            continue;
        end

        % Gaze positions
        xL = dataET.trial{t}(strcmp(dataET.label, 'L-GAZE-X'), idx);
        yL = dataET.trial{t}(strcmp(dataET.label, 'L-GAZE-Y'), idx);

        % Screen bounds check
        validIdx = xL >= 0 & xL <= 800 & yL >= 0 & yL <= 600;
        xL = xL(validIdx);
        yL = yL(validIdx);

        if isempty(xL)
            trialsToKeep(t) = false;
            invalidTrials(end+1) = t;
            continue;
        end

        % Compute Euclidean distance from screen centre
        distL{t} = sqrt((xL - screenCentreX).^2 + (yL - screenCentreY).^2);

        % Proportion of valid samples within allowed distance
        okSamples = (distL{t} <= distOK);
        if mean(okSamples) < fixThresh
            trialsToKeep(t) = false;
        end
    end

    % Indices of trials to exclude: both poor fixation and invalid gaze
    excludedTrialIdx = find(~trialsToKeep);
end
