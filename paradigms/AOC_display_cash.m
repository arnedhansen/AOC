%% AOC_display_cash.m

%% Calculate cash payment from time and accuracy

% Load timing
cd(DATA_PATH)
for block = {0, 6}
    % Load Sternberg accuracy data
    tsk = 'AOC_Sternberg';
    if block{1, 1} == 0
        fileName = [subjectID '_', tsk, '_block', num2str(block{1, 1}), '_training.mat'];
    elseif block{1, 1} == 6
        fileName = [subjectID '_', tsk, '_block', num2str(block{1, 1}), '_task.mat'];
    end
    dataTimeSternberg = load(['/home/methlab/Desktop/AOC_data/', subjectID, '/', fileName]);
    if block{1, 1} == 0
        startTimeSternberg = dataTimeSternberg.saves.timing.startTime;
    elseif block{1, 1} == 6
        endTimeSternberg = dataTimeSternberg.saves.timing.endTime;
    end

    % Load Nback accuracy data
    tsk = 'AOC_NBack';
    dataDir = ['/home/methlab/Desktop/AOC_data/', subjectID, '/'];
    if block{1, 1} == 0
        filePattern = [subjectID, '_', tsk, '_block', num2str(block{1, 1}), '_training.mat'];
    elseif block{1, 1} == 6
        filePattern = [subjectID, '_', tsk, '_block', num2str(block{1, 1}), '_*_task.mat'];
    end
    files = dir(fullfile(dataDir, filePattern));
    fileName = files.name;
    dataTimeNback = load(['/home/methlab/Desktop/AOC_data/', subjectID, '/', fileName]);
    if block{1, 1} == 0
        startTimeNback = dataTimeNback.saves.timing.startTime;
    elseif block{1, 1} == 6
        endTimeNback = dataTimeNback.saves.timing.endTime;
    end
end

% Convert times to datetime format
startTimes = [datetime(startTimeNback, 'InputFormat', 'dd/MM/yy-HH:mm:ss'), datetime(startTimeSternberg, 'InputFormat', 'dd/MM/yy-HH:mm:ss')];
endTimes = [datetime(endTimeNback, 'InputFormat', 'dd/MM/yy-HH:mm:ss'), datetime(endTimeSternberg, 'InputFormat', 'dd/MM/yy-HH:mm:ss')];

% Find earliest start time and latest end time
startTime = min(startTimes);
endTime = max(endTimes);

% Adjust start time if it's later than 12:00
if hour(startTime) > 12
    startTime = dateshift(startTime, 'start', 'day') + hours(13.5); % 13:30
else
    startTime = dateshift(startTime, 'start', 'day') + hours(9); % 9:00
end

% Calculate duration in minutes
durationMinutes = minutes(endTime - startTime) + 10; % Add 10 minutes for showering and hair drying

% Calculate bonus payment for accuracy
bonusPayment = 0;
threshSternberg = 80; % Minimum accuracy for percCorr(:, 1)
threshNback = 90; % Minimum accuracy for percCorr(:, 2)
bonusPerBlock = 15 / 12; % Max bonus split across 12 blocks
for i = 1:size(percCorr, 1)
    if percCorr(i, 1) >= threshSternberg
        blockBonus1 = ((percCorr(i, 1) - threshSternberg) / (100 - threshSternberg)) * bonusPerBlock;
    else
        blockBonus1 = 0;
    end
    if percCorr(i, 2) >= threshNback
        blockBonus2 = ((percCorr(i, 2) - threshNback) / (100 - threshNback)) * bonusPerBlock;
    else
        blockBonus2 = 0;
    end
    bonusPayment = bonusPayment + blockBonus1 + blockBonus2;
end
bonusPayment = min(bonusPayment, 15); % Cap the bonus at 15 CHF

% Total payment
basePayment = durationMinutes * 0.4167; % Assuming 0.4167 CHF per minute or 25 CHF per hour
cashPayment = basePayment + bonusPayment;
cashPaymentRounded = ceil(cashPayment / 5) * 5;

% Display results in a popup
resultMessage = sprintf(['Start Time: %s\n', ...
    'End Time: %s\n', ...
    'Duration (incl. 10 mins for showering): %.2f minutes\n', ...
    'Base Payment: %.2f CHF\n', ...
    'Bonus Payment: %.2f CHF\n', ...
    'Total Cash Payment: %.2f CHF\n', ...
    'Recommended Payment: %.2f CHF'], ...
    datestr(startTime), datestr(max(endTimes)),  ...
    durationMinutes, basePayment, bonusPayment, cashPayment, cashPaymentRounded);

msgbox(resultMessage, 'Cash Payment Calculation');