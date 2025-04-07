% Define participants
subjects = {'subj1', 'subj2', 'subj3', ... }; % <-- replace with your actual subject IDs
excludedCounts = zeros(1, numel(subjects));
keptCounts = zeros(1, numel(subjects));
totalCounts = zeros(1, numel(subjects));

% Loop through each subject and load data
for subj = 1:numel(subjects)
    % Set the correct path depending on OS
    if ispc
        loadPath = fullfile('W:\Students\Arne\AOC\data\controls\preStimFixation\', subjects{subj});
    else
        loadPath = fullfile('/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/', subjects{subj});
    end

    % Load the file
    load(fullfile(loadPath, ['AOC_preStimFixation_', subjects{subj}, '_nback.mat']), 'preStimFixInfo');

    % Store counts
    excludedCounts(subj) = numel(preStimFixInfo.excludedTrials);
    keptCounts(subj) = numel(preStimFixInfo.keptTrials);
    totalCounts(subj) = preStimFixInfo.totalTrials;
end

% Visualise
figure;
bar(categorical(subjects), [excludedCounts; keptCounts]', 'stacked');
ylabel('Number of Trials');
title('Excluded vs Kept Trials per Participant');
legend({'Excluded', 'Kept'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
xtickangle(45);
grid on;
