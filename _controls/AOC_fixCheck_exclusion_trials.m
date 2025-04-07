% Define participants
startup
clear
dirs = dir('/Volumes/methlab/Students/Arne/AOC/data/merged/');
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%%
close all
excludedCounts = zeros(2, numel(subjects));  % Row 1 = nback, Row 2 = sternberg
keptCounts = zeros(2, numel(subjects));
totalCounts = zeros(2, numel(subjects));

% Loop through each subject and load data
for subj = 1:numel(subjects)
    try
        % Load nback data
        loadPathNback = fullfile('/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/', subjects{subj});
        load(fullfile(loadPathNback, ['AOC_preStimFixation_', subjects{subj}, '_nback.mat']), 'preStimFixInfo');
        excludedCounts(1, subj) = numel(preStimFixInfo.excludedTrials);
        keptCounts(1, subj) = numel(preStimFixInfo.keptTrials);
        totalCounts(1, subj) = preStimFixInfo.totalTrials;

        % Load sternberg data
        load(fullfile(loadPathNback, ['AOC_preStimFixation_', subjects{subj}, '_sternberg.mat']), 'preStimFixInfo');
        excludedCounts(2, subj) = numel(preStimFixInfo.excludedTrials);
        keptCounts(2, subj) = numel(preStimFixInfo.keptTrials);
        totalCounts(2, subj) = preStimFixInfo.totalTrials;
    end
end

% Plotting
tasks = {'N-back', 'Sternberg'};

for t = 1:2
    figure;
    b = bar(categorical(subjects), [excludedCounts(t, :); keptCounts(t, :)]', 'stacked');
    b(1).FaceColor = [1 0 0];   % Red for excluded
    b(2).FaceColor = [0 1 0];   % Green for kept
    ylabel('Number of Trials');
    title([tasks{t}, ': Excluded vs Kept Trials per Participant']);
    legend({'Excluded', 'Kept'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
    xtickangle(45);
    grid on;
end
