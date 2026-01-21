%% AOC Pre-Stimulus Fixation Check — Exclusion Trials
% Loads AOC_preStimFixation_*_nback/sternberg.mat, aggregates excluded/kept/invalid trial counts per subject and task. Plots bar or summary. No files saved.
%
% Key outputs:
%   Bar figures (excluded/kept/invalid counts per task)

startup
clear
dirs = dir('/Volumes/methlab/Students/Arne/AOC/data/merged/');
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%%
close all
excludedCounts = zeros(2, numel(subjects));
keptCounts = zeros(2, numel(subjects));
invalidCounts = zeros(2, numel(subjects));
totalCounts = zeros(2, numel(subjects));

% Loop through each subject and load data
for subj = 1:numel(subjects)
    try
        % Load nback data
        loadPathNback = fullfile('/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/', subjects{subj});
        load(fullfile(loadPathNback, ['AOC_preStimFixation_', subjects{subj}, '_nback.mat']), 'preStimFixInfo');
        excludedCounts(1, subj) = numel(preStimFixInfo.excludedTrials);
        keptCounts(1, subj) = numel(preStimFixInfo.keptTrials);
        invalidCounts(1, subj) = numel(preStimFixInfo.invalidTrials);
        totalCounts(1, subj) = preStimFixInfo.totalTrials;

        % Load sternberg data
        load(fullfile(loadPathNback, ['AOC_preStimFixation_', subjects{subj}, '_sternberg.mat']), 'preStimFixInfo');
        excludedCounts(2, subj) = numel(preStimFixInfo.excludedTrials);
        keptCounts(2, subj) = numel(preStimFixInfo.keptTrials);
        invalidCounts(2, subj) = numel(preStimFixInfo.invalidTrials);
        totalCounts(2, subj) = preStimFixInfo.totalTrials;
    end
end

% Plotting
tasks = {'N-back', 'Sternberg'};

for t = 1:2
    figure;
    if t == 1
        set(gcf, "Position", [0 0 750 1000])
    else
        set(gcf, "Position", [750 0 750 1000])
    end

    % Combine data for plotting
    dataToPlot = [excludedCounts(t, :); keptCounts(t, :); invalidCounts(t, :)];

    b = bar(categorical(subjects), dataToPlot', 'stacked');
    b(1).FaceColor = [1 0 0];   % Red for excluded
    b(2).FaceColor = [0 1 0];   % Green for kept
    b(3).FaceColor = [0 0 1];   % Blue for invalid

    ylabel('Number of Trials');
    title([tasks{t}, ': Excluded vs Kept vs Invalid Trials'], 'FontSize', 20);
    legend({'Excluded', 'Kept', 'Invalid'}, 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 15);
    xtickangle(90);
    
    if t == 1
        ylim([0 305])
        saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/AOC_preStimFixation_nback.png')
    else
        saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/data/controls/preStimFixation/AOC_preStimFixation_sternberg.png')
    end
end

