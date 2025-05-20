%% AOC Alpha Power Sternberg (trial-by-trial)

%% Load data struct
clear
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat');

%% Plot trial-by-trial ALPHA POWER GRAND AVERAGE
%% Grand-average α-power across subjects (trial-by-trial)

% Load data (if not already in workspace)
% load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat');

% Identify unique subjects
subjIDs = unique([eeg_data_sternberg_trials.ID]);

% Determine blocks & trials per block from first subject
sel0 = [eeg_data_sternberg_trials.ID] == subjIDs(1);
firstData = eeg_data_sternberg_trials(sel0);
blocks0    = unique([firstData.Block]);
nBlocks    = numel(blocks0);
nTrials    = numel(firstData) / nBlocks;   % assumes equal trials per block
totalTrials = nBlocks * nTrials;

% Preallocate: rows = subjects, cols = trials
alphaMatrix = nan(numel(subjIDs), totalTrials);

for si = 1:numel(subjIDs)
    sid = subjIDs(si);
    sel = [eeg_data_sternberg_trials.ID] == sid;
    subj = eeg_data_sternberg_trials(sel);
    % assume subj is already ordered by block then trial
    alphaMatrix(si, :) = [subj.AlphaPower];
end

% Compute grand mean
grandMean = nanmean(alphaMatrix, 1);

% Plot in same style
figure;
set(gcf,'Position',[0 0 1800 1200],'Colour','w');
hold on;
bar(grandMean, 1, 'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');

% Add block separators
for b = 1:(nBlocks-1)
    xline(b * nTrials, '--k', 'LineWidth', 1);
end

% Labels and formatting
xlabel('Trial');
ylabel('Alpha Power [µV^2 / Hz]');
title('Grand-average Sternberg α-Power over Trials');
set(gca,'FontSize',16);

%% Plot trial-by-trial ALPHA POWER for each INDIVIDUAL SUBJECT
for s = 1:numel(subjIDs)
    close all
    subjectID = num2str(subjIDs(s));
    
    % Extract just this subject’s trials
    sel = [eeg_data_sternberg_trials.ID] == subjIDs(s);
    subjData = eeg_data_sternberg_trials(sel);
    blocks = repelem(1:6,50)';
    blocks = num2cell(blocks);
    [subjData.Block] = blocks{:};
    
    % Determine how many blocks you have
    blocks   = unique([subjData.Block]);
    nBlocks  = numel(blocks);
    
    % Get colours for conditions 2,4,6
    colors    = color_def('AOC');
    condVals  = [2 4 6];
    condIdx   = arrayfun(@(x) find(condVals==x), condVals);
    condCols  = colors(condIdx, :);
    
    % Create the figure
    figure;
    set(gcf, 'Position', [0 0 1800 1200], 'Color', 'w');
    hold on;
    xPos = 0;
    xTickPositions = [];
    xTickLabels = [];
    
    for b = 1:nBlocks
        thisBlock = blocks(b);
        blkData   = subjData([subjData.Block] == thisBlock);
        nTrials(b)   = numel(blkData);
        
        % trial‐by‐trial values
        alphaVals = [blkData.AlphaPower];
        setSizes  = [blkData.Condition];  % 2,4 or 6
        
        % x positions for this block
        xB = xPos + (1:nTrials(b));
        
        % bar each trial
        for t = 1:nTrials(b)
            cs = setSizes(t);
            ci = find(condVals==cs);
            bar(xB(t), alphaVals(t), 0.8, ...
                'FaceColor', condCols(ci,:), 'EdgeColor','none');
        end

        % Add block infos
        xTickPositions = [xTickPositions, xPos + nTrials(b)/2];
        xTickLabels{end+1} = sprintf('Block %d', b);

        % draw vertical separator after block (except last)
        if b < nBlocks
            xline(xB(end), '--k','LineWidth',1);
        end

        % advance xPos: end of this block + gap of 2
        xPos = xB(end);

    end
    
    % Finish axes
    xlabel('Trial');
    ylabel('Alpha Power [µV^2 / Hz]');
    title(sprintf('Sternberg Alpha Power over Trials — Subject %s', subjectID));
    set(gca, 'XTick', xTickPositions, 'XTickLabel', xTickLabels, 'FontSize', 16);
    
    % Legend
    h = gobjects(3,1);
    for i = 1:3
        h(i) = bar(nan, nan, 0.8, 'FaceColor', condCols(i,:), 'EdgeColor','none');
    end
    legend(h, {'WM load 2','WM load 4','WM load 6'}, 'Location','northeast');
    
    % Save
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_over_trials/AOC_alpha_power_over_trials_sternberg_sub' subjectID '.png']);
    fprintf('Plot trial-by-trial ALPHA POWER for Subject %s done... \n', subjectID);
end


