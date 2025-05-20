%% AOC Alpha Power Sternberg (trial-by-trial)

%% Load data struct
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat');

%% Plot trial-by-trial ALPHA POWER for each INDIVIDUAL SUBJECT
% Unique subject IDs
subjIDs = unique([eeg_data_sternberg_trials.ID]);

for s = 1:numel(subjIDs)
    subjectID = num2str(subjIDs(s));
    
    % Extract just this subject’s trials
    sel = [eeg_data_sternberg_trials.ID] == subjIDs(s);
    subjData = eeg_data_sternberg_trials(sel);
    
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
    xTicks   = zeros(1,nBlocks);
    xLabels  = cell(1,nBlocks);
    
    for b = 1:nBlocks
        thisBlock = blocks(b);
        blkData   = subjData([subjData.Block] == thisBlock);
        nTrials   = numel(blkData);
        
        % trial‐by‐trial values
        alphaVals = [blkData.AlphaPower];
        setSizes  = [blkData.Condition];  % 2,4 or 6
        
        % x positions for this block
        xB = xPos + (1:nTrials);
        
        % bar each trial
        for t = 1:nTrials
            cs = setSizes(t);
            ci = find(condVals==cs);
            bar(xB(t), alphaVals(t), 0.8, ...
                'FaceColor', condCols(ci,:), 'EdgeColor','none');
        end
        
        % mark centre for x‐tick
        xTicks(b)  = xPos + nTrials/2;
        xLabels{b} = sprintf('Block %d', thisBlock);
        
        % draw vertical separator after block (except last)
        if b < nBlocks
            xline(xB(end)+1, '--k','LineWidth',1);
        end
        
        % advance xPos: end of this block + gap of 2
        xPos = xB(end) + 2;
    end
    
    % Finish axes
    xlabel('Trials (blocks separated by gap)');
    ylabel('Alpha Power (µV^2)');
    title(sprintf('Sternberg α-Power over Trials — Subject %s', subjectID));
    set(gca, 'XTick', xTicks, 'XTickLabel', xLabels, 'FontSize', 16);
    
    % Legend
    h = gobjects(3,1);
    for i = 1:3
        h(i) = bar(nan, nan, 0.8, 'FaceColor', condCols(i,:), 'EdgeColor','none');
    end
    legend(h, {'WM load 2','WM load 4','WM load 6'}, 'Location','northeast');
    
    % Save
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_over_trials/AOC_alpha_power_over_trials_sternberg_sub' subjectID '.png']);
    close;
    fprintf('Plot trial-by-trial ALPHA POWER for Subject %s done... \n', subjectID);
end
