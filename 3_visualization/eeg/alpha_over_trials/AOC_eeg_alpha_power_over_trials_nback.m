%% Analyse Alpha Power Changes Over Trials Within Each Block (All Conditions Mixed)
% This script loads each block's EEG data for the N‑back task, extracts epochs
% for all conditions (1‑back: trigger '21', 2‑back: '22', 3‑back: '23'), combines
% them for each block (while retaining the trial's condition label), computes the
% alpha (8–14 Hz) power for each trial averaged over occipital channels, and finally
% plots a bar graph with a bar for each trial (grouped by block) and bars coloured
% by condition.

startup;
close all; clear; clc;
addEEGLab

%% Set data path based on operating system
if ispc
    basePath = 'W:\Students\Arne\AOC\data\merged\';
else
    basePath = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
end

%% List subjects and choose one subject for analysis
dirs = dir(basePath);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');  % Exclude unwanted subjects

for subs = 1:length(subjects)
    close all
    subjectID = subjects{subs};  % Change index if you wish to analyse a different subject

    %% Define conditions and corresponding triggers
    conds = {'21', '22', '23'};   % Triggers for 1‑back, 2‑back, 3‑back respectively
    condLabels = [1, 2, 3];        % Numeric labels (can later be mapped to condition names)

    nBlocks = 6;  % Number of blocks

    % Preallocate a cell array to store results for each block
    blockAlpha = cell(1, nBlocks);

    %% Loop over each block
    for block = 1:nBlocks
        clc
        fprintf('Processing Block %d for Subject %s...\n', block, subjectID);
        datapath = fullfile(basePath, subjectID);
        cd(datapath);

        % Load the block's EEG file
        try
            fileName = sprintf('%s_EEG_ET_Nback_block%d_merged.mat', subjectID, block);
            load(fileName);  % Loads variable EEG
            EEG_block = EEG;
        catch ME
            warning('Could not load block %d: %s', block, ME.message);
            continue;
        end

        % Container for FieldTrip data from each condition in this block
        ftDataAll = {};

        % Loop over each condition (trigger)
        for c = 1:length(conds)
            trig = conds{c};
            epoch_window = [-1.5 2.5];  % Time window (in seconds)
            try
                EEG_epoch = pop_epoch(EEG_block, {trig}, epoch_window);
                % Exclude trials with motor responses (trigger '4')
                matching_trials = find(strcmp({EEG_epoch.event.type}, '4'));
                if ~isempty(matching_trials)
                    exclude_epochs = unique([EEG_epoch.event(matching_trials).epoch]);
                    EEG_epoch = pop_select(EEG_epoch, 'notrial', exclude_epochs);
                end
            catch ME
                warning('Epoching error for condition %s in block %d: %s', trig, block, ME.message);
                continue;
            end

            % Convert to FieldTrip format
            try
                ftData = eeglab2fieldtrip(EEG_epoch, 'raw');
            catch ME
                warning('Conversion to FieldTrip failed for condition %s in block %d: %s', trig, block, ME.message);
                continue;
            end

            % Save the condition label (1 for 1‑back, 2 for 2‑back, 3 for 3‑back)
            nTrials = numel(ftData.trial);
            ftData.trialinfo = repmat(condLabels(c), nTrials, 1);

            ftDataAll{end+1} = ftData;
        end

        % If no data was collected for this block, skip it
        if isempty(ftDataAll)
            warning('No valid data for block %d.', block);
            continue;
        end

        % Combine data from all conditions for this block
        cfg = [];
        cfg.keepsampleinfo = 'no';
        try
            dataCombined = ft_appenddata(cfg, ftDataAll{:});
        catch ME
            warning('Data combination failed for block %d: %s', block, ME.message);
            continue;
        end

        % Perform frequency analysis to compute power in the alpha band (8–14 Hz)
        cfg = [];
        cfg.output     = 'pow';
        cfg.method     = 'mtmfft';
        cfg.taper      = 'dpss';
        cfg.tapsmofrq  = 1;          % Smoothing frequency (Hz)
        cfg.foilim     = [8 14];     % Alpha band range
        cfg.keeptrials = 'yes';
        cfg.pad        = 10;
        try
            freqData = ft_freqanalysis(cfg, dataCombined);
        catch ME
            warning('Frequency analysis failed for block %d: %s', block, ME.message);
            continue;
        end

        % Identify occipital channels (assumes labels contain 'O' or 'I')
        occIdx = find(~cellfun(@isempty, regexp(freqData.label, 'O|I', 'once')));
        if isempty(occIdx)
            warning('No occipital channels found in block %d.', block);
            continue;
        end

        % Compute mean alpha power for each trial (averaging over occipital channels and frequency bins)
        % Dimensions of freqData.powspctrm: trials x channels x frequency bins
        trialAlpha = squeeze(mean(mean(freqData.powspctrm(:, occIdx, :), 3), 2));

        % Save the results for this block: alpha power and the trial's condition label
        blockAlpha{block}.alpha   = trialAlpha;
        blockAlpha{block}.cond    = dataCombined.trialinfo;  % Numeric condition labels (1, 2, 3)
        blockAlpha{block}.nTrials = numel(trialAlpha);

        fprintf('Block %d processed: %d trials.\n', block, numel(trialAlpha));
    end

    %% Plotting the results
    % Create a bar plot of alpha power over trials (grouped by block), colouring bars
    % according to the condition.
    figure;
    set(gcf, 'Position', [0 0 1800 1200], 'Color', 'W')
    hold on;
    % Define colours for the conditions (blue for 1‑back, green for 2‑back, red for 3‑back)
    colors = color_def('AOC');
    condColors = [colors(1, :); colors(2, :); colors(3, :)];
    xTickPositions = [];
    xTickLabels = {};
    xPos = 0;  % Initialise x position

    for block = 1:nBlocks
        if isempty(blockAlpha{block})
            continue;
        end
        nTrials = blockAlpha{block}.nTrials;
        alphaVals = blockAlpha{block}.alpha;
        condVals  = blockAlpha{block}.cond;  % Condition for each trial (1, 2 or 3)
        % Determine x positions for the current block's trials
        xBlock = xPos + (1:nTrials);

        % Plot each trial as an individual bar
        for t = 1:nTrials
            thisCond = condVals(t);  % Numeric condition label
            bar(xBlock(t), alphaVals(t), 0.8, 'FaceColor', condColors(thisCond, :), 'EdgeColor', 'none');
        end

        % Store x-tick position (centre of the current block)
        xTickPositions = [xTickPositions, xPos + nTrials/2];
        xTickLabels{end+1} = sprintf('Block %d', block);
        % Advance xPos for the next block (adding a gap of 2 units)
        xPos = xBlock(end) + 2;
    end

    xlabel('Trial & Block');
    ylabel('Alpha Power');
    title(sprintf('N-back Alpha Power over Trials for Subject %s', subjectID));
    set(gca, 'XTick', xTickPositions, 'XTickLabel', xTickLabels);

    % Create a legend for the condition colours
    h1 = bar(nan, nan, 0.8, 'FaceColor', condColors(1, :));
    h2 = bar(nan, nan, 0.8, 'FaceColor', condColors(2, :));
    h3 = bar(nan, nan, 0.8, 'FaceColor', condColors(3, :));
    legend([h1, h2, h3], {'1-back', '2-back', '3-back'}, 'Location', 'northeast');
    set(gca, 'FontSize', 25)

    % Save
    allSubjBlockAlpha{subs} = blockAlpha;
    if ispc == 1
        saveas(gcf, ['W:\Students\Arne\AOC\figures\eeg\alpha_over_trials\AOC_alpha_power_over_trials_nback_sub', subjectID, '.png'])
    else
        saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_over_trials/AOC_alpha_power_over_trials_nback_sub', subjectID, '.png'])
    end
    fprintf('Subject %s done. \n', subjectID);
end

%% Compute and Visualise Grand Average Across Subjects
% Each blockAlpha{block} is a structure with fields:
%   .alpha   : vector of alpha power values for each trial in that block
%   .cond    : numeric condition labels (1 for 1‐back, 2 for 2‐back, 3 for 3‐back)
%   .nTrials : number of trials in the block

nSubjects = length(allSubjBlockAlpha);  % Total number of subjects processed
nBlocks   = 6;    % Number of blocks (as in your original script)
nConds    = 3;    % Three conditions: 1‑back, 2‑back, 3‑back

% Preallocate a 3D array to store the mean alpha power per subject, block and condition.
% Dimensions: subjects x blocks x conditions
meanAlphaPerSubj = nan(nSubjects, nBlocks, nConds);

% Loop over subjects and blocks to compute the mean alpha power per condition
for s = 1:nSubjects
    blockAlpha = allSubjBlockAlpha{s};
    for block = 1:nBlocks
        if ~isempty(blockAlpha{block})
            alphaVals = blockAlpha{block}.alpha;
            condVals  = blockAlpha{block}.cond;
            for cond = 1:nConds
                idx = find(condVals == cond);
                if ~isempty(idx)
                    meanAlphaPerSubj(s, block, cond) = mean(alphaVals(idx));
                end
            end
        end
    end
end

% Compute the grand average and standard error (SEM) across subjects for each block and condition
grandMean = squeeze(nanmean(meanAlphaPerSubj, 1));  % Size: nBlocks x nConds
grandSEM  = nan(nBlocks, nConds);
for block = 1:nBlocks
    for cond = 1:nConds
        data = squeeze(meanAlphaPerSubj(:, block, cond));
        data = data(~isnan(data));  % Exclude subjects missing data for this block/condition
        if ~isempty(data)
            grandSEM(block, cond) = std(data) / sqrt(length(data));
        end
    end
end

% Visualise Grand Average
figure;
b = bar(grandMean);  % Grouped bar plot (each group = a block, each bar = a condition)
hold on;

% Determine the number of groups (blocks) and number of bars per group (conditions)
ngroups = size(grandMean, 1);
nbars   = size(grandMean, 2);
groupwidth = min(0.8, nbars/(nbars+1.5));

% Add error bars to each bar
for i = 1:nbars
    % Calculate centre of each bar within the group
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, grandMean(:, i), grandSEM(:, i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
end

set(gca, 'XTick', 1:nBlocks, 'XTickLabel', arrayfun(@(x) sprintf('Block %d', x), 1:nBlocks, 'UniformOutput', false));
xlabel('Block');
ylabel('Mean Alpha Power');
title('Grand Average Alpha Power by Block and Condition');
legend({'1-back', '2-back', '3-back'}, 'Location', 'northeast');
set(gca, 'FontSize', 14);
if ispc == 1
    saveas(gcf, 'W:\Students\Arne\AOC\figures\eeg\alpha_over_trials\AOC_alpha_power_over_trials_allsubs_nback.png')
else
    saveas(gcf, '/Volumes/methlab/Students/Arne/AOC/figures/eeg/alpha_over_trials/AOC_alpha_power_over_trials_allsubs_nback.png')
end
