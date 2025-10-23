%% AOC Sternberg Fivos Blink Visualization
% Aggregates blink times relative to stimulus onsets across subjects/blocks

%% Setup
startup
clear
[~, ~, ~ , ~] = setup('AOC');
path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/merged/';

dirs     = dir(path);
folders  = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');
% subjects = subjects(1:20)

%% Parameters
twin          = [-1 5];          % window around stimulus [s]
binWidth      = 0.05;            % histogram bin width [s]
edges         = twin(1):binWidth:twin(2);
centres       = edges(1:end-1) + diff(edges)/2;
 
blinkTypeStr  = 'L_blink';       % blink event label as stored in EEG.event.type
stimTypesNum  = [22 24 26];      % stimulus onsets
probeTypesNum = [4 5];           % probe onsets (occurring ~+3 s)
probeLatencyS = 3.0;             % nominal probe at +3 s after stimulus
 
colors        = color_def('AOC');
fontSize      = 25;

%% Read data and collect blink times relative to stimulus
tic
allBlinkRelTimes      = [];
perSubjHistCounts     = nan(numel(subjects), numel(centres));
perSubjTrialCounts    = nan(numel(subjects), 1);
perSubjHasData        = false(numel(subjects),1);

for subj = 1:length(subjects)
    subjID  = subjects{subj};
    datapath = strcat(path, subjID);
    cd(datapath)
    clc; disp(['Loading ET data for Subject AOC ', subjID])

    subj_relTimes = [];
    subj_stimN    = 0;

    for block = 1:6
        try
            load(strcat(subjID, '_EEG_ET_Sternberg_block', num2str(block), '_merged.mat'), 'EEG');
        catch
            continue
        end

        if ~exist('EEG','var') || ~isfield(EEG,'event') || isempty(EEG.event)
            clear EEG
            continue
        end

        srate = EEG.srate;

        % Extract event types and latencies
        ev = EEG.event;
        nEv = numel(ev);
        evTypeStr = cell(nEv,1);
        evTypeNum = nan(nEv,1);
        evLatSamp = nan(nEv,1);

        for e = 1:nEv
            t = [];
            if isfield(ev(e),'type')
                t = ev(e).type;
            end
            % capture both string and numeric encodings
            if isstring(t) || ischar(t)
                evTypeStr{e} = char(t);
                numTry = str2double(evTypeStr{e});
                if ~isnan(numTry)
                    evTypeNum(e) = numTry;
                end
            elseif isnumeric(t)
                evTypeNum(e) = double(t);
                evTypeStr{e} = num2str(evTypeNum(e));
            end

            if isfield(ev(e),'latency') && ~isempty(ev(e).latency)
                evLatSamp(e) = double(ev(e).latency);
            end
        end

        % Find stimulus and blink events
        stimMask  = ismember(evTypeNum, stimTypesNum);
        stimLat   = evLatSamp(stimMask);
        stimLat   = stimLat(isfinite(stimLat));

        blinkMask = strcmp(evTypeStr, blinkTypeStr);
        blinkLat  = evLatSamp(blinkMask);
        blinkLat  = blinkLat(isfinite(blinkLat));

        if isempty(stimLat) || isempty(blinkLat)
            clear EEG
            continue
        end

        preSamp  = round(twin(1) * srate);
        postSamp = round(twin(2) * srate);

        % Iterate stimulus events and collect blink times in window
        for si = 1:numel(stimLat)
            s0 = stimLat(si);
            subj_stimN = subj_stimN + 1;

            inWin = blinkLat >= (s0 + preSamp) & blinkLat <= (s0 + postSamp);
            theseBlinks = blinkLat(inWin);

            if ~isempty(theseBlinks)
                relTimes = (theseBlinks - s0) / srate;   % seconds
                relTimes = relTimes(:)';                 % force 1×N row vector
                relTimes(relTimes < twin(1) | relTimes > twin(2)) = [];
                subj_relTimes = [subj_relTimes, relTimes]; %#ok<AGROW>
            end
        end


        clear EEG
    end

    if ~isempty(subj_relTimes)
        perSubjHasData(subj) = true;
        allBlinkRelTimes     = [allBlinkRelTimes, subj_relTimes]; %#ok<AGROW>
        perSubjHistCounts(subj, :) = histcounts(subj_relTimes, edges);
        perSubjTrialCounts(subj)   = subj_stimN;
    end
end
toc

%% Aggregate across subjects
validSubs   = find(perSubjHasData);
nValid      = numel(validSubs);

aggCounts   = histcounts(allBlinkRelTimes, edges);

% Per-subject normalisation: counts per trial per bin
countsPerTrial = perSubjHistCounts ./ perSubjTrialCounts;
countsPerTrial(~isfinite(countsPerTrial)) = 0;

meanPerSubj  = nanmean(countsPerTrial(validSubs, :), 1);
semPerSubj   = nanstd(countsPerTrial(validSubs, :), [], 1) ./ sqrt(nValid);

%% Visualisation
close all
figure('Color','w','Position',[0 0 2000 1200])

% Aggregate counts
subplot(2,1,1); hold on
bar(centres, aggCounts, 'BarWidth', 1.0, 'FaceColor', colors(2,:), 'EdgeColor', 'none')
xline(0,'k--','LineWidth',2)          % stimulus
xline(probeLatencyS,'k--','LineWidth',2)  % nominal probe at +3 s
xlim(twin)
text(-0.15, -50, 'STIMULUS', 'FontSize', 15, 'FontWeight', 'bold')
text(2.9, -50, 'PROBE', 'FontSize', 15, 'FontWeight', 'bold')
text(3.25, 330, ['N = ', num2str(nValid), ' Subjects'], 'FontSize', 15, 'FontWeight', 'bold')
xlabel('Time [s]')
ylabel('Pooled Blink Count')
title('Aggregated Blinks')
set(gca,'FontSize',fontSize)

% Per-subject mean ± SEM (counts/trial/bin)
subplot(2,1,2); hold on
bar(centres, meanPerSubj, 'BarWidth', 1.0, 'FaceColor', colors(2,:), 'EdgeColor', 'none')
plot(centres, meanPerSubj + semPerSubj, '-', 'LineWidth', 1.5, 'Color', 'k')
plot(centres, meanPerSubj - semPerSubj, '-', 'LineWidth', 1.5, 'Color', 'k')
xline(0,'k--','LineWidth',2)
xline(probeLatencyS,'k--','LineWidth',2)
xlim(twin)
text(-0.15, -0.0125, 'STIMULUS', 'FontSize', 15, 'FontWeight', 'bold')
text(2.9, -0.0125, 'PROBE', 'FontSize', 15, 'FontWeight', 'bold')
xlabel('Time [s]')
ylabel('Averaged Binned Blinks')
title('Averaged Blinks')
set(gca,'FontSize',fontSize)

% Save
saveas(gcf, '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/AOC_sternberg_blink_hist.png');