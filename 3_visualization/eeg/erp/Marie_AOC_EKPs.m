%% Marie AOC EKPs
%
%  Visualization of AOC N-back and Sternberg EKPs
%
%  P300 (P3b)
%  Function: Reflects stimulus evaluation and context updating.
%  Modulation: Reduced amplitude with increasing load, as working memory demands increase.

%% Setup
clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/merged';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};
subjects = exclude_subjects(subjects, 'AOC');
subjects = subjects(1:5);

%% Initialize EEGLab
addpath('eeglab2021.1')
eeglab
close all
clc

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(path, filesep, subjects{subj});
    cd(datapath)

    %% Read blocks
    for block = 1:6
        try % Do not load emtpy blocks
            load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
            alleeg{block} = EEG;
            clear EEG
            fprintf('Subject %.1s: Block %.1d loaded \n', subjects{subj}, block)
        catch ME
            ME.message
            disp(['ERROR loading Block ' num2str(block) '!'])
        end
    end

    %% Segment data into epochs -2s before and 3.5s after stim onset and
    %  convert to Fieldtrip data structure
    epoch_window = [-2 3.5];
    for block = 1:6
        % PRESENTATION2 = 22 (Trigger for stimuli presentation (setSize = 2))
        % PRESENTATION4 = 24 (Trigger for stimuli presentation (setSize = 4))
        % PRESENTATION6 = 26 (Trigger for stimuli presentation (setSize = 6))

        %% Segment WM load 2 data
        try
            EEGload2 = pop_epoch(alleeg{block}, {'22'}, epoch_window);
            data2{block} = eeglab2fieldtrip(EEGload2, 'raw');
            fprintf('WM load 2 data processed for block %d.\n', block);
        catch ME
            fprintf('ERROR processing WM load 2 for block %d: %s\n', block, ME.message);
            data2{block} = struct;
        end

        %% Segment WM load 4 data
        try
            EEGload4 = pop_epoch(alleeg{block}, {'24'}, epoch_window);
            data4{block} = eeglab2fieldtrip(EEGload4, 'raw');
            fprintf('WM load 4 data processed for block %d.\n', block);
        catch ME
            fprintf('ERROR processing WM load 4 for block %d: %s\n', block, ME.message);
            data4{block} = struct;
        end

        %% Segment WM load 6 data
        try
            EEGload6 = pop_epoch(alleeg{block}, {'26'}, epoch_window);
            data6{block} = eeglab2fieldtrip(EEGload6, 'raw');
            fprintf('WM load 6 data processed for block %d.\n', block);
        catch ME
            fprintf('ERROR processing WM load 6 for block %d: %s\n', block, ME.message);
            data6{block} = struct;
        end
    end
end

%% Iterate over all subjects
clear ERP1 ERP2 ERP3 ERP1_65 ERP2_65 ERP3_65
for sub = 1 : size(d, 1) 
    
    disp(sub)
    
    % load EEG data
    load(fullfile(d(sub).folder, d(sub).name))

    % epoch
    % OUTEEG = pop_epoch( EEG, events, timelimits);
    EEGC1 = pop_epoch(EEG, {5, 6, 7}, [-0.2 0.8]); % familiar
    EEGC2 = pop_epoch(EEG, {13, 14, 15}, [-0.2 0.8]); % unfamiliar
    EEGC3 = pop_epoch(EEG, {17, 18, 19}, [-0.2 0.8]); % scrambled

    % basline correction
    EEGC1 = pop_rmbase(EEGC1, [-200 0]);
    EEGC2 = pop_rmbase(EEGC2, [-200 0]);
    EEGC3 = pop_rmbase(EEGC3, [-200 0]);

    % compute means: cond 1
    ERP1(sub, :, :) = mean(EEGC1.data(:, :, :), 3);
    ERP2(sub, :, :) = mean(EEGC2.data(:, :, :), 3);
    ERP3(sub, :, :) = mean(EEGC3.data(:, :, :), 3);
    
    % electrode 65
    i65 = find(ismember({EEGC1.chanlocs.labels}, 'EEG065'));

    ERP1_65(sub, :) = mean(EEGC1.data(i65, :, :), 3);
    ERP2_65(sub, :) = mean(EEGC2.data(i65, :, :), 3);
    ERP3_65(sub, :) = mean(EEGC3.data(i65, :, :), 3);

end
toc

%% Plotten
% variability between subjects
figure;
subplot(3, 1, 1)
plot(EEGC1.times, ERP1_65)
title('Familiar')
subplot(3, 1, 2)
plot(EEGC1.times, ERP2_65)
title('Unfamiliar')
subplot(3, 1, 3)
plot(EEGC1.times, ERP3_65)
title('Scrambled')
xlabel('Time (ms)')
sgtitle('Variability between subjects across conditions')

% grand average across conditions
GM1 = squeeze(mean(ERP1, 1));
GM2 = squeeze(mean(ERP2, 1));
GM3 = squeeze(mean(ERP3, 1));

% EKP plotten using GM
figure;
plot(EEGC1.times, GM1(61, :), 'LineWidth', 1, 'Color', 'red') % try setting LineWidth to 2
hold on
plot(EEGC1.times, GM2(61, :), 'Color', 'blue')
plot(EEGC1.times, GM3(61, :))
xlabel('Time (ms)')
ylabel('Voltage')
set(gca, 'FontSize', 16)
set(gcf, 'Color', 'w')
legend('familiar', 'unfamiliar', 'scrambled')

% EKP plotten EEG065
figure;
plot(EEGC1.times, mean(ERP1_65, 1), 'LineWidth', 2, 'Color', 'red')
hold on
plot(EEGC1.times, mean(ERP2_65, 1), 'Color', 'blue')
plot(EEGC1.times, mean(ERP3_65, 1))
xlabel('Time (ms)')
ylabel('Voltage')
set(gca, 'FontSize', 16)
legend('familiar', 'unfamiliar', 'scrambled')

%% Topography
figure; 
topoplot([],EEGC1.chanlocs, 'style','blank','electrodes','labels'); % blank topo

a = find(EEGC1.times == 152, 1) % find time points between 152ms 
b = find(EEGC1.times == 212, 1) % and 212 ms (time window around N170 peak)

figure; 
topoplot(mean(GM1(:, a:b), 2),EEGC1.chanlocs,'electrodes','on'); % mean(GM1(:, a:b), 2) => average over time window = [152 : 212]
colorbar;
% für extra punkte: colorbar zu red-white-blue ändern

%% shaded error bar
figure;
% Input: 1. time, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
shadedErrorBar(EEGC1.times, mean(ERP1_65, 1), std(ERP1_65)/sqrt(size(ERP1_65, 1)),'lineprops','b', 'patchSaturation', 0.05);
hold on
shadedErrorBar(EEGC2.times, mean(ERP2_65, 1), std(ERP2_65)/sqrt(size(ERP2_65, 1)), 'lineprops','r', 'patchSaturation', 0.05);
shadedErrorBar(EEGC3.times, mean(ERP3_65, 1), std(ERP3_65)/sqrt(size(ERP3_65, 1)), 'lineprops','g', 'patchSaturation', 0.05);
alpha(0.3)
legend('familiar', 'unfamiliar', 'scrambled')
xlabel('Time (ms)')
ylabel('Voltage')
set(gcf, 'color', 'w')
set(gca, 'FontSize', 16)

%% FUNCTIONS - DO NOT CHANGE
% Function to exclude subjects from rtf file
function filteredSubjects = exclude_subjects(subjects, project)
if strcmp(project, 'AOC')
    % Open and read the RTF file
    if ispc == 1
        fileID = fopen('W:\Students\Arne\AOC\data\controls\AOC_exclusion_participants.rtf', 'r');
    else
        fileID = fopen('/Volumes/methlab/Students/Arne/AOC/data/controls/AOC_exclusion_participants.rtf', 'r');
    end

    % Read the content as a raw string
    rawText = fread(fileID, '*char')';

    % Close the file
    fclose(fileID);

    % Extract numbers from the text (assuming numbers are separated by spaces or new lines)
    exclusionSubjects = regexp(rawText, '\d+', 'match'); % Extracts numeric strings
    exclusionSubjects = str2double(exclusionSubjects); % Convert to numeric array
    subjectsNumeric = str2double(subjects); % Convert cell array to numeric

    % Exclude subjects
    filteredSubjectsList = setdiff(subjectsNumeric, exclusionSubjects);

    % Convert back to cell array of strings
    filteredSubjects = cellstr(string(filteredSubjectsList));

    % Display the updated subjects list
    disp('Filtered Subjects:');
    disp(filteredSubjectsList');
end
end