%% AOC EEG Lateralization — Sternberg (By-Trial Scatter)
% Loads eeg_matrix_sternberg_trials, plots lateralization vs trial per subject. Saves PNG per subject.
%
% Key outputs:
%   AOC_eeg_lateralization_by_trial_sternberg_subj_*.png

clear
clc
close all
path = '/Volumes/methlab/Students/Arne/AOC/data/automagic_nohp';
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

%% Load data
load('/Volumes/methlab/Students/Arne/AOC/data/features/eeg_matrix_sternberg_trials.mat');

%% Visualize Lateralization
for subs = 1:length(subjects)
    close all
    figure;
    set(gcf, "Position", [0 0 2000 1400], "Color", "W")
    subNum = str2double(subjects(subs));
    subDat = eeg_data_sternberg_trials([eeg_data_sternberg_trials.ID] == subNum);
    scatter([subDat.Lateralisation], [subDat.Trial], 150, 'k', 'filled')
    ylabel('Trial')
    xlabel('Alpha Power Lateralization [%]')
    ylim([0 300])
    xlim([-1 1])
    xline(0, '--', 'LineWidth', 2)
    text(-0.95, 275, 'LEFT', 'FontSize', 40, 'FontWeight', 'bold')
    text(0.8, 275, 'RIGHT', 'FontSize', 40, 'FontWeight', 'bold')
    set(gca, 'FontSize', 25)
    title(['Subj ',  num2str(subNum), ' Sternberg Alpha Power Lateralization over Trials | Avg: ', num2str(round(nanmean([subDat.Lateralisation]), 3))], 'FontSize', 35)
    saveas(gcf, ['/Volumes/methlab/Students/Arne/AOC/figures/eeg/lateralization/AOC_eeg_lateralization_by_trial_sternberg_subj_', num2str(subNum),'.png'])
end