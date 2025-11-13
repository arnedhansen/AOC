%% AOC Sternberg — Split participants by PupilSize

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');
%subjects = subjects(1:10)

% Figure config
color_map = cbrewer('seq', 'Reds', 64);
colors = color_def('AOC');
fontSize = 30;

%% Load data
load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_sternberg_trials.mat')

