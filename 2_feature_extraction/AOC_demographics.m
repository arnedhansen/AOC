% AOC Demographics
clear

% demos = readtable('/Volumes/methlab_vp/OCC/AOC/AOC_VPs.xlsx');
% demos = demos(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'Ausbildung_1_6_Andere_', 'OcularDominance'});
% demos = table2struct(demos(1:120, :));
% histogram([demos.Alter])

load('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback.mat');
dat = merged_data_nback;

%% Basic participant characteristics

% Extract vectors
ages       = [dat.Age];
genders    = {dat.Gender};
handedness = {dat.Handedness};

% Mean and SD age
mean_age = mean(ages, 'omitnan');
sd_age   = std(ages, 'omitnan');

% Percentage female (assuming 'W' indicates female)
n_female = sum(strcmp(genders, 'W'));
perc_female = (n_female / numel(genders)) * 100;

% Percentage right-handed (assuming 'R' indicates right-handed)
n_right = sum(strcmp(handedness, 'R'));
perc_right = (n_right / numel(handedness)) * 100;

% Display
clc
fprintf('N: %.f participants \n', numel(ages));
fprintf('Mean age: %.2f years\n', mean_age);
fprintf('SD age: %.2f years\n', sd_age);
fprintf('Female: %.1f%%\n', perc_female);
fprintf('Right-handed: %.1f%%\n', perc_right);
