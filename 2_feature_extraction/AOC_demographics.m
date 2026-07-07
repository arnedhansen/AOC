%% AOC Demographics
% Compute participant summary statistics from merged N back demographics.
% Output: Console summary of N, age mean and SD, female percent, and right handed percent.

startup
[~, paths, ~, ~] = setup('AOC', 0);
featPath = paths.features;
load(fullfile(featPath, 'AOC_merged_data_nback.mat'));
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

% Display summary
clc
fprintf('N: %.f participants \n', numel(ages));
fprintf('Mean age: %.2f years\n', mean_age);
fprintf('SD age: %.2f years\n', sd_age);
fprintf('Female: %.1f%%\n', perc_female);
fprintf('Right-handed: %.1f%%\n', perc_right);
