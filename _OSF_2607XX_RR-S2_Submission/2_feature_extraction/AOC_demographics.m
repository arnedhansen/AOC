%% AOC Demographics
% Computes participant summary statistics from anonymized participant metadata.
% Output: Console summary of N, age mean and SD, female percent, and right-handed percent.

startup
[~, paths, ~, ~] = setup('AOC', 0);

if ~isfile(paths.vp_table)
    error('AOC_demographics:MissingVPTable', 'Anonymized participant table not found: %s', paths.vp_table);
end

dat = readtable(paths.vp_table);
requiredCols = {'ID', 'Sex', 'Age', 'Handedness'};
missingCols = requiredCols(~ismember(requiredCols, dat.Properties.VariableNames));
if ~isempty(missingCols)
    error('AOC_demographics:MissingColumns', 'Missing required columns in anonymized table: %s', strjoin(missingCols, ', '));
end

ages = dat.Age;
sex = dat.Sex;
handedness = dat.Handedness;

mean_age = mean(ages, 'omitnan');
sd_age = std(ages, 'omitnan');

n_female = sum(strcmpi(string(sex), 'W'));
perc_female = (n_female / numel(genders)) * 100;

n_right = sum(strcmpi(string(handedness), 'R'));
perc_right = (n_right / numel(handedness)) * 100;

fprintf('N: %.f participants\n', numel(ages));
fprintf('Mean age: %.2f years\n', mean_age);
fprintf('SD age: %.2f years\n', sd_age);
fprintf('Female: %.1f%%\n', perc_female);
fprintf('Right-handed: %.1f%%\n', perc_right);
