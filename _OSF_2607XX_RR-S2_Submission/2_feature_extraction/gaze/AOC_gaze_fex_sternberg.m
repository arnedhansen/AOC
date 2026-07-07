%% AOC Gaze Feature Extraction — Sternberg
%
% Extracted features:

%% Setup 
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;

gaze_data_sternberg = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeDeviationFullBL', {}, 'GazeDeviationEarlyBL', {}, 'GazeDeviationLateBL', {}, ...
    'MSRateFullBL', {}, 'MSRateEarlyBL', {}, 'MSRateLateBL', {}, ...

%% Load all eye movements
for subj = 1:length(subjects)
    clc; fprintf('[GAZE FEX - STERNBERG] Gaze feature extraction for Subject %d / %d \n', subj, length(subjects))
    datapath = fullfile(path, subjects{subj}, 'gaze');
    load([datapath, filesep, 'dataET_sternberg'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = size(dataETlong.trialinfo, 1);
    condition = [];
    gazeDev = [];
    microsaccadeRate = [];
    scanPathLength = [];
    % Baselined (% change) gaze metrics
    gdevEarlyBL = []; gdevLateBL = []; gdevFullBL = [];
    msEarlyBL = [];   msLateBL = [];   msFullBL = [];

    %% Get trial-by-trial gaze data
    for trl = 1:size(dataETlong.trialinfo, 1)
        close all
        data = dataETlong.trial{trl};

        %% Define windows
        t = dataETlong.time{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices);
        data = valid_data;
        t = t(valid_data_indices);
        idx_base   = t >= -1.5 & t <= -0.5;
        idx_early  = t >= 0    & t <= 1;
        idx_late   = t >= 1    & t <= 2;
        idx_full   = t >= 0    & t <= 2;
        data(2, :) = 600 - data(2, :); % Invert Y-axis

        %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        data_late = data(:, idx_late);
        gaze_x{subj, trl} = data_late(1, :);
        gaze_y{subj, trl} = data_late(2, :);

        %% Compute gaze deviation as euclidean distances from the center
        x_coords = gaze_x{subj, trl};
        y_coords = gaze_y{subj, trl};

        % Calculate Euclidean distances
        dx = x_coords - 400; % Distance from center of the screen
        dy = y_coords - 300; % Distance from center of the screen
        gaze_euclidean_dev = sqrt(dx.^2 + dy.^2);

        % Calculate the mean Euclidean distance
        mean_euclidean_distance = nanmean(gaze_euclidean_dev);
        gaze_standard_deviation_x = nanstd(x_coords);
        gaze_standard_deviation_y = nanstd(y_coords);

        %% Compute microsaccades
        fsample = 500; % Sample rate of 500 Hz
        velData = [gaze_x{subj, trl}; gaze_y{subj, trl}]; % Concatenate x and y gaze coordinates to compute the velocity of eye movements in a 2D space
        trlLength = size(velData, 2);
        [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength);
        ms_data(trl) = microsaccade_details;

        %% Compute Scan Path Length
        dxf_s = diff(x_coords);
        dyf_s = diff(y_coords);

        if numel(x_bv) >= 10
            sx_b = std(x_bv); sy_b = std(y_bv);
            rho_b = corr(x_bv(:), y_bv(:));
            k95 = -log(1 - 0.95);
        else
        end

        %% Baselined metrics from per-trial baseline window
        fsample = 500;

        xb = data(1, idx_base); yb = data(2, idx_base); pb = data(3, idx_base);
        gd_base = mean(sqrt((xb-400).^2 + (yb-300).^2), 'omitnan');
        Tb = sum(isfinite(xb) & isfinite(yb)) / fsample;
        [~, msb] = detect_microsaccades(fsample, [xb; yb], numel(xb));
        ms_base = numel(msb.Onset) / Tb;
        if ~isfinite(ms_base) || ms_base <= 0, ms_base = NaN; end
        vb = isfinite(xb) & isfinite(yb);
        xb2 = double(xb(vb)); yb2 = double(yb(vb));
        if numel(xb2) >= 10
            k95 = -log(1 - 0.95);
        else
        end

        winDefs = {idx_early, idx_late, idx_full};
        gd_bl = nan(1,3); ms_bl = nan(1,3);
        for wi = 1:3
            xw = data(1, winDefs{wi}); yw = data(2, winDefs{wi}); pw = data(3, winDefs{wi});
            gd_w = mean(sqrt((xw-400).^2 + (yw-300).^2), 'omitnan');
            Tw = sum(isfinite(xw) & isfinite(yw)) / fsample;
            [~, msw] = detect_microsaccades(fsample, [xw; yw], numel(xw));
            ms_w = numel(msw.Onset) / Tw;
            if ~isfinite(ms_w), ms_w = NaN; end

            vw = isfinite(xw) & isfinite(yw);
            xw2 = double(xw(vw)); yw2 = double(yw(vw));
            if numel(xw2) >= 10
                k95 = -log(1 - 0.95);
            else
            end

            if isfinite(gd_w)   && isfinite(gd_base)   && gd_base > 0, gd_bl(wi)   = 100*(gd_w-gd_base)/gd_base; end
            if isfinite(ms_w)   && isfinite(ms_base)   && ms_base > 0, ms_bl(wi)   = 100*(ms_w-ms_base)/ms_base; end
        end

        %% Append data for this trial
        subject_id = [subject_id; str2double(subjects{subj})];
        trial_num = [trial_num; dataETlong.trialinfo(trl, 2)];
        condition = [condition; dataETlong.trialinfo(trl, 1)-20];
        gazeDev = [gazeDev; mean_euclidean_distance];
        microsaccadeRate = [microsaccadeRate; microsaccade_rate];
        gdevEarlyBL = [gdevEarlyBL; gd_bl(1)]; gdevLateBL = [gdevLateBL; gd_bl(2)]; gdevFullBL = [gdevFullBL; gd_bl(3)];
        msEarlyBL = [msEarlyBL; ms_bl(1)];     msLateBL = [msLateBL; ms_bl(2)];     msFullBL = [msFullBL; ms_bl(3)];

    end

    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), ...
        'Trial', num2cell(trial_num), ...
        'Condition', num2cell(condition), ...
        'GazeDeviation', num2cell(gazeDev), ...
        'MSRate', num2cell(microsaccadeRate), ...

    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');

    l4 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 4);
    l4gdev = mean([l4.GazeDeviation], 'omitnan');
    l4msrate = mean([l4.MSRate], 'omitnan');

    l6 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 6);
    l6gdev = mean([l6.GazeDeviation], 'omitnan');
    l6msrate = mean([l6.MSRate], 'omitnan');

    %% Load gaze metrics (extracted in GCP_preprocessing.m)
    load([datapath, filesep, 'gaze_metrics_sternberg'])

    % Aggregate baselined trial-level values per condition
    l2gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');

    %% Create across condition structure
    subID = str2double(subjects{subj});
    subj_data_gaze = struct('ID', num2cell([subID; subID; subID]), ...
        'Condition', num2cell([2; 4; 6]), ...
        'GazeDeviation', num2cell([l2gdev; l4gdev; l6gdev]), ...
        'MSRate', num2cell([l2msrate; l4msrate; l6msrate]), ...
        'Blinks', num2cell([blinks_l2; blinks_l4; blinks_l6]), ...
        'Fixations', num2cell([fixations_l2; fixations_l4; fixations_l6]), ...
        'Saccades', num2cell([saccades_l2; saccades_l4; saccades_l6]), ...
        'GazeDeviationFullBL', num2cell([l2gdevBLf; l4gdevBLf; l6gdevBLf]), ...
        'GazeDeviationEarlyBL', num2cell([l2gdevBLe; l4gdevBLe; l6gdevBLe]), ...
        'GazeDeviationLateBL', num2cell([l2gdevBLl; l4gdevBLl; l6gdevBLl]), ...
        'MSRateFullBL', num2cell([l2msBLf; l4msBLf; l6msBLf]), ...
        'MSRateEarlyBL', num2cell([l2msBLe; l4msBLe; l6msBLe]), ...
        'MSRateLateBL', num2cell([l2msBLl; l4msBLl; l6msBLl]), ...

    %% Save data
    savepath = fullfile(paths.features, subjects{subj}, 'gaze');
    if ~isfolder(savepath)
        mkdir(savepath)
    end
    cd(savepath)
    save gaze_matrix_sternberg_trial subj_data_gaze_trial
    save gaze_matrix_sternberg subj_data_gaze
    save gaze_dev_sternberg l2gdev l4gdev l6gdev
    save gaze_std_sternberg l2gSDx l4gSDx l6gSDx l2gSDy l4gSDy l6gSDy
    save ms_rate_sternberg l2msrate l4msrate l6msrate

    % Append to the final structure array
    gaze_data_sternberg = [gaze_data_sternberg; subj_data_gaze];
end
trialinfo = dataETlong.trialinfo';

save(fullfile(paths.features, 'AOC_gaze_sternberg.mat'), 'gaze_x', 'gaze_y', 'trialinfo')
save(fullfile(paths.features, 'AOC_gaze_matrix_sternberg.mat'), 'gaze_data_sternberg')
