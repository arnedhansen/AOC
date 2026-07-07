%% AOC Gaze Feature Extraction — N-Back
%
% Extracted features:

%% Setup 
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;

gaze_data_nback = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeDeviationFullBL', {}, 'GazeDeviationEarlyBL', {}, 'GazeDeviationLateBL', {}, ...
    'MSRateFullBL', {}, 'MSRateEarlyBL', {}, 'MSRateLateBL', {}, ...

%% Load all eye movements
for subj = 1:length(subjects)
    clc; fprintf('[GAZE FEX - NBACK] Gaze feature extraction for Subject %d / %d \n', subj, length(subjects))
    datapath = fullfile(path, subjects{subj}, 'gaze');
    load([datapath, filesep, 'dataET_nback'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = size(dataETlong.trialinfo, 1);
    condition = [];
    gazeDev = [];
    microsaccadeRate = [];
    scanPathLength = [];
    % Baselined (% change or delta) gaze metrics
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
        idx_base  = t >= -1.5 & t <= -0.5;
        idx_early = t >= 0    & t <= 1;
        idx_late  = t >= 1    & t <= 2;
        idx_full  = t >= 0    & t <= 2;
        data(2, :) = 600 - data(2, :); % Invert Y-axis

        %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        gaze_x{subj, trl} = data(1, :);
        gaze_y{subj, trl} = data(2, :);

        %% Compute gaze deviation as euclidean distances from the center
        x_coords = gaze_x{subj, trl};
        y_coords = gaze_y{subj, trl};
        x_full = data(1, idx_full);
        y_full = data(2, idx_full);

        % Calculate Euclidean distances
        dx = x_full - 400; % Distance from center of the screen
        dy = y_full - 300; % Distance from center of the screen
        gaze_euclidean_dev = sqrt(dx.^2 + dy.^2);

        % Calculate the mean Euclidean distance
        mean_euclidean_distance = nanmean(gaze_euclidean_dev);
        gaze_standard_deviation_x = nanstd(x_coords);
        gaze_standard_deviation_y = nanstd(y_coords);

        %% Compute microsaccades
        fsample = 500; % Sample rate of 500 Hz
        velData = [x_full; y_full]; % Compute raw MS rate in the same 0-2 s window
        trlLength = size(velData, 2);
        [microsaccade_rate, microsaccade_details] = detect_microsaccades(fsample, velData, trlLength);
        ms_data(trl) = microsaccade_details;

        %% Compute Scan Path Length
        dxf_s = diff(x_full);
        dyf_s = diff(y_full);

        if numel(x_bv) >= 10
            sx_b = std(x_bv); sy_b = std(y_bv);
            rho_b = corr(x_bv(:), y_bv(:));
            k95 = -log(1 - 0.95);
        else
        end

        %% Baselined metrics from per-trial baseline window
        fsample = 500;

        % baseline
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

        % helper inline for each active window
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

    l1 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 1);
    l1gdev = mean([l1.GazeDeviation], 'omitnan');
    l1msrate = mean([l1.MSRate], 'omitnan');

    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');

    l3 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 3);
    l3gdev = mean([l3.GazeDeviation], 'omitnan');
    l3msrate = mean([l3.MSRate], 'omitnan');

    %% Load gaze metrics (extracted in GCP_preprocessing.m)
    load([datapath, filesep, 'gaze_metrics_nback'])

    % Aggregate baselined trial-level values per condition
    l1gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3gdevBLf = mean(gdevFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3gdevBLe = mean(gdevEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3gdevBLl = mean(gdevLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');

    %% Create across condition structure
    subID = str2double(subjects{subj});
    subj_data_gaze = struct('ID', num2cell([subID; subID; subID]), ...
        'Condition', num2cell([1; 2; 3]), ...
        'GazeDeviation', num2cell([l1gdev; l2gdev; l3gdev]), ...
        'MSRate', num2cell([l1msrate; l2msrate; l3msrate]), ...
        'Blinks', num2cell([blinks_1back; blinks_2back; blinks_3back]), ...
        'Fixations', num2cell([fixations_1back; fixations_2back; fixations_3back]), ...
        'Saccades', num2cell([saccades_1back; saccades_2back; saccades_3back]), ...
        'GazeDeviationFullBL', num2cell([l1gdevBLf; l2gdevBLf; l3gdevBLf]), ...
        'GazeDeviationEarlyBL', num2cell([l1gdevBLe; l2gdevBLe; l3gdevBLe]), ...
        'GazeDeviationLateBL', num2cell([l1gdevBLl; l2gdevBLl; l3gdevBLl]), ...
        'MSRateFullBL', num2cell([l1msBLf; l2msBLf; l3msBLf]), ...
        'MSRateEarlyBL', num2cell([l1msBLe; l2msBLe; l3msBLe]), ...
        'MSRateLateBL', num2cell([l1msBLl; l2msBLl; l3msBLl]), ...

    %% Save
    savepath = fullfile(paths.features, subjects{subj}, 'gaze');
    if ~isfolder(savepath)
        mkdir(savepath)
    end
    cd(savepath)
    save gaze_matrix_nback_trial subj_data_gaze_trial
    save gaze_matrix_nback subj_data_gaze
    save gaze_dev_nback l1gdev l2gdev l3gdev
    save gaze_std_nback l1gSDx l2gSDx l3gSDx l1gSDy l2gSDy l3gSDy
    save ms_rate_nback l1msrate l2msrate l3msrate
    
    % Append to the final structure array
    gaze_data_nback = [gaze_data_nback; subj_data_gaze];
end

save(fullfile(paths.features, 'AOC_gaze_nback.mat'), 'gaze_x', 'gaze_y')
save(fullfile(paths.features, 'AOC_gaze_matrix_nback.mat'), 'gaze_data_nback')
