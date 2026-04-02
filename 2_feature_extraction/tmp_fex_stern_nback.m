%% AOC Gaze Feature Extraction — N-Back
% Loads dataET_nback, computes gaze deviation, std, pupil, microsaccade rate, and scan-path length per trial and condition. Saccades, blinks and fixations come from AOC_preprocessing_nback. Saves gaze_matrix_nback.mat.
%
% Extracted features:
%   Gaze deviation (Euclidean), GazeStdX/Y, PupilSize, MSRate
%   Blinks, Fixations, Saccades, ScanPathLength (from preprocessing)

%% Setup 
startup
setup('AOC');
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\features\';
else
    path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

gaze_data_nback = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeStdX', {}, 'GazeStdY', {}, 'PupilSize', {}, 'MSRate', {}, ...
    'Blinks', {}, 'Fixations', {}, 'Saccades', {}, 'ScanPathLength', {}, ...
    'BCEA', {}, 'BCEALateralization', {}, ...
    'GazeDeviationFullBL', {}, 'GazeDeviationEarlyBL', {}, 'GazeDeviationLateBL', {}, ...
    'ScanPathLengthFullBL', {}, 'ScanPathLengthEarlyBL', {}, 'ScanPathLengthLateBL', {}, ...
    'PupilSizeFullBL', {}, 'PupilSizeEarlyBL', {}, 'PupilSizeLateBL', {}, ...
    'MSRateFullBL', {}, 'MSRateEarlyBL', {}, 'MSRateLateBL', {}, ...
    'BCEAFullBL', {}, 'BCEAEarlyBL', {}, 'BCEALateBL', {}, ...
    'BCEALatFullBL', {}, 'BCEALatEarlyBL', {}, 'BCEALatLateBL', {});

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, filesep, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET_nback'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = size(dataet.trialinfo, 1);
    condition = [];
    gazeDev = [];
    gazeSDx = [];
    gazeSDy = [];
    pupilSize = [];
    microsaccadeRate = [];
    scanPathLength = [];
    bcea95Area = [];
    bceaLateralization = [];
    % Baselined (% change or delta) gaze metrics
    gdevEarlyBL = []; gdevLateBL = []; gdevFullBL = [];
    splEarlyBL = [];  splLateBL = [];  splFullBL = [];
    pupEarlyBL = [];  pupLateBL = [];  pupFullBL = [];
    msEarlyBL = [];   msLateBL = [];   msFullBL = [];
    bceaEarlyBL = []; bceaLateBL = []; bceaFullBL = [];
    blatEarlyBL = []; blatLateBL = []; blatFullBL = [];

    %% Get trial-by-trial gaze data
    for trl = 1:size(dataet.trialinfo, 1)
        close all
        data = dataet.trial{trl};

        %% Define windows
        t = dataet.time{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices);
        data = valid_data;
        t = t(valid_data_indices);
        idx_base  = t >= -0.5 & t <= -0.25;
        idx_early = t >= 0    & t <= 1;
        idx_late  = t >= 1    & t <= 2;
        idx_full  = t >= 0    & t <= 2;
        data(2, :) = 600 - data(2, :); % Invert Y-axis

        %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        %% Extract gaze data and pupil size
        gaze_x{subj, trl} = data(1, :);
        gaze_y{subj, trl} = data(2, :);
        pupil_size{subj, trl} = mean(data(3, :), 'omitnan') / 1000;
        pups = pupil_size{subj, trl};

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
        spl = sum(sqrt(dxf_s.^2 + dyf_s.^2), 'omitnan');

        %% Compute BCEA (95%) and Lateralization
        valid_bcea = isfinite(x_coords) & isfinite(y_coords);
        x_bv = double(x_coords(valid_bcea)); y_bv = double(y_coords(valid_bcea));
        if numel(x_bv) >= 10
            sx_b = std(x_bv); sy_b = std(y_bv);
            rho_b = corr(x_bv(:), y_bv(:));
            k95 = -log(1 - 0.95);
            bcea = 2 * k95 * pi * sx_b * sy_b * sqrt(1 - rho_b^2);
            bcea_lat = (mean(x_bv) - 400) / 400; % -1=left, 0=centre, +1=right
        else
            bcea = NaN; bcea_lat = NaN;
        end

        %% Baselined metrics from per-trial baseline window
        fsample = 500;
        min_valid = 100;

        % baseline
        xb = data(1, idx_base); yb = data(2, idx_base); pb = data(3, idx_base);
        ok_b = sum(isfinite(xb) & isfinite(yb)) >= min_valid;
        if ok_b
            gd_base = mean(sqrt((xb-400).^2 + (yb-300).^2), 'omitnan');
            spl_base = sum(sqrt(diff(xb).^2 + diff(yb).^2), 'omitnan');
            pup_base = mean(pb, 'omitnan');
            Tb = sum(isfinite(xb) & isfinite(yb)) / fsample;
            [~, msb] = detect_microsaccades(fsample, [xb; yb], numel(xb));
            ms_base = numel(msb.Onset) / Tb;
            if ~isfinite(ms_base) || ms_base <= 0, ms_base = NaN; end
            vb = isfinite(xb) & isfinite(yb);
            xb2 = double(xb(vb)); yb2 = double(yb(vb));
            if numel(xb2) >= 10
                k95 = -log(1 - 0.95);
                bcea_base = 2*k95*pi*std(xb2)*std(yb2)*sqrt(1-corr(xb2(:), yb2(:))^2);
                blat_base = (mean(xb2) - 400) / 400;
            else
                bcea_base = NaN; blat_base = NaN;
            end
        else
            gd_base = NaN; spl_base = NaN; pup_base = NaN; ms_base = NaN; bcea_base = NaN; blat_base = NaN;
        end

        % helper inline for each active window
        winDefs = {idx_early, idx_late, idx_full};
        gd_bl = nan(1,3); spl_bl = nan(1,3); pup_bl = nan(1,3); ms_bl = nan(1,3); bcea_bl = nan(1,3); blat_bl = nan(1,3);
        for wi = 1:3
            xw = data(1, winDefs{wi}); yw = data(2, winDefs{wi}); pw = data(3, winDefs{wi});
            ok_w = sum(isfinite(xw) & isfinite(yw)) >= min_valid;
            if ~ok_w, continue; end

            gd_w = mean(sqrt((xw-400).^2 + (yw-300).^2), 'omitnan');
            spl_w = sum(sqrt(diff(xw).^2 + diff(yw).^2), 'omitnan');
            pup_w = mean(pw, 'omitnan');
            Tw = sum(isfinite(xw) & isfinite(yw)) / fsample;
            [~, msw] = detect_microsaccades(fsample, [xw; yw], numel(xw));
            ms_w = numel(msw.Onset) / Tw;
            if ~isfinite(ms_w), ms_w = NaN; end

            vw = isfinite(xw) & isfinite(yw);
            xw2 = double(xw(vw)); yw2 = double(yw(vw));
            if numel(xw2) >= 10
                k95 = -log(1 - 0.95);
                bcea_w = 2*k95*pi*std(xw2)*std(yw2)*sqrt(1-corr(xw2(:), yw2(:))^2);
                blat_w = (mean(xw2) - 400) / 400;
            else
                bcea_w = NaN; blat_w = NaN;
            end

            if isfinite(gd_w)   && isfinite(gd_base)   && gd_base > 0, gd_bl(wi)   = 100*(gd_w-gd_base)/gd_base; end
            if isfinite(spl_w)  && isfinite(spl_base)  && spl_base > 0, spl_bl(wi)  = 100*(spl_w-spl_base)/spl_base; end
            if isfinite(pup_w)  && isfinite(pup_base)  && pup_base ~= 0, pup_bl(wi)  = 100*(pup_w-pup_base)/pup_base; end
            if isfinite(ms_w)   && isfinite(ms_base)   && ms_base > 0, ms_bl(wi)   = 100*(ms_w-ms_base)/ms_base; end
            if isfinite(bcea_w) && isfinite(bcea_base) && bcea_base > 0, bcea_bl(wi) = 100*(bcea_w-bcea_base)/bcea_base; end
            if isfinite(blat_w) && isfinite(blat_base), blat_bl(wi) = blat_w - blat_base; end
        end

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; dataet.trialinfo(trl, 2)];
        condition = [condition; dataet.trialinfo(trl, 1)-20];
        gazeDev = [gazeDev; mean_euclidean_distance];
        gazeSDx = [gazeSDx; gaze_standard_deviation_x];
        gazeSDy = [gazeSDy; gaze_standard_deviation_y];
        pupilSize = [pupilSize; pups];
        microsaccadeRate = [microsaccadeRate; microsaccade_rate];
        scanPathLength = [scanPathLength; spl];
        bcea95Area = [bcea95Area; bcea];
        bceaLateralization = [bceaLateralization; bcea_lat];
        gdevEarlyBL = [gdevEarlyBL; gd_bl(1)]; gdevLateBL = [gdevLateBL; gd_bl(2)]; gdevFullBL = [gdevFullBL; gd_bl(3)];
        splEarlyBL = [splEarlyBL; spl_bl(1)];  splLateBL = [splLateBL; spl_bl(2)];  splFullBL = [splFullBL; spl_bl(3)];
        pupEarlyBL = [pupEarlyBL; pup_bl(1)];  pupLateBL = [pupLateBL; pup_bl(2)];  pupFullBL = [pupFullBL; pup_bl(3)];
        msEarlyBL = [msEarlyBL; ms_bl(1)];     msLateBL = [msLateBL; ms_bl(2)];     msFullBL = [msFullBL; ms_bl(3)];
        bceaEarlyBL = [bceaEarlyBL; bcea_bl(1)]; bceaLateBL = [bceaLateBL; bcea_bl(2)]; bceaFullBL = [bceaFullBL; bcea_bl(3)];
        blatEarlyBL = [blatEarlyBL; blat_bl(1)]; blatLateBL = [blatLateBL; blat_bl(2)]; blatFullBL = [blatFullBL; blat_bl(3)];

    end

    %% Check data by visualizing raw gaze data
    % close all
    % % Preallocate arrays for averaged gaze data
    % num_samples = 1000;
    % mean_gaze_x = nan(num_samples, 1);
    % mean_gaze_y = nan(num_samples, 1);
    % 
    % % Stack trials into matrices for averaging
    % gaze_x_matrix = cell2mat(cellfun(@(x) x(:), gaze_x(subj, :), 'UniformOutput', false)');
    % gaze_y_matrix = cell2mat(cellfun(@(y) y(:), gaze_y(subj, :), 'UniformOutput', false)');
    % 
    % % Calculate the mean over trials for each sample
    % mean_gaze_x = nanmean(gaze_x_matrix, 2);
    % mean_gaze_y = nanmean(gaze_y_matrix, 2);
    % 
    % % Plot the averaged gaze data
    % figure;
    % set(gcf, "Position", [200, 200, 1000, 600]);
    % plot(mean_gaze_x, mean_gaze_y, 'o');
    % hold on;
    % plot(400, 300, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Centre point
    % title('Averaged Gaze Data Distribution Across Samples');
    % xlabel('X Coordinates');
    % ylabel('Y Coordinates');
    % xlim([0 800]);
    % ylim([0 600]);

    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), ...
        'Trial', num2cell(trial_num), ...
        'Condition', num2cell(condition), ...
        'GazeDeviation', num2cell(gazeDev), ...
        'GazeStdX', num2cell(gazeSDx), ...
        'GazeStdY', num2cell(gazeSDy), ...
        'PupilSize', num2cell(pupilSize), ...
        'MSRate', num2cell(microsaccadeRate), ...
        'ScanPathLength', num2cell(scanPathLength), ...
        'BCEA', num2cell(bcea95Area), ...
        'BCEALateralization', num2cell(bceaLateralization));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    l1 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 1);
    l1gdev = mean([l1.GazeDeviation], 'omitnan');
    l1gSDx = mean([l1.GazeStdX], 'omitnan');
    l1gSDy = mean([l1.GazeStdY], 'omitnan');
    l1pups = mean([l1.PupilSize], 'omitnan');
    l1msrate = mean([l1.MSRate], 'omitnan');
    l1spl = mean([l1.ScanPathLength], 'omitnan');
    l1bcea = mean([l1.BCEA], 'omitnan');
    l1blat = mean([l1.BCEALateralization], 'omitnan');

    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2gSDx = mean([l2.GazeStdX], 'omitnan');
    l2gSDy = mean([l2.GazeStdY], 'omitnan');
    l2pups = mean([l2.PupilSize], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');
    l2spl = mean([l2.ScanPathLength], 'omitnan');
    l2bcea = mean([l2.BCEA], 'omitnan');
    l2blat = mean([l2.BCEALateralization], 'omitnan');

    l3 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 3);
    l3gdev = mean([l3.GazeDeviation], 'omitnan');
    l3gSDx = mean([l3.GazeStdX], 'omitnan');
    l3gSDy = mean([l3.GazeStdY], 'omitnan');
    l3pups = mean([l3.PupilSize], 'omitnan');
    l3msrate = mean([l3.MSRate], 'omitnan');
    l3spl = mean([l3.ScanPathLength], 'omitnan');
    l3bcea = mean([l3.BCEA], 'omitnan');
    l3blat = mean([l3.BCEALateralization], 'omitnan');

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
    l1splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');
    l1blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 1), 'omitnan');
    l2blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l3blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 3), 'omitnan');

    %% Create across condition structure
    subID = str2double(subjects{subj});
    subj_data_gaze = struct('ID', num2cell([subID; subID; subID]), ...
        'Condition', num2cell([1; 2; 3]), ...
        'GazeDeviation', num2cell([l1gdev; l2gdev; l3gdev]), ...
        'GazeStdX', num2cell([l1gSDx; l2gSDx; l3gSDx]), ...
        'GazeStdY', num2cell([l1gSDy; l2gSDy; l3gSDy]), ...
        'PupilSize', num2cell([l1pups; l2pups; l3pups]), ...
        'MSRate', num2cell([l1msrate; l2msrate; l3msrate]), ...
        'Blinks', num2cell([blinks_1back; blinks_2back; blinks_3back]), ...
        'Fixations', num2cell([fixations_1back; fixations_2back; fixations_3back]), ...
        'Saccades', num2cell([saccades_1back; saccades_2back; saccades_3back]), ...
        'ScanPathLength', num2cell([l1spl; l2spl; l3spl]), ...
        'BCEA', num2cell([l1bcea; l2bcea; l3bcea]), ...
        'BCEALateralization', num2cell([l1blat; l2blat; l3blat]), ...
        'GazeDeviationFullBL', num2cell([l1gdevBLf; l2gdevBLf; l3gdevBLf]), ...
        'GazeDeviationEarlyBL', num2cell([l1gdevBLe; l2gdevBLe; l3gdevBLe]), ...
        'GazeDeviationLateBL', num2cell([l1gdevBLl; l2gdevBLl; l3gdevBLl]), ...
        'ScanPathLengthFullBL', num2cell([l1splBLf; l2splBLf; l3splBLf]), ...
        'ScanPathLengthEarlyBL', num2cell([l1splBLe; l2splBLe; l3splBLe]), ...
        'ScanPathLengthLateBL', num2cell([l1splBLl; l2splBLl; l3splBLl]), ...
        'PupilSizeFullBL', num2cell([l1pupBLf; l2pupBLf; l3pupBLf]), ...
        'PupilSizeEarlyBL', num2cell([l1pupBLe; l2pupBLe; l3pupBLe]), ...
        'PupilSizeLateBL', num2cell([l1pupBLl; l2pupBLl; l3pupBLl]), ...
        'MSRateFullBL', num2cell([l1msBLf; l2msBLf; l3msBLf]), ...
        'MSRateEarlyBL', num2cell([l1msBLe; l2msBLe; l3msBLe]), ...
        'MSRateLateBL', num2cell([l1msBLl; l2msBLl; l3msBLl]), ...
        'BCEAFullBL', num2cell([l1bceaBLf; l2bceaBLf; l3bceaBLf]), ...
        'BCEAEarlyBL', num2cell([l1bceaBLe; l2bceaBLe; l3bceaBLe]), ...
        'BCEALateBL', num2cell([l1bceaBLl; l2bceaBLl; l3bceaBLl]), ...
        'BCEALatFullBL', num2cell([l1blatBLf; l2blatBLf; l3blatBLf]), ...
        'BCEALatEarlyBL', num2cell([l1blatBLe; l2blatBLe; l3blatBLe]), ...
        'BCEALatLateBL', num2cell([l1blatBLl; l2blatBLl; l3blatBLl]));

    %% Save
    if ispc == 1
        savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\gaze\');
    else
        savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
    end
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_nback_trial subj_data_gaze_trial
    save gaze_matrix_nback subj_data_gaze
    save gaze_dev_nback l1gdev l2gdev l3gdev
    save gaze_std_nback l1gSDx l2gSDx l3gSDx l1gSDy l2gSDy l3gSDy
    save pupil_size_nback l1pups l2pups l3pups
    save ms_rate_nback l1msrate l2msrate l3msrate
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_nback = [gaze_data_nback; subj_data_gaze];
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_gaze_nback gaze_x gaze_y
    save W:\Students\Arne\AOC\data\features\AOC_gaze_matrix_nback gaze_data_nback
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_gaze_nback gaze_x gaze_y
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_gaze_matrix_nback gaze_data_nback
end

%% AOC Gaze Feature Extraction — Sternberg
% Loads dataET_sternberg, computes gaze deviation, std, pupil, microsaccade rate, and scan-path length per trial and condition. Saccades, blinks and fixations come from AOC_preprocessing_sternberg. Saves gaze_matrix_sternberg.mat.
%
% Extracted features:
%   Gaze deviation (Euclidean), GazeStdX/Y, PupilSize, MSRate
%   Blinks, Fixations, Saccades, ScanPathLength (from preprocessing)

%% Setup
startup
setup('AOC');
if ispc == 1
    path = 'W:\Students\Arne\AOC\data\features\';
else
    path = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/';
end
dirs = dir(path);
folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
subjects = {folders.name};

gaze_data_sternberg = struct('ID', {}, 'Condition', {}, 'GazeDeviation', {}, ...
    'GazeStdX', {}, 'GazeStdY', {}, 'PupilSize', {}, 'MSRate', {}, ...
    'Blinks', {}, 'Fixations', {}, 'Saccades', {}, 'ScanPathLength', {}, ...
    'BCEA', {}, 'BCEALateralization', {}, ...
    'GazeDeviationFullBL', {}, 'GazeDeviationEarlyBL', {}, 'GazeDeviationLateBL', {}, ...
    'ScanPathLengthFullBL', {}, 'ScanPathLengthEarlyBL', {}, 'ScanPathLengthLateBL', {}, ...
    'PupilSizeFullBL', {}, 'PupilSizeEarlyBL', {}, 'PupilSizeLateBL', {}, ...
    'MSRateFullBL', {}, 'MSRateEarlyBL', {}, 'MSRateLateBL', {}, ...
    'BCEAFullBL', {}, 'BCEAEarlyBL', {}, 'BCEALateBL', {}, ...
    'BCEALatFullBL', {}, 'BCEALatEarlyBL', {}, 'BCEALatLateBL', {});

%% Load all eye movements
for subj = 1:length(subjects)
    datapath = strcat(path, filesep, subjects{subj}, '/gaze');
    load([datapath, filesep, 'dataET_sternberg'])

    %% Initialize arrays
    subject_id = [];
    trial_num = [];
    num_trials = size(dataet.trialinfo, 1);
    condition = [];
    gazeDev = [];
    gazeSDx = [];
    gazeSDy = [];
    pupilSize = [];
    microsaccadeRate = [];
    scanPathLength = [];
    bcea95Area = [];
    bceaLateralization = [];
    % Baselined (% change or delta) gaze metrics
    gdevEarlyBL = []; gdevLateBL = []; gdevFullBL = [];
    splEarlyBL = [];  splLateBL = [];  splFullBL = [];
    pupEarlyBL = [];  pupLateBL = [];  pupFullBL = [];
    msEarlyBL = [];   msLateBL = [];   msFullBL = [];
    bceaEarlyBL = []; bceaLateBL = []; bceaFullBL = [];
    blatEarlyBL = []; blatLateBL = []; blatFullBL = [];

    %% Get trial-by-trial gaze data
    for trl = 1:size(dataet.trialinfo, 1)
        close all
        data = dataet.trial{trl};

        %% Define windows
        t = dataet.time{trl};

        %% Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(1:3, valid_data_indices);
        data = valid_data;
        t = t(valid_data_indices);
        idx_base   = t >= -0.5 & t <= -0.25;
        idx_early  = t >= 0    & t <= 1;
        idx_late   = t >= 1    & t <= 2;
        idx_full   = t >= 0    & t <= 2;
        idx_legacy = idx_late; % keep legacy non-baselined metrics on [1 2]s
        data(2, :) = 600 - data(2, :); % Invert Y-axis

        %% Remove blinks with a window of 100ms (= 50 samples/timepoints)
        win_size = 50;
        data = remove_blinks(data, win_size);

        %% Extract gaze data and pupil size (legacy non-baselined on [1 2]s)
        data_legacy = data(:, idx_legacy);
        gaze_x{subj, trl} = data_legacy(1, :);
        gaze_y{subj, trl} = data_legacy(2, :);
        pupil_size{subj, trl} = mean(data_legacy(3, :), 'omitnan') / 1000;
        pups = pupil_size{subj, trl};

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
        spl = sum(sqrt(dxf_s.^2 + dyf_s.^2), 'omitnan');

        %% Compute BCEA (95%) and Lateralization
        valid_bcea = isfinite(x_coords) & isfinite(y_coords);
        x_bv = double(x_coords(valid_bcea)); y_bv = double(y_coords(valid_bcea));
        if numel(x_bv) >= 10
            sx_b = std(x_bv); sy_b = std(y_bv);
            rho_b = corr(x_bv(:), y_bv(:));
            k95 = -log(1 - 0.95);
            bcea = 2 * k95 * pi * sx_b * sy_b * sqrt(1 - rho_b^2);
            bcea_lat = (mean(x_bv) - 400) / 400; % -1=left, 0=centre, +1=right
        else
            bcea = NaN; bcea_lat = NaN;
        end

        %% Baselined metrics from per-trial baseline window
        fsample = 500;
        min_valid = 100;

        xb = data(1, idx_base); yb = data(2, idx_base); pb = data(3, idx_base);
        ok_b = sum(isfinite(xb) & isfinite(yb)) >= min_valid;
        if ok_b
            gd_base = mean(sqrt((xb-400).^2 + (yb-300).^2), 'omitnan');
            spl_base = sum(sqrt(diff(xb).^2 + diff(yb).^2), 'omitnan');
            pup_base = mean(pb, 'omitnan');
            Tb = sum(isfinite(xb) & isfinite(yb)) / fsample;
            [~, msb] = detect_microsaccades(fsample, [xb; yb], numel(xb));
            ms_base = numel(msb.Onset) / Tb;
            if ~isfinite(ms_base) || ms_base <= 0, ms_base = NaN; end
            vb = isfinite(xb) & isfinite(yb);
            xb2 = double(xb(vb)); yb2 = double(yb(vb));
            if numel(xb2) >= 10
                k95 = -log(1 - 0.95);
                bcea_base = 2*k95*pi*std(xb2)*std(yb2)*sqrt(1-corr(xb2(:), yb2(:))^2);
                blat_base = (mean(xb2) - 400) / 400;
            else
                bcea_base = NaN; blat_base = NaN;
            end
        else
            gd_base = NaN; spl_base = NaN; pup_base = NaN; ms_base = NaN; bcea_base = NaN; blat_base = NaN;
        end

        winDefs = {idx_early, idx_late, idx_full};
        gd_bl = nan(1,3); spl_bl = nan(1,3); pup_bl = nan(1,3); ms_bl = nan(1,3); bcea_bl = nan(1,3); blat_bl = nan(1,3);
        for wi = 1:3
            xw = data(1, winDefs{wi}); yw = data(2, winDefs{wi}); pw = data(3, winDefs{wi});
            ok_w = sum(isfinite(xw) & isfinite(yw)) >= min_valid;
            if ~ok_w, continue; end

            gd_w = mean(sqrt((xw-400).^2 + (yw-300).^2), 'omitnan');
            spl_w = sum(sqrt(diff(xw).^2 + diff(yw).^2), 'omitnan');
            pup_w = mean(pw, 'omitnan');
            Tw = sum(isfinite(xw) & isfinite(yw)) / fsample;
            [~, msw] = detect_microsaccades(fsample, [xw; yw], numel(xw));
            ms_w = numel(msw.Onset) / Tw;
            if ~isfinite(ms_w), ms_w = NaN; end

            vw = isfinite(xw) & isfinite(yw);
            xw2 = double(xw(vw)); yw2 = double(yw(vw));
            if numel(xw2) >= 10
                k95 = -log(1 - 0.95);
                bcea_w = 2*k95*pi*std(xw2)*std(yw2)*sqrt(1-corr(xw2(:), yw2(:))^2);
                blat_w = (mean(xw2) - 400) / 400;
            else
                bcea_w = NaN; blat_w = NaN;
            end

            if isfinite(gd_w)   && isfinite(gd_base)   && gd_base > 0, gd_bl(wi)   = 100*(gd_w-gd_base)/gd_base; end
            if isfinite(spl_w)  && isfinite(spl_base)  && spl_base > 0, spl_bl(wi)  = 100*(spl_w-spl_base)/spl_base; end
            if isfinite(pup_w)  && isfinite(pup_base)  && pup_base ~= 0, pup_bl(wi)  = 100*(pup_w-pup_base)/pup_base; end
            if isfinite(ms_w)   && isfinite(ms_base)   && ms_base > 0, ms_bl(wi)   = 100*(ms_w-ms_base)/ms_base; end
            if isfinite(bcea_w) && isfinite(bcea_base) && bcea_base > 0, bcea_bl(wi) = 100*(bcea_w-bcea_base)/bcea_base; end
            if isfinite(blat_w) && isfinite(blat_base), blat_bl(wi) = blat_w - blat_base; end
        end

        %% Append data for this trial
        subject_id = [subject_id; str2num(subjects{subj})];
        trial_num = [trial_num; dataet.trialinfo(trl, 2)];
        condition = [condition; dataet.trialinfo(trl, 1)-20];
        gazeDev = [gazeDev; mean_euclidean_distance];
        gazeSDx = [gazeSDx; gaze_standard_deviation_x];
        gazeSDy = [gazeSDy; gaze_standard_deviation_y];
        pupilSize = [pupilSize; pups];
        microsaccadeRate = [microsaccadeRate; microsaccade_rate];
        scanPathLength = [scanPathLength; spl];
        bcea95Area = [bcea95Area; bcea];
        bceaLateralization = [bceaLateralization; bcea_lat];
        gdevEarlyBL = [gdevEarlyBL; gd_bl(1)]; gdevLateBL = [gdevLateBL; gd_bl(2)]; gdevFullBL = [gdevFullBL; gd_bl(3)];
        splEarlyBL = [splEarlyBL; spl_bl(1)];  splLateBL = [splLateBL; spl_bl(2)];  splFullBL = [splFullBL; spl_bl(3)];
        pupEarlyBL = [pupEarlyBL; pup_bl(1)];  pupLateBL = [pupLateBL; pup_bl(2)];  pupFullBL = [pupFullBL; pup_bl(3)];
        msEarlyBL = [msEarlyBL; ms_bl(1)];     msLateBL = [msLateBL; ms_bl(2)];     msFullBL = [msFullBL; ms_bl(3)];
        bceaEarlyBL = [bceaEarlyBL; bcea_bl(1)]; bceaLateBL = [bceaLateBL; bcea_bl(2)]; bceaFullBL = [bceaFullBL; bcea_bl(3)];
        blatEarlyBL = [blatEarlyBL; blat_bl(1)]; blatLateBL = [blatLateBL; blat_bl(2)]; blatFullBL = [blatFullBL; blat_bl(3)];

    end

    %% Check data by visualizing raw gaze data
    % close all
    % % Preallocate arrays for averaged gaze data
    % num_samples = 1000;
    % mean_gaze_x = nan(num_samples, 1);
    % mean_gaze_y = nan(num_samples, 1);
    %
    % % Stack trials into matrices for averaging
    % gaze_x_matrix = cell2mat(cellfun(@(x) x(:), gaze_x(subj, :), 'UniformOutput', false)');
    % gaze_y_matrix = cell2mat(cellfun(@(y) y(:), gaze_y(subj, :), 'UniformOutput', false)');
    %
    % % Calculate the mean over trials for each sample
    % mean_gaze_x = nanmean(gaze_x_matrix, 2);
    % mean_gaze_y = nanmean(gaze_y_matrix, 2);
    %
    % % Plot the averaged gaze data
    % figure;
    % set(gcf, "Position", [200, 200, 1000, 600]);
    % plot(mean_gaze_x, mean_gaze_y, 'o');
    % hold on;
    % plot(400, 300, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Centre point
    % title('Averaged Gaze Data Distribution Across Samples');
    % xlabel('X Coordinates');
    % ylabel('Y Coordinates');
    % xlim([0 800]);
    % ylim([0 600]);

    %% Create a trial-by-trial structure array for this subject
    subj_data_gaze_trial = struct('ID', num2cell(subject_id), ...
        'Trial', num2cell(trial_num), ...
        'Condition', num2cell(condition), ...
        'GazeDeviation', num2cell(gazeDev), ...
        'GazeStdX', num2cell(gazeSDx), ...
        'GazeStdY', num2cell(gazeSDy), ...
        'PupilSize', num2cell(pupilSize), ...
        'MSRate', num2cell(microsaccadeRate), ...
        'ScanPathLength', num2cell(scanPathLength), ...
        'BCEA', num2cell(bcea95Area), ...
        'BCEALateralization', num2cell(bceaLateralization));

    %% Calculate subject-specific data by condition (GazeDev, PupilSize, MSRate)
    l2 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 2);
    l2gdev = mean([l2.GazeDeviation], 'omitnan');
    l2gSDx = mean([l2.GazeStdX], 'omitnan');
    l2gSDy = mean([l2.GazeStdY], 'omitnan');
    l2pups = mean([l2.PupilSize], 'omitnan');
    l2msrate = mean([l2.MSRate], 'omitnan');
    l2spl = mean([l2.ScanPathLength], 'omitnan');
    l2bcea = mean([l2.BCEA], 'omitnan');
    l2blat = mean([l2.BCEALateralization], 'omitnan');

    l4 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 4);
    l4gdev = mean([l4.GazeDeviation], 'omitnan');
    l4gSDx = mean([l4.GazeStdX], 'omitnan');
    l4gSDy = mean([l4.GazeStdY], 'omitnan');
    l4pups = mean([l4.PupilSize], 'omitnan');
    l4msrate = mean([l4.MSRate], 'omitnan');
    l4spl = mean([l4.ScanPathLength], 'omitnan');
    l4bcea = mean([l4.BCEA], 'omitnan');
    l4blat = mean([l4.BCEALateralization], 'omitnan');

    l6 = subj_data_gaze_trial([subj_data_gaze_trial.Condition] == 6);
    l6gdev = mean([l6.GazeDeviation], 'omitnan');
    l6gSDx = mean([l6.GazeStdX], 'omitnan');
    l6gSDy = mean([l6.GazeStdY], 'omitnan');
    l6pups = mean([l6.PupilSize], 'omitnan');
    l6msrate = mean([l6.MSRate], 'omitnan');
    l6spl = mean([l6.ScanPathLength], 'omitnan');
    l6bcea = mean([l6.BCEA], 'omitnan');
    l6blat = mean([l6.BCEALateralization], 'omitnan');

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
    l2splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6splBLf = mean(splFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6splBLe = mean(splEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6splBLl = mean(splLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6pupBLf = mean(pupFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6pupBLe = mean(pupEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6pupBLl = mean(pupLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLf = mean(msFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLe = mean(msEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6msBLl = mean(msLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6bceaBLf = mean(bceaFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6bceaBLe = mean(bceaEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6bceaBLl = mean(bceaLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6blatBLf = mean(blatFullBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6blatBLe = mean(blatEarlyBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');
    l2blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 2), 'omitnan');
    l4blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 4), 'omitnan');
    l6blatBLl = mean(blatLateBL([subj_data_gaze_trial.Condition] == 6), 'omitnan');

    %% Create across condition structure
    subID = str2double(subjects{subj});
    subj_data_gaze = struct('ID', num2cell([subID; subID; subID]), ...
        'Condition', num2cell([2; 4; 6]), ...
        'GazeDeviation', num2cell([l2gdev; l4gdev; l6gdev]), ...
        'GazeStdX', num2cell([l2gSDx; l4gSDx; l6gSDx]), ...
        'GazeStdY', num2cell([l2gSDy; l4gSDy; l6gSDy]), ...
        'PupilSize', num2cell([l2pups; l4pups; l6pups]), ...
        'MSRate', num2cell([l2msrate; l4msrate; l6msrate]), ...
        'Blinks', num2cell([blinks_l2; blinks_l4; blinks_l6]), ...
        'Fixations', num2cell([fixations_l2; fixations_l4; fixations_l6]), ...
        'Saccades', num2cell([saccades_l2; saccades_l4; saccades_l6]), ...
        'ScanPathLength', num2cell([l2spl; l4spl; l6spl]), ...
        'BCEA', num2cell([l2bcea; l4bcea; l6bcea]), ...
        'BCEALateralization', num2cell([l2blat; l4blat; l6blat]), ...
        'GazeDeviationFullBL', num2cell([l2gdevBLf; l4gdevBLf; l6gdevBLf]), ...
        'GazeDeviationEarlyBL', num2cell([l2gdevBLe; l4gdevBLe; l6gdevBLe]), ...
        'GazeDeviationLateBL', num2cell([l2gdevBLl; l4gdevBLl; l6gdevBLl]), ...
        'ScanPathLengthFullBL', num2cell([l2splBLf; l4splBLf; l6splBLf]), ...
        'ScanPathLengthEarlyBL', num2cell([l2splBLe; l4splBLe; l6splBLe]), ...
        'ScanPathLengthLateBL', num2cell([l2splBLl; l4splBLl; l6splBLl]), ...
        'PupilSizeFullBL', num2cell([l2pupBLf; l4pupBLf; l6pupBLf]), ...
        'PupilSizeEarlyBL', num2cell([l2pupBLe; l4pupBLe; l6pupBLe]), ...
        'PupilSizeLateBL', num2cell([l2pupBLl; l4pupBLl; l6pupBLl]), ...
        'MSRateFullBL', num2cell([l2msBLf; l4msBLf; l6msBLf]), ...
        'MSRateEarlyBL', num2cell([l2msBLe; l4msBLe; l6msBLe]), ...
        'MSRateLateBL', num2cell([l2msBLl; l4msBLl; l6msBLl]), ...
        'BCEAFullBL', num2cell([l2bceaBLf; l4bceaBLf; l6bceaBLf]), ...
        'BCEAEarlyBL', num2cell([l2bceaBLe; l4bceaBLe; l6bceaBLe]), ...
        'BCEALateBL', num2cell([l2bceaBLl; l4bceaBLl; l6bceaBLl]), ...
        'BCEALatFullBL', num2cell([l2blatBLf; l4blatBLf; l6blatBLf]), ...
        'BCEALatEarlyBL', num2cell([l2blatBLe; l4blatBLe; l6blatBLe]), ...
        'BCEALatLateBL', num2cell([l2blatBLl; l4blatBLl; l6blatBLl]));

    %% Save data
    if ispc == 1
        savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\gaze\');
    else
        savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/gaze/');
    end
    mkdir(savepath)
    cd(savepath)
    save gaze_matrix_sternberg_trial subj_data_gaze_trial
    save gaze_matrix_sternberg subj_data_gaze
    save gaze_dev_sternberg l2gdev l4gdev l6gdev
    save gaze_std_sternberg l2gSDx l4gSDx l6gSDx l2gSDy l4gSDy l6gSDy
    save pupil_size_sternberg l2pups l4pups l6pups
    save ms_rate_sternberg l2msrate l4msrate l6msrate
    clc
    disp(['Subject ' num2str(subj) '/' num2str(length(subjects)) ' done.'])

    % Append to the final structure array
    gaze_data_sternberg = [gaze_data_sternberg; subj_data_gaze];
end
trialinfo = dataet.trialinfo';
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_gaze_sternberg gaze_x gaze_y trialinfo
    save W:\Students\Arne\AOC\data\features\AOC_gaze_matrix_sternberg gaze_data_sternberg
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_gaze_sternberg gaze_x gaze_y trialinfo
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_gaze_matrix_sternberg gaze_data_sternberg
end



%% AOC EEG Feature Extraction — Sternberg
% Computes subject-level power spectra in early/late/full windows (raw + baselined),
% IAF, lateralization, and FOOOF-based alpha (raw + baselined) from preprocessed
% Sternberg EEG. Saves per-subject and grand results.
%
% Extracted features:
%   Power Spectrum (Early [0 1], Late [1 2], Full [0 2]) + Baseline [-0.5 -0.25]
%   Baselined spectra (dB) for each window
%   IAF (subject-level, pooled across conditions, full window)
%   Alpha power in IAF band (early/late/full; raw + dB)
%   Lateralization index (late baselined)
%   FOOOF-based spectra (model fit - aperiodic) for each window + baselined (absolute, log space)
%   FOOOF alpha power (8-14 Hz) full + baselined (full/early/late)

%% POWSPCTRM (Baseline + Early/Late/Full) + FOOOF window spectra
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Setup logging
if ispc == 1
    logDir = 'W:\Students\Arne\AOC\data\controls\logs';
else
    logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
end
scriptName = 'AOC_eeg_fex_sternberg';

for subj = 1:length(subjects)
        clc
        disp(['Processing windowed POWSPCTRM + FOOOF (subject-level) for Subject AOC ', num2str(subjects{subj})])
        % Load data
        datapath = strcat(path, filesep, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_sternberg   % time-domain FieldTrip data struct (used for tf transforms)

        % Identify indices of trials belonging to conditions
        ind2 = find(dataTFR.trialinfo(:, 1) == 22); % WM load 2
        ind4 = find(dataTFR.trialinfo(:, 1) == 24); % WM load 4
        ind6 = find(dataTFR.trialinfo(:, 1) == 26); % WM load 6

        % Windows
        t_base = [-0.5 -0.25];
        t_early = [0 1];
        t_late  = [1 2];
        t_full  = [0 2];

        % ----------------------
        % Time-frequency transform (raw) per condition (keeps time for window selection)
        % ----------------------
        cfg            = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 3:1:30;
        cfg.t_ftimwin  = ones(size(cfg.foi))*1;     % 1 s windows
        cfg.toi        = -1.5:0.05:3;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR); % chan_freq_time
        cfg.trials = ind4; tfr4 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind6; tfr6 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = t_base;
        cfgb.baselinetype = 'db';
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr4_bl = ft_freqbaseline(cfgb, tfr4);
        tfr6_bl = ft_freqbaseline(cfgb, tfr6);

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [3 30];
        pow2_early     = remove_time_dimension(select_data(t_early, freq_range, tfr2));
        pow4_early     = remove_time_dimension(select_data(t_early, freq_range, tfr4));
        pow6_early     = remove_time_dimension(select_data(t_early, freq_range, tfr6));
        pow2_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr2_bl));
        pow4_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr4_bl));
        pow6_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr6_bl));

        pow2_late      = remove_time_dimension(select_data(t_late, freq_range, tfr2));
        pow4_late      = remove_time_dimension(select_data(t_late, freq_range, tfr4));
        pow6_late      = remove_time_dimension(select_data(t_late, freq_range, tfr6));
        pow2_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr2_bl));
        pow4_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr4_bl));
        pow6_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr6_bl));

        pow2_full      = remove_time_dimension(select_data(t_full, freq_range, tfr2));
        pow4_full      = remove_time_dimension(select_data(t_full, freq_range, tfr4));
        pow6_full      = remove_time_dimension(select_data(t_full, freq_range, tfr6));
        pow2_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr2_bl));
        pow4_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr4_bl));
        pow6_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr6_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_stern_windows.mat', ...
            'pow2_early','pow4_early','pow6_early', ...
            'pow2_late','pow4_late','pow6_late', ...
            'pow2_full','pow4_full','pow6_full', ...
            'pow2_early_bl','pow4_early_bl','pow6_early_bl', ...
            'pow2_late_bl','pow4_late_bl','pow6_late_bl', ...
            'pow2_full_bl','pow4_full_bl','pow6_full_bl')

        % ----------------------
        % FOOOF spectra per window (subject-level, no sliding)
        % Output convention: pow*_fooof = (model fit - aperiodic) in log space
        % Baselining: absolute difference (window - baseline) in log space
        % ----------------------
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = freq_range;
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';

        % Baseline window per condition (for baselining in log space)
        dat2_base = ft_selectdata(struct('latency', t_base, 'trials', ind2), dataTFR);
        dat4_base = ft_selectdata(struct('latency', t_base, 'trials', ind4), dataTFR);
        dat6_base = ft_selectdata(struct('latency', t_base, 'trials', ind6), dataTFR);

        fooof2_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_base);
        fooof4_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_base);
        fooof6_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_base);

        % Full window (0-2) non-baselined
        dat2_full = ft_selectdata(struct('latency', t_full, 'trials', ind2), dataTFR);
        dat4_full = ft_selectdata(struct('latency', t_full, 'trials', ind4), dataTFR);
        dat6_full = ft_selectdata(struct('latency', t_full, 'trials', ind6), dataTFR);
        pow2_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_full);
        pow4_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_full);
        pow6_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_full);

        % Baselined full/early/late (absolute in log space)
        pow2_fooof_bl = pow2_fooof; pow2_fooof_bl.powspctrm = pow2_fooof.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl = pow4_fooof; pow4_fooof_bl.powspctrm = pow4_fooof.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl = pow6_fooof; pow6_fooof_bl.powspctrm = pow6_fooof.powspctrm - fooof6_base.powspctrm;

        dat2_early = ft_selectdata(struct('latency', t_early, 'trials', ind2), dataTFR);
        dat4_early = ft_selectdata(struct('latency', t_early, 'trials', ind4), dataTFR);
        dat6_early = ft_selectdata(struct('latency', t_early, 'trials', ind6), dataTFR);
        fooof2_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_early);
        fooof4_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_early);
        fooof6_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_early);
        pow2_fooof_bl_early = fooof2_early; pow2_fooof_bl_early.powspctrm = fooof2_early.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl_early = fooof4_early; pow4_fooof_bl_early.powspctrm = fooof4_early.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl_early = fooof6_early; pow6_fooof_bl_early.powspctrm = fooof6_early.powspctrm - fooof6_base.powspctrm;

        dat2_late = ft_selectdata(struct('latency', t_late, 'trials', ind2), dataTFR);
        dat4_late = ft_selectdata(struct('latency', t_late, 'trials', ind4), dataTFR);
        dat6_late = ft_selectdata(struct('latency', t_late, 'trials', ind6), dataTFR);
        fooof2_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_late);
        fooof4_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat4_late);
        fooof6_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat6_late);
        pow2_fooof_bl_late = fooof2_late; pow2_fooof_bl_late.powspctrm = fooof2_late.powspctrm - fooof2_base.powspctrm;
        pow4_fooof_bl_late = fooof4_late; pow4_fooof_bl_late.powspctrm = fooof4_late.powspctrm - fooof4_base.powspctrm;
        pow6_fooof_bl_late = fooof6_late; pow6_fooof_bl_late.powspctrm = fooof6_late.powspctrm - fooof6_base.powspctrm;

        save('power_stern_fooof.mat', ...
            'pow2_fooof','pow4_fooof','pow6_fooof', ...
            'pow2_fooof_bl','pow4_fooof_bl','pow6_fooof_bl', ...
            'pow2_fooof_bl_early','pow4_fooof_bl_early','pow6_fooof_bl_early', ...
            'pow2_fooof_bl_late','pow4_fooof_bl_late','pow6_fooof_bl_late')
end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Define channels
subj = 1;
datapath = strcat(path, filesep, subjects{subj}, filesep, 'eeg');
cd(datapath);
load('power_stern_windows.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(pow2_full.label)
    label = pow2_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;
% Left and right channels
left_channels = {};
right_channels = {};
for i = 1:length(channels)
    try
        ch = channels{i};
        % Find the first numeric part in the channel name
        numStr = regexp(ch, '\d+', 'match');
        % Convert the first numerical token to a number
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch;
        else
            right_channels{end+1} = ch;
        end
    catch ME
        ME.message
        disp(['Midline channel: ', ch])
    end
end

% Load data and calculate alpha power, IAF and lateralization index
alphaRange = [8 14];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
IAF_results = struct();
eeg_data_sternberg = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {}, ...
    'AlphaPower_FOOOF', {}, 'AlphaPower_FOOOF_bl', {}, 'AlphaPower_FOOOF_bl_early', {}, 'AlphaPower_FOOOF_bl_late', {});

for subj = 1:length(subjects)
        clc
        disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])
        datapath = strcat(path, filesep, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_stern_windows.mat');
        load('power_stern_fooof.mat');

        % Channel selection
        channelIdx = find(ismember(pow2_full.label, channels));

        % Extract FULL-window power spectra for selected channels (used for IAF)
        powspctrm2 = mean(pow2_full.powspctrm(channelIdx, :), 1);
        powspctrm4 = mean(pow4_full.powspctrm(channelIdx, :), 1);
        powspctrm6 = mean(pow6_full.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(pow2_full.freq >= alphaRange(1) & pow2_full.freq <= alphaRange(2));

        % Calculate IAF for WM load 2
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = pow2_full.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(pow2_full.freq > (IAF2-4) & pow2_full.freq < (IAF2+2));
            powerIAF2 = mean(powspctrm2(IAF_range2));
        end

        % Calculate IAF for WM load 4
        alphaPower4 = powspctrm4(alphaIndices);
        [pks4,locs] = findpeaks(alphaPower4);
        if isempty(pks4)
            IAF4 = NaN;
            IAF_range4 = NaN;
            powerIAF4 = NaN;
        else
            [~, ind] = max(pks4);
            IAF4 = pow4_full.freq(alphaIndices(locs(ind)));
            IAF_range4 = find(pow4_full.freq > (IAF4-4) & pow4_full.freq < (IAF4+2));
            powerIAF4 = mean(powspctrm4(IAF_range4));
        end

        % Calculate IAF for WM load 6
        alphaPower6 = powspctrm6(alphaIndices);
        [pks6,locs] = findpeaks(alphaPower6);
        if isempty(pks6)
            IAF6 = NaN;
            IAF_range6 = NaN;
            powerIAF6 = NaN;
        else
            [~, ind] = max(pks6);
            IAF6 = pow6_full.freq(alphaIndices(locs(ind)));
            IAF_range6 = find(pow6_full.freq > (IAF6-4) & pow6_full.freq < (IAF6+2));
            powerIAF6 = mean(powspctrm6(IAF_range6));
        end

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF4 == alphaRange(1) || IAF4 == alphaRange(2)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if IAF6 == alphaRange(1) || IAF6 == alphaRange(2)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF4 > max(pks4)
            powerIAF4 = NaN;
            IAF4 = NaN;
        end
        if powerIAF6 > max(pks6)
            powerIAF6 = NaN;
            IAF6 = NaN;
        end

        % Compute lateralization index on LATE BASELINED spectra (dB)
        powloads = {pow2_late_bl, pow4_late_bl, pow6_late_bl};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx2 = LatIdx(1);
        LatIdx4 = LatIdx(2);
        LatIdx6 = LatIdx(3);

        % Alpha power in IAF band (raw + baselined) for early/late/full
        IAF_band2 = [IAF2-4 IAF2+2]; if any(isnan(IAF_band2)); IAF_band2 = alphaRange; end
        IAF_band4 = [IAF4-4 IAF4+2]; if any(isnan(IAF_band4)); IAF_band4 = alphaRange; end
        IAF_band6 = [IAF6-4 IAF6+2]; if any(isnan(IAF_band6)); IAF_band6 = alphaRange; end

        % helper anonymous to average ROI x freq
        roi_pow = @(S, band) mean(mean(S.powspctrm(channelIdx, S.freq>=band(1) & S.freq<=band(2)), 2, 'omitnan'), 1, 'omitnan');

        AlphaPowerEarly    = [roi_pow(pow2_early, IAF_band2);    roi_pow(pow4_early, IAF_band4);    roi_pow(pow6_early, IAF_band6)];
        AlphaPowerLate     = [roi_pow(pow2_late,  IAF_band2);    roi_pow(pow4_late,  IAF_band4);    roi_pow(pow6_late,  IAF_band6)];
        AlphaPowerFull     = [roi_pow(pow2_full,  IAF_band2);    roi_pow(pow4_full,  IAF_band4);    roi_pow(pow6_full,  IAF_band6)];
        AlphaPowerEarlyBL  = [roi_pow(pow2_early_bl, IAF_band2); roi_pow(pow4_early_bl, IAF_band4); roi_pow(pow6_early_bl, IAF_band6)];
        AlphaPowerLateBL   = [roi_pow(pow2_late_bl,  IAF_band2); roi_pow(pow4_late_bl,  IAF_band4); roi_pow(pow6_late_bl,  IAF_band6)];
        AlphaPowerFullBL   = [roi_pow(pow2_full_bl,  IAF_band2); roi_pow(pow4_full_bl,  IAF_band4); roi_pow(pow6_full_bl,  IAF_band6)];
        AlphaPower_FOOOF         = [roi_pow(pow2_fooof,          alphaRange); roi_pow(pow4_fooof,          alphaRange); roi_pow(pow6_fooof,          alphaRange)];
        AlphaPower_FOOOF_bl      = [roi_pow(pow2_fooof_bl,       alphaRange); roi_pow(pow4_fooof_bl,       alphaRange); roi_pow(pow6_fooof_bl,       alphaRange)];
        AlphaPower_FOOOF_bl_early= [roi_pow(pow2_fooof_bl_early, alphaRange); roi_pow(pow4_fooof_bl_early, alphaRange); roi_pow(pow6_fooof_bl_early, alphaRange)];
        AlphaPower_FOOOF_bl_late = [roi_pow(pow2_fooof_bl_late,  alphaRange); roi_pow(pow4_fooof_bl_late,  alphaRange); roi_pow(pow6_fooof_bl_late,  alphaRange)];

        % Keep legacy AlphaPower as LATE raw (registered retention)
        AlphaPower = AlphaPowerLate;

        % Create a structure array for this subject (legacy fields preserved)
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([2; 4; 6]), ...
            'AlphaPower', num2cell(AlphaPower), 'IAF', num2cell([IAF2; IAF4; IAF6]), ...
            'Lateralization', num2cell([LatIdx2; LatIdx4; LatIdx6]));

        % Add windowed alpha power fields (raw + baselined, IAF band)
        tmp = num2cell(AlphaPowerEarly);   [subj_data_eeg.AlphaPowerEarly]   = tmp{:};
        tmp = num2cell(AlphaPowerLate);    [subj_data_eeg.AlphaPowerLate]    = tmp{:};
        tmp = num2cell(AlphaPowerFull);    [subj_data_eeg.AlphaPowerFull]    = tmp{:};
        tmp = num2cell(AlphaPowerEarlyBL); [subj_data_eeg.AlphaPowerEarlyBL] = tmp{:};
        tmp = num2cell(AlphaPowerLateBL);  [subj_data_eeg.AlphaPowerLateBL]  = tmp{:};
        tmp = num2cell(AlphaPowerFullBL);  [subj_data_eeg.AlphaPowerFullBL]  = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF);          [subj_data_eeg.AlphaPower_FOOOF]          = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl);       [subj_data_eeg.AlphaPower_FOOOF_bl]       = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_early); [subj_data_eeg.AlphaPower_FOOOF_bl_early] = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_late);  [subj_data_eeg.AlphaPower_FOOOF_bl_late]  = tmp{:};

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_sternberg_subj subj_data_eeg
        save alpha_power_sternberg powerIAF2 powerIAF4 powerIAF6
        save IAF_sternberg IAF2 IAF4 IAF6
        save lateralization_sternberg LatIdx2 LatIdx4 LatIdx6
        eeg_data_sternberg = [eeg_data_sternberg; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: WM2: %f Hz (Power: %f), WM4: %f Hz (Power: %f), ' ...
            'WM6: %f Hz (Power: %f) | Lateralization: %f %f %f \n'], subjects{subj}, ...
            IAF2, powerIAF2, IAF4, powerIAF4, IAF6, powerIAF6, LatIdx2, LatIdx4, LatIdx6);
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_sternberg eeg_data_sternberg
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_sternberg eeg_data_sternberg
end

%% AOC EEG Feature Extraction — N-Back
% Computes subject-level power spectra in early/late/full windows (raw + baselined),
% IAF, lateralization, and FOOOF-based alpha (raw + baselined) from preprocessed
% N-back EEG. Saves per-subject and grand results.
%
% Extracted features:
%   Power Spectrum
%   IAF, Power at IAF, and Lateralization Index
%   POWER TRIAL-BY-TRIAL
%   TFR (Raw, FOOOF and Baselined) and FOOOFed POWSPCTRM

%% POWSPCTRM (Baseline + Early/Late/Full) + FOOOF window spectra (subject-level)
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Setup logging
if ispc == 1
    logDir = 'W:\Students\Arne\AOC\data\controls\logs';
else
    logDir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/logs';
end
scriptName = 'AOC_eeg_fex_nback';

for subj = 1:length(subjects)
        clc
        disp(['Processing windowed POWSPCTRM + FOOOF (subject-level) for Subject AOC ', num2str(subjects{subj})])
        % Load data
        datapath = strcat(path, filesep, subjects{subj}, filesep, 'eeg');
        cd(datapath)
        close all
        load dataEEG_TFR_nback

        % Identify indices of trials belonging to conditions
        ind1 = find(dataTFR.trialinfo(:, 1) == 21);
        ind2 = find(dataTFR.trialinfo(:, 1) == 22);
        ind3 = find(dataTFR.trialinfo(:, 1) == 23);

        % Windows
        t_base = [-0.5 -0.25];
        t_early = [0 1];
        t_late  = [1 2];
        t_full  = [0 2];

        % ----------------------
        % Time-frequency transform (raw) per condition (keeps time for window selection)
        % ----------------------
        cfg            = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'pow';
        cfg.taper      = 'hanning';
        cfg.foi        = 3:1:30;
        cfg.t_ftimwin  = ones(size(cfg.foi))*1;     % 1 s windows
        cfg.toi        = -1.25:0.05:2.25;
        cfg.pad        = 'maxperlen';
        cfg.keeptrials = 'no';

        cfg.trials = ind1; tfr1 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind2; tfr2 = ft_freqanalysis(cfg, dataTFR);
        cfg.trials = ind3; tfr3 = ft_freqanalysis(cfg, dataTFR);

        % Baselined TFR (dB)
        cfgb              = [];
        cfgb.baseline     = t_base;
        cfgb.baselinetype = 'db';
        tfr1_bl = ft_freqbaseline(cfgb, tfr1);
        tfr2_bl = ft_freqbaseline(cfgb, tfr2);
        tfr3_bl = ft_freqbaseline(cfgb, tfr3);

        % Convert to window-collapsed POWSPCTRM (chan x freq)
        freq_range = [3 30];
        pow1_early     = remove_time_dimension(select_data(t_early, freq_range, tfr1));
        pow2_early     = remove_time_dimension(select_data(t_early, freq_range, tfr2));
        pow3_early     = remove_time_dimension(select_data(t_early, freq_range, tfr3));
        pow1_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr1_bl));
        pow2_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr2_bl));
        pow3_early_bl  = remove_time_dimension(select_data(t_early, freq_range, tfr3_bl));

        pow1_late      = remove_time_dimension(select_data(t_late, freq_range, tfr1));
        pow2_late      = remove_time_dimension(select_data(t_late, freq_range, tfr2));
        pow3_late      = remove_time_dimension(select_data(t_late, freq_range, tfr3));
        pow1_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr1_bl));
        pow2_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr2_bl));
        pow3_late_bl   = remove_time_dimension(select_data(t_late, freq_range, tfr3_bl));

        pow1_full      = remove_time_dimension(select_data(t_full, freq_range, tfr1));
        pow2_full      = remove_time_dimension(select_data(t_full, freq_range, tfr2));
        pow3_full      = remove_time_dimension(select_data(t_full, freq_range, tfr3));
        pow1_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr1_bl));
        pow2_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr2_bl));
        pow3_full_bl   = remove_time_dimension(select_data(t_full, freq_range, tfr3_bl));

        % Save windowed spectra
        cd(datapath)
        save('power_nback_windows.mat', ...
            'pow1_early','pow2_early','pow3_early', ...
            'pow1_late','pow2_late','pow3_late', ...
            'pow1_full','pow2_full','pow3_full', ...
            'pow1_early_bl','pow2_early_bl','pow3_early_bl', ...
            'pow1_late_bl','pow2_late_bl','pow3_late_bl', ...
            'pow1_full_bl','pow2_full_bl','pow3_full_bl')

        % ----------------------
        % FOOOF spectra per window (subject-level, no sliding)
        % Baselining: absolute difference (window - baseline) in log space
        % ----------------------
        cfg_fooof            = [];
        cfg_fooof.method     = 'mtmfft';
        cfg_fooof.taper      = 'hanning';
        cfg_fooof.foilim     = freq_range;
        cfg_fooof.pad        = 5;
        cfg_fooof.output     = 'fooof';
        cfg_fooof.keeptrials = 'no';

        dat1_base = ft_selectdata(struct('latency', t_base, 'trials', ind1), dataTFR);
        dat2_base = ft_selectdata(struct('latency', t_base, 'trials', ind2), dataTFR);
        dat3_base = ft_selectdata(struct('latency', t_base, 'trials', ind3), dataTFR);
        fooof1_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat1_base);
        fooof2_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_base);
        fooof3_base = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat3_base);

        dat1_full = ft_selectdata(struct('latency', t_full, 'trials', ind1), dataTFR);
        dat2_full = ft_selectdata(struct('latency', t_full, 'trials', ind2), dataTFR);
        dat3_full = ft_selectdata(struct('latency', t_full, 'trials', ind3), dataTFR);
        pow1_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat1_full);
        pow2_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_full);
        pow3_fooof = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat3_full);

        pow1_fooof_bl = pow1_fooof; pow1_fooof_bl.powspctrm = pow1_fooof.powspctrm - fooof1_base.powspctrm;
        pow2_fooof_bl = pow2_fooof; pow2_fooof_bl.powspctrm = pow2_fooof.powspctrm - fooof2_base.powspctrm;
        pow3_fooof_bl = pow3_fooof; pow3_fooof_bl.powspctrm = pow3_fooof.powspctrm - fooof3_base.powspctrm;

        dat1_early = ft_selectdata(struct('latency', t_early, 'trials', ind1), dataTFR);
        dat2_early = ft_selectdata(struct('latency', t_early, 'trials', ind2), dataTFR);
        dat3_early = ft_selectdata(struct('latency', t_early, 'trials', ind3), dataTFR);
        fooof1_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat1_early);
        fooof2_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_early);
        fooof3_early = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat3_early);
        pow1_fooof_bl_early = fooof1_early; pow1_fooof_bl_early.powspctrm = fooof1_early.powspctrm - fooof1_base.powspctrm;
        pow2_fooof_bl_early = fooof2_early; pow2_fooof_bl_early.powspctrm = fooof2_early.powspctrm - fooof2_base.powspctrm;
        pow3_fooof_bl_early = fooof3_early; pow3_fooof_bl_early.powspctrm = fooof3_early.powspctrm - fooof3_base.powspctrm;

        dat1_late = ft_selectdata(struct('latency', t_late, 'trials', ind1), dataTFR);
        dat2_late = ft_selectdata(struct('latency', t_late, 'trials', ind2), dataTFR);
        dat3_late = ft_selectdata(struct('latency', t_late, 'trials', ind3), dataTFR);
        fooof1_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat1_late);
        fooof2_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat2_late);
        fooof3_late = ft_freqanalysis_Arne_FOOOF(cfg_fooof, dat3_late);
        pow1_fooof_bl_late = fooof1_late; pow1_fooof_bl_late.powspctrm = fooof1_late.powspctrm - fooof1_base.powspctrm;
        pow2_fooof_bl_late = fooof2_late; pow2_fooof_bl_late.powspctrm = fooof2_late.powspctrm - fooof2_base.powspctrm;
        pow3_fooof_bl_late = fooof3_late; pow3_fooof_bl_late.powspctrm = fooof3_late.powspctrm - fooof3_base.powspctrm;

        save('power_nback_fooof.mat', ...
            'pow1_fooof','pow2_fooof','pow3_fooof', ...
            'pow1_fooof_bl','pow2_fooof_bl','pow3_fooof_bl', ...
            'pow1_fooof_bl_early','pow2_fooof_bl_early','pow3_fooof_bl_early', ...
            'pow1_fooof_bl_late','pow2_fooof_bl_late','pow3_fooof_bl_late')
end

%% ALPHA POWER, IAF and LATERALIZATION INDEX
% Setup
startup
[subjects, paths, ~ , ~] = setup('AOC');
path = paths.features;

% Define channels
datapath = strcat(path, filesep, subjects{1}, filesep, 'eeg');
cd(datapath);
load('power_nback_windows.mat');
% Occipital channels
occ_channels = {};
for i = 1:length(pow2_full.label)
    label = pow2_full.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        occ_channels{end+1} = label;
    end
end
channels = occ_channels;

% Left and right channels
left_channels = {};
right_channels = {};
for i = 1:length(channels)
    try
        ch = channels{i};
        % Find the first numeric part in the channel name
        numStr = regexp(ch, '\d+', 'match');
        % Convert the first numerical token to a number
        numVal = str2double(numStr{1});
        if mod(numVal, 2) == 1
            left_channels{end+1} = ch;
        else
            right_channels{end+1} = ch;
        end
    catch ME
        ME.message
        disp(['Midline channel: ', ch])
    end
end

% Load data and calculate power, IAF and lateralization index
close all
alphaRange = [8 14];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();
eeg_data_nback = struct('ID', {}, 'Condition', {}, 'AlphaPower', {}, 'IAF', {}, 'Lateralization', {}, ...
    'AlphaPower_FOOOF', {}, 'AlphaPower_FOOOF_bl', {}, 'AlphaPower_FOOOF_bl_early', {}, 'AlphaPower_FOOOF_bl_late', {});

for subj = 1:length(subjects)
        clc
        disp(['Processing Alpha Power, IAF and Lateralization for Subject AOC ', num2str(subjects{subj})])
        datapath = strcat(path, filesep, subjects{subj}, filesep, 'eeg');
        cd(datapath);
        load('power_nback_windows.mat');
        load('power_nback_fooof.mat');

        % Channels selection based on CBPT
        channelIdx = find(ismember(pow1_full.label, channels));

        % FULL-window power spectra for IAF
        powspctrm1 = mean(pow1_full.powspctrm(channelIdx, :), 1);
        powspctrm2 = mean(pow2_full.powspctrm(channelIdx, :), 1);
        powspctrm3 = mean(pow3_full.powspctrm(channelIdx, :), 1);

        % Find the indices corresponding to the alpha range
        alphaIndices = find(pow1_full.freq >= alphaRange(1) & pow1_full.freq <= alphaRange(2));

        % Calculate IAF for 1-back
        alphaPower1 = powspctrm1(alphaIndices);
        [pks1,locs] = findpeaks(alphaPower1);
        if isempty(pks1)
            IAF1 = NaN;
            IAF_range1 = NaN;
            powerIAF1 = NaN;
        else
            [~, ind] = max(pks1);
            IAF1 = pow1_full.freq(alphaIndices(locs(ind)));
            IAF_range1 = find(pow1_full.freq > (IAF1-4) & pow1_full.freq < (IAF1+2));
            powerIAF1 = mean(powspctrm1(IAF_range1));
        end

        % Calculate IAF for 2-back
        alphaPower2 = powspctrm2(alphaIndices);
        [pks2,locs] = findpeaks(alphaPower2);
        if isempty(pks2)
            IAF2 = NaN;
            IAF_range2 = NaN;
            powerIAF2 = NaN;
        else
            [~, ind] = max(pks2);
            IAF2 = pow2_full.freq(alphaIndices(locs(ind)));
            IAF_range2 = find(pow2_full.freq > (IAF2-4) & pow2_full.freq < (IAF2+2));
            powerIAF2 = mean(powspctrm2(IAF_range2));
        end

        % Calculate IAF for 3-back
        alphaPower3 = powspctrm3(alphaIndices);
        [pks3,locs] = findpeaks(alphaPower3);
        if isempty(pks3)
            IAF3 = NaN;
            IAF_range3 = NaN;
            powerIAF3 = NaN;
        else
            [~, ind] = max(pks3);
            IAF3 = pow3_full.freq(alphaIndices(locs(ind)));
            IAF_range3 = find(pow3_full.freq > (IAF3-4) & pow3_full.freq < (IAF3+2));
            powerIAF3 = mean(powspctrm3(IAF_range3));
        end

        % Do not extract alpha peak if there is no clear peak
        % Check if any IAF is 8 or 14 and set the corresponding power to NaN
        if IAF1 == alphaRange(1) || IAF1 == alphaRange(2)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if IAF2 == alphaRange(1) || IAF2 == alphaRange(2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if IAF3 == alphaRange(1) || IAF3 == alphaRange(2)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Check if the averaged power at IAF -4/+2 Hz is more than the peak
        % of power. If so, set power to NaN
        if powerIAF1 > max(pks1)
            powerIAF1 = NaN;
            IAF1 = NaN;
        end
        if powerIAF2 > max(pks2)
            powerIAF2 = NaN;
            IAF2 = NaN;
        end
        if powerIAF3 > max(pks3)
            powerIAF3 = NaN;
            IAF3 = NaN;
        end

        % Compute lateralization index on LATE BASELINED spectra (dB)
        powloads = {pow1_late_bl, pow2_late_bl, pow3_late_bl};
        for i = 1:3
            curr_load = powloads{i};
            left_idx = find(ismember(curr_load.label, left_channels));
            right_idx = find(ismember(curr_load.label, right_channels));
            alpha_idx = find(curr_load.freq >= alphaRange(1) & curr_load.freq <= alphaRange(2));
            left_alpha_power = mean(mean(curr_load.powspctrm(left_idx, alpha_idx), 2));
            right_alpha_power = mean(mean(curr_load.powspctrm(right_idx, alpha_idx), 2));
            LatIdx(i) = (right_alpha_power - left_alpha_power) / (right_alpha_power + left_alpha_power);
        end
        LatIdx1 = LatIdx(1);
        LatIdx2 = LatIdx(2);
        LatIdx3 = LatIdx(3);

        % Alpha power in IAF band (raw + baselined) for early/late/full
        IAF_band1 = [IAF1-4 IAF1+2]; if any(isnan(IAF_band1)); IAF_band1 = alphaRange; end
        IAF_band2 = [IAF2-4 IAF2+2]; if any(isnan(IAF_band2)); IAF_band2 = alphaRange; end
        IAF_band3 = [IAF3-4 IAF3+2]; if any(isnan(IAF_band3)); IAF_band3 = alphaRange; end

        roi_pow = @(S, band) mean(mean(S.powspctrm(channelIdx, S.freq>=band(1) & S.freq<=band(2)), 2, 'omitnan'), 1, 'omitnan');

        AlphaPowerEarly    = [roi_pow(pow1_early, IAF_band1);    roi_pow(pow2_early, IAF_band2);    roi_pow(pow3_early, IAF_band3)];
        AlphaPowerLate     = [roi_pow(pow1_late,  IAF_band1);    roi_pow(pow2_late,  IAF_band2);    roi_pow(pow3_late,  IAF_band3)];
        AlphaPowerFull     = [roi_pow(pow1_full,  IAF_band1);    roi_pow(pow2_full,  IAF_band2);    roi_pow(pow3_full,  IAF_band3)];
        AlphaPowerEarlyBL  = [roi_pow(pow1_early_bl, IAF_band1); roi_pow(pow2_early_bl, IAF_band2); roi_pow(pow3_early_bl, IAF_band3)];
        AlphaPowerLateBL   = [roi_pow(pow1_late_bl,  IAF_band1); roi_pow(pow2_late_bl,  IAF_band2); roi_pow(pow3_late_bl,  IAF_band3)];
        AlphaPowerFullBL   = [roi_pow(pow1_full_bl,  IAF_band1); roi_pow(pow2_full_bl,  IAF_band2); roi_pow(pow3_full_bl,  IAF_band3)];
        AlphaPower_FOOOF         = [roi_pow(pow1_fooof,          alphaRange); roi_pow(pow2_fooof,          alphaRange); roi_pow(pow3_fooof,          alphaRange)];
        AlphaPower_FOOOF_bl      = [roi_pow(pow1_fooof_bl,       alphaRange); roi_pow(pow2_fooof_bl,       alphaRange); roi_pow(pow3_fooof_bl,       alphaRange)];
        AlphaPower_FOOOF_bl_early= [roi_pow(pow1_fooof_bl_early, alphaRange); roi_pow(pow2_fooof_bl_early, alphaRange); roi_pow(pow3_fooof_bl_early, alphaRange)];
        AlphaPower_FOOOF_bl_late = [roi_pow(pow1_fooof_bl_late,  alphaRange); roi_pow(pow2_fooof_bl_late,  alphaRange); roi_pow(pow3_fooof_bl_late,  alphaRange)];

        % Keep legacy AlphaPower as LATE raw (registered retention)
        AlphaPower = AlphaPowerLate;

        % Create a structure array for this subject (legacy fields preserved)
        subID = str2num(subjects{subj});
        subj_data_eeg = struct('ID', num2cell([subID; subID; subID]), 'Condition', num2cell([1; 2; 3]), ...
            'AlphaPower', num2cell(AlphaPower), 'IAF', num2cell([IAF1; IAF2; IAF3]), 'Lateralization', num2cell([LatIdx1; LatIdx2; LatIdx3]));

        tmp = num2cell(AlphaPowerEarly);   [subj_data_eeg.AlphaPowerEarly]   = tmp{:};
        tmp = num2cell(AlphaPowerLate);    [subj_data_eeg.AlphaPowerLate]    = tmp{:};
        tmp = num2cell(AlphaPowerFull);    [subj_data_eeg.AlphaPowerFull]    = tmp{:};
        tmp = num2cell(AlphaPowerEarlyBL); [subj_data_eeg.AlphaPowerEarlyBL] = tmp{:};
        tmp = num2cell(AlphaPowerLateBL);  [subj_data_eeg.AlphaPowerLateBL]  = tmp{:};
        tmp = num2cell(AlphaPowerFullBL);  [subj_data_eeg.AlphaPowerFullBL]  = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF);          [subj_data_eeg.AlphaPower_FOOOF]          = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl);       [subj_data_eeg.AlphaPower_FOOOF_bl]       = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_early); [subj_data_eeg.AlphaPower_FOOOF_bl_early] = tmp{:};
        tmp = num2cell(AlphaPower_FOOOF_bl_late);  [subj_data_eeg.AlphaPower_FOOOF_bl_late]  = tmp{:};

        % Save
        if ispc == 1
            savepath = strcat('W:\Students\Arne\AOC\data\features\', subjects{subj}, '\eeg\');
        else
            savepath = strcat('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, '/eeg/');
        end
        mkdir(savepath)
        cd(savepath)
        save eeg_matrix_nback_subj subj_data_eeg
        save alpha_power_nback powerIAF1 powerIAF2 powerIAF3
        save IAF_nback IAF1 IAF2 IAF3
        save lateralization_nback LatIdx1 LatIdx2 LatIdx3
        eeg_data_nback = [eeg_data_nback; subj_data_eeg];
        clc
        fprintf(['Subject %s IAF: 1-back: %f Hz (Power: %f), 2-back: %f Hz (Power: %f), ' ...
            '3-back: %f Hz (Power: %f) |Lateralization: %f %f %f \n'], subjects{subj}, IAF1, ...
            powerIAF1, IAF2, powerIAF2, IAF3, powerIAF3, LatIdx1, LatIdx2, LatIdx3);
end
if ispc == 1
    save W:\Students\Arne\AOC\data\features\AOC_eeg_matrix_nback eeg_data_nback
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/AOC_eeg_matrix_nback eeg_data_nback
end
