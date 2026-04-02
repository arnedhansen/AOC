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
    datapath = fullfile(path, subjects{subj}, 'gaze');
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
        savepath = fullfile('W:\Students\Arne\AOC\data\features\', subjects{subj}, 'gaze');
    else
        savepath = fullfile('/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/', subjects{subj}, 'gaze');
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
