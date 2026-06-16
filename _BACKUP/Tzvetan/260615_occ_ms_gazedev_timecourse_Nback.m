%% N-back gaze deviation and microsaccade-rate time courses

clear all; close all; clc

load('/Volumes/Homestore/OCC/arne/subjects.mat');

base_dir = '/Volumes/Homestore/OCC/arne/merged';
addpath('/Volumes/Homestore/OCC/arne/funcs');

%% parameters
fsample = 500;

% blink cleaning settings
thresh = 20;
pad_ms = 150;

% EyeLink screen centre
screen_cx = 400;
screen_cy = 300;

% time courses
latency_win = [-1 2];

% microsaccade-rate smoothing window, controls the temporal rate estimate, not detection
ms_win_ms = 100;

%% output containers
allDev = cell(1,3);
allMS  = cell(1,3);

%% loop subjects
for s = 1:length(subjects)

    subj = subjects{s};
    fprintf('\nSubject %s (%d/%d)\n', subj, s, length(subjects));

    cd(fullfile(base_dir, subj));
    load([subj '_Nback_all.mat']);

    % select eye channels
    cfg = [];
    cfg.channel = {'L-GAZE-X','L-GAZE-Y'};
    dataet = ft_selectdata(cfg, dataNback);

    nTrials = numel(dataet.trial);

    % blink cleaning
    dataetnan    = dataet;
    dataetinterp = dataet;

    dataetinterp.trialinfo_blink = false(nTrials,1);
    dataetinterp.blink_mask      = cell(nTrials,1);

    valid_trials = true(nTrials,1);

    for i = 1:nTrials

        x = dataet.trial{i}(1,:);
        y = dataet.trial{i}(2,:);
        t = dataet.time{i};

        [x_nan, y_nan, x_interp, y_interp, blink_mask, is_valid] = ...
            removeAndInterpolateBlinks_checktrials(x, y, t, fsample, thresh, pad_ms);

        if ~is_valid
            valid_trials(i) = false;
            continue
        end

        dataetnan.trial{i}(1,:) = x_nan;
        dataetnan.trial{i}(2,:) = y_nan;

        dataetinterp.trial{i}(1,:) = x_interp;
        dataetinterp.trial{i}(2,:) = y_interp;

        dataetinterp.blink_mask{i} = blink_mask;
        dataetinterp.trialinfo_blink(i) = any(blink_mask);
    end

    % remove invalid trials
    cfg = [];
    cfg.trials = find(valid_trials);

    dataetnan    = ft_selectdata(cfg, dataetnan);
    dataetinterp = ft_selectdata(cfg, dataetinterp);
    dataNback    = ft_selectdata(cfg, dataNback);

    % create gaze deviation FieldTrip structure

    dataGazeDev = [];
    dataGazeDev.label     = {'GAZEDEV'};
    dataGazeDev.fsample   = dataetinterp.fsample;
    dataGazeDev.trialinfo = dataetinterp.trialinfo;

    if isfield(dataetinterp,'sampleinfo')
        dataGazeDev.sampleinfo = dataetinterp.sampleinfo;
    end

    % create microsaccade-rate FieldTrip structure

    dataMSrate = [];
    dataMSrate.label     = {'MSRATE'};
    dataMSrate.fsample   = dataetinterp.fsample;
    dataMSrate.trialinfo = dataetinterp.trialinfo;

    if isfield(dataetinterp,'sampleinfo')
        dataMSrate.sampleinfo = dataetinterp.sampleinfo;
    end

    % compute trial-wise time courses

    scalar_dev = nan(numel(dataetinterp.trial),1);
    scalar_ms  = nan(numel(dataetinterp.trial),1);

    for tr = 1:numel(dataetinterp.trial)

        t = dataetinterp.time{tr};

        %  Euclidean deviation from screen centre in pixels
  
        x_interp = dataetinterp.trial{tr}(1,:);
        y_interp = dataetinterp.trial{tr}(2,:);

        dev = sqrt((x_interp - screen_cx).^2 + ...
                   (y_interp - screen_cy).^2);

        dataGazeDev.trial{tr} = dev;
        dataGazeDev.time{tr}  = t;

        scalar_dev(tr) = mean(dev,'omitnan');

  
   
        %  MS rate Engbert & Kliegl kernel, threshold 6, mindur 6 samples
        % uses NaN blink-cleaned gaze as detector input
        
        velData = dataetnan.trial{tr};      % 2 x time
        trlLength = size(velData,2);

        [msrate, microsaccades, ms_scalar] = ...
            detect_microsaccades_timecourse_student( ...
                dataetinterp.fsample, velData, trlLength, ms_win_ms);

        dataMSrate.trial{tr} = msrate;
        dataMSrate.time{tr}  = t;

        scalar_ms(tr) = ms_scalar;
    end

    % select conditions

    trl1 = find(dataGazeDev.trialinfo == 21); % 1-back
    trl2 = find(dataGazeDev.trialinfo == 22); % 2-back
    trl3 = find(dataGazeDev.trialinfo == 23); % 3-back

    % subject-level condition timelocks

    for cond = 1:3

        eval(sprintf('trl = trl%d;',cond));

        cfg = [];
        cfg.trials  = trl;
        cfg.latency = latency_win;

        dev_cond = ft_selectdata(cfg, dataGazeDev);
        ms_cond  = ft_selectdata(cfg, dataMSrate);

        cfg = [];
        tl_dev = ft_timelockanalysis(cfg, dev_cond);
        tl_ms  = ft_timelockanalysis(cfg, ms_cond);

        allDev{cond}{s} = tl_dev;
        allMS{cond}{s}  = tl_ms;
    end

    %  scalar checks per subject
    l1gdev_check = mean(scalar_dev(trl1),'omitnan');
    l2gdev_check = mean(scalar_dev(trl2),'omitnan');
    l3gdev_check = mean(scalar_dev(trl3),'omitnan');

    l1msrate_check = mean(scalar_ms(trl1),'omitnan');
    l2msrate_check = mean(scalar_ms(trl2),'omitnan');
    l3msrate_check = mean(scalar_ms(trl3),'omitnan');

    fprintf('GazeDev px: %.3f %.3f %.3f\n', ...
        l1gdev_check, l2gdev_check, l3gdev_check);

    fprintf('MS rate Hz: %.3f %.3f %.3f\n', ...
        l1msrate_check, l2msrate_check, l3msrate_check);
end
%% plot single subj 
cfg = [];
cfg.channel = 'GAZEDEV';
ft_singleplotER(cfg, allDev{1}{1}, allDev{2}{1}, allDev{3}{1});
legend({'1-back','2-back','3-back'});
ylabel('Euclidean gaze deviation [px]');
xlabel('Time [s]');
title('N-back gaze deviation');
%%

cfg = [];
cfg.channel = 'MSRATE';

figure('Color','w');
ft_singleplotER(cfg, allMS{1}{1}, allMS{2}{1}, allMS{3}{1});
legend({'1-back','2-back','3-back'});
ylabel('Microsaccade rate [s^{-1}]');
xlabel('Time [s]');
title('N-back microsaccade rate');
%% grand average across subjects

for cond = 1:3

    cfg = [];
    cfg.keepindividual = 'yes';

    GA_dev{cond} = ft_timelockgrandaverage(cfg, allDev{cond}{:});
    GA_ms{cond}  = ft_timelockgrandaverage(cfg, allMS{cond}{:});
end

%% plot gaze deviation

cfg = [];
cfg.channel = 'GAZEDEV';

figure('Color','w');
ft_singleplotER(cfg, GA_dev{1}, GA_dev{2}, GA_dev{3});
legend({'1-back','2-back','3-back'});
ylabel('Euclidean gaze deviation [px]');
xlabel('Time [s]');
title('N-back gaze deviation');

%% plot microsaccade rate

cfg = [];
cfg.channel = 'MSRATE';

figure('Color','w');
ft_singleplotER(cfg, GA_ms{1}, GA_ms{2}, GA_ms{3});
legend({'1-back','2-back','3-back'});
ylabel('Microsaccade rate [s^{-1}]');
xlabel('Time [s]');
title('N-back microsaccade rate');

%% helper function

function [msrate, microsaccades, microsaccade_rate_scalar] = ...
    detect_microsaccades_timecourse_student(fsample, velData, trlLength, win_ms)

    if nargin < 4
        win_ms = 100;
    end

    % Convolution kernel as per Engbert et al., 2003
    kernel = [1 1 0 -1 -1] .* (fsample/6);

    % parameters
    velthres = 6;
    mindur   = 6;

    % padding and convolution
    n = size(kernel,2);
    pad = ceil(n/2);

    dat = ft_preproc_padding(velData, 'localmean', pad);
    vel = convn(dat, kernel, 'same');
    vel = ft_preproc_padding(vel, 'remove', pad);

    % compute velocity thresholds
    medianstd = sqrt(median(vel.^2, 2, 'omitnan') - ...
                     (median(vel, 2, 'omitnan')).^2);

    radius = velthres * medianstd;

    % protection against bad trials
    if any(radius == 0) || any(isnan(radius))
        msrate = zeros(1,trlLength);
        microsaccades = [];
        microsaccade_rate_scalar = NaN;
        return
    end

    % microsaccade detection based on threshold crossing
    test = sum((vel ./ radius(:, ones(1, size(vel,2)))).^2, 1);
    sacsmp = find(test > 1);

    microsaccades = [];

    if ~isempty(sacsmp)

        j = find(diff(sacsmp) == 1);

        if ~isempty(j)

            j1 = [j; j + 1];
            com = intersect(j, j + 1);
            cut = ~ismember(j1, com);

            sacidx = reshape(j1(cut), 2, []);

            for k = 1:size(sacidx,2)

                duration = sacidx(1,k):sacidx(2,k);

                if length(duration) >= mindur

                    onset  = sacsmp(duration(1));
                    offset = sacsmp(duration(end));
                    peak   = sacsmp(duration(round(length(duration)/2)));

                    microsaccades = [microsaccades; onset, offset, peak];
                end
            end
        end
    end

    % scalar rate
    microsaccade_rate_scalar = ...
        size(microsaccades,1) / ((trlLength-1) / fsample);

    % event train
    ms_binary = zeros(1,trlLength);

    if ~isempty(microsaccades)
        onsets = microsaccades(:,1);
        onsets = onsets(onsets >= 1 & onsets <= trlLength);
        ms_binary(onsets) = 1;
    end

    % time-resolved rate estimate in Hz
    win_samp = round(win_ms/1000 * fsample);

    msrate = smoothdata(ms_binary, 'movmean', win_samp) * fsample;
end