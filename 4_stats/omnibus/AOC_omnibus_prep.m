%% AOC Omnibus Data Preparation
%
% Load subject list and task-specific TFR data (Sternberg: loads 2/4/6; N-back: loads 1/2/3)
% Apply baseline correction to single-subject TFRs and compute high–low load contrasts per task
% Compute grand average TFRs (per task and per load) and visualise time–frequency/topography patterns

%% Setup
tic
startup
[subjects, path, colors, headmodel] = setup('AOC');

%% Load N-back TFR FOOOF data and apply baseline
for subj = 1:length(subjects)
    clc
    disp(['Loading N-back TFR FOOOF data for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_nback
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr1_fooof);
    load1nb{subj} = tfr;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr2_fooof);
    load2nb{subj} = tfr;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr3_fooof);
    load3nb{subj} = tfr;
end

% check data
meanspctrm = squeeze(mean(tfr.powspctrm, 1));   % freq x time
meanspctrm = squeeze(mean(meanspctrm, 2));      % freq x 1

figure('Color','w'); hold on
plot(tfr.freq, meanspctrm, 'LineWidth', 3)
xlabel('Frequency (Hz)')
ylabel('Power (FOOOF space)')
set(gca, 'FontSize', 15)
xlim([tfr.freq(1) tfr.freq(end)])
%set(gca, 'YScale', 'log')

%% Compute diff nback
clc
disp('Computing N-back Diff')
for subj = 1:length(subjects)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    nb_high_low{subj} = ft_math(cfg,load3nb{subj},load1nb{subj});
    cfg = [];
    cfg.latency = [-.5 2];
    nb_high_low{subj} = ft_selectdata(cfg,nb_high_low{subj});
end

%% Load Sternberg TFR FOOOF data and apply baseline
for subj = 1:length(subjects)
    clc
    disp(['Loading Sternberg TFR FOOOF data for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_stern
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr2_fooof);
    load2{subj} = tfr;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr4_fooof);
    load4{subj} = tfr;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr6_fooof);
    load6{subj} = tfr;
end

%% Compute diff stern
clc
disp('Computing Sternberg Diff')
for subj = 1:length(subjects)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    sb_high_low{subj} = ft_math(cfg,load6{subj},load2{subj});
    cfg = [];
    cfg.latency = [-.5 2];
    sb_high_low{subj} = ft_selectdata(cfg,sb_high_low{subj});
end

%% Grand average of differences
clc
disp('Computing Grand Averages')
cfg = [];
ga_nb = ft_freqgrandaverage(cfg,nb_high_low{:});
ga_sb = ft_freqgrandaverage(cfg,sb_high_low{:});

%% Save variables
disp('Saving...')
if ispc
    save W:\Students\Arne\AOC\data\features\omnibus_data.mat ...
        load1nb load2nb load3nb nb_high_low ga_nb ...
        load2 load4 load6 sb_high_low ga_sb
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data.mat ...
        load1nb load2nb load3nb nb_high_low ga_nb ...
        load2 load4 load6 sb_high_low ga_sb
end 
toc