%% AOC Omnibus Data Preparation
% Loads TFR FOOOF (Sternberg 2/4/6, N-back 1/2/3), applies baseline, 
% computes high–low contrasts and grand averages. Also builds a no-baseline version.
%
% Key outputs:
%   omnibus_data.mat (baselined TFR, load diffs, GAs)
%   omnibus_data_NOBASELINE.mat (per-condition TFR only, no baseline; no diffs; _nob suffix)
%
% run time approx. 30 minutes

%% Setup
tic
startup
[subjects, path, colors, headmodel] = setup('AOC');

%% Load N-back TFR FOOOF data; store no-baseline and apply baseline
for subj = 1:length(subjects)
    clc
    disp(['Loading N-back TFR FOOOF data for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_nback
    load1nb_nob{subj} = tfr1_fooof;
    load2nb_nob{subj} = tfr2_fooof;
    load3nb_nob{subj} = tfr3_fooof;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr1_fooof);
    load1nb{subj} = tfr;
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr  = ft_freqbaseline(cfg,tfr2_fooof);
    load2nb{subj} = tfr;
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr  = ft_freqbaseline(cfg,tfr3_fooof);
    load3nb{subj} = tfr;
end

% check data (for last subject)
% meanspctrm = squeeze(mean(tfr.powspctrm, 1));   % freq x time
% meanspctrm = squeeze(mean(meanspctrm, 2));      % freq x 1
% 
% disp(['Loaded N-back data for Subject ', subjects{end}, ' with dimensions: ', num2str(size(load1nb{end}.powspctrm))]);
% fprintf('powspctrm: min %.3e | max %.3e\n', min(tfr.powspctrm(:)), max(tfr.powspctrm(:)))
% fprintf('meanspctrm: min %.3e | max %.3e\n', min(meanspctrm(:)), max(meanspctrm(:)))
% 
% figure('Color','w'); hold on
% plot(tfr.freq, meanspctrm, 'LineWidth', 3)
% xlabel('Frequency (Hz)')
% ylabel('Power (FOOOF space)')
% set(gca, 'FontSize', 15)
% xlim([tfr.freq(1) tfr.freq(end)])

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

%% Load Sternberg TFR FOOOF data; store no-baseline and apply baseline
for subj = 1:length(subjects)
    clc
    disp(['Loading Sternberg TFR FOOOF data for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_stern
    load2_nob{subj} = tfr2_fooof;
    load4_nob{subj} = tfr4_fooof;
    load6_nob{subj} = tfr6_fooof;
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute'; % FOOOF already in log space
    tfr  = ft_freqbaseline(cfg,tfr2_fooof);
    load2{subj} = tfr;
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr  = ft_freqbaseline(cfg,tfr4_fooof);
    load4{subj} = tfr;
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr  = ft_freqbaseline(cfg,tfr6_fooof);
    load6{subj} = tfr;
end

% check data (for last subject)
% meanspctrm = squeeze(mean(tfr.powspctrm, 1));   % freq x time
% meanspctrm = squeeze(mean(meanspctrm, 2));      % freq x 1
% 
% disp(['Loaded Sternberg data for Subject ', subjects{end}, ' with dimensions: ', num2str(size(load2{end}.powspctrm))]);
% fprintf('powspctrm: min %.3e | max %.3e\n', min(tfr.powspctrm(:)), max(tfr.powspctrm(:)))
% fprintf('meanspctrm: min %.3e | max %.3e\n', min(meanspctrm(:)), max(meanspctrm(:)))
% 
% figure('Color','w'); hold on
% plot(tfr.freq, meanspctrm, 'LineWidth', 3)
% xlabel('Frequency (Hz)')
% ylabel('Power (FOOOF space)')
% set(gca, 'FontSize', 15)
% xlim([tfr.freq(1) tfr.freq(end)])

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
    save W:\Students\Arne\AOC\data\features\omnibus_data_NOBASELINE.mat ...
        load1nb_nob load2nb_nob load3nb_nob load2_nob load4_nob load6_nob
else
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data.mat ...
        load1nb load2nb load3nb nb_high_low ga_nb ...
        load2 load4 load6 sb_high_low ga_sb
    save /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data_NOBASELINE.mat ...
        load1nb_nob load2nb_nob load3nb_nob load2_nob load4_nob load6_nob
end
disp('Saved omnibus_data.mat and omnibus_data_NOBASELINE.mat')
toc