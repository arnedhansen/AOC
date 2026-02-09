%% AOC Omnibus Data Preparation (plain: raw TFR, no FOOOF)
% Loads raw TFR (Sternberg 2/4/6, N-back 1/2/3) from same subject files,
% restricts to 3–30 Hz, applies baseline (dB), computes high–low contrasts and grand averages.
%
% Key outputs:
%   omnibus_data.mat (baselined raw TFR, load diffs, GAs)
%
% run time approx. 30 minutes

%% Setup
tic
startup
[subjects, path, colors, headmodel] = setup('AOC');

cfg_freq = [];
cfg_freq.frequency = [3 30];

%% Load N-back raw TFR; restrict 3–30 Hz; apply baseline (dB)
for subj = 1:length(subjects)
    clc
    disp(['Loading N-back raw TFR for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_nback
    tfr1_sel = ft_selectdata(cfg_freq, tfr1);
    tfr2_sel = ft_selectdata(cfg_freq, tfr2);
    tfr3_sel = ft_selectdata(cfg_freq, tfr3);
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'db';
    load1nb{subj} = ft_freqbaseline(cfg, tfr1_sel);
    load2nb{subj} = ft_freqbaseline(cfg, tfr2_sel);
    load3nb{subj} = ft_freqbaseline(cfg, tfr3_sel);
end

%% Compute diff nback
clc
disp('Computing N-back Diff')
for subj = 1:length(subjects)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    nb_high_low{subj} = ft_math(cfg, load3nb{subj}, load1nb{subj});
    cfg = [];
    cfg.latency = [-.5 2];
    nb_high_low{subj} = ft_selectdata(cfg, nb_high_low{subj});
end

%% Load Sternberg raw TFR; restrict 3–30 Hz; apply baseline (dB)
for subj = 1:length(subjects)
    clc
    disp(['Loading Sternberg raw TFR for Subject ', subjects{subj}])
    datapath = strcat(path, subjects{subj}, filesep, 'eeg');
    cd(datapath)
    load tfr_stern
    tfr2_sel = ft_selectdata(cfg_freq, tfr2);
    tfr4_sel = ft_selectdata(cfg_freq, tfr4);
    tfr6_sel = ft_selectdata(cfg_freq, tfr6);
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'db';
    load2{subj} = ft_freqbaseline(cfg, tfr2_sel);
    load4{subj} = ft_freqbaseline(cfg, tfr4_sel);
    load6{subj} = ft_freqbaseline(cfg, tfr6_sel);
end

%% Compute diff stern
clc
disp('Computing Sternberg Diff')
for subj = 1:length(subjects)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    sb_high_low{subj} = ft_math(cfg, load6{subj}, load2{subj});
    cfg = [];
    cfg.latency = [-.5 2];
    sb_high_low{subj} = ft_selectdata(cfg, sb_high_low{subj});
end

%% Grand average of differences
clc
disp('Computing Grand Averages')
cfg = [];
ga_nb = ft_freqgrandaverage(cfg, nb_high_low{:});
ga_sb = ft_freqgrandaverage(cfg, sb_high_low{:});

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
disp('Saved omnibus_data.mat')
toc
