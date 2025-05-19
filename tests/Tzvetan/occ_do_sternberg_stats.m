clear all
close all

% Subject IDs
subjects = {'301','302','304','309','310','312','313','314','315', ...
    '321','325','327','328','329','330','331','335', ...
    '336','339','341','342','343','344','345','346','348', ...
    '349','352','355','356','359','361','362','364','365', ...
    '368','372','373','374','375','377','379','386','388', ...
    '389','390','391','394','397','398','401','403','404','406', ...
   '413','416','419'};%,'306' '319', '320',,'347' ,'358','378'  '412',

base_dir = '/Volumes/TOURO/arne/merged';
%% Loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
%     load(strcat(subjects{s},'_Sternberg_cond52_tfr.mat'));
    load(strcat(subjects{s},'_Sternberg_cond52_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2{s}=tfr_fooof;
%     load(strcat(subjects{s},'_Sternberg_cond54_tfr.mat'));
    load(strcat(subjects{s},'_Sternberg_cond54_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load4{s}=tfr_fooof;
%     load(strcat(subjects{s},'_Sternberg_cond56_tfr.mat'));
    load(strcat(subjects{s},'_Sternberg_cond56_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load6{s} = tfr_fooof;
end
%%
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/layANThead.mat');
cfg = [];
ga2 = ft_freqgrandaverage(cfg,load2{:});
ga4 = ft_freqgrandaverage(cfg,load4{:});
ga6 = ft_freqgrandaverage(cfg,load6{:});
%%
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-2 2];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, ga2);
%%
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
diff6m2 = ft_math(cfg,ga6,ga2);
%%
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-2 2];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, diff6m2);
%%
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/elec_aligned.mat');% adapt the path according to your setup
cfg =[];
cfg.method ='distance';
cfg.elec = elec_aligned;
cfg.feedback      = 'yes' ;
cfg.neighbourdist=40;
neighbours = ft_prepare_neighbours(cfg);
%%
cfg = [];
% cfg.latency          = [1.5 3];
% cfg.frequency        = [5 30];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;
subj = numel(load6);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, load6{:}, load2{:});
%%
n = numel(load2);  
cohens_d = stat.stat ./ sqrt(n);
stat.effectsize = cohens_d;
%%
cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'effectsize';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
figure; ft_multiplotTFR(cfg,stat);