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
    load(strcat(subjects{s},'_Sternberg_cond52_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2{s}=tfr_fooof;
    load(strcat(subjects{s},'_Sternberg_cond54_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load4{s}=tfr_fooof;
    load(strcat(subjects{s},'_Sternberg_cond56_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load6{s} = tfr_fooof;
end
%% compute diff stern
for s = 1:length(load6)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    sb_high_low{s} = ft_math(cfg,load6{s},load2{s});
    % select retention
    cfg = [];
%     cfg.latency = [1.5 3];
    cfg.latency = [0 3];
    sb_high_low{s} = ft_selectdata(cfg,sb_high_low{s});
end
%% handle nback loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
    load(strcat(subjects{s},'_Nback_cond21_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load1{s}=tfr_fooof;
    load(strcat(subjects{s},'_Nback_cond22_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2nb{s}=tfr_fooof;
    load(strcat(subjects{s},'_Nback_cond23_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load3{s} = tfr_fooof;
end
%% compute diff nback
for s = 1:length(load3)
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'powspctrm';
    nb_high_low{s} = ft_math(cfg,load3{s},load1{s});
    % select retention
    cfg = [];
%     cfg.latency = [.5 2];
    cfg.latency = [0 3];
    nb_high_low{s} = ft_selectdata(cfg,nb_high_low{s});
    % equalize time dim
%     nb_high_low{s}.time = sb_high_low{s}.time;
end
%%
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/layANThead.mat');
cfg = [];
ga_nb = ft_freqgrandaverage(cfg,nb_high_low{:});
ga_sb = ft_freqgrandaverage(cfg,sb_high_low{:});
%%
% close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-2 2];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, ga_nb);
%%
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
omnibus = ft_math(cfg,ga_sb,ga_nb);
%%
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, omnibus);
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
cfg.latency          = [0 3];
% cfg.frequency        = [5 30];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
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

[stat] = ft_freqstatistics(cfg, sb_high_low{:}, nb_high_low{:});
%%
n = numel(load3);  % number of subjects (paired samples)
cohens_d = stat.stat ./ sqrt(n);
stat.effectsize = cohens_d;
% make 0 be the onset of retention for both tasks, i.e. 1.5 sec after
% letter presentation in Nback and the last 1.5 sec of retention in
% sternberg
% stat.time=stat.time-1.5;
%%

close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'effectsize';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
cfg.zlim = [-.8 .8];
figure; ft_multiplotTFR(cfg,stat);
%%
stattfr=stat;
stattfr.stat= stat.effectsize;
figure;
cfg = [];
% cfg.channel = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan = 'yes';
cfg.frequency = [5 30];
% cfg.latency   = [0 1.5];
freq = ft_selectdata(cfg,stattfr);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(5, 30, 500);
mask_interp = linspace(5, 30, 500);
% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, meanmask,...
    tim_grid_interp, freq_grid_interp, 'spline');
figure;
% subplot(2,1,2);ft_plot_matrix(flip(pow_interp));
subplot(2,1,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))));

yticks([1 100 200 300 400 500 ]);
yticklabels({ '30','25','20','15','10','5'});
xticks([1 250 500])
xticklabels({num2str(freq.time(1)),'1.5',  num2str(freq.time(end))});
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.8 .8]);
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.8 0 .8];
title(c,'Effect size \it d')

cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
cfg.parameter = 'effectsize';
 cfg.xlim = [0.195536 0.595536];
 cfg.ylim = [10.288321 16.740876];
cfg.zlim = [-.8 .8];
cfg.marker             = 'off';
cfg.highlight          = 'on';
% cfg.highlightchannel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.highlightchannel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

cfg.highlightsymbol    = '.';
 cfg.highlightsize      = 14;
 cfg.comment = 'no';
% figure; 
subplot(2,1,2);ft_topoplotTFR(cfg,stat);
%% extract values sternberg
for s = 1:length(load6)
    % select retention
    cfg = [];
%     cfg.latency = [0.2 0.6]+1.5;
%     cfg.frequency = [10 16];
        cfg.latency = [0.5 3];
    cfg.frequency = [8 14];
    cfg.channel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
   cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

    cfg.avgoverfreq = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    val6{s} = ft_selectdata(cfg,load6{s});
    sb6(s) = val6{s}.powspctrm;
    
    val4{s} = ft_selectdata(cfg,load4{s});
    sb4(s) = val4{s}.powspctrm;
    
    val2{s} = ft_selectdata(cfg,load2{s});
    sb2(s) = val2{s}.powspctrm;
end
%%
figure(30); clf;
subplot(1,2,1);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(sb2);
fill(positions(1) + f_sb2*0.2, xi_sb2, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(sb4);
fill(positions(2) + f_sb4*0.2, xi_sb4, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(sb6);
fill(positions(3) + f_sb6*0.2, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2', 'load 4', 'load 6'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(sb2))-0.5)*jitter, sb2, 'b.');
scatter(positions(2) + (rand(size(sb4))-0.5)*jitter, sb4, 'k.');
scatter(positions(3) + (rand(size(sb6))-0.5)*jitter, sb6, 'r.');

% Axis & Label
% ylabel('\muV');
ylabel('\alpha power [a.u.] ');
xlim([0 1]);
title('Sternberg');
box on;
set(gcf,'color','w');
set(gca,'Fontsize',20);

% Significance tests
[~, p_24] = ttest(sb2, sb4);
[~, p_46] = ttest(sb4, sb6);
[~, p_26] = ttest(sb2, sb6);

% Significance annotations
y_max = max([sb2(:); sb4(:); sb6(:)]) + 0.5;  % smaller margin above data overall height above
y_step = 0.1;  % tighter vertical spacing between lines

sig_label = getSigLabel(p_24);
if ~isempty(sig_label)
    line([positions(1), positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(2)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_46);
if ~isempty(sig_label)
    line([positions(2), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(2), positions(3)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_26);
if ~isempty(sig_label)
    line([positions(1), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(3)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
end
ylim([min([sb2(:); sb4(:); sb6(:)]) - 0.5, y_max + 0.1]);
xlim([0 1.3])

% % Significance annotations
% y_max = max([sb2(:); sb4(:); sb6(:)]) + 5;
% 
% sig_label = getSigLabel(p_24);
% if ~isempty(sig_label)
%     line([positions(1), positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', .5);
%     text(mean([positions(1), positions(2)]), y_max, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'Center');
%     y_max = y_max + 2;
% end
% 
% sig_label = getSigLabel(p_46);
% if ~isempty(sig_label)
%     line([positions(2), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', .5);
%     text(mean([positions(2), positions(3)]), y_max, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'Center');
%     y_max = y_max + 2;
% end
% 
% sig_label = getSigLabel(p_26);
% if ~isempty(sig_label)
%     line([positions(1), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', .5);
%     text(mean([positions(1), positions(3)]), y_max, sig_label, 'FontSize', 14, 'HorizontalAlignment', 'Center');
% end
% 
% ylim([min([sb2(:); sb4(:); sb6(:)]) - 5, y_max + 2]);
%% extract values nback
for s = 1:length(load6)
    % select retention
    cfg = [];
%     cfg.latency = [0.2 0.6]+.5;
%     cfg.frequency = [10 16];
            cfg.latency = [0.5 3];
    cfg.frequency = [8 14];
    cfg.channel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

    cfg.avgoverfreq = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    val3{s} = ft_selectdata(cfg,load3{s});
    nb3(s) = val3{s}.powspctrm;
    
    val2nb{s} = ft_selectdata(cfg,load2nb{s});
    nb2(s) = val2nb{s}.powspctrm;
    
    val1{s} = ft_selectdata(cfg,load1{s});
    nb1(s) = val1{s}.powspctrm;
end
%%
% figure(31); clf;
subplot(1,2,2);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(nb1);
fill(positions(1) + f_sb2*0.2, xi_sb2, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(nb2);
fill(positions(2) + f_sb4*0.2, xi_sb4, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(nb3);
fill(positions(3) + f_sb6*0.2, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 1', 'load 2', 'load 3'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(sb2))-0.5)*jitter, sb2, 'b.');
scatter(positions(2) + (rand(size(sb4))-0.5)*jitter, sb4, 'k.');
scatter(positions(3) + (rand(size(sb6))-0.5)*jitter, sb6, 'r.');

% Axis & Label
% ylabel('\muV');
ylabel('\alpha power [a.u.] ');
xlim([0 1]);
title('N back');
box on;
set(gcf,'color','w');
set(gca,'Fontsize',20);

% Significance tests
[~, p_24] = ttest(nb1, nb2);
[~, p_46] = ttest(nb2, nb3);
[~, p_26] = ttest(nb1, nb3);

% Significance annotations
y_max = max([sb2(:); sb4(:); sb6(:)]) + 0.5;  % smaller margin above data overall height above
y_step = 0.1;  % tighter vertical spacing between lines

sig_label = getSigLabel(p_24);
if ~isempty(sig_label)
    line([positions(1), positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(2)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_46);
if ~isempty(sig_label)
    line([positions(2), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(2), positions(3)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_26);
if ~isempty(sig_label)
    line([positions(1), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(3)]), y_max + 0.02, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
end
ylim([min([sb2(:); sb4(:); sb6(:)]) - 0.5, y_max + 0.1]);
xlim([0 1.3])

%% Helper function
function sig_label = getSigLabel(p)
    if p < 0.001
        sig_label = '***';
    elseif p < 0.01
        sig_label = '**';
    elseif p < 0.05
        sig_label = '*';
    else
        sig_label = '';
    end
end
