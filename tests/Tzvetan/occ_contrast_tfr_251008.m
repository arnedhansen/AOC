clear all
close all

% Subject IDs
load('/Volumes/Homestore/OCC/arne/subjects.mat');

base_dir = '/Volumes/Homestore/OCC/arne/merged';

%% Loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
    load(strcat(subjects{s},'_Sternberg_cond22_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond24_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load4{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond26_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load6{s} = tfr;
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
    cfg.latency = [-1 3];
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
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load1{s}=tfr;
    load(strcat(subjects{s},'_Nback_cond22_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2nb{s}=tfr;
    load(strcat(subjects{s},'_Nback_cond23_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-Inf -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load3{s} = tfr;
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
    cfg.latency = [-1 3];
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
%% sternberg per condition
cfg = [];
ga_sb_2 = ft_freqgrandaverage(cfg,load2{:});
ga_sb_4 = ft_freqgrandaverage(cfg,load4{:});
ga_sb_6 = ft_freqgrandaverage(cfg,load6{:});
ga_sb = ft_freqgrandaverage(cfg,load2{:},load4{:},load6{:});
ga_nb_2 = ft_freqgrandaverage(cfg,load2nb{:});
ga_nb_1 = ft_freqgrandaverage(cfg,load1{:});
ga_nb_3 = ft_freqgrandaverage(cfg,load3{:});
ga_nb = ft_freqgrandaverage(cfg,load1{:},load2nb{:},load3{:});

%%
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, ga_nb);
%% single 
cfg = [];
cfg.figure = 'gcf';
cfg.zlim = [-.2 .2];
cfg.xlim = [-1 2];
cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};% sb
cfg.channel = {'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb
cfg.layout = layANThead;
figure; ft_singleplotTFR(cfg, ga_nb);
%% plot SB tfr and topo
figure;
cfg = [];
cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};
cfg.avgoverchan = 'yes';
cfg.frequency = [1 40];
cfg.latency   = [-1 3];
freq = ft_selectdata(cfg,ga_sb);
meanpow = squeeze(mean(freq.powspctrm, 1));

% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(1, 40, 500);
% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');

figure;
subplot(2,1,2);ft_plot_matrix(flip(pow_interp));
% subplot(2,1,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))));
ax = gca; hold(ax,'on');

% map 0 sec to the interpolated column index
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);

% ticks (still index-based)
xticks( round(interp1(tim_interp, 1:numel(tim_interp), [-1 0 1 2 3])) );
xticklabels({'-1','0','1','2','3'});
yticks([1 125 250 375 ]); % positions in the interpolated grid
yticklabels({'40','30','20','10'}); % corresponding freq values

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.2 .2]);
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.2 0 .2];
title(c,"Power change \newline from baseline")
xline(0,'k--','LineWidth',2); % black dashed line at 0 sec
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
 cfg.xlim = [1.5 3];
 cfg.ylim = [8 12];
cfg.zlim = [-.1 .1];
cfg.marker             = 'off';
cfg.highlight          = 'on';
cfg.highlightchannel = {'Pz', 'POz', 'P2', 'PPO2'};

cfg.highlightsymbol    = '.';
 cfg.highlightsize      = 14;
 cfg.comment = 'no';
% figure; 
subplot(2,1,1);ft_topoplotTFR(cfg,ga_sb);
%% plot NB tfr and topo
figure;
cfg = [];
cfg.channel = {'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb
cfg.avgoverchan = 'yes';
cfg.frequency = [1 40];
cfg.latency   = [-1 2];
freq = ft_selectdata(cfg,ga_nb);
meanpow = squeeze(mean(freq.powspctrm, 1));

% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(1, 40, 500);
% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');

figure;
subplot(2,1,2);ft_plot_matrix(flip(pow_interp));
% subplot(2,1,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))));
ax = gca; hold(ax,'on');

% map 0 sec to the interpolated column index
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);

% ticks (still index-based)
xticks( round(interp1(tim_interp, 1:numel(tim_interp), [-1 0 1 2])) );
xticklabels({'-1','0','1','2'});
yticks([1 125 250 375 ]); % positions in the interpolated grid
yticklabels({'40','30','20','10'}); % corresponding freq values

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.2 .2]);
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.2 0 .2];
title(c,"Power change \newline from baseline")
xline(0,'k--','LineWidth',2); % black dashed line at 0 sec
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
 cfg.xlim = [0 2];
 cfg.ylim = [8 12];
cfg.zlim = [-.15 .15];
cfg.marker             = 'off';
cfg.highlight          = 'on';
% cfg.highlightchannel = {'Pz', 'POz', 'P2', 'PPO2'};
cfg.highlightchannel = {'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb


cfg.highlightsymbol    = '.';
 cfg.highlightsize      = 14;
 cfg.comment = 'no';
% figure; 
subplot(2,1,1);ft_topoplotTFR(cfg,ga_nb);
%% plot SB posterior
cfg = [];
cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};
% cfg.channel = {'P8', 'PO4', 'PO8', 'PPO6h', 'POO4h'};
cfg.figure = 'gcf';
cfg.ylim = [3 40];
cfg.zlim = [-.2 .2];
cfg.xlim = [-.5 3];
figure; 
subplot(3,1,1); ft_singleplotTFR(cfg,ga_sb_2);
subplot(3,1,2); ft_singleplotTFR(cfg,ga_sb_4);
subplot(3,1,3); ft_singleplotTFR(cfg,ga_sb_6);
%% plot NB posterior
cfg = [];
% cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};
cfg.channel = {'POz', 'PO3', 'PO4', 'PPO1', 'PPO6h'};
cfg.figure = 'gcf';
cfg.ylim = [3 40];
cfg.zlim = [-.2 .2];
cfg.xlim = [-.5 2];
figure; 
subplot(3,1,1); ft_singleplotTFR(cfg,ga_nb_1);
subplot(3,1,2); ft_singleplotTFR(cfg,ga_nb_2);
subplot(3,1,3); ft_singleplotTFR(cfg,ga_nb_3);
%% extract power spectra SB
cfg = [];
cfg.latency = [1 3];
cfg.avgovertime = 'yes';
for s = 1 :length(subjects)
load2_pow{s}= ft_selectdata(cfg,load2{s});
load4_pow{s}= ft_selectdata(cfg,load4{s});
load6_pow{s}= ft_selectdata(cfg,load6{s});
load2_pow{s}.dimord = 'chan_freq';
load2_pow{s} = rmfield(load2_pow{s},'time');
load4_pow{s}.dimord = 'chan_freq';
load4_pow{s} = rmfield(load4_pow{s},'time');
load6_pow{s}.dimord = 'chan_freq';
load6_pow{s} = rmfield(load6_pow{s},'time');
end
cfg.latency = [.5 2];
for s = 1 :length(subjects)
load1_pow{s}= ft_selectdata(cfg,load1{s});
load2nb_pow{s}= ft_selectdata(cfg,load2nb{s});
load3_pow{s}= ft_selectdata(cfg,load3{s});
load1_pow{s}.dimord = 'chan_freq';
load1_pow{s} = rmfield(load1_pow{s},'time');
load2nb_pow{s}.dimord = 'chan_freq';
load2nb_pow{s} = rmfield(load2nb_pow{s},'time');
load3_pow{s}.dimord = 'chan_freq';
load3_pow{s} = rmfield(load3_pow{s},'time');
end
%% sternberg pow per condition
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2pow = ft_freqgrandaverage(cfg,load2_pow{:});
ga_sb_4pow = ft_freqgrandaverage(cfg,load4_pow{:});
ga_sb_6pow = ft_freqgrandaverage(cfg,load6_pow{:});
ga_nb_2pow = ft_freqgrandaverage(cfg,load2nb_pow{:});
ga_nb_1pow = ft_freqgrandaverage(cfg,load1_pow{:});
ga_nb_3pow = ft_freqgrandaverage(cfg,load3_pow{:});

%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotER(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
%% plot with SE sternberg
% close all
figure;
subplot(2,1,1);
cfg = [];
cfg.channel ={'Pz', 'POz', 'P2', 'PPO2'};
cfg.channel = {'P5', 'PPO5h'};% based on F stat
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga_sb_2pow);
tlk4_ind        = ft_selectdata(cfg,ga_sb_4pow);
tlk6_ind        = ft_selectdata(cfg,ga_sb_6pow);
% plot load 6
x = tlk6_ind.freq'; % x-axis def
y = mean(squeeze(tlk6_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk6_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound


hp1 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r', 'HandleVisibility', 'off');
hold on;
hl1 = line(x, y);
set(hp1, 'facecolor', [0.97, 0.26, 0.26], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl1, 'color', [0.97, 0.26, 0.26], 'linewidth', 2);

% plot load 4
x = tlk4_ind.freq'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% plot load 2
x = tlk2_ind.freq'; % x-axis def
y = mean(squeeze(tlk2_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk2_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp3 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'k', 'HandleVisibility', 'off');
hl3 = line(x, y);
set(hp3, 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl3, 'color', 'k', 'linewidth', 2);

% Label the axes
set(gca, 'FontSize', 22);
title('');
xlabel('Frequency [Hz]');
ylabel("Power change \newline from baseline");
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',20);
%% plot with SE nback
% close all
% figure;
subplot(2,1,2);
cfg = [];
cfg.channel ={'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb
cfg.channel = {'P7', 'PPO9h'};% based on Fstat
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga_nb_1pow);
tlk4_ind        = ft_selectdata(cfg,ga_nb_2pow);
tlk6_ind        = ft_selectdata(cfg,ga_nb_3pow);
% plot load 3
x = tlk6_ind.freq'; % x-axis def
y = mean(squeeze(tlk6_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk6_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound


hp1 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r', 'HandleVisibility', 'off');
hold on;
hl1 = line(x, y);
set(hp1, 'facecolor', [0.97, 0.26, 0.26], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl1, 'color', [0.97, 0.26, 0.26], 'linewidth', 2);

% plot load 2
x = tlk4_ind.freq'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% plot load 1
x = tlk2_ind.freq'; % x-axis def
y = mean(squeeze(tlk2_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk2_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp3 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'k', 'HandleVisibility', 'off');
hl3 = line(x, y);
set(hp3, 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl3, 'color', 'k', 'linewidth', 2);

% Label the axes
set(gca, 'FontSize', 22);
title('');
xlabel('Frequency [Hz]');
ylabel("Power change \newline from baseline");
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','southeast','Fontsize',20);
%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotER(cfg, ga_nb_1pow,ga_nb_2pow,ga_nb_3pow);
%%
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/elec_aligned.mat');% adapt the path according to your setup
cfg =[];
cfg.method ='distance';
cfg.elec = elec_aligned;
cfg.feedback      = 'yes' ;
cfg.neighbourdist=40;
neighbours = ft_prepare_neighbours(cfg);
%%
cfg                  = [];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesFunivariate'; % use the dependent samples F-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                                              % permutation distribution.
cfg.minnbchan        = 2;          % minimum number of neighborhood channels that is
                                                  %    required for a selected sample to be included
                                                %   in the clustering algorithm (default=0).
cfg.neighbours       = neighbours; % see below
cfg.tail             = 1;          % 1 as the F distribution is skewed
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.025;      % alpha level of the permutation test
cfg.numrandomization = 1000;        % number of draws from the permutation distribution

n_U  = numel(subjects);
n_P  = numel(subjects);
n_N = numel(subjects);

clear design
design = zeros(2,3*n_U);
cfg.design(1,:)           = [ones(1,n_U), ones(1,n_P)*2,ones(1,n_N)*3]; % design matrix
cfg.design(2,:)           = [1:n_U,1:n_P, 1:n_N]; 
cfg.ivar                  = 1; % number or list with indices indicating the independent variable(s)
cfg.uvar                  = 2;% units-of-observation (subjects or trials
[statFnb] = ft_freqstatistics(cfg, ga_nb_1pow,ga_nb_2pow,ga_nb_3pow);
% [statFsb] = ft_freqstatistics(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
%%
close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
% cfg.maskstyle = 'outline';
% cfg.zlim = [-.8 .8];
figure; ft_multiplotER(cfg,statFsb);
%% 

cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
cfg.linecolor = 'k';
cfg.linewidth = 2;
cfg.comment = 'no';
figure; ft_multiplotER(cfg,statFnb);
set(gcf,'color','w');
%%
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
cfg.marker             = 'off';
cfg.highlight          = 'on';
cfg.colormap = 'YlOrRd';
% cfg.highlightchannel = {'P7', 'P3', 'P5', 'PO3', 'TP7', 'PO7', 'TPP9h', 'PO9', 'P9', 'TPP7h', 'PPO9h', 'PPO5h', 'POO9h'};% sb
cfg.highlightchannel = {'P7', 'P3', 'O1', 'P5', 'P1', 'PO3', 'TP7', 'PO7', 'TPP9h', 'PO9', 'P9', 'CPP5h', 'CPP3h', 'PPO1', 'I1', 'TPP7h', 'PPO9h', 'PPO5h', 'POO9h', 'POO3h', 'OI1h'};% sb
cfg.highlightsymbol    = '.';
 cfg.highlightsize      = 14;
 cfg.comment = 'no';
cfg.parameter = 'stat';
% cfg.xlim = [7.267857 9.963710];% sb
cfg.xlim = [8.975230 17.601959];% nb
% cfg.zlim = [0 7];%sb
cfg.zlim = [0 20];%sb

figure; ft_topoplotER(cfg,statFnb);
set(gcf,'color','w');
caxis([0 20]);
c = colorbar;
c.LineWidth = 1;
c.FontSize = 30;
c.Ticks = [0 20];
title(c,'F-values')
%% do omnibus again and test
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
figure; ft_multiplotTFR(cfg, ga_sb);
%%
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
omnibus = ft_math(cfg,ga_sb,ga_nb);
%%
% close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotTFR(cfg, omnibus);
%% compute omnibus stat test freq time and elec
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
n = numel(load3);  % number of subjects (paired samples)
cohens_d = stat.stat ./ sqrt(n);
stat.effectsize = cohens_d;
%% now compute stat but only for the electrodes per pre registration
cfg = [];
cfg.latency          = [0 3];
cfg.channel = {'M1', 'M2', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan ='yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
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

[statprereg] = ft_freqstatistics(cfg, sb_high_low{:}, nb_high_low{:});
n = numel(load3);  % number of subjects (paired samples)
cohens_d = statprereg.stat ./ sqrt(n);
statprereg.effectsize = cohens_d;
%%

% close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'effectsize';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
% cfg.zlim = [-.8 .8];
figure; ft_multiplotTFR(cfg,stat);
%%
stattfr=statprereg;
stattfr.stat= statprereg.effectsize;
% figure;
cfg = [];
% cfg.channel = {'CP2', 'Pz','P2', 'CPP4h', 'CPP2h', 'CPz'};
cfg.avgoverchan = 'yes';
% cfg.frequency = [2 40];
% cfg.latency   = [0 1.5];
freq = ft_selectdata(cfg,stattfr);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(2, 40, 500);
mask_interp = linspace(2, 40, 500);
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
freq_flipped = fliplr(freq_interp);

target_freqs = [10 20 30 40];
ytick_idx = round(interp1(freq_flipped, 1:length(freq_flipped), target_freqs));
yticks(fliplr(ytick_idx));
yticklabels({'40','30','20','10'});

xticks([1 250 500])
xticklabels({num2str(freq.time(1)),'1.5',  num2str(freq.time(end))});
set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.5 .5]);
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [-.5 0 .5];
title(c,'Effect size \it d')

stattfr=stat;
stattfr.stat= stat.effectsize;
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
cfg.parameter = 'effectsize';
cfg.xlim = [1.308497 1.716100];
cfg.ylim = [8 11];
cfg.zlim = [-.5 .5];
cfg.marker             = 'off';
cfg.highlight          = 'on';
% cfg.highlightchannel = {'CP2', 'Pz','P2', 'CPP4h', 'CPP2h', 'CPz'};
cfg.highlightchannel = {'M1', 'M2', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.highlightsymbol    = '.';
 cfg.highlightsize      = 14;
 cfg.comment = 'no';
% figure; 
subplot(2,1,2);ft_topoplotTFR(cfg,stat);
%% extract power spectra omnibus
cfg = [];
cfg.latency = [1 2];
cfg.avgovertime = 'yes';
for s = 1 :length(subjects)
sb_hl_pow{s}= ft_selectdata(cfg,sb_high_low{s});
nb_hl_pow{s}= ft_selectdata(cfg,nb_high_low{s});
sb_hl_pow{s}.dimord = 'chan_freq';
sb_hl_pow{s} = rmfield(sb_hl_pow{s},'time');
nb_hl_pow{s}.dimord = 'chan_freq';
nb_hl_pow{s} = rmfield(nb_hl_pow{s},'time');
end
%% sternberg pow per condition
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_hl_pow = ft_freqgrandaverage(cfg,sb_hl_pow{:});
ga_nb_hl_pow = ft_freqgrandaverage(cfg,nb_hl_pow{:});
%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotER(cfg, ga_sb_hl_pow,ga_nb_hl_pow);
%% plot with SE sternberg
% close all
figure;
subplot(2,1,1);
cfg = [];
cfg.channel ={'CP2', 'Pz','P2', 'CPP4h', 'CPP2h', 'CPz'};
% cfg.channel = {'M1', 'M2', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};

cfg.avgoverchan = 'yes';
tlk_sb_ind        = ft_selectdata(cfg,ga_sb_hl_pow);
tlk_nb_ind        = ft_selectdata(cfg,ga_nb_hl_pow);

% plot load sb
x = tlk_sb_ind.freq'; % x-axis def
y = mean(squeeze(tlk_sb_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk_sb_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound


hp1 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r', 'HandleVisibility', 'off');
hold on;
hl1 = line(x, y);
set(hp1, 'facecolor', [0.97, 0.26, 0.26], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl1, 'color', [0.97, 0.26, 0.26], 'linewidth', 2);

% plot load nb
x = tlk_nb_ind.freq'; % x-axis def
y = mean(squeeze(tlk_nb_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk_nb_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);


% Label the axes
set(gca, 'FontSize', 22);
title('');
xlabel('Frequency [Hz]');
ylabel("Power change \newline from baseline");
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Sternberg high-low','N-back high-low'}, 'Location','southeast','Fontsize',20);

%% extract values sternberg
for s = 1:length(load6)
    % select retention
    cfg = [];
%     cfg.latency = [0.2 0.6]+1.5;
%     cfg.frequency = [10 16];
        cfg.latency = [0 3];
    cfg.frequency = [8 14];
    cfg.channel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
   cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.channel = {'CP5', 'P5', 'TP7', 'CPP5h', 'TPP7h'};
cfg.channel = {'P5', 'PPO5h'};% based on F stat
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
fill(positions(1) + f_sb2*0.05, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(sb4);
fill(positions(2) + f_sb4*0.05, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(sb6);
fill(positions(3) + f_sb6*0.05, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2', 'load 4', 'load 6'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2) + (rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3) + (rand(size(sb6))-0.5)*jitter, sb6, 'r.');

% Axis & Label
% ylabel('\muV');
ylabel('\alpha power [change from bl] ');
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
cfg.channel = {'CP5', 'P5', 'TP7', 'CPP5h', 'TPP7h'};
cfg.channel = {'P7', 'PPO9h'};% based on Fstat
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
fill(positions(1) + f_sb2*0.05, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(nb2);
fill(positions(2) + f_sb4*0.05, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(nb3);
fill(positions(3) + f_sb6*0.05, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 1', 'load 2', 'load 3'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2) + (rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3) + (rand(size(sb6))-0.5)*jitter, sb6, 'r.');

% Axis & Label
% ylabel('\muV');
ylabel('\alpha power [change from bl] ');
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
