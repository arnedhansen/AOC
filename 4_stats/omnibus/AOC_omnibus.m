%% AOC Omnibus — Cluster Stats and Alpha Plots
% Loads omnibus_data.mat. Runs cluster-based permutation (N-back load, Sternberg vs N-back), extracts ROI alpha, raincloud/box/scatter plots, paired t-tests. Produces figures and stats.
%
% Key outputs:
%   Cluster permutation results; TFR/alpha figures; rainclouds with significance; paired t-tests

%% Setup
startup
[subjects, path, colors, headmodel] = setup('AOC');

% Set up save directories
control_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/omnibus';
figures_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/omnibus';
if ~exist(control_dir, 'dir'), mkdir(control_dir); end
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end 

%% Load variables
tic
disp('Loading omnibus data...')
if ispc
    load W:\Students\Arne\AOC\data\features\omnibus_data.mat
else
    load /Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/omnibus_data.mat
end 
toc

%% Visualize N-back
% close all
% cfg = [];
% cfg.figure = 'gcf';
% %cfg.ylim = [3 40];
% %cfg.zlim = [-2 2];
% cfg.layout = headmodel.layANThead;
% figure; ft_multiplotTFR(cfg, ga_nb);

%% Compute omnibus GA Sternberg vs. GA N-back
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
%omnibus = ft_math(cfg,ga_sb,ga_nb);

%% Visualize Omnibus
% close all
% cfg = [];
% cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% % cfg.zlim = [-3 3];
% cfg.layout = headmodel.layANThead;
% figure; ft_multiplotTFR(cfg, omnibus);

%% Sternberg per condition
cfg = [];
ga_sb_2 = ft_freqgrandaverage(cfg,load2{:});
ga_sb_4 = ft_freqgrandaverage(cfg,load4{:});
ga_sb_6 = ft_freqgrandaverage(cfg,load6{:});
ga_sb   = ft_freqgrandaverage(cfg,load2{:},load4{:},load6{:});
ga_nb_1 = ft_freqgrandaverage(cfg,load1nb{:});
ga_nb_2 = ft_freqgrandaverage(cfg,load2nb{:});
ga_nb_3 = ft_freqgrandaverage(cfg,load3nb{:});
ga_nb   = ft_freqgrandaverage(cfg,load1nb{:},load2nb{:},load3nb{:});

% %%
% close all
% cfg = [];
% cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% % cfg.zlim = [-3 3];
% cfg.layout = headmodel.layANThead;
% figure; ft_multiplotTFR(cfg, ga_nb);
% 
% %% single 
% cfg = [];
% cfg.figure = 'gcf';
% cfg.zlim = [-.2 .2];
% cfg.xlim = [-1 2];
% cfg.channel = {'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb
% cfg.layout = headmodel.layANThead;
% figure; ft_singleplotTFR(cfg, ga_nb);

%% plot SB tfr and topo
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
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

set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.2 .2]);
color_map = flipud(cbrewer('div', 'RdBu', 64)); % Red-Blue diverging color map
colormap(color_map);
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 18;
cb.Ticks = [-.2 0 .2];
title(cb,"Power change \newline from baseline")
xline(0,'k--','LineWidth',2); % black dashed line at 0 sec
cfg = [];
cfg.layout = headmodel.layANThead;
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
colormap(color_map); % Use same colormap for topo
% Save figure
sgtitle('Sternberg TFR and Topography', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_TFR_topo.png'));

%% plot NB tfr and topo
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
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

subplot(2,1,2);ft_plot_matrix(flip(pow_interp));
% subplot(2,1,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))));
ax = gca; hold(ax,'on');

% map 0 sec to the interpolated column index
x0 = interp1(tim_interp, 1:numel(tim_interp), 0, 'linear', 'extrap');
xline(ax, x0, 'k-', 'LineWidth', 1);

% ticks (still index-based)
%xticks( round(interp1(tim_interp, 1:numel(tim_interp), [-1 0 1 2])) );
xticklabels({'-1','0','1','2'});
yticks([1 125 250 375 ]); % positions in the interpolated grid
yticklabels({'40','30','20','10'}); % corresponding freq values

set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.2 .2]);
color_map = flipud(cbrewer('div', 'RdBu', 64)); % Red-Blue diverging color map
colormap(color_map);
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 18;
cb.Ticks = [-.2 0 .2];
title(cb,"Power change \newline from baseline")
xline(0,'k--','LineWidth',2); % black dashed line at 0 sec
cfg = [];
cfg.layout = headmodel.layANThead;
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
colormap(color_map); % Use same colormap for topo
% Save figure
sgtitle('N-back TFR and Topography', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_TFR_topo.png'));

%% plot SB posterior
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
cfg = [];
cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};
% cfg.channel = {'P8', 'PO4', 'PO8', 'PPO6h', 'POO4h'};
cfg.figure = 'gcf';
cfg.ylim = [3 40];
cfg.zlim = [-.2 .2];
cfg.xlim = [-.5 3];
subplot(3,1,1); ft_singleplotTFR(cfg,ga_sb_2);
subplot(3,1,2); ft_singleplotTFR(cfg,ga_sb_4);
subplot(3,1,3); ft_singleplotTFR(cfg,ga_sb_6);
color_map = flipud(cbrewer('div', 'RdBu', 64));
colormap(color_map);
% Save figure
sgtitle('Sternberg Posterior TFR by Load', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_posterior_TFR.png'));

%% plot NB posterior
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
cfg = [];
% cfg.channel = {'Pz', 'POz', 'P2', 'PPO2'};
cfg.channel = {'POz', 'PO3', 'PO4', 'PPO1', 'PPO6h'};
cfg.figure = 'gcf';
cfg.ylim = [3 40];
cfg.zlim = [-.2 .2];
cfg.xlim = [-.5 2];
subplot(3,1,1); ft_singleplotTFR(cfg,ga_nb_1);
subplot(3,1,2); ft_singleplotTFR(cfg,ga_nb_2);
subplot(3,1,3); ft_singleplotTFR(cfg,ga_nb_3);
color_map = flipud(cbrewer('div', 'RdBu', 64));
colormap(color_map);
% Save figure
sgtitle('N-back Posterior TFR by Load', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_posterior_TFR.png'));

%% Control: Check baseline stability and data quality
disp('Running data quality controls...');
% Check baseline values across conditions
cfg_base = [];
cfg_base.latency = [-0.5 -0.25];
cfg_base.avgovertime = 'yes';
cfg_base.avgoverchan = 'yes';
cfg_base.frequency = [8 14];  % alpha band
cfg_base.avgoverfreq = 'yes';

baseline_sb2 = zeros(length(subjects), 1);
baseline_sb4 = zeros(length(subjects), 1);
baseline_sb6 = zeros(length(subjects), 1);
baseline_nb1 = zeros(length(subjects), 1);
baseline_nb2 = zeros(length(subjects), 1);
baseline_nb3 = zeros(length(subjects), 1);

for s = 1:length(subjects)
    tmp = ft_selectdata(cfg_base, load2{s});
    baseline_sb2(s) = tmp.powspctrm;
    tmp = ft_selectdata(cfg_base, load4{s});
    baseline_sb4(s) = tmp.powspctrm;
    tmp = ft_selectdata(cfg_base, load6{s});
    baseline_sb6(s) = tmp.powspctrm;
    tmp = ft_selectdata(cfg_base, load1nb{s});
    baseline_nb1(s) = tmp.powspctrm;
    tmp = ft_selectdata(cfg_base, load2nb{s});
    baseline_nb2(s) = tmp.powspctrm;
    tmp = ft_selectdata(cfg_base, load3nb{s});
    baseline_nb3(s) = tmp.powspctrm;
end

% Plot baseline stability
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
subplot(2,1,1);
boxplot([baseline_sb2, baseline_sb4, baseline_sb6], 'Labels', {'Load 2', 'Load 4', 'Load 6'});
ylabel('Baseline Power (alpha band)');
title('Sternberg Baseline Stability');
grid on;

subplot(2,1,2);
boxplot([baseline_nb1, baseline_nb2, baseline_nb3], 'Labels', {'Load 1', 'Load 2', 'Load 3'});
ylabel('Baseline Power (alpha band)');
title('N-back Baseline Stability');
grid on;
sgtitle('Baseline Stability Check', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(control_dir, 'AOC_omnibus_baseline_stability.png'));

% Check for baseline differences (should be ~0 after baseline correction)
fprintf('Baseline means (should be ~0):\n');
fprintf('  SB Load 2: %.4f, Load 4: %.4f, Load 6: %.4f\n', mean(baseline_sb2), mean(baseline_sb4), mean(baseline_sb6));
fprintf('  NB Load 1: %.4f, Load 2: %.4f, Load 3: %.4f\n', mean(baseline_nb1), mean(baseline_nb2), mean(baseline_nb3));

%% extract power spectra SB
cfg = [];
cfg.latency = [1 3];
cfg.avgovertime = 'yes';
for subj = 1 :length(subjects)
load2_pow{subj}= ft_selectdata(cfg,load2{subj});
load4_pow{subj}= ft_selectdata(cfg,load4{subj});
load6_pow{subj}= ft_selectdata(cfg,load6{subj});
load2_pow{subj}.dimord = 'chan_freq';
load2_pow{subj} = rmfield(load2_pow{subj},'time');
load4_pow{subj}.dimord = 'chan_freq';
load4_pow{subj} = rmfield(load4_pow{subj},'time');
load6_pow{subj}.dimord = 'chan_freq';
load6_pow{subj} = rmfield(load6_pow{subj},'time');
end
cfg.latency = [.5 2];
for subj = 1 :length(subjects)
load1nb_pow{subj}= ft_selectdata(cfg,load1nb{subj});
load2nb_pow{subj}= ft_selectdata(cfg,load2nb{subj});
load3nb_pow{subj}= ft_selectdata(cfg,load3nb{subj});
load1nb_pow{subj}.dimord = 'chan_freq';
load1nb_pow{subj} = rmfield(load1nb_pow{subj},'time');
load2nb_pow{subj}.dimord = 'chan_freq';
load2nb_pow{subj} = rmfield(load2nb_pow{subj},'time');
load3nb_pow{subj}.dimord = 'chan_freq';
load3nb_pow{subj} = rmfield(load3nb_pow{subj},'time');
end

%% sternberg pow per condition
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2pow = ft_freqgrandaverage(cfg,load2_pow{:});
ga_sb_4pow = ft_freqgrandaverage(cfg,load4_pow{:});
ga_sb_6pow = ft_freqgrandaverage(cfg,load6_pow{:});
ga_nb_2pow = ft_freqgrandaverage(cfg,load2nb_pow{:});
ga_nb_1pow = ft_freqgrandaverage(cfg,load1nb_pow{:});
ga_nb_3pow = ft_freqgrandaverage(cfg,load3nb_pow{:});

%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = headmodel.layANThead;
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotER(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_power_spectra_multiplot.fig'));

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
e = std(squeeze(tlk6_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound


% plot load 6 (highest - red)
hp1 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r', 'HandleVisibility', 'off');
hold on;
hl1 = line(x, y);
set(hp1, 'facecolor', [0.97, 0.26, 0.26], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl1, 'color', [0.97, 0.26, 0.26], 'linewidth', 2);

% plot load 4 (middle - green)
x = tlk4_ind.freq'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'g', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.2, 0.8, 0.2], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.2, 0.8, 0.2], 'linewidth', 2);

% plot load 2 (lowest - blue)
x = tlk2_ind.freq'; % x-axis def
y = mean(squeeze(tlk2_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk2_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound

hp3 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl3 = line(x, y);
set(hp3, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl3, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% Label the axes
set(gca, 'FontSize', 22);
title('');
xlabel('Frequency [Hz]');
ylabel("Power change \newline from baseline");
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',20);
% Save figure
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_power_spectra_SE.png'));
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
% plot load 3 (highest - red)
x = tlk6_ind.freq'; % x-axis def
y = mean(squeeze(tlk6_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk6_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound

hp1 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r', 'HandleVisibility', 'off');
hold on;
hl1 = line(x, y);
set(hp1, 'facecolor', [0.97, 0.26, 0.26], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl1, 'color', [0.97, 0.26, 0.26], 'linewidth', 2);

% plot load 2 (middle - green)
x = tlk4_ind.freq'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'g', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.2, 0.8, 0.2], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.2, 0.8, 0.2], 'linewidth', 2);

% plot load 1 (lowest - blue)
x = tlk2_ind.freq'; % x-axis def
y = mean(squeeze(tlk2_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk2_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
low = y - e; % lower bound
high = y + e; % upper bound

hp3 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl3 = line(x, y);
set(hp3, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl3, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% Label the axes
set(gca, 'FontSize', 22);
title('');
xlabel('Frequency [Hz]');
ylabel("Power change \newline from baseline");
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','southeast','Fontsize',20);
% Save figure
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_power_spectra_SE.png'));
%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = headmodel.layANThead;
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotER(cfg, ga_nb_1pow,ga_nb_2pow,ga_nb_3pow);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_power_spectra_multiplot.fig'));
%%
load('/Volumes/g_psyplafor_methlab$/Students/Arne/toolboxes/headmodel/elec_aligned.mat');
cfg =[];
cfg.method ='distance';
cfg.elec = elec_aligned;
cfg.layout = headmodel.layANThead;
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
cfg.alpha            = 0.05;       % alpha level of the permutation test (two-tailed)
cfg.numrandomization = 1000;        % number of draws from the permutation distribution

n_U  = numel(subjects);
n_P  = numel(subjects);
n_N = numel(subjects);

clear design
design = zeros(2,3*n_U);
cfg.design(1,:)           = [ones(1,n_U), ones(1,n_P)*2,ones(1,n_N)*3]; % design matrix
cfg.design(2,:)           = [1:n_U,1:n_P, 1:n_N]; 
cfg.ivar                  = 1; % number or list with indices indicating the independent variable(subj)
cfg.uvar                  = 2;% units-of-observation (subjects or trials
[statFnb] = ft_freqstatistics(cfg, ga_nb_1pow,ga_nb_2pow,ga_nb_3pow);
[statFsb] = ft_freqstatistics(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
%%
close all
cfg = [];
cfg.layout = headmodel.layANThead;
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
% cfg.maskstyle = 'outline';
% cfg.zlim = [-.8 .8];
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotER(cfg,statFsb);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_Fstat_multiplot.fig'));
%% 

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
cfg.linecolor = 'k';
cfg.linewidth = 2;
cfg.comment = 'no';
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotER(cfg,statFnb);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_Fstat_multiplot.fig'));
%%
cfg = [];
cfg.layout = headmodel.layANThead;
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

figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_topoplotER(cfg,statFnb);
caxis([0 20]);
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 30;
cb.Ticks = [0 20];
title(cb,'F-values')
% Save figure
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_nback_Fstat_topoplot.png'));
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
cfg.layout = headmodel.layANThead;
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotTFR(cfg, ga_sb);
color_map = flipud(cbrewer('div', 'RdBu', 64));
colormap(color_map);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_stern_highlow_TFR_multiplot.fig'));
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
cfg.layout = headmodel.layANThead;
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotTFR(cfg, omnibus);
color_map = flipud(cbrewer('div', 'RdBu', 64));
colormap(color_map);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternvsnback_TFR_multiplot.fig'));
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
cfg.alpha            = 0.05;       % alpha level of the permutation test (two-tailed)
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
n_subj = numel(load6);
design = zeros(2,2*n_subj);
for i = 1:n_subj
  design(1,i) = i;
end
for i = 1:n_subj
  design(1,n_subj+i) = i;
end
design(2,1:n_subj)        = 1;
design(2,n_subj+1:2*n_subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, sb_high_low{:}, nb_high_low{:});
% Calculate Cohen's d properly for paired samples
% Extract data values for effect size calculation
n = numel(load6);  % number of subjects (paired samples)
% Calculate mean difference and SD of differences for Cohen's d
% For each time-frequency point, calculate d = mean_diff / SD_diff
n_freq = numel(stat.freq);
n_time = numel(stat.time);
sb_data = zeros(n, n_freq, n_time);
nb_data = zeros(n, n_freq, n_time);
for s = 1:n
    % Select matching time and frequency ranges (use ranges, not vectors)
    cfg_sel = [];
    cfg_sel.latency = [min(stat.time) max(stat.time)];
    cfg_sel.frequency = [min(stat.freq) max(stat.freq)];
    cfg_sel.avgoverchan = 'yes';
    sb_sel = ft_selectdata(cfg_sel, sb_high_low{s});
    nb_sel = ft_selectdata(cfg_sel, nb_high_low{s});
    
    % Extract data and match to stat dimensions
    sb_tmp = squeeze(sb_sel.powspctrm);
    nb_tmp = squeeze(nb_sel.powspctrm);
    
    % If we got a scalar (averaged everything), we need to extract differently
    if isscalar(sb_tmp) || (ndims(sb_tmp) == 2 && (size(sb_tmp, 1) ~= n_freq || size(sb_tmp, 2) ~= n_time))
        % Extract data channel by channel and average, matching stat dimensions
        % Get one channel to check structure
        cfg_check = [];
        cfg_check.latency = [min(stat.time) max(stat.time)];
        cfg_check.frequency = [min(stat.freq) max(stat.freq)];
        cfg_check.channel = sb_high_low{s}.label(1);
        tmp_check = ft_selectdata(cfg_check, sb_high_low{s});
        if ndims(tmp_check.powspctrm) == 3 && size(tmp_check.powspctrm, 1) == 1
            % Extract all channels and average
            all_sb = zeros(length(sb_high_low{s}.label), length(tmp_check.freq), length(tmp_check.time));
            all_nb = zeros(length(nb_high_low{s}.label), length(tmp_check.freq), length(tmp_check.time));
            for ch = 1:length(sb_high_low{s}.label)
                cfg_ch = [];
                cfg_ch.latency = [min(stat.time) max(stat.time)];
                cfg_ch.frequency = [min(stat.freq) max(stat.freq)];
                cfg_ch.channel = sb_high_low{s}.label(ch);
                tmp_sb = ft_selectdata(cfg_ch, sb_high_low{s});
                tmp_nb = ft_selectdata(cfg_ch, nb_high_low{s});
                all_sb(ch, :, :) = squeeze(tmp_sb.powspctrm);
                all_nb(ch, :, :) = squeeze(tmp_nb.powspctrm);
            end
            sb_tmp = squeeze(mean(all_sb, 1));  % average over channels
            nb_tmp = squeeze(mean(all_nb, 1));  % average over channels
        end
    end
    
    % Match dimensions to stat structure
    if ndims(sb_tmp) == 2
        if size(sb_tmp, 1) == n_freq && size(sb_tmp, 2) == n_time
            % Exact match
            sb_data(s, :, :) = sb_tmp;
            nb_data(s, :, :) = nb_tmp;
        else
            % Need to interpolate to match stat dimensions
            [freq_orig, time_orig] = meshgrid(sb_sel.freq, sb_sel.time);
            [freq_stat, time_stat] = meshgrid(stat.freq, stat.time);
            sb_data(s, :, :) = interp2(freq_orig, time_orig, sb_tmp', freq_stat, time_stat, 'linear', NaN)';
            nb_data(s, :, :) = interp2(freq_orig, time_orig, nb_tmp', freq_stat, time_stat, 'linear', NaN)';
        end
    else
        error('Unexpected dimensions in powspctrm');
    end
end
% Calculate differences (Sternberg - N-back)
diff_data = sb_data - nb_data;
mean_diff = squeeze(mean(diff_data, 1, 'omitnan'));  % mean across subjects [freq x time]
sd_diff = squeeze(std(diff_data, 0, 1, 'omitnan'));   % sample SD across subjects [freq x time]
% Cohen's d = mean_diff / sd_diff (for paired samples)
cohens_d = mean_diff ./ (sd_diff + eps);  % add eps to avoid division by zero
% Expand effectsize to match stat dimensions [chan×freq×time]
% Since effectsize was calculated channel-averaged, replicate across all channels
% cohens_d is [freq×time] = [136×41], need [chan×freq×time] = [125×136×41]
cohens_d_3d = reshape(cohens_d, [1, size(cohens_d, 1), size(cohens_d, 2)]); % [1×136×41]
stat.effectsize = repmat(cohens_d_3d, [length(stat.label), 1, 1]); % [125×136×41]
%% now compute stat but only for the electrodes per pre registration
cfg                  = [];
cfg.latency          = [0 2];
cfg.channel          = {'M1', 'M2', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;       % alpha level of the permutation test (two-tailed)
cfg.numrandomization = 1000;
n_subj = numel(load6);
design = zeros(2,2*n_subj);
for i = 1:n_subj
  design(1,i) = i;
end
for i = 1:n_subj
  design(1,n_subj+i) = i;
end
design(2,1:n_subj)        = 1;
design(2,n_subj+1:2*n_subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[statprereg] = ft_freqstatistics(cfg, sb_high_low{:}, nb_high_low{:});
% Calculate Cohen's d properly for paired samples
n = numel(load6);  % number of subjects (paired samples)
% Extract data values for effect size calculation
n_freq = numel(statprereg.freq);
n_time = numel(statprereg.time);
sb_data = zeros(n, n_freq, n_time);
nb_data = zeros(n, n_freq, n_time);
for s = 1:n
    % Select matching time and frequency ranges (use ranges, not vectors)
    cfg_sel = [];
    cfg_sel.latency = [min(statprereg.time) max(statprereg.time)];
    cfg_sel.frequency = [min(statprereg.freq) max(statprereg.freq)];
    cfg_sel.avgoverchan = 'yes';
    sb_sel = ft_selectdata(cfg_sel, sb_high_low{s});
    nb_sel = ft_selectdata(cfg_sel, nb_high_low{s});
    
    % Extract data and match to stat dimensions
    sb_tmp = squeeze(sb_sel.powspctrm);
    nb_tmp = squeeze(nb_sel.powspctrm);
    
    % If we got a scalar (averaged everything), we need to extract differently
    if isscalar(sb_tmp) || (ndims(sb_tmp) == 2 && (size(sb_tmp, 1) ~= n_freq || size(sb_tmp, 2) ~= n_time))
        % Extract data channel by channel and average, matching stat dimensions
        % Get one channel to check structure
        cfg_check = [];
        cfg_check.latency = [min(statprereg.time) max(statprereg.time)];
        cfg_check.frequency = [min(statprereg.freq) max(statprereg.freq)];
        cfg_check.channel = sb_high_low{s}.label(1);
        tmp_check = ft_selectdata(cfg_check, sb_high_low{s});
        if ndims(tmp_check.powspctrm) == 3 && size(tmp_check.powspctrm, 1) == 1
            % Extract all channels and average
            all_sb = zeros(length(sb_high_low{s}.label), length(tmp_check.freq), length(tmp_check.time));
            all_nb = zeros(length(nb_high_low{s}.label), length(tmp_check.freq), length(tmp_check.time));
            for ch = 1:length(sb_high_low{s}.label)
                cfg_ch = [];
                cfg_ch.latency = [min(statprereg.time) max(statprereg.time)];
                cfg_ch.frequency = [min(statprereg.freq) max(statprereg.freq)];
                cfg_ch.channel = sb_high_low{s}.label(ch);
                tmp_sb = ft_selectdata(cfg_ch, sb_high_low{s});
                tmp_nb = ft_selectdata(cfg_ch, nb_high_low{s});
                all_sb(ch, :, :) = squeeze(tmp_sb.powspctrm);
                all_nb(ch, :, :) = squeeze(tmp_nb.powspctrm);
            end
            sb_tmp = squeeze(mean(all_sb, 1));  % average over channels
            nb_tmp = squeeze(mean(all_nb, 1));  % average over channels
        end
    end
    
    % Match dimensions to stat structure
    if ndims(sb_tmp) == 2
        if size(sb_tmp, 1) == n_freq && size(sb_tmp, 2) == n_time
            % Exact match
            sb_data(s, :, :) = sb_tmp;
            nb_data(s, :, :) = nb_tmp;
        else
            % Need to interpolate to match stat dimensions
            [freq_orig, time_orig] = meshgrid(sb_sel.freq, sb_sel.time);
            [freq_stat, time_stat] = meshgrid(statprereg.freq, statprereg.time);
            sb_data(s, :, :) = interp2(freq_orig, time_orig, sb_tmp', freq_stat, time_stat, 'linear', NaN)';
            nb_data(s, :, :) = interp2(freq_orig, time_orig, nb_tmp', freq_stat, time_stat, 'linear', NaN)';
        end
    else
        error('Unexpected dimensions in powspctrm');
    end
end
% Calculate differences (Sternberg - N-back)
diff_data = sb_data - nb_data;
mean_diff = squeeze(mean(diff_data, 1, 'omitnan'));  % mean across subjects [freq x time]
sd_diff = squeeze(std(diff_data, 0, 1, 'omitnan'));   % sample SD across subjects [freq x time]
% Cohen's d = mean_diff / sd_diff (for paired samples)
cohens_d = mean_diff ./ (sd_diff + eps);  % add eps to avoid division by zero
% Expand effectsize to match statprereg dimensions [chan×freq×time]
% Since effectsize was calculated channel-averaged, replicate across all channels
% cohens_d is [freq×time], need [chan×freq×time]
cohens_d_3d = reshape(cohens_d, [1, size(cohens_d, 1), size(cohens_d, 2)]); % [1×freq×time]
statprereg.effectsize = repmat(cohens_d_3d, [length(statprereg.label), 1, 1]); % [chan×freq×time]
%%

% close all
cfg = [];
cfg.layout = headmodel.layANThead;
cfg.parameter = 'effectsize';
cfg.maskparameter ='mask';
cfg.maskstyle = 'outline';
% cfg.zlim = [-.8 .8];
figure; ft_multiplotTFR(cfg,stat);

%% do stats
stattfr=statprereg;
stattfr.stat= statprereg.effectsize;
% figure;
cfg = [];
cfg.channel = {'CP2', 'Pz','P2', 'CPP4h', 'CPP2h', 'CPz'};
cfg.avgoverchan = 'yes';
%cfg.frequency = [2 40];
%cfg.latency   = [0 2];

% cfg.latency          = [0 3];
% cfg.channel          = {'M1', 'M2', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.05;
% cfg.numrandomization = 1000;
% cfg.neighbours       = neighbours;
% cfg.minnbchan        = 2;

freq = ft_selectdata(cfg,stattfr);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 500);
freq_interp = linspace(2, 40, 500);
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
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
% subplot(2,1,2);ft_plot_matrix(flip(pow_interp));
subplot(2,1,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))));
freq_flipped = fliplr(freq_interp);

target_freqs = [10 20 30 40];
ytick_idx = round(interp1(freq_flipped, 1:length(freq_flipped), target_freqs));
yticks(fliplr(ytick_idx));
yticklabels({'40','30','20','10'});

xticks([1 250 500])
xticklabels({num2str(freq.time(1)),'1.5',  num2str(freq.time(end))});
set(gca,'Fontsize',20);
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
caxis([-.5 .5]);
color_map = flipud(cbrewer('div', 'RdBu', 64));
colormap(color_map);
cb = colorbar;
cb.LineWidth = 1;
cb.FontSize = 18;
cb.Ticks = [-.5 0 .5];
title(cb,'Effect size \it d')
% Save figure
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternvsnback_effectsize_TFR.png'));

stattfr=stat;
stattfr.stat= stat.effectsize;
cfg = [];
cfg.layout = headmodel.layANThead;
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
colormap(color_map); % Use same colormap for topo
% Save figure
sgtitle('Sternberg vs N-back Effect Size', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_sternvsnback_effectsize_topo.png'));
%% extract power spectra omnibus
cfg = [];
cfg.latency = [1 2];
cfg.avgovertime = 'yes';
for subj = 1 :length(subjects)
sb_hl_pow{subj}= ft_selectdata(cfg,sb_high_low{subj});
nb_hl_pow{subj}= ft_selectdata(cfg,nb_high_low{subj});
sb_hl_pow{subj}.dimord = 'chan_freq';
sb_hl_pow{subj} = rmfield(sb_hl_pow{subj},'time');
nb_hl_pow{subj}.dimord = 'chan_freq';
nb_hl_pow{subj} = rmfield(nb_hl_pow{subj},'time');
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
cfg.layout = headmodel.layANThead;
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
ft_multiplotER(cfg, ga_sb_hl_pow,ga_nb_hl_pow);
% Save figure (multiplot: .fig only)
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_highlow_power_spectra_multiplot.fig'));
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
e = std(squeeze(tlk_sb_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
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
e = std(squeeze(tlk_nb_ind.powspctrm), 0)' ./ sqrt(numel(subjects)); % SEM (using sample SD)
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
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Sternberg high-low','N-back high-low'}, 'Location','southeast','Fontsize',20);
% Save figure
saveas(gcf, fullfile(figures_dir, 'AOC_omnibus_highlow_power_spectra_SE.png'));

%% timewins
tsb = [1 2];
tnb = [0 2];

%% extract values sternberg
clc
for subj = 1:length(load6)
    % select retention
    cfg = [];
    cfg.latency = tsb;
    cfg.frequency = [8 14];
    cfg.channel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'}; elecsName = 'prereg';
    % cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
    % cfg.channel = {'CP5', 'P5', 'TP7', 'CPP5h', 'TPP7h'};
    % cfg.channel = {'P5', 'PPO5h'};% based on F stat
    cfg.avgoverfreq = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    val6{subj} = ft_selectdata(cfg,load6{subj});
    sb6(subj) = val6{subj}.powspctrm;
    
    val4{subj} = ft_selectdata(cfg,load4{subj});
    sb4(subj) = val4{subj}.powspctrm;
    
    val2{subj} = ft_selectdata(cfg,load2{subj});
    sb2(subj) = val2{subj}.powspctrm;
end

%% exclude outliers
% exclude outliers from sb6, sb4 and sb2 using z-score
sb2 = sb2(:);
sb4 = sb4(:);
sb6 = sb6(:);

% Control: Visualize outliers before exclusion
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
subplot(2,2,1);
scatter(1:length(sb2), sb2, 'filled', 'k'); hold on;
scatter(1:length(sb4), sb4, 'filled', 'b');
scatter(1:length(sb6), sb6, 'filled', 'r');
xlabel('Subject'); ylabel('Alpha Power');
title('Sternberg: Before Outlier Exclusion');
legend({'Load 2', 'Load 4', 'Load 6'}, 'Location', 'best');
grid on;

Z = zscore([sb2 sb4 sb6], 0, 1); % z per condition across subjects
keepIdx = all(abs(Z) < 2, 2); % keep subjects that are not outliers in any condition

% Control: Show which subjects are excluded
n_excluded = sum(~keepIdx);
fprintf('Sternberg: Excluding %d subjects (%.1f%%) as outliers\n', n_excluded, 100*n_excluded/length(keepIdx));
if n_excluded > 0
    fprintf('  Excluded subject indices: %s\n', num2str(find(~keepIdx)));
end

subplot(2,2,2);
imagesc(Z); colorbar;
xlabel('Condition'); ylabel('Subject');
title('Z-scores (Sternberg)');
xticklabels({'Load 2', 'Load 4', 'Load 6'});
hold on;
for i = 1:length(keepIdx)
    if ~keepIdx(i)
        plot([0.5 3.5], [i i], 'r-', 'LineWidth', 2);
    end
end

sb2 = sb2(keepIdx);
sb4 = sb4(keepIdx);
sb6 = sb6(keepIdx);
% Note: After outlier exclusion, all arrays should have the same length

% Control: Visualize after exclusion
subplot(2,2,3);
scatter(1:length(sb2), sb2, 'filled', 'k'); hold on;
scatter(1:length(sb4), sb4, 'filled', 'b');
scatter(1:length(sb6), sb6, 'filled', 'r');
xlabel('Subject'); ylabel('Alpha Power');
title('Sternberg: After Outlier Exclusion');
legend({'Load 2', 'Load 4', 'Load 6'}, 'Location', 'best');
grid on;

subplot(2,2,4);
Z_after = zscore([sb2 sb4 sb6], 0, 1);
imagesc(Z_after); colorbar;
xlabel('Condition'); ylabel('Subject');
title('Z-scores After Exclusion');
xticklabels({'Load 2', 'Load 4', 'Load 6'});
sgtitle('Sternberg Outlier Exclusion Control', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(control_dir, 'AOC_omnibus_stern_outlier_exclusion.png'));

%%
close all
figure(30);
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
clf;
subplot(1,2,2);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(sb2);
fill(positions(1) + f_sb2*0.3, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(sb4);
fill(positions(2) + f_sb4*0.3, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(sb6);
fill(positions(3) + f_sb6*0.3, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

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
tsb = cfg.latency;
title(['Sternberg (time: ', num2str(tsb), ', elecs ',elecsName, ')']);
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
% ylim([min([sb2(:); sb4(:); sb6(:)]) - 0.5, y_max + 0.1]);
ylim([-y_max y_max])
xlim([0 1.3])

% 
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
for subj = 1:length(load3nb)
    % select retention
    cfg = [];
    cfg.latency = tnb;
    cfg.frequency = [8 14];
    cfg.channel   = {'O1', 'O2', 'PO7', 'PO8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'}; elecsName = 'prereg';
    % cfg.channel = {'M2', 'CP5', 'P7', 'P3', 'Pz', 'P4', 'P8', 'POz', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO3', 'PO4', 'TP7', 'TP8', 'PO7', 'PO8', 'TPP9h', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'CPP5h', 'CPP3h', 'CPP4h', 'CPP6h', 'PPO1', 'PPO2', 'I1', 'Iz', 'I2', 'TPP7h', 'CPP1h', 'CPP2h', 'TPP8h', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
    % cfg.channel = {'CP5', 'P5', 'TP7', 'CPP5h', 'TPP7h'};
    % cfg.channel = {'P7', 'PPO9h'};% based on Fstat
    cfg.avgoverfreq = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    val3{subj} = ft_selectdata(cfg,load3nb{subj});
    nb3(subj) = val3{subj}.powspctrm;
    
    val2nb{subj} = ft_selectdata(cfg,load2nb{subj});
    nb2(subj) = val2nb{subj}.powspctrm;
    
    val1{subj} = ft_selectdata(cfg,load1nb{subj});
    nb1(subj) = val1{subj}.powspctrm;
end

%% exclude outliers
nb1 = nb1(:);
nb2 = nb2(:);
nb3 = nb3(:);

% Control: Visualize outliers before exclusion
figure;
set(gcf, 'Position', [0, 0, 1512, 982], 'Color', 'w');
subplot(2,2,1);
scatter(1:length(nb1), nb1, 'filled', 'k'); hold on;
scatter(1:length(nb2), nb2, 'filled', 'b');
scatter(1:length(nb3), nb3, 'filled', 'r');
xlabel('Subject'); ylabel('Alpha Power');
title('N-back: Before Outlier Exclusion');
legend({'Load 1', 'Load 2', 'Load 3'}, 'Location', 'best');
grid on;

Z = zscore([nb1 nb2 nb3], 0, 1); % z per condition across subjects
keepIdx = all(abs(Z) < 2, 2); % keep subjects that are not outliers in any condition

% Control: Show which subjects are excluded
n_excluded = sum(~keepIdx);
fprintf('N-back: Excluding %d subjects (%.1f%%) as outliers\n', n_excluded, 100*n_excluded/length(keepIdx));
if n_excluded > 0
    fprintf('  Excluded subject indices: %s\n', num2str(find(~keepIdx)));
end

subplot(2,2,2);
imagesc(Z); colorbar;
xlabel('Condition'); ylabel('Subject');
title('Z-scores (N-back)');
xticklabels({'Load 1', 'Load 2', 'Load 3'});
hold on;
for i = 1:length(keepIdx)
    if ~keepIdx(i)
        plot([0.5 3.5], [i i], 'r-', 'LineWidth', 2);
    end
end

nb1 = nb1(keepIdx);
nb2 = nb2(keepIdx);
nb3 = nb3(keepIdx);
% Note: After outlier exclusion, all arrays should have the same length

% Control: Visualize after exclusion
subplot(2,2,3);
scatter(1:length(nb1), nb1, 'filled', 'k'); hold on;
scatter(1:length(nb2), nb2, 'filled', 'b');
scatter(1:length(nb3), nb3, 'filled', 'r');
xlabel('Subject'); ylabel('Alpha Power');
title('N-back: After Outlier Exclusion');
legend({'Load 1', 'Load 2', 'Load 3'}, 'Location', 'best');
grid on;

subplot(2,2,4);
Z_after = zscore([nb1 nb2 nb3], 0, 1);
imagesc(Z_after); colorbar;
xlabel('Condition'); ylabel('Subject');
title('Z-scores After Exclusion');
xticklabels({'Load 1', 'Load 2', 'Load 3'});
sgtitle('N-back Outlier Exclusion Control', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(control_dir, 'AOC_omnibus_nback_outlier_exclusion.png'));

%%
% figure(31); clf;
subplot(1,2,1);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(nb1);
fill(positions(1) + f_sb2*0.3, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(nb2);
fill(positions(2) + f_sb4*0.3, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(nb3);
fill(positions(3) + f_sb6*0.3, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 1', 'load 2', 'load 3'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(nb1))-0.5)*jitter, nb1, 'k.');
scatter(positions(2) + (rand(size(nb2))-0.5)*jitter, nb2, 'b.');
scatter(positions(3) + (rand(size(nb3))-0.5)*jitter, nb3, 'r.');

% Axis & Label
% ylabel('\muV');
ylabel('\alpha power [change from bl] ');
xlim([0 1]);
title(['N-back (time: ', num2str(tnb), ', elecs ',elecsName, ')']);
box on;
set(gcf,'color','w');
set(gca,'Fontsize',20);

% Significance tests
[~, p_24] = ttest(nb1, nb2);
[~, p_46] = ttest(nb2, nb3);
[~, p_26] = ttest(nb1, nb3);

% Significance annotations
y_max = max([nb1(:); nb2(:); nb3(:)]) + 0.5;  % smaller margin above data overall height above
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
% ylim([min([sb2(:); sb4(:); sb6(:)]) - 0.5, y_max + 0.1]);
ylim([-y_max y_max])
xlim([0 1.3])

% Save
tnb_str = sprintf('%d', tnb);
tsb_str = sprintf('%d', tsb);
saveas(gcf, fullfile(figures_dir, ['AOC_omnibus_rainclouds_tnb', num2str(tnb_str), '_tsb', num2str(tsb_str), '.png']));

%% Control: Summary statistics and data quality report
disp('Generating summary control report...');
summary_file = fullfile(control_dir, 'AOC_omnibus_summary_report.txt');
fid = fopen(summary_file, 'w');
fprintf(fid, 'AOC Omnibus Analysis Summary Report\n');
fprintf(fid, '====================================\n\n');
fprintf(fid, 'Analysis Date: %s\n\n', datestr(now));

fprintf(fid, 'Sample Information:\n');
fprintf(fid, '  Total subjects: %d\n', length(subjects));
fprintf(fid, '  Sternberg after outlier exclusion: %d\n', length(sb2));
fprintf(fid, '  N-back after outlier exclusion: %d\n\n', length(nb1));

fprintf(fid, 'Baseline Stability (alpha band, should be ~0):\n');
fprintf(fid, '  Sternberg Load 2: %.4f (SD: %.4f)\n', mean(baseline_sb2), std(baseline_sb2));
fprintf(fid, '  Sternberg Load 4: %.4f (SD: %.4f)\n', mean(baseline_sb4), std(baseline_sb4));
fprintf(fid, '  Sternberg Load 6: %.4f (SD: %.4f)\n', mean(baseline_sb6), std(baseline_sb6));
fprintf(fid, '  N-back Load 1: %.4f (SD: %.4f)\n', mean(baseline_nb1), std(baseline_nb1));
fprintf(fid, '  N-back Load 2: %.4f (SD: %.4f)\n', mean(baseline_nb2), std(baseline_nb2));
fprintf(fid, '  N-back Load 3: %.4f (SD: %.4f)\n\n', mean(baseline_nb3), std(baseline_nb3));

fprintf(fid, 'Alpha Power Statistics (after outlier exclusion):\n');
fprintf(fid, '  Sternberg Load 2: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', mean(sb2), std(sb2), min(sb2), max(sb2));
fprintf(fid, '  Sternberg Load 4: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', mean(sb4), std(sb4), min(sb4), max(sb4));
fprintf(fid, '  Sternberg Load 6: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', mean(sb6), std(sb6), min(sb6), max(sb6));
fprintf(fid, '  N-back Load 1: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', mean(nb1), std(nb1), min(nb1), max(nb1));
fprintf(fid, '  N-back Load 2: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n', mean(nb2), std(nb2), min(nb2), max(nb2));
fprintf(fid, '  N-back Load 3: M=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n\n', mean(nb3), std(nb3), min(nb3), max(nb3));

fprintf(fid, 'Statistical Tests (paired t-tests, uncorrected):\n');
[~, p_24_sb] = ttest(sb2, sb4);
[~, p_46_sb] = ttest(sb4, sb6);
[~, p_26_sb] = ttest(sb2, sb6);
[~, p_12_nb] = ttest(nb1, nb2);
[~, p_23_nb] = ttest(nb2, nb3);
[~, p_13_nb] = ttest(nb1, nb3);
fprintf(fid, '  Sternberg Load 2 vs 4: p=%.4f\n', p_24_sb);
fprintf(fid, '  Sternberg Load 4 vs 6: p=%.4f\n', p_46_sb);
fprintf(fid, '  Sternberg Load 2 vs 6: p=%.4f\n', p_26_sb);
fprintf(fid, '  N-back Load 1 vs 2: p=%.4f\n', p_12_nb);
fprintf(fid, '  N-back Load 2 vs 3: p=%.4f\n', p_23_nb);
fprintf(fid, '  N-back Load 1 vs 3: p=%.4f\n\n', p_13_nb);

fprintf(fid, 'Analysis Parameters:\n');
fprintf(fid, '  Baseline window: [-0.5 -0.25] s\n');
fprintf(fid, '  Alpha band: 8-14 Hz\n');
fprintf(fid, '  Sternberg time window: [%.1f %.1f] s\n', tsb(1), tsb(2));
fprintf(fid, '  N-back time window: [%.1f %.1f] s\n', tnb(1), tnb(2));
fprintf(fid, '  Cluster permutation alpha: 0.05 (two-tailed)\n');
fprintf(fid, '  Outlier threshold: |z| > 2\n\n');

fprintf(fid, 'Files Saved:\n');
fprintf(fid, '  Controls: %s\n', control_dir);
fprintf(fid, '  Figures: %s\n', figures_dir);
fclose(fid);
fprintf('Summary report saved to: %s\n', summary_file);

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
