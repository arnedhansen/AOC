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
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load2{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond24_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load4{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond26_fooof.mat'));
    cfg = [];
    cfg.baseline     = [-.5 -.25];
    cfg.baselinetype = 'absolute';
    tfr = ft_freqbaseline(cfg,tfr_fooof);
    load6{s} = tfr;
end
%% compute grand average
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/layANThead.mat');
cfg = [];
cfg.keepindividual = 'yes';
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
figure; ft_multiplotTFR(cfg, ga6);
%% select alpha power
cfg = [];
cfg.frequency = [8 14];
cfg.latency = [1 3];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
% cfg.channel = {'Pz', 'P4', 'POz', 'P1', 'P2', 'P6','PO3', 'PO4', 'PPO1', 'PPO2', 'PPO6h', 'POO3h', 'POO4h'};% raw pow
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
val2 = ft_selectdata(cfg,ga2);
val4 = ft_selectdata(cfg,ga4);
val6 = ft_selectdata(cfg,ga6);
alpha2 = val2.powspctrm;
alpha4 = val4.powspctrm;
alpha6 = val6.powspctrm;
%% individual slope across loads
% X = [2 4 6]; % load levels

% for s = 1:length(alpha2)
%     y = [alpha2(s), alpha4(s), alpha6(s)];
%     
%     p = polyfit(X, y, 1); % linear slope
%     slope(s) = p(1);
% end
%% use regression to identify slopes
for s = 1:length(alpha2)
    y = [alpha2(s), alpha4(s), alpha6(s)];
    Xmat = [ones(3,1), [2;4;6]];
    
    b = Xmat \ y';
    slope(s) = b(2);
end
%% plot slopes
figure
histogram(slope, 20)
xlabel('Alpha-load slope')
ylabel('Count')
%% split
thr = 0.015; 

idx_jensen = slope > thr;
idx_nback  = slope < -thr;
idx_flat   = abs(slope) <= thr;
sum(idx_jensen)
sum(idx_nback)
sum(idx_flat)
%% plot included
close all
figure; hold on

% counts
n_j = sum(idx_jensen);
n_n = sum(idx_nback);
n_f = sum(idx_flat);

histogram(slope(idx_jensen), 20, 'FaceColor', [0.8 0 0], 'FaceAlpha', 0.6);
histogram(slope(idx_nback),  20, 'FaceColor', [0 0 0.8], 'FaceAlpha', 0.6);
histogram(slope(idx_flat),   20, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);

xline(thr,  'k--', 'LineWidth', 2);
xline(-thr, 'k--', 'LineWidth', 2);

xlabel('\alpha power slope')
ylabel('# participants')
% title('Distribution of \alpha power modulation across subjects')
title('Subject-wise linear slope of \alpha power across memory load (2, 4, 6 items)')
legend({sprintf('\\alpha increase with load (N=%d)', n_j), ...
        sprintf('\\alpha decrease with load (N=%d)', n_n), ...
        sprintf('intermediate (N=%d)', n_f), ...
        'threshold'})
    box on
%% plot
cfg = [];
cfg.keepindividual = 'yes';
ga2jensen = ft_freqgrandaverage(cfg,load2{idx_jensen});
ga4jensen = ft_freqgrandaverage(cfg,load4{idx_jensen});
ga6jensen = ft_freqgrandaverage(cfg,load6{idx_jensen});
ga2nback = ft_freqgrandaverage(cfg,load2{idx_nback});
ga4nback = ft_freqgrandaverage(cfg,load4{idx_nback});
ga6nback = ft_freqgrandaverage(cfg,load6{idx_nback});
%% plot jensen
close all
cfg = [];
cfg.figure = 'gcf';
cfg.ylim = [3 40];
cfg.zlim = [-.25 .25];
% cfg.zlim = [0 1];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
% cfg.channel = {'Pz', 'P4', 'POz', 'P1', 'P2', 'P6','PO3', 'PO4', 'PPO1', 'PPO2', 'PPO6h', 'POO3h', 'POO4h'};% raw pow
cfg.layout = layANThead;
figure; 
subplot(3,2,1);ft_singleplotTFR(cfg, ga2jensen);
subplot(3,2,3);ft_singleplotTFR(cfg, ga4jensen);
subplot(3,2,5);ft_singleplotTFR(cfg, ga6jensen);

subplot(3,2,2);ft_singleplotTFR(cfg, ga2nback);
subplot(3,2,4);ft_singleplotTFR(cfg, ga4nback);
subplot(3,2,6);ft_singleplotTFR(cfg, ga6nback);
%% select powspctrm
cfg = [];
cfg.latency = [1 3];
cfg.avgovertime = 'yes';
ga2jensen_spctrm = ft_selectdata(cfg,ga2jensen);
ga2jensen_spctrm = rmfield(ga2jensen_spctrm,'time');

ga4jensen_spctrm = ft_selectdata(cfg,ga4jensen);
ga4jensen_spctrm = rmfield(ga4jensen_spctrm,'time');

ga6jensen_spctrm = ft_selectdata(cfg,ga6jensen);
ga6jensen_spctrm = rmfield(ga6jensen_spctrm,'time');

ga2nback_spctrm = ft_selectdata(cfg,ga2nback);
ga2nback_spctrm = rmfield(ga2nback_spctrm,'time');

ga4nback_spctrm = ft_selectdata(cfg,ga4nback);
ga4nback_spctrm = rmfield(ga4nback_spctrm,'time');

ga6nback_spctrm = ft_selectdata(cfg,ga6nback);
ga6nback_spctrm = rmfield(ga6nback_spctrm,'time');
%% plot spectra
close all
cfg = [];
cfg.figure = 'gcf';
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
% cfg.channel = {'Pz', 'P4', 'POz', 'P1', 'P2', 'P6','PO3', 'PO4', 'PPO1', 'PPO2', 'PPO6h', 'POO3h', 'POO4h'};% raw pow
cfg.layout = layANThead;
figure; 
subplot(2,1,1);ft_singleplotER(cfg, ga2jensen_spctrm,ga4jensen_spctrm,ga6jensen_spctrm);
subplot(2,1,2);ft_singleplotER(cfg, ga2nback_spctrm,ga4nback_spctrm,ga6nback_spctrm);
%% do nicer figure
close all
figure;
cfg = [];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
cfg.avgoverchan = 'yes';
cfg.frequency = [1 40];
cfg.latency   = [-1 3];
freq = ft_selectdata(cfg,ga2jensen);
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
subplot(3,2,1);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
% load 4
freq = ft_selectdata(cfg,ga4jensen);
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


subplot(3,2,3);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
% load 6
freq = ft_selectdata(cfg,ga6jensen);
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


subplot(3,2,5);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
%% decrease with load
cfg = [];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
cfg.avgoverchan = 'yes';
cfg.frequency = [1 40];
cfg.latency   = [-1 3];
freq = ft_selectdata(cfg,ga2nback);
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


subplot(3,2,2);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
% load 4
freq = ft_selectdata(cfg,ga4nback);
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


subplot(3,2,4);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
% load 6
freq = ft_selectdata(cfg,ga6nback);
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


subplot(3,2,6);ft_plot_matrix(flip(pow_interp));
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
% title(c,"Power change \newline from baseline")
title(c,'dB')
%% plot with SE sternberg
cfg = [];
cfg.figure = 'gcf';
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
% cfg.channel = {'Pz', 'P4', 'POz', 'P1', 'P2', 'P6','PO3', 'PO4', 'PPO1', 'PPO2', 'PPO6h', 'POO3h', 'POO4h'};% raw pow
cfg.layout = layANThead;
figure; 
subplot(2,1,1);ft_singleplotER(cfg, ga2jensen_spctrm,ga4jensen_spctrm,ga6jensen_spctrm);
subplot(2,1,2);ft_singleplotER(cfg, ga2nback_spctrm,ga4nback_spctrm,ga6nback_spctrm);
%%
% close all
figure;
subplot(2,1,1);
cfg = [];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga2jensen_spctrm);
tlk4_ind        = ft_selectdata(cfg,ga4jensen_spctrm);
tlk6_ind        = ft_selectdata(cfg,ga6jensen_spctrm);
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
ylabel("Change from \newline baseline [dB]");
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',20);
% plot with SE nback
subplot(2,1,2);
cfg = [];
cfg.channel = {'P7', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO3h', 'POO4h'};% posterior chan
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga2nback_spctrm);
tlk4_ind        = ft_selectdata(cfg,ga4nback_spctrm);
tlk6_ind        = ft_selectdata(cfg,ga6nback_spctrm);
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
ylabel("Change from \newline baseline [dB]");
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','northeast','Fontsize',20);
%% now handle gaze
load('sterngaze_norm.mat')
load('sterngaze.mat')

%% baseline gaze
% for subj=1:length(allgazebase6_norm)
%     load2{subj}=allgazebase2_norm{subj};
%     load4{subj}=allgazebase4_norm{subj};
%     load6{subj}=allgazebase6_norm{subj};
%     load2{subj}.powspctrm=allgazebase2_norm{subj}.powspctrm-allgazetasklate2_norm{subj}.powspctrm;
%     load4{subj}.powspctrm=allgazebase4_norm{subj}.powspctrm-allgazetasklate4_norm{subj}.powspctrm;
%     load6{subj}.powspctrm=allgazebase6_norm{subj}.powspctrm-allgazetasklate6_norm{subj}.powspctrm;
% end
for subj=1:length(allgazebase6)
    load2{subj}=allgazebase2{subj};
    load4{subj}=allgazebase4{subj};
    load6{subj}=allgazebase6{subj};
    load2{subj}.powspctrm=allgazebase2{subj}.powspctrm-allgazetasklate2{subj}.powspctrm;
    load4{subj}.powspctrm=allgazebase4{subj}.powspctrm-allgazetasklate4{subj}.powspctrm;
    load6{subj}.powspctrm=allgazebase6{subj}.powspctrm-allgazetasklate6{subj}.powspctrm;
end
%% plot
cfg = [];
cfg.keepindividual = 'yes';
ga2jensen_gaze = ft_freqgrandaverage(cfg,load2{idx_jensen});
ga4jensen_gaze = ft_freqgrandaverage(cfg,load4{idx_jensen});
ga6jensen_gaze = ft_freqgrandaverage(cfg,load6{idx_jensen});
ga2nback_gaze = ft_freqgrandaverage(cfg,load2{idx_nback});
ga4nback_gaze = ft_freqgrandaverage(cfg,load4{idx_nback});
ga6nback_gaze = ft_freqgrandaverage(cfg,load6{idx_nback});
%% plot jensen
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
cfg.zlim = [-.05 .05];
figure; 
subplot(3,2,1);ft_singleplotTFR(cfg, ga2jensen_gaze);
subplot(3,2,3);ft_singleplotTFR(cfg, ga4jensen_gaze);
subplot(3,2,5);ft_singleplotTFR(cfg, ga6jensen_gaze);

subplot(3,2,2);ft_singleplotTFR(cfg, ga2nback_gaze);
subplot(3,2,4);ft_singleplotTFR(cfg, ga4nback_gaze);
subplot(3,2,6);ft_singleplotTFR(cfg, ga6nback_gaze);
%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.latency =[300 500];
% cfg.frequency =[200 400];
cfg.statistic        = 'ft_statfun_diff';
cfg.clusterthreshold ='nonparametric_common';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
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

% [stat_inc_2] = ft_freqstatistics(cfg, allgazetasklate2_norm{idx_jensen},allgazebase2_norm{idx_jensen});
% [stat_inc_4] = ft_freqstatistics(cfg, allgazetasklate4_norm{idx_jensen},allgazebase4_norm{idx_jensen});
% [stat_inc_6] = ft_freqstatistics(cfg, allgazetasklate6_norm{idx_jensen},allgazebase6_norm{idx_jensen});
[stat_inc_2] = ft_freqstatistics(cfg, allgazetasklate2{idx_jensen},allgazebase2{idx_jensen});
[stat_inc_4] = ft_freqstatistics(cfg, allgazetasklate4{idx_jensen},allgazebase4{idx_jensen});
[stat_inc_6] = ft_freqstatistics(cfg, allgazetasklate6{idx_jensen},allgazebase6{idx_jensen});

subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
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

% [stat_inc_n_2] = ft_freqstatistics(cfg, allgazetasklate2_norm{idx_nback},allgazebase2_norm{idx_nback});
% [stat_inc_n_4] = ft_freqstatistics(cfg, allgazetasklate4_norm{idx_nback},allgazebase4_norm{idx_nback});
% [stat_inc_n_6] = ft_freqstatistics(cfg, allgazetasklate6_norm{idx_nback},allgazebase6_norm{idx_nback});
[stat_inc_n_2] = ft_freqstatistics(cfg, allgazetasklate2{idx_nback},allgazebase2{idx_nback});
[stat_inc_n_4] = ft_freqstatistics(cfg, allgazetasklate4{idx_nback},allgazebase4{idx_nback});
[stat_inc_n_6] = ft_freqstatistics(cfg, allgazetasklate6{idx_nback},allgazebase6{idx_nback});


% cohensd=((stat_inc_2.stat)./sqrt(subj));
% stat_inc_2.stat=cohensd;
% stat_inc_2.stat(stat_inc_2.mask==0)=0;% set everything not relevant to zero
% 
% cohensd=((stat_inc_4.stat)./sqrt(subj));
% stat_inc_4.stat=cohensd;
% stat_inc_4.stat(stat_inc_4.mask==0)=0;% set everything not relevant to zero
% 
% cohensd=((stat_inc_6.stat)./sqrt(subj));
% stat_inc_6.stat=cohensd;
% stat_inc_6.stat(stat_inc_6.mask==0)=0;% set everything not relevant to zero

%%
stat_inc_n_2.cfg=[];
stat_inc_n_4.cfg=[];
stat_inc_n_6.cfg=[];
% close all
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';  
% cfg.zlim = [-10e-5 10e-5];
cfg.zlim = [-.001 .001];
% cfg.xlim =[300 500];
% cfg.ylim =[200 400];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,stat_inc_2);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')

subplot(3,2,3);
ft_singleplotTFR(cfg,stat_inc_4);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')

subplot(3,2,5);
ft_singleplotTFR(cfg,stat_inc_6);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')

% plot decreaseing
subplot(3,2,2);
ft_singleplotTFR(cfg,stat_inc_n_2);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')

subplot(3,2,4);
ft_singleplotTFR(cfg,stat_inc_n_4);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')

subplot(3,2,6);
ft_singleplotTFR(cfg,stat_inc_n_6);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
% title(c,'Gaze density [a.u.]');%\it d
c.Label.String = 'Gaze density [a.u.]';
c.Label.FontSize = 18;   % optional
title('')
%%

cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_depsamplesFunivariate'; 
                              
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       
                                  
cfg.clusterstatistic = 'maxsum';  
                                                             
cfg.tail             = 1;          
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        

subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3]; 
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statF_gaze_norm] = ft_freqstatistics(cfg, allgazetasklate2_norm{idx_jensen},allgazetasklate4_norm{idx_jensen},allgazetasklate6_norm{idx_jensen});
[statF_gaze] = ft_freqstatistics(cfg, load2{idx_jensen},load4{idx_jensen},load6{idx_jensen});
statF_gaze.stat(statF_gaze.mask==0)=0;% set everything not relevant to zero

cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_depsamplesFunivariate'; 
                              
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       
                                  
cfg.clusterstatistic = 'maxsum';  
                                                             
cfg.tail             = 1;          
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        

subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3]; 
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statF_gaze_n] = ft_freqstatistics(cfg, load2{idx_nback},load4{idx_nback},load6{idx_nback});
statF_gaze_n.stat(statF_gaze_n.mask==0)=0;% set everything not relevant to zero
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';  
cfg.colormap = 'YlOrRd';
cfg.zlim = [0 10];
% cfg.xlim =[300 500];
% cfg.ylim =[200 400];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,statF_gaze);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [0 10];
title(c,'F-values')
title('')

subplot(3,2,2);
ft_singleplotTFR(cfg,statF_gaze_n);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [0 10];
title(c,'F-values')
title('')

% zoomed version
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';  
cfg.colormap = 'YlOrRd';
cfg.zlim = [0 10];
cfg.xlim =[300 500];
cfg.ylim =[200 400];
cfg.figure = 'gcf';

subplot(3,2,3);
ft_singleplotTFR(cfg,statF_gaze);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [0 10];
title(c,'F-values')
title('')

subplot(3,2,4);
ft_singleplotTFR(cfg,statF_gaze_n);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [0 10];
title(c,'F-values')
title('')
%% test linear trend
addpath('/Volumes/Homestore/OCC/arne/funcs');
cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_loadtrend'; % positive gaze increase with load at that location, negative gaze decrease with load at that location
% cfg.statistic        = 'ft_statfun_loadquadratic'; % positive U shape lowest at the middle load, negative inverted U shape peaks at middle load
                                    
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       
                                  
cfg.clusterstatistic = 'maxsum';  
                                                             
cfg.tail             = 0;          
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        

subj = numel(ga2jensen_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3]; 
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statF_gaze] = ft_freqstatistics(cfg, load2{idx_jensen},load4{idx_jensen},load6{idx_jensen});
statF_gaze.stat(statF_gaze.mask==0)=0;% set everything not relevant to zero


cfg                  = [];
cfg.method           = 'montecarlo'; 
cfg.statistic        = 'ft_statfun_loadtrend'; % positive gaze increase with load at that location, negative gaze decrease with load at that location
% cfg.statistic        = 'ft_statfun_loadquadratic'; % positive U shape lowest at the middle load, negative inverted U shape peaks at middle load
                                    
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       
                                  
cfg.clusterstatistic = 'maxsum';  
                                                             
cfg.tail             = 0;          
cfg.clustertail      = cfg.tail;  
cfg.alpha            = 0.05;      
cfg.numrandomization = 1000;        
subj = numel(ga2nback_gaze.powspctrm(:,1,1,1));
n_2  = subj;
n_4  = subj;
n_6 =  subj;

cfg.design(1,:)           = [ones(1,n_2), ones(1,n_4)*2,ones(1,n_6)*3]; 
cfg.design(2,:)           = [1:n_2,1:n_4, 1:n_6]; 
cfg.ivar                  = 1; 
cfg.uvar                  = 2;
[statF_gaze_n] = ft_freqstatistics(cfg, load2{idx_nback},load4{idx_nback},load6{idx_nback});
statF_gaze_n.stat(statF_gaze_n.mask==0)=0;% set everything not relevant to zero
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';  
% cfg.colormap = 'YlOrRd';
cfg.zlim = [-5 5];
% cfg.xlim =[300 500];
% cfg.ylim =[200 400];
cfg.figure = 'gcf';
figure;
subplot(3,2,1);
ft_singleplotTFR(cfg,statF_gaze);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
title(c,'t-values')
title('')

subplot(3,2,2);
ft_singleplotTFR(cfg,statF_gaze_n);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
title(c,'t-values')
title('')
%% zoom
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';  
% cfg.colormap = 'YlOrRd';
cfg.zlim = [-5 5];
cfg.xlim =[300 500];
cfg.ylim =[200 400];
cfg.figure = 'gcf';

subplot(3,2,3);
ft_singleplotTFR(cfg,statF_gaze);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
title(c,'t-values')
title('')

subplot(3,2,4);
ft_singleplotTFR(cfg,statF_gaze_n);

set(gcf,'color','w');
set(gca,'Fontsize',20);
xlabel('x [px]');
ylabel('y [px]');
% title('high- low during baseline')
grid on
c = colorbar;
c.LineWidth = 1;
c.FontSize = 18;
c.Ticks = [cfg.zlim(1) 0 cfg.zlim(2)];
title(c,'t-values')
title('')
%% 
load('/Volumes/Homestore/OCC/arne/behavioral_matrix_sternberg.mat');
% Load subjects (must match slope order!)
load('/Volumes/Homestore/OCC/arne/subjects.mat');

nSubj = length(subjects);

% Preallocate
RT2 = nan(nSubj,1); RT4 = nan(nSubj,1); RT6 = nan(nSubj,1);
ACC2 = nan(nSubj,1); ACC4 = nan(nSubj,1); ACC6 = nan(nSubj,1);

for i = 1:nSubj
    
    subj_id = str2double(subjects{i}); % assumes subjects like '301'
    
    idx = [behav_data_sternberg.ID] == subj_id;
    
    cond = [behav_data_sternberg(idx).Condition];
    RT   = [behav_data_sternberg(idx).ReactionTime];
    ACC  = [behav_data_sternberg(idx).Accuracy];
    
    RT2(i) = RT(cond==2);
    RT4(i) = RT(cond==4);
    RT6(i) = RT(cond==6);
    
    ACC2(i) = ACC(cond==2);
    ACC4(i) = ACC(cond==4);
    ACC6(i) = ACC(cond==6);
end
%%
% Reaction Time example
sb2 = RT2(idx_jensen);
sb4 = RT4(idx_jensen);
sb6 = RT6(idx_jensen);

nb1 = RT2(idx_nback);
nb2 = RT4(idx_nback);
nb3 = RT6(idx_nback);
%%
figure(30); clf;

positions = [.3, .6, .9];
jitter = 0.05;

% LEFT: alpha increase group
subplot(2,2,1); hold on

[f,x] = ksdensity(sb2); fill(positions(1)+f*0.05,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(sb4); fill(positions(2)+f*0.05,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(sb6); fill(positions(3)+f*0.05,x,'r','FaceAlpha',0.5);

boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2)+(rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3)+(rand(size(sb6))-0.5)*jitter, sb6, 'r.');

ylabel('Reaction Time [s]')
title(sprintf('\\alpha increase with load (N=%d)', sum(idx_jensen)))
set(gca,'FontSize',18); box on;
ylim([0 1.5])

% RIGHT: alpha decrease group
subplot(2,2,2); hold on

[f,x] = ksdensity(nb1); fill(positions(1)+f*0.03,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(nb2); fill(positions(2)+f*0.03,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(nb3); fill(positions(3)+f*0.03,x,'r','FaceAlpha',0.5);

boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(nb1))-0.5)*jitter, nb1, 'k.');
scatter(positions(2)+(rand(size(nb2))-0.5)*jitter, nb2, 'b.');
scatter(positions(3)+(rand(size(nb3))-0.5)*jitter, nb3, 'r.');

ylabel('Reaction Time [s]')
title(sprintf('\\alpha decrease with load (N=%d)', sum(idx_nback)))
set(gca,'FontSize',18); box on;
ylim([0 1.5])
set(gcf,'color','w');

%% Accuracy example
sb2 = ACC2(idx_jensen);
sb4 = ACC4(idx_jensen);
sb6 = ACC6(idx_jensen);

nb1 = ACC2(idx_nback);
nb2 = ACC4(idx_nback);
nb3 = ACC6(idx_nback);

% remove NaNs (important!)
sb2 = sb2(~isnan(sb2)); sb4 = sb4(~isnan(sb4)); sb6 = sb6(~isnan(sb6));
nb1 = nb1(~isnan(nb1)); nb2 = nb2(~isnan(nb2)); nb3 = nb3(~isnan(nb3));

%%
% figure(31); clf;

positions = [.3, .6, .9];
jitter = 0.05;

% ================= LEFT: alpha increase =================
subplot(2,2,3); hold on

[f,x] = ksdensity(sb2); fill(positions(1)+f*3,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(sb4); fill(positions(2)+f*2,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(sb6); fill(positions(3)+f*2,x,'r','FaceAlpha',0.5);

boxplot([sb2(:), sb4(:), sb6(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(sb2))-0.5)*jitter, sb2, 'k.');
scatter(positions(2)+(rand(size(sb4))-0.5)*jitter, sb4, 'b.');
scatter(positions(3)+(rand(size(sb6))-0.5)*jitter, sb6, 'r.');

ylabel('Accuracy [%]')
title(sprintf('\\alpha increase with load (N=%d)', sum(idx_jensen)))
set(gca,'FontSize',18); box on;

ylim([50 110])   % adjust if needed


% ================= RIGHT: alpha decrease =================
subplot(2,2,4); hold on

[f,x] = ksdensity(nb1); fill(positions(1)+f*3,x,'k','FaceAlpha',0.5);
[f,x] = ksdensity(nb2); fill(positions(2)+f*3,x,'b','FaceAlpha',0.5);
[f,x] = ksdensity(nb3); fill(positions(3)+f*3,x,'r','FaceAlpha',0.5);

boxplot([nb1(:), nb2(:), nb3(:)], ...
    'Labels', {'load 2','load 4','load 6'}, ...
    'Positions', positions, 'Widths', 0.05);

scatter(positions(1)+(rand(size(nb1))-0.5)*jitter, nb1, 'k.');
scatter(positions(2)+(rand(size(nb2))-0.5)*jitter, nb2, 'b.');
scatter(positions(3)+(rand(size(nb3))-0.5)*jitter, nb3, 'r.');

ylabel('Accuracy [%]')
title(sprintf('\\alpha decrease with load (N=%d)', sum(idx_nback)))
set(gca,'FontSize',18); box on;

ylim([50 110])   % same scale for comparison

set(gcf,'color','w');
%% test effects 

load('/Volumes/Homestore/OCC/arne/subjects.mat');

nSubj = length(subjects);


RT2 = nan(nSubj,1); RT4 = nan(nSubj,1); RT6 = nan(nSubj,1);
ACC2 = nan(nSubj,1); ACC4 = nan(nSubj,1); ACC6 = nan(nSubj,1);

for i = 1:nSubj
    
    subj_id = str2double(subjects{i});
    
    idx = [behav_data_sternberg.ID] == subj_id;
    
    if sum(idx)==0
        continue
    end
    
    cond = [behav_data_sternberg(idx).Condition];
    RT   = [behav_data_sternberg(idx).ReactionTime];
    ACC  = [behav_data_sternberg(idx).Accuracy];
    
    if sum(cond==2)==0 || sum(cond==4)==0 || sum(cond==6)==0
        continue
    end
    
    RT2(i) = RT(cond==2);
    RT4(i) = RT(cond==4);
    RT6(i) = RT(cond==6);
    
    ACC2(i) = ACC(cond==2);
    ACC4(i) = ACC(cond==4);
    ACC6(i) = ACC(cond==6);
end

%% ================= GROUPING =================
thr = 0.015;

idx_jensen = slope > thr;
idx_nback  = slope < -thr;

%% ================= BUILD TABLE =================
Subject = [];
Load = [];
Group = [];
RT = [];
ACC = [];

for i = 1:nSubj
    
    if idx_jensen(i)
        g = 1;
    elseif idx_nback(i)
        g = -1;
    else
        continue
    end
    
    if any(isnan([RT2(i), RT4(i), RT6(i), ACC2(i), ACC4(i), ACC6(i)]))
        continue
    end
    
    Subject = [Subject; i; i; i];
    Load = [Load; 2; 4; 6];
    Group = [Group; g; g; g];
    
    RT = [RT; RT2(i); RT4(i); RT6(i)];
    ACC = [ACC; ACC2(i); ACC4(i); ACC6(i)];
end

tbl = table(Subject, Load, Group, RT, ACC);

tbl.Subject = categorical(tbl.Subject);
tbl.Group   = categorical(tbl.Group);

%% ================= LME ANALYSIS =================
disp('--- RT model ---')
lme_RT = fitlme(tbl, 'RT ~ Load * Group + (1|Subject)');
disp(anova(lme_RT))

disp('--- ACC model ---')
lme_ACC = fitlme(tbl, 'ACC ~ Load * Group + (1|Subject)');
disp(anova(lme_ACC))


%% ================= TOST ANALYSIS =================
delta_RT = 0.05;   % 50 ms
delta_ACC = 5;     % 5%

RT_jensen = [RT2(idx_jensen); RT4(idx_jensen); RT6(idx_jensen)];
RT_nback  = [RT2(idx_nback);  RT4(idx_nback);  RT6(idx_nback)];

ACC_jensen = [ACC2(idx_jensen); ACC4(idx_jensen); ACC6(idx_jensen)];
ACC_nback  = [ACC2(idx_nback);  ACC4(idx_nback);  ACC6(idx_nback)];

fprintf('\n--- TOST (strict) ---\n')

[p1, p2, eq_RT] = tost_welch(RT_jensen, RT_nback, delta_RT);
fprintf('RT equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_RT);

[p1, p2, eq_ACC] = tost_welch(ACC_jensen, ACC_nback, delta_ACC);
fprintf('ACC equivalence: p1=%.4f, p2=%.4f, equivalent=%d\n', p1, p2, eq_ACC);

%% ================= RT SENSITIVITY =================
fprintf('\n--- RT equivalence sensitivity ---\n')

tost_welch(RT_jensen, RT_nback, 0.05)
tost_welch(RT_jensen, RT_nback, 0.1)
tost_welch(RT_jensen, RT_nback, 0.15)

%% ================= ACC ROBUSTNESS (optional) =================
fprintf('\n--- ACC robustness ---\n')

tost_welch(ACC_jensen, ACC_nback, 3)
tost_welch(ACC_jensen, ACC_nback, 5)
tost_welch(ACC_jensen, ACC_nback, 10)
%% Helper functions
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
function [p1, p2, equivalent] = tost(x1, x2, delta, alpha)

if nargin < 4
    alpha = 0.05;
end

% remove NaNs
x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

n1 = length(x1);
n2 = length(x2);

m1 = mean(x1);
m2 = mean(x2);

s1 = var(x1);
s2 = var(x2);

% pooled SE
SE = sqrt(s1/n1 + s2/n2);

df = n1 + n2 - 2;

% TOST tests
t1 = (m1 - m2 + delta) / SE; % lower bound
t2 = (m1 - m2 - delta) / SE; % upper bound

p1 = 1 - tcdf(t1, df);
p2 = tcdf(t2, df);

equivalent = (p1 < alpha) && (p2 < alpha);
end

%% ================= TOST FUNCTION =================
function [p1, p2, equivalent] = tost_welch(x1, x2, delta, alpha)

if nargin < 4
    alpha = 0.05;
end

x1 = x1(~isnan(x1));
x2 = x2(~isnan(x2));

n1 = length(x1);
n2 = length(x2);

m1 = mean(x1);
m2 = mean(x2);

v1 = var(x1);
v2 = var(x2);

SE = sqrt(v1/n1 + v2/n2);

df = (v1/n1 + v2/n2)^2 / ...
    ((v1^2)/(n1^2*(n1-1)) + (v2^2)/(n2^2*(n2-1)));

t1 = (m1 - m2 + delta) / SE;
t2 = (m1 - m2 - delta) / SE;

p1 = 1 - tcdf(t1, df);
p2 = tcdf(t2, df);

equivalent = (p1 < alpha) && (p2 < alpha);
end