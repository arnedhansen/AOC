%% check baseline per load per task
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
    load2{s}=tfr_fooof;
    load(strcat(subjects{s},'_Sternberg_cond24_fooof.mat'));
    load4{s}=tfr_fooof;
    load(strcat(subjects{s},'_Sternberg_cond26_fooof.mat'));
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
%     cfg.latency = [-1 3];
    sb_high_low{s} = ft_selectdata(cfg,sb_high_low{s});
end
%% handle nback loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
    load(strcat(subjects{s},'_Nback_cond21_fooof.mat'));
    load1{s}=tfr_fooof;
    load(strcat(subjects{s},'_Nback_cond22_fooof.mat'));
    load2nb{s}=tfr_fooof;
    load(strcat(subjects{s},'_Nback_cond23_fooof.mat'));
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
%     cfg.latency = [-1 3];
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
figure; ft_multiplotTFR(cfg, ga_sb);
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
%% extract power spectra SB
cfg = [];
cfg.latency = [-1.5 -.5];
cfg.avgovertime = 'yes';
for s = 1 :length(subjects)
load2_powbl{s}= ft_selectdata(cfg,load2{s});
load4_powbl{s}= ft_selectdata(cfg,load4{s});
load6_powbl{s}= ft_selectdata(cfg,load6{s});
load2_powbl{s}.dimord = 'chan_freq';
load2_powbl{s} = rmfield(load2_powbl{s},'time');
load4_powbl{s}.dimord = 'chan_freq';
load4_powbl{s} = rmfield(load4_powbl{s},'time');
load6_powbl{s}.dimord = 'chan_freq';
load6_powbl{s} = rmfield(load6_powbl{s},'time');
end
% sb task power
cfg = [];
cfg.latency = [1 2];
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

cfg.latency = [-1 0];
for s = 1 :length(subjects)
load1_powbl{s}= ft_selectdata(cfg,load1{s});
load2nb_powbl{s}= ft_selectdata(cfg,load2nb{s});
load3_powbl{s}= ft_selectdata(cfg,load3{s});
load1_powbl{s}.dimord = 'chan_freq';
load1_powbl{s} = rmfield(load1_powbl{s},'time');
load2nb_powbl{s}.dimord = 'chan_freq';
load2nb_powbl{s} = rmfield(load2nb_powbl{s},'time');
load3_powbl{s}.dimord = 'chan_freq';
load3_powbl{s} = rmfield(load3_powbl{s},'time');
end
% do nback task power
cfg.latency = [.5 1.5];
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
ga_sb_2powbl = ft_freqgrandaverage(cfg,load2_powbl{:});
ga_sb_4powbl = ft_freqgrandaverage(cfg,load4_powbl{:});
ga_sb_6powbl = ft_freqgrandaverage(cfg,load6_powbl{:});

ga_nb_2pow = ft_freqgrandaverage(cfg,load2nb_pow{:});
ga_nb_1pow = ft_freqgrandaverage(cfg,load1_pow{:});
ga_nb_3pow = ft_freqgrandaverage(cfg,load3_pow{:});
ga_nb_2powbl = ft_freqgrandaverage(cfg,load2nb_powbl{:});
ga_nb_1powbl = ft_freqgrandaverage(cfg,load1_powbl{:});
ga_nb_3powbl = ft_freqgrandaverage(cfg,load3_powbl{:});

%%
close all
cfg = [];
cfg.figure = 'gcf';
% cfg.ylim = [3 40];
% cfg.zlim = [-3 3];
cfg.layout = layANThead;
figure; ft_multiplotER(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
%% plot with SE sternberg
close all
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
ylabel('Power [a.u.]');
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
ylabel('Power [a.u.]');
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','southeast','Fontsize',20);
%% compare nback power baseline and task
figure;
subplot(2,1,1);
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
ylabel('Power [a.u.]');
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','northeast','Fontsize',20);
title('task .5 to 1.5 sec')

subplot(2,1,2);
cfg = [];
cfg.channel ={'P3', 'P4', 'POz', 'PO3', 'PO4', 'PPO1', 'PPO2', 'PPO5h', 'PPO6h'};% nb
cfg.channel = {'P7', 'PPO9h'};% based on Fstat
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga_nb_1powbl);
tlk4_ind        = ft_selectdata(cfg,ga_nb_2powbl);
tlk6_ind        = ft_selectdata(cfg,ga_nb_3powbl);
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
ylabel('Power [a.u.]');
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 3','Load 2','Load 1'}, 'Location','northeast','Fontsize',20);
title('baseline -1.5 to -.5')
%% compare sb task and baseline power
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

xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',20);
title('task 1 to 2 sec');

subplot(2,1,2);
cfg = [];
cfg.channel ={'Pz', 'POz', 'P2', 'PPO2'};
cfg.channel = {'P5', 'PPO5h'};% based on F stat
cfg.avgoverchan = 'yes';
tlk2_ind        = ft_selectdata(cfg,ga_sb_2powbl);
tlk4_ind        = ft_selectdata(cfg,ga_sb_4powbl);
tlk6_ind        = ft_selectdata(cfg,ga_sb_6powbl);
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

xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(gcf,'color','w');
box on;
grid on;

%     xticks([-1 0 .5 1]);
% xticklabels({'o' '500' '1000'})
xlim([1  40]);
% ylim([-1.5 2.5]);
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',20);
title('baseline -1.5 to -.5 sec');
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
% [statFnb] = ft_freqstatistics(cfg, ga_nb_1pow,ga_nb_2pow,ga_nb_3pow);
[statFsb] = ft_freqstatistics(cfg, ga_sb_2pow,ga_sb_4pow,ga_sb_6pow);
%%
% close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter = 'stat';
cfg.maskparameter ='mask';
% cfg.maskstyle = 'outline';
% cfg.zlim = [-.8 .8];
figure; ft_multiplotER(cfg,statFsb);
set(gcf,'color','w');