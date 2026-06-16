clear all
close all
load('/Volumes/Homestore/OCC/arne/subjects.mat');
base_dir = '/Volumes/Homestore/OCC/arne/merged';
%% loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
    load(strcat(subjects{s},'_Sternberg_cond22_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load2{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond24_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load4{s}=tfr;
    load(strcat(subjects{s},'_Sternberg_cond26_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load6{s} = tfr;
end
%% handle nback loop over subjects
for s = 1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    cd(subj_dir)
    load(strcat(subjects{s},'_Nback_cond21_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load1{s}=tfr;
    load(strcat(subjects{s},'_Nback_cond22_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load2nb{s}=tfr;
    load(strcat(subjects{s},'_Nback_cond23_tfr.mat'));
    cfg = [];
    cfg.baseline     = [-1.5 -.5];
    cfg.baselinetype = 'db';
    tfr = ft_freqbaseline(cfg,tfr);
    load3{s} = tfr;
end
%% grand average 
load('/Users/tpopov/Documents/DATA4FT/DeepEye/headmodel_ant/layANThead.mat');

cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2 = ft_freqgrandaverage(cfg,load2{:});
ga_sb_4 = ft_freqgrandaverage(cfg,load4{:});
ga_sb_6 = ft_freqgrandaverage(cfg,load6{:});

ga_nb_2 = ft_freqgrandaverage(cfg,load2nb{:});
ga_nb_1 = ft_freqgrandaverage(cfg,load1{:});
ga_nb_3 = ft_freqgrandaverage(cfg,load3{:});
%% select ERD time course
cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_nb_1_erd = ft_selectdata(cfg,ga_nb_1);
ga_nb_2_erd = ft_selectdata(cfg,ga_nb_2);
ga_nb_3_erd = ft_selectdata(cfg,ga_nb_3);

ga_sb_2_erd = ft_selectdata(cfg,ga_sb_2);
ga_sb_4_erd = ft_selectdata(cfg,ga_sb_4);
ga_sb_6_erd = ft_selectdata(cfg,ga_sb_6);
%% plot with SEM
% close all
figure;
subplot(2,2,1);
cfg = [];
cfg.channel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan = 'yes';
cfg.latency = [-.5 2];
tlk2_ind        = ft_selectdata(cfg,ga_nb_1_erd);
tlk4_ind        = ft_selectdata(cfg,ga_nb_2_erd);
tlk6_ind        = ft_selectdata(cfg,ga_nb_3_erd);
% plot load 6
x = tlk6_ind.time'; % x-axis def
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
x = tlk4_ind.time'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% plot load 2
x = tlk2_ind.time'; % x-axis def
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
title('N-back task');
xlabel('Time [sec]');
% ylabel("Power change \newline from baseline");
ylabel("Power change [dB]");
set(gcf,'color','w');
box on;
grid on;
legend({'Load 3','Load 2','Load 1'}, 'Location','northeast','Fontsize',18);
xlim([-.5 2])
ylim([-3 3]);
% plot sternberg
subplot(2,2,2);
cfg = [];
cfg.channel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan = 'yes';
cfg.latency = [-.5 3];
tlk2_ind        = ft_selectdata(cfg,ga_sb_2_erd);
tlk4_ind        = ft_selectdata(cfg,ga_sb_4_erd);
tlk6_ind        = ft_selectdata(cfg,ga_sb_6_erd);
% plot load 6
x = tlk6_ind.time'; % x-axis def
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
x = tlk4_ind.time'; % x-axis def
y = mean(squeeze(tlk4_ind.powspctrm), 1)'; % y-axis def
e = std(squeeze(tlk4_ind.powspctrm), 1)' ./ sqrt(numel(subjects)); % sme
low = y - e; % lower bound
high = y + e; % upper bound

hp2 = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'b', 'HandleVisibility', 'off');
hl2 = line(x, y);
set(hp2, 'facecolor', [0.30, 0.75, 0.93], 'edgecolor', 'none', 'facealpha', 0.2);
set(hl2, 'color', [0.30, 0.75, 0.93], 'linewidth', 2);

% plot load 2
x = tlk2_ind.time'; % x-axis def
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
title('Sternberg task');
xlabel('Time [sec]');
% ylabel("Power change \newline from baseline");
ylabel("Power change [dB]");
set(gcf,'color','w');
box on;
grid on;
legend({'Load 6','Load 4','Load 2'}, 'Location','northeast','Fontsize',18);
xlim([-.5 2])
ylim([-3 3]);
%% plot nback topos
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [0 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.highlightsymbol    = '.';
cfg.highlightsize      = 20;
cfg.comment = 'no';
figure;
subplot(1,3,1);
ft_topoplotER(cfg,ga_nb_1_erd);
subplot(1,3,2);
ft_topoplotER(cfg,ga_nb_2_erd);
subplot(1,3,3);
ft_topoplotER(cfg,ga_nb_3_erd);
%% plot sternback topos
cfg = [];
cfg.layout = layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [1 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.highlightsymbol    = '.';
cfg.highlightsize      = 20;
cfg.comment = 'no';
figure;
subplot(1,3,1);
ft_topoplotER(cfg,ga_sb_2_erd);
subplot(1,3,2);
ft_topoplotER(cfg,ga_sb_4_erd);
subplot(1,3,3);
ft_topoplotER(cfg,ga_sb_6_erd);
%% extract alpha values
cfg =[];
cfg.channel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan = 'yes';
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
cfg.latency = [0 2];
cfg.avgovertime = 'yes';

val_nb_1 = ft_selectdata(cfg,ga_nb_1);
val_nb_2 = ft_selectdata(cfg,ga_nb_2);
val_nb_3 = ft_selectdata(cfg,ga_nb_3);
val_nb_1 = val_nb_1.powspctrm;
val_nb_2 = val_nb_2.powspctrm;
val_nb_3 = val_nb_3.powspctrm;

cfg =[];
cfg.channel = { 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'PPO1', 'PPO2', 'PPO9h', 'PPO5h', 'PPO6h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h'};
cfg.avgoverchan = 'yes';
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
cfg.latency = [1 2];
cfg.avgovertime = 'yes';

val_sb_2 = ft_selectdata(cfg,ga_sb_2);
val_sb_4 = ft_selectdata(cfg,ga_sb_4);
val_sb_6 = ft_selectdata(cfg,ga_sb_6);
val_sb_2 = val_sb_2.powspctrm;
val_sb_4 = val_sb_4.powspctrm;
val_sb_6 = val_sb_6.powspctrm;
%% plot distributions nback
figure; clf;
subplot(1,2,1);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(val_nb_1);
fill(positions(1) + f_sb2*0.4, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(val_nb_2);
fill(positions(2) + f_sb4*0.4, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(val_nb_3);
fill(positions(3) + f_sb6*0.4, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([val_nb_1(:), val_nb_2(:), val_nb_3(:)], ...
    'Labels', {'load 1', 'load 2', 'load 3'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(val_nb_1))-0.5)*jitter, val_nb_1, 'k.');
scatter(positions(2) + (rand(size(val_nb_2))-0.5)*jitter, val_nb_2, 'b.');
scatter(positions(3) + (rand(size(val_nb_3))-0.5)*jitter, val_nb_3, 'r.');


ylabel('\alpha power [dB] ');
xlim([0 1]);
title('');
box on;
set(gcf,'color','w');
set(gca,'Fontsize',20);

% Significance tests
[~, p_24] = ttest(val_nb_1, val_nb_2);
[~, p_46] = ttest(val_nb_2, val_nb_3);
[~, p_26] = ttest(val_nb_1, val_nb_3);

% Significance annotations
y_max = max([val_nb_1(:); val_nb_2(:); val_nb_3(:)]) + 2;  % smaller margin above data overall height above
y_step = 0.5;  % tighter vertical spacing between lines

sig_label = getSigLabel(p_24);
if ~isempty(sig_label)
    line([positions(1), positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(2)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_46);
if ~isempty(sig_label)
    line([positions(2), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(2), positions(3)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_26);
if ~isempty(sig_label)
    line([positions(1), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(3)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
end
% ylim([min([val_nb_1(:); val_nb_2(:); val_nb_3(:)]) - 2, y_max + 2]);
ylim([-8 10])
xlim([.1 1.1])
%% compute linear mixed model nback

Y = [val_nb_1(:); val_nb_2(:); val_nb_3(:)];

Load = [ones(numel(val_nb_1),1); ...
        2*ones(numel(val_nb_2),1); ...
        3*ones(numel(val_nb_3),1)];

Subj = [(1:numel(val_nb_1))'; ...
        (1:numel(val_nb_2))'; ...
        (1:numel(val_nb_3))'];

Tlong = table(Y, categorical(Load), categorical(Subj), ...
    'VariableNames', {'AlphaPower','Load','Subj'});


lme = fitlme(Tlong, ...
    'AlphaPower ~ Load + (1|Subj)');


disp(lme)
%% simple t-test
[H,P,C,stat]=ttest(val_nb_1, val_nb_2)
[H,P,C,stat]=ttest(val_nb_2, val_nb_3)
[H,P,C,stat]=ttest(val_nb_1, val_nb_3)
%% plot distributions sternberg
% figure; clf;
subplot(1,2,2);
positions = [.3, .6, .9];  % Adjusted x-positions

% Kernel Density Plots
[f_sb2, xi_sb2] = ksdensity(val_sb_2);
fill(positions(1) + f_sb2*0.4, xi_sb2, 'k', 'FaceAlpha', 0.5); hold on;

[f_sb4, xi_sb4] = ksdensity(val_sb_4);
fill(positions(2) + f_sb4*0.4, xi_sb4, 'b', 'FaceAlpha', 0.5); hold on;

[f_sb6, xi_sb6] = ksdensity(val_sb_6);
fill(positions(3) + f_sb6*0.4, xi_sb6, 'r', 'FaceAlpha', 0.5); hold on;

% Boxplots
box_h = boxplot([val_sb_2(:), val_sb_4(:), val_sb_6(:)], ...
    'Labels', {'load 2', 'load 4', 'load 6'}, ...
    'Widths', 0.05, ...
    'Positions', positions);
set(box_h, 'LineWidth', 2);

% Raindrop scatter
jitter = 0.05;
scatter(positions(1) + (rand(size(val_nb_1))-0.5)*jitter, val_sb_2, 'k.');
scatter(positions(2) + (rand(size(val_nb_2))-0.5)*jitter, val_sb_4, 'b.');
scatter(positions(3) + (rand(size(val_nb_3))-0.5)*jitter, val_sb_6, 'r.');


ylabel('\alpha power [dB] ');
xlim([0 1]);
title('');
box on;
set(gcf,'color','w');
set(gca,'Fontsize',20);

% Significance tests
[~, p_24] = ttest(val_sb_2, val_sb_4);
[~, p_46] = ttest(val_sb_4, val_sb_6);
[~, p_26] = ttest(val_sb_2, val_sb_6);

% Significance annotations
y_max = max([val_sb_2(:); val_sb_4(:); val_sb_6(:)]) + 2;  % smaller margin above data overall height above
y_step = 0.5;  % tighter vertical spacing between lines

sig_label = getSigLabel(p_24);
if ~isempty(sig_label)
    line([positions(1), positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(2)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_46);
if ~isempty(sig_label)
    line([positions(2), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(2), positions(3)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    y_max = y_max + y_step;
end

sig_label = getSigLabel(p_26);
if ~isempty(sig_label)
    line([positions(1), positions(3)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.2);
    text(mean([positions(1), positions(3)]), y_max + 0.2, sig_label, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
end
% ylim([min([val_nb_1(:); val_nb_2(:); val_nb_3(:)]) - 2, y_max + 2]);
ylim([-8 10])
xlim([.1 1.1])
%% compute linear mixed model nback

Y = [val_sb_2(:); val_sb_4(:); val_sb_6(:)];

Load = [ones(numel(val_sb_2),1); ...
        2*ones(numel(val_sb_4),1); ...
        3*ones(numel(val_sb_6),1)];

Subj = [(1:numel(val_sb_2))'; ...
        (1:numel(val_sb_4))'; ...
        (1:numel(val_sb_6))'];

Tlong = table(Y, categorical(Load), categorical(Subj), ...
    'VariableNames', {'AlphaPower','Load','Subj'});


lme = fitlme(Tlong, ...
    'AlphaPower ~ Load + (1|Subj)');
disp(lme)
%% helper function
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
