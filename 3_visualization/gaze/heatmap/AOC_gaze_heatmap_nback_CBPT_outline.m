%% AOC Gaze Heatmap - N-Back CBPT Outline (Effect Size, Full Grid)
% Builds raw smoothed-count gaze heatmaps (0 to 2 s) as FieldTrip freq
% structures, runs cluster-based permutation testing (dependent-samples t,
% 3-back vs 1-back) on the full grid without occupancy masking, then plots
% Cohen's d (observed t divided by sqrt(n subjects)) with significant cluster
% outlines (cfg.maskstyle outline, cfg.interactivecolor [0 0 0]).
% Gaze preprocessing matches AOC_gaze_fex_nback.m: in-bounds mask on raw
% coordinates, Y -> screen coordinates via (600 - y), remove_blinks (50 samples).

%% Setup
startup
[subjects, paths, ~, ~] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'gaze', 'heatmap');
if ~isfolder(figDir), mkdir(figDir); end

latencyWindow = [0 2];
num_bins = 200;
smoothing_factor =[2.604 4.623];
x_edges = linspace(0, 800, num_bins);
y_edges = linspace(0, 600, num_bins);
x_centers = x_edges(2:end);
y_centers = y_edges(2:end);

blink_win = 50; % 100 ms at 500 Hz; matches AOC_gaze_fex_nback.m
screenH = 600;

%% Subject-level heatmaps and FieldTrip freq structs
nSub = length(subjects);
dataStim1Allsubs = cell(1, nSub);
dataStim3Allsubs = cell(1, nSub);

for subj = 1:nSub
    datapath = fullfile(path, subjects{subj}, 'gaze');
    load([datapath, filesep 'dataET_nback']);
    clc
    disp(upper(['Loading ET data for subject ' num2str(subj) '/' num2str(nSub) '...']))

    trialinfo = dataETlong.trialinfo;
    if size(trialinfo, 2) > 1
        trialinfo = trialinfo(:, 1);
    end

    condCodes = [21 22 23];
    cfg = [];
    cfg.avgovertime = 'no';
    cfg.keeptrials = 'yes';
    cfg.latency = latencyWindow;

    for condIdx = [1 3]
        cfg.trials = find(trialinfo == condCodes(condIdx));
        dataSel = ft_selectdata(cfg, dataETlong);

        parts_x = cell(1, numel(dataSel.trial));
        parts_y = cell(1, numel(dataSel.trial));
        for tr = 1:numel(dataSel.trial)
            data = double(dataSel.trial{tr});
            if size(data, 1) < 3 || isempty(data)
                continue
            end
            valid_tr = data(1, :) >= 0 & data(1, :) <= 800 & ...
                       data(2, :) >= 0 & data(2, :) <= screenH;
            data = data(1:3, valid_tr);
            if isempty(data)
                continue
            end
            data(2, :) = screenH - data(2, :);
            data = remove_blinks(data, blink_win);
            x_positions = data(1, :);
            y_positions = data(2, :);
            fin = isfinite(x_positions) & isfinite(y_positions);
            x_positions = x_positions(fin);
            y_positions = y_positions(fin);
            ok = ~(x_positions == 0 & y_positions == 0);
            parts_x{tr} = x_positions(ok);
            parts_y{tr} = y_positions(ok);
        end
        x_positions = [parts_x{:}];
        y_positions = [parts_y{:}];

        binned_data = histcounts2(x_positions, y_positions, x_edges, y_edges);
        hm = imgaussfilt(binned_data, smoothing_factor);
        %hm = binned_data; %%%%%

        freq = [];
        freq.label = {'et'};
        freq.dimord = 'chan_freq_time';
        freq.time = x_centers;
        freq.freq = y_centers;
        freq.powspctrm(1, :, :) = hm.';

        switch condIdx
            case 1
                dataStim1Allsubs{subj} = freq;
            case 3
                dataStim3Allsubs{subj} = freq;
        end
    end
end

%% CBPT: 3-back vs 1-back, full grid
cfg = [];
cfg.spmversion = 'spm12';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;
cfg.neighbours = [];
cfg.minnbchan = 0;

design = zeros(2, 2 * nSub);
for i = 1:nSub
    design(1, i) = i;
end
for i = 1:nSub
    design(1, nSub + i) = i;
end
design(2, 1:nSub) = 1;
design(2, nSub + 1:2 * nSub) = 2;

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

clc
disp('COMPUTING CBPT (3-back vs 1-back, full grid)...')
statDIFF_full = ft_freqstatistics(cfg, dataStim3Allsubs{:}, dataStim1Allsubs{:});

%% Cohen's d from observed t
statD = statDIFF_full;
statD.stat = statDIFF_full.stat ./ sqrt(nSub);
if isfield(statDIFF_full, 'mask') && ~isempty(statDIFF_full.mask)
    maskL = logical(statDIFF_full.mask);
    statD.stat(~maskL) = 0;
end

%% Plot
close all
overallFontSize = 35;
centerX = 400;
centerY = 300;
diffColMap = customcolormap_preset('red-white-blue');
plotFtEffectFullWithOutline(statD, diffColMap, overallFontSize, centerX, centerY, ...
    '', ...
    fullfile(figDir, 'AOC_gaze_heatmap_nback_CBPT_outline_stat_full_d.png'));

%%
function lim = prctileFinite(x, p)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    lim = 1;
elseif numel(x) < 5
    lim = max(x);
else
    lim = prctile(x, p);
end
end

function plotFtEffectFullWithOutline(statData, cmap, fontSize, centerX, centerY, panelTitle, outPath)
dLim = prctileFinite(abs(statData.stat(:)), 99.5);
dLim = max(dLim, 0.01);
dLim = min(dLim, 1);
figure('Position', [0 0 1512 982], 'Color', 'w');
cfg = [];
cfg.figure = 'gcf';
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.interactivecolor = [0 0 0];
cfg.zlim = [-dLim dLim];
cfg.colormap = cmap;
ft_singleplotTFR(cfg, statData);
xlim([0 800]);
ylim([0 600]);
yticks([0 150 300 450 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar;
ylabel(cb, 'Effect size (Cohen''s d)', 'FontSize', fontSize);
hold on
plot(centerX, centerY, '+', 'MarkerSize', 15, 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', fontSize);
title(panelTitle, 'FontSize', 30);
saveas(gcf, outPath);
end
