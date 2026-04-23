%% AOC Topos — Raw EEG (N-back + Sternberg)
% Creates one figure with:
%   - rows   = tasks (N-back, Sternberg)
%   - cols   = conditions (3 per task)
% Uses raw EEG power sources and the same colormap style as the raw EEG topo scripts.

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('AOC');
pathAOC = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'topos');
if ~isfolder(figDir), mkdir(figDir); end

fontSize = 22;

%% ---------- N-back: load raw power windows ----------
subjects_nback = setdiff(subjects, {'361'}); % keep exclusion used in existing raw n-back script
powl1 = {};
powl2 = {};
powl3 = {};

for subj = 1:numel(subjects_nback)
    datapath = fullfile(pathAOC, subjects_nback{subj}, 'eeg');
    if ~isfolder(datapath), continue, end
    cd(datapath);
    if ~isfile('power_nback_windows.mat'), continue, end
    S = load('power_nback_windows.mat', 'pow1_full', 'pow2_full', 'pow3_full');
    if isfield(S, 'pow1_full'), powl1{end+1} = S.pow1_full; end
    if isfield(S, 'pow2_full'), powl2{end+1} = S.pow2_full; end
    if isfield(S, 'pow3_full'), powl3{end+1} = S.pow3_full; end
end

assert(~isempty(powl1) && ~isempty(powl2) && ~isempty(powl3), ...
    'Missing N-back power windows for at least one condition.');

ga_nback = cell(1, 3);
ga_nback{1} = ft_freqgrandaverage([], powl1{:});
ga_nback{2} = ft_freqgrandaverage([], powl2{:});
ga_nback{3} = ft_freqgrandaverage([], powl3{:});
nback_labels = {'1-back', '2-back', '3-back'};

%% ---------- Sternberg: load raw power windows ----------
pows2 = {};
pows4 = {};
pows6 = {};

for subj = 1:numel(subjects)
    datapath = fullfile(pathAOC, subjects{subj}, 'eeg');
    if ~isfolder(datapath), continue, end
    cd(datapath);
    if ~isfile('power_stern_windows.mat'), continue, end
    S = load('power_stern_windows.mat', 'pow2_late', 'pow4_late', 'pow6_late');
    if isfield(S, 'pow2_late'), pows2{end+1} = S.pow2_late; end
    if isfield(S, 'pow4_late'), pows4{end+1} = S.pow4_late; end
    if isfield(S, 'pow6_late'), pows6{end+1} = S.pow6_late; end
end

assert(~isempty(pows2) && ~isempty(pows4) && ~isempty(pows6), ...
    'Missing Sternberg power windows for at least one condition.');

ga_stern = cell(1, 3);
ga_stern{1} = ft_freqgrandaverage([], pows2{:});
ga_stern{2} = ft_freqgrandaverage([], pows4{:});
ga_stern{3} = ft_freqgrandaverage([], pows6{:});
stern_labels = {'WM load 2', 'WM load 4', 'WM load 6'};

%% ---------- Channels to highlight (occipital + infero-occipital style) ----------
channels = occ_channels_from_labels(ga_nback{1}.label);

%% ---------- Shared plotting config ----------
cfg_base = [];
cfg_base.layout = headmodel.layANThead;
allchannels = cfg_base.layout.label;
cfg_base.channel = allchannels(1:end-2);
cfg_base.channel = cfg_base.channel(~strcmp(cfg_base.channel, 'M2'));
cfg_base.channel = cfg_base.channel(~strcmp(cfg_base.channel, 'M1'));
cfg_base.marker = 'off';
cfg_base.highlight = 'on';
cfg_base.highlightchannel = channels;
cfg_base.highlightsymbol = '.';
cfg_base.highlightsize = 10;
cfg_base.gridscale = 300;
cfg_base.comment = 'no';
cfg_base.xlim = [8 14];

% Match raw EEG topo scripts' colormap
cfg_base.colormap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);

%% ---------- Compute per-task shared zlim across conditions ----------
zlim_nback = compute_task_zlim(ga_nback);
zlim_stern = compute_task_zlim(ga_stern);

%% ---------- Plot ----------
figure('Position', [0 0 1512 982], 'Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Row 1: N-back
for c = 1:3
    ax = nexttile(c);
    cfg = cfg_base;
    cfg.figure = ax;
    cfg.zlim = zlim_nback;
    ft_topoplotER(cfg, ga_nback{c});
    cb = colorbar(ax);
    cb.Label.String = 'Power [\muV^2/Hz]';
    set(ax, 'FontSize', fontSize);
    title(ax, nback_labels{c}, 'Interpreter', 'none');
end

% Row 2: Sternberg
for c = 1:3
    ax = nexttile(3 + c);
    cfg = cfg_base;
    cfg.figure = ax;
    cfg.zlim = zlim_stern;
    ft_topoplotER(cfg, ga_stern{c});
    cb = colorbar(ax);
    cb.Label.String = 'Power [\muV^2/Hz]';
    set(ax, 'FontSize', fontSize);
    title(ax, stern_labels{c}, 'Interpreter', 'none');
end

saveas(gcf, fullfile(figDir, 'AOC_eeg_topos_raw_all_tasks.png'));

%% ---------- Local functions ----------
function zlim = compute_task_zlim(ga_cells)
% Shared zlim for all conditions in one task (alpha range 8-14 Hz).
all_alpha = [];
for i = 1:numel(ga_cells)
    g = ga_cells{i};
    fidx = g.freq >= 8 & g.freq <= 14;
    if ~any(fidx), continue, end
    Ai = mean(g.powspctrm(:, fidx), 2, 'omitnan');
    all_alpha = [all_alpha; Ai(:)]; %#ok<AGROW>
end
all_alpha = all_alpha(isfinite(all_alpha));
if isempty(all_alpha)
    zlim = [0 1];
    return
end
mx = prctile(all_alpha, 99);
if ~isfinite(mx) || mx <= 0
    mx = max(all_alpha);
end
if ~isfinite(mx) || mx <= 0
    mx = 1;
end
zlim = [0 mx];
end

function ch = occ_channels_from_labels(labels)
ch = {};
for i = 1:numel(labels)
    lab = labels{i};
    if contains(lab, {'O'}) || contains(lab, {'I'})
        ch{end+1} = lab; %#ok<AGROW>
    end
end
if isempty(ch)
    ch = labels;
end
end
