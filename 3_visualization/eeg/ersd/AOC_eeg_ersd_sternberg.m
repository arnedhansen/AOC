%% AOC EEG ERD/ERS Sternberg time course plus topos
% Time course is loaded from saved ERSD outputs produced during feature extraction.
% No ERSD recomputation is performed here.

%% Setup
startup
[subjects, paths, ~, headmodel] = setup('AOC');
path = paths.features;
figDir = fullfile(paths.figures, 'eeg', 'ersd');
if ~isfolder(figDir), mkdir(figDir); end
subjects = setdiff(subjects, {'361'});

%% Load subject ERSD timecourses and cached TFR for topoplots
tc2 = [];
tc4 = [];
tc6 = [];
timeVec = [];
tfr2_all = {};
tfr4_all = {};
tfr6_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tc = fullfile(dp, 'ersd_sternberg_timecourse.mat');
    f_tfr = fullfile(dp, 'tfr_stern.mat');
    if ~isfile(f_tc) || ~isfile(f_tfr)
        warning('Missing ERSD or TFR file for subject %s', subjects{subj});
        continue
    end

    S = load(f_tc, 'ersd_timecourse');
    if ~isfield(S, 'ersd_timecourse')
        warning('Invalid ERSD timecourse file for subject %s', subjects{subj});
        continue
    end
    E = S.ersd_timecourse;
    conds = E.condition(:);
    if ~all(ismember([2; 4; 6], conds))
        warning('Missing Sternberg conditions in ERSD timecourse for subject %s', subjects{subj});
        continue
    end

    if isempty(timeVec)
        timeVec = E.time(:)';
    elseif numel(timeVec) ~= numel(E.time) || any(abs(timeVec - E.time(:)') > 1e-12)
        warning('Time vector mismatch in subject %s', subjects{subj});
        continue
    end

    tcMat = E.ersd_occ_8_14_db;
    tc2 = [tc2; tcMat(conds == 2, :)];
    tc4 = [tc4; tcMat(conds == 4, :)];
    tc6 = [tc6; tcMat(conds == 6, :)];

    T = load(f_tfr, 'tfr2_bl', 'tfr4_bl', 'tfr6_bl');
    tfr2_all{end + 1} = T.tfr2_bl;
    tfr4_all{end + 1} = T.tfr4_bl;
    tfr6_all{end + 1} = T.tfr6_bl;
    fprintf('Loaded Sternberg ERSD timecourse %d / %d\n', subj, numel(subjects));
end

nSubj = size(tc2, 1);
if nSubj == 0
    error('No ERSD timecourse files found. Run AOC_eeg_fex_sternberg.m first.');
end

ref = tfr2_all{1};
chPlot = {};
for i = 1:numel(ref.label)
    label = ref.label{i};
    if contains(label, {'O'}) || contains(label, {'I'})
        chPlot{end + 1} = label;
    end
end
if isempty(chPlot)
    error('No occipital O or I channels in TFR labels.');
end

%% Timecourse figure from cached ERSD timecourses
figure('Position', [0 0 1512 982], 'Color', 'w');
fontSize = 22;
mask = timeVec >= -0.5 & timeVec <= 3;
x = timeVec(mask);

y6 = mean(tc6(:, mask), 1, 'omitnan');
e6 = std(tc6(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y4 = mean(tc4(:, mask), 1, 'omitnan');
e4 = std(tc4(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y2 = mean(tc2(:, mask), 1, 'omitnan');
e2 = std(tc2(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);

h6 = shadedErrorBars(x, y6, e6, [0.97, 0.26, 0.26]); hold on;
h4 = shadedErrorBars(x, y4, e4, [0.30, 0.75, 0.93]);
h2 = shadedErrorBars(x, y2, e2, [0, 0, 0]);

set(gca, 'FontSize', fontSize);
title('Sternberg task');
xlabel('Time [s]');
ylabel('Power change [dB]');
xlim([-0.5 2]);
ylim([-3 3]);
grid on;
box on;
legend([h6.line, h4.line, h2.line], {'WM 6', 'WM 4', 'WM 2'}, 'Location', 'northeast', 'FontSize', 18);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_timecourse.png'));

%% Topoplots from cached baselined TFR
cfg = [];
cfg.keepindividual = 'yes';
ga_sb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_sb_4 = ft_freqgrandaverage(cfg, tfr4_all{:});
ga_sb_6 = ft_freqgrandaverage(cfg, tfr6_all{:});

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_sb_2_erd = ft_selectdata(cfg, ga_sb_2);
ga_sb_4_erd = ft_selectdata(cfg, ga_sb_4);
ga_sb_6_erd = ft_selectdata(cfg, ga_sb_6);

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [1 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';

figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1, 3, 1); ft_topoplotER(cfg, ga_sb_2_erd);
subplot(1, 3, 2); ft_topoplotER(cfg, ga_sb_4_erd);
subplot(1, 3, 3); ft_topoplotER(cfg, ga_sb_6_erd);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_sternberg_topos.png'));

function h = shadedErrorBars(x, y, e, col)
x = x(:);
y = y(:);
e = e(:);
low = y - e;
high = y + e;
h.patch = patch([x; flipud(x)], [low; flipud(high)], col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on
h.line = plot(x, y, 'Color', col, 'LineWidth', 2);
end
