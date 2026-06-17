%% AOC EEG ERD/ERS N Back time course plus topos
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
tc1 = [];
tc2 = [];
tc3 = [];
timeVec = [];
tfr1_all = {};
tfr2_all = {};
tfr3_all = {};

for subj = 1:numel(subjects)
    dp = fullfile(path, subjects{subj}, 'eeg');
    f_tc = fullfile(dp, 'ersd_nback_timecourse.mat');
    f_tfr = fullfile(dp, 'tfr_nback.mat');
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
    if ~all(ismember([1; 2; 3], conds))
        warning('Missing n back conditions in ERSD timecourse for subject %s', subjects{subj});
        continue
    end

    if isempty(timeVec)
        timeVec = E.time(:)';
    elseif numel(timeVec) ~= numel(E.time) || any(abs(timeVec - E.time(:)') > 1e-12)
        warning('Time vector mismatch in subject %s', subjects{subj});
        continue
    end

    tcMat = E.ersd_occ_8_14_db;
    tc1 = [tc1; tcMat(conds == 1, :)];
    tc2 = [tc2; tcMat(conds == 2, :)];
    tc3 = [tc3; tcMat(conds == 3, :)];

    T = load(f_tfr, 'tfr1_bl', 'tfr2_bl', 'tfr3_bl');
    tfr1_all{end + 1} = T.tfr1_bl;
    tfr2_all{end + 1} = T.tfr2_bl;
    tfr3_all{end + 1} = T.tfr3_bl;
    fprintf('Loaded n back ERSD timecourse %d / %d\n', subj, numel(subjects));
end

nSubj = size(tc1, 1);
if nSubj == 0
    error('No ERSD timecourse files found. Run AOC_eeg_fex_nback.m first.');
end

ref = tfr1_all{1};
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
mask = timeVec >= -0.5 & timeVec <= 2;
x = timeVec(mask);

y3 = mean(tc3(:, mask), 1, 'omitnan');
e3 = std(tc3(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y2 = mean(tc2(:, mask), 1, 'omitnan');
e2 = std(tc2(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);
y1 = mean(tc1(:, mask), 1, 'omitnan');
e1 = std(tc1(:, mask), 0, 1, 'omitnan') ./ sqrt(nSubj);

h3 = shadedErrorBars(x, y3, e3, [0.97, 0.26, 0.26]); hold on;
h2 = shadedErrorBars(x, y2, e2, [0.30, 0.75, 0.93]);
h1 = shadedErrorBars(x, y1, e1, [0, 0, 0]);

set(gca, 'FontSize', fontSize);
title('N back task');
xlabel('Time [s]');
ylabel('Power change [dB]');
xlim([-0.5 2]);
ylim([-3 3]);
grid on;
box on;
legend([h3.line, h2.line, h1.line], {'3 back', '2 back', '1 back'}, 'Location', 'northeast', 'FontSize', 18);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_timecourse.png'));

%% Topoplots from cached baselined TFR
cfg = [];
cfg.keepindividual = 'yes';
ga_nb_1 = ft_freqgrandaverage(cfg, tfr1_all{:});
ga_nb_2 = ft_freqgrandaverage(cfg, tfr2_all{:});
ga_nb_3 = ft_freqgrandaverage(cfg, tfr3_all{:});

cfg = [];
cfg.frequency = [8 14];
cfg.avgoverfreq = 'yes';
ga_nb_1_erd = ft_selectdata(cfg, ga_nb_1);
ga_nb_2_erd = ft_selectdata(cfg, ga_nb_2);
ga_nb_3_erd = ft_selectdata(cfg, ga_nb_3);

cfg = [];
cfg.layout = headmodel.layANThead;
cfg.figure = 'gcf';
cfg.zlim = [-2 2];
cfg.xlim = [0 2];
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = chPlot;
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.comment = 'no';

figure('Position', [0 0 1512 982], 'Color', 'w');
subplot(1, 3, 1); ft_topoplotER(cfg, ga_nb_1_erd);
subplot(1, 3, 2); ft_topoplotER(cfg, ga_nb_2_erd);
subplot(1, 3, 3); ft_topoplotER(cfg, ga_nb_3_erd);
saveas(gcf, fullfile(figDir, 'AOC_eeg_ersd_nback_topos.png'));

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
