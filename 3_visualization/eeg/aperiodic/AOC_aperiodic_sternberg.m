%% AOC Aperiodic Exponent & Offset — Sternberg (Sliding-Window Time Course)
%
% Fits FOOOF at each time step (sliding 500 ms Hanning window, 50 ms steps)
% across occipital channels, per subject and condition. Grand-averages and
% plots aperiodic exponent and offset as a function of time.
%
% Outputs:
%   - Per-condition line plots (exponent, offset) with SEM shading
%   - Difference plot (WM load 6 minus WM load 2)
%
% Requires: FieldTrip, ft_freqanalysis_Arne_FOOOF, startup/setup('AOC')

%% Setup
startup
[subjects, path, ~, ~] = setup('AOC');

fig_dir = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/eeg/aperiodic';
if ~isfolder(fig_dir), mkdir(fig_dir); end

%% Parameters
toi      = -0.5:0.05:2;
win_half = 0.25;
n_toi    = length(toi);
cond_codes = [52 54 56];
cond_vals  = [2 4 6];
n_cond   = length(cond_codes);

cfg_fooof        = [];
cfg_fooof.method = 'mtmfft';
cfg_fooof.taper  = 'hanning';
cfg_fooof.foilim = [2 40];
cfg_fooof.pad    = 5;
cfg_fooof.output = 'fooof';
cfg_fooof.keeptrials = 'no';

%% Resolve occipital channels from first available subject
occ_channels = {};
for s = 1:length(subjects)
    eeg_file = fullfile(path, subjects{s}, 'eeg', 'dataEEG_sternberg.mat');
    if isfile(eeg_file)
        tmp = load(eeg_file);
        for i = 1:length(tmp.dataEEG.label)
            lb = tmp.dataEEG.label{i};
            if contains(lb, 'O') || contains(lb, 'I')
                occ_channels{end+1} = lb;
            end
        end
        occ_channels = unique(occ_channels, 'stable');
        break
    end
end
fprintf('Occipital channels (%d): %s\n', length(occ_channels), strjoin(occ_channels, ', '));

%% Loop over subjects
exp_all = nan(length(subjects), n_cond, n_toi);
off_all = nan(length(subjects), n_cond, n_toi);

for s = 1:length(subjects)
    sid = subjects{s};
    eeg_file = fullfile(path, sid, 'eeg', 'dataEEG_sternberg.mat');
    if ~isfile(eeg_file)
        fprintf('  Skipping %s — no EEG file\n', sid);
        continue
    end
    fprintf('Subject %s (%d/%d)\n', sid, s, length(subjects));
    load(eeg_file, 'dataEEG');
    data = dataEEG;

    for c = 1:n_cond
        trial_idx = find(data.trialinfo == cond_codes(c));
        if isempty(trial_idx), continue; end

        cfg_sel = [];
        cfg_sel.trials = trial_idx;
        cfg_sel.channel = occ_channels;
        d_cond = ft_selectdata(cfg_sel, data);

        cfg_avg = [];
        cfg_avg.keeptrials = 'no';
        d_cond = ft_timelockanalysis(cfg_avg, d_cond);

        cfg_ch = [];
        cfg_ch.avgoverchan = 'yes';
        d_cond = ft_selectdata(cfg_ch, d_cond);
        d_cond.label = {'ROI'};

        for t = 1:n_toi
            t_start = toi(t) - win_half;
            t_end   = toi(t) + win_half;
            if t_start < d_cond.time(1) || t_end > d_cond.time(end)
                continue
            end

            cfg_lat = [];
            cfg_lat.latency = [t_start t_end];
            d_win = ft_selectdata(cfg_lat, d_cond);

            d_raw = [];
            d_raw.label   = d_win.label;
            d_raw.fsample = 1 / mean(diff(d_win.time));
            d_raw.trial   = {d_win.avg};
            d_raw.time    = {d_win.time};

            try
                out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, d_raw);
                rep = out.fooofparams;
                if iscell(rep), rep = rep{1}; end
                if isfield(rep, 'aperiodic_params') && numel(rep.aperiodic_params) >= 2
                    exp_all(s, c, t) = rep.aperiodic_params(2);
                    off_all(s, c, t) = rep.aperiodic_params(1);
                end
            catch
            end
        end
    end
end

%% Grand average
exp_mean = squeeze(nanmean(exp_all, 1));
exp_sem  = squeeze(nanstd(exp_all, 0, 1)) ./ sqrt(sum(~isnan(exp_all(:,1,1))));
off_mean = squeeze(nanmean(off_all, 1));
off_sem  = squeeze(nanstd(off_all, 0, 1)) ./ sqrt(sum(~isnan(off_all(:,1,1))));

%% Colors
colors = [0.2 0.6 1.0;   % WM load 2 — blue
          0.5 0.5 0.5;   % WM load 4 — grey
          0.9 0.2 0.2];  % WM load 6 — red
cond_labels = {'WM load 2', 'WM load 4', 'WM load 6'};
fontSize = 24;

%% Plot: Exponent per condition
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:n_cond
    m = squeeze(exp_mean(c, :));
    se = squeeze(exp_sem(c, :));
    fill([toi fliplr(toi)], [m+se fliplr(m-se)], colors(c,:), ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
for c = 1:n_cond
    plot(toi, squeeze(exp_mean(c,:)), 'Color', colors(c,:), 'LineWidth', 2.5);
end
xline(0, '--k', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', fontSize);
ylabel('Aperiodic Exponent', 'FontSize', fontSize);
title('Sternberg — Aperiodic Exponent', 'FontSize', fontSize);
legend(cond_labels, 'FontSize', 16, 'Location', 'best');
set(gca, 'FontSize', fontSize);
xlim([toi(1) toi(end)]);
hold off
saveas(gcf, fullfile(fig_dir, 'AOC_aperiodic_exponent_sternberg.png'));
fprintf('Saved: AOC_aperiodic_exponent_sternberg.png\n');

%% Plot: Offset per condition
figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
for c = 1:n_cond
    m = squeeze(off_mean(c, :));
    se = squeeze(off_sem(c, :));
    fill([toi fliplr(toi)], [m+se fliplr(m-se)], colors(c,:), ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
for c = 1:n_cond
    plot(toi, squeeze(off_mean(c,:)), 'Color', colors(c,:), 'LineWidth', 2.5);
end
xline(0, '--k', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', fontSize);
ylabel('Aperiodic Offset', 'FontSize', fontSize);
title('Sternberg — Aperiodic Offset', 'FontSize', fontSize);
legend(cond_labels, 'FontSize', 16, 'Location', 'best');
set(gca, 'FontSize', fontSize);
xlim([toi(1) toi(end)]);
hold off
saveas(gcf, fullfile(fig_dir, 'AOC_aperiodic_offset_sternberg.png'));
fprintf('Saved: AOC_aperiodic_offset_sternberg.png\n');

%% Plot: Exponent difference (WM load 6 minus WM load 2)
diff_exp = squeeze(exp_mean(3, :) - exp_mean(1, :));
diff_exp_sem = sqrt(squeeze(exp_sem(3,:)).^2 + squeeze(exp_sem(1,:)).^2);

figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
fill([toi fliplr(toi)], [diff_exp+diff_exp_sem fliplr(diff_exp-diff_exp_sem)], ...
    [0.6 0.3 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(toi, diff_exp, 'Color', [0.6 0.3 0.7], 'LineWidth', 2.5);
yline(0, '--k', 'LineWidth', 1);
xline(0, '--k', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', fontSize);
ylabel('\Delta Exponent (WM6 \minus WM2)', 'FontSize', fontSize);
title('Sternberg — Aperiodic Exponent Difference', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([toi(1) toi(end)]);
hold off
saveas(gcf, fullfile(fig_dir, 'AOC_aperiodic_exponent_diff_sternberg.png'));
fprintf('Saved: AOC_aperiodic_exponent_diff_sternberg.png\n');

%% Plot: Offset difference (WM load 6 minus WM load 2)
diff_off = squeeze(off_mean(3, :) - off_mean(1, :));
diff_off_sem = sqrt(squeeze(off_sem(3,:)).^2 + squeeze(off_sem(1,:)).^2);

figure('Position', [0 0 1512 982], 'Color', 'w');
hold on
fill([toi fliplr(toi)], [diff_off+diff_off_sem fliplr(diff_off-diff_off_sem)], ...
    [0.6 0.3 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(toi, diff_off, 'Color', [0.6 0.3 0.7], 'LineWidth', 2.5);
yline(0, '--k', 'LineWidth', 1);
xline(0, '--k', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', fontSize);
ylabel('\Delta Offset (WM6 \minus WM2)', 'FontSize', fontSize);
title('Sternberg — Aperiodic Offset Difference', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xlim([toi(1) toi(end)]);
hold off
saveas(gcf, fullfile(fig_dir, 'AOC_aperiodic_offset_diff_sternberg.png'));
fprintf('Saved: AOC_aperiodic_offset_diff_sternberg.png\n');

fprintf('=== Sternberg aperiodic time course DONE ===\n');
