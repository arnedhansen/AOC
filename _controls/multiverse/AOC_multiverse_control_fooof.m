%% AOC Multiverse — FOOOF R² Quality Control
% Loads per-trial (and per-subject) FOOOF R² CSVs produced by the
% multiverse prep scripts and generates diagnostic scatterplots.
%
% Figures per task (4 × 2 = 8 figures total):
%   1. R² per subject — jittered scatterplot, faceted by latency, colored by electrode set
%   2. R² vs aperiodic exponent — faceted by latency
%   3. R² vs aperiodic offset — faceted by latency
%   4. R² distribution histogram — with threshold line at 0.90
%
% Input:  /Volumes/.../data/controls/multiverse/fooof_r2_{task}.csv
%         /Volumes/.../data/controls/multiverse/fooof_r2_{task}_subject.csv
% Output: /Volumes/.../figures/controls/multiverse/FOOOF/

%% Paths
base_data = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
if ispc
    base_data = 'W:\Students\Arne\AOC';
end
csv_dir = fullfile(base_data, 'data', 'controls', 'multiverse');
fig_dir = fullfile(base_data, 'figures', 'controls', 'multiverse', 'FOOOF');
if ~isfolder(fig_dir), mkdir(fig_dir); end

r2_thresh = 0.90;

tasks = {'sternberg', 'nback'};
levels = {'trial', 'subject'};

for t = 1:length(tasks)
    task = tasks{t};
    for lv = 1:length(levels)
        level = levels{lv};
        if strcmp(level, 'trial')
            csv_file = fullfile(csv_dir, ['fooof_r2_' task '.csv']);
            suffix = '';
        else
            csv_file = fullfile(csv_dir, ['fooof_r2_' task '_subject.csv']);
            suffix = '_subject';
        end
        if ~isfile(csv_file)
            fprintf('Skipping %s %s-level: %s not found.\n', task, level, csv_file);
            continue
        end

        T = readtable(csv_file);
        fprintf('Loaded %s %s-level: %d FOOOF fits.\n', task, level, height(T));

        %% --- Exclude extreme aperiodic parameter values (IQR-based) ---
        iqr_factor = 3;
        n_before = height(T);

        valid_off = isfinite(T.aperiodic_offset);
        Q1_off = prctile(T.aperiodic_offset(valid_off), 25);
        Q3_off = prctile(T.aperiodic_offset(valid_off), 75);
        IQR_off = Q3_off - Q1_off;
        extreme_off = valid_off & ...
            (T.aperiodic_offset < Q1_off - iqr_factor * IQR_off | ...
             T.aperiodic_offset > Q3_off + iqr_factor * IQR_off);

        valid_exp = isfinite(T.aperiodic_exponent);
        Q1_exp = prctile(T.aperiodic_exponent(valid_exp), 25);
        Q3_exp = prctile(T.aperiodic_exponent(valid_exp), 75);
        IQR_exp = Q3_exp - Q1_exp;
        extreme_exp = valid_exp & ...
            (T.aperiodic_exponent < Q1_exp - iqr_factor * IQR_exp | ...
             T.aperiodic_exponent > Q3_exp + iqr_factor * IQR_exp);

        extreme_mask = extreme_off | extreme_exp;
        n_extreme = sum(extreme_mask);
        fprintf('  Extreme aperiodic values (%.0f× IQR): %d / %d (%.1f%%)\n', ...
            iqr_factor, n_extreme, n_before, 100 * n_extreme / max(n_before, 1));
        fprintf('    Offset  bounds: [%.2f, %.2f]  (excluded %d)\n', ...
            Q1_off - iqr_factor * IQR_off, Q3_off + iqr_factor * IQR_off, sum(extreme_off));
        fprintf('    Exponent bounds: [%.2f, %.2f]  (excluded %d)\n', ...
            Q1_exp - iqr_factor * IQR_exp, Q3_exp + iqr_factor * IQR_exp, sum(extreme_exp));
        T(extreme_mask, :) = [];

        n_below = sum(T.r_squared < r2_thresh & isfinite(T.r_squared));
        n_valid = sum(isfinite(T.r_squared));
        fprintf('  After exclusion: %d fits remaining.\n', height(T));
        fprintf('  R² < %.2f: %d / %d (%.1f%%)\n', r2_thresh, n_below, n_valid, ...
            100 * n_below / max(n_valid, 1));

        lat_labels = unique(T.latency_ms, 'stable');
        elec_labels = unique(T.electrodes, 'stable');
        subj_ids = unique(T.subjectID);
        n_subj = length(subj_ids);
        subj_map = containers.Map(subj_ids, 1:n_subj);

        colors_elec = [0.2 0.4 0.8; 0.8 0.3 0.2];  % posterior = blue, occipital = red

        %% --- Figure 1: R² per subject (faceted by latency, colored by electrode) ---
        fig1 = figure('Position', [0 0 1512 982], 'Color', 'w');
        task_title = task; task_title(1) = upper(task_title(1));
        sgtitle(sprintf('FOOOF R^2 per subject — %s (%s-level)', task_title, level), ...
            'FontSize', 16, 'FontWeight', 'bold');
        for il = 1:length(lat_labels)
            subplot(2, 2, il); hold on;
            mask = strcmp(T.latency_ms, lat_labels{il});
            for ie = 1:length(elec_labels)
                sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                if ~any(sub_mask), continue, end
                subj_x = arrayfun(@(sid) subj_map(sid), T.subjectID(sub_mask));
                jitter = (ie - 1.5) * 0.2 + 0.05 * randn(sum(sub_mask), 1);
                scatter(subj_x + jitter, T.r_squared(sub_mask), 8, ...
                    colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.4, ...
                    'DisplayName', elec_labels{ie});
            end
            yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                'Label', sprintf('R^2 = %.2f', r2_thresh), 'HandleVisibility', 'off');
            xticks(1:n_subj);
            xticklabels(arrayfun(@num2str, subj_ids, 'UniformOutput', false));
            xlabel('Subject'); ylabel('R^2');
            title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
            if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
            ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
            hold off;
        end
        fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_per_subject_%s%s.png', task, suffix));
        exportgraphics(fig1, fname, 'Resolution', 300);
        fprintf('  Saved: %s\n', fname);
        close(fig1);

        %% --- Figure 2: R² vs aperiodic exponent (faceted by latency) ---
        fig2 = figure('Position', [0 0 1512 982], 'Color', 'w');
        sgtitle(sprintf('FOOOF R^2 vs Aperiodic Exponent — %s (%s-level)', task_title, level), ...
            'FontSize', 16, 'FontWeight', 'bold');
        for il = 1:length(lat_labels)
            subplot(2, 2, il); hold on;
            mask = strcmp(T.latency_ms, lat_labels{il});
            for ie = 1:length(elec_labels)
                sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                if ~any(sub_mask), continue, end
                scatter(T.aperiodic_exponent(sub_mask), T.r_squared(sub_mask), 10, ...
                    colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.35, ...
                    'DisplayName', elec_labels{ie});
            end
            yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
            xlabel('Aperiodic Exponent'); ylabel('R^2');
            title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
            if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
            ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
            hold off;
        end
        fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_vs_exponent_%s%s.png', task, suffix));
        exportgraphics(fig2, fname, 'Resolution', 300);
        fprintf('  Saved: %s\n', fname);
        close(fig2);

        %% --- Figure 3: R² vs aperiodic offset (faceted by latency) ---
        fig3 = figure('Position', [0 0 1512 982], 'Color', 'w');
        sgtitle(sprintf('FOOOF R^2 vs Aperiodic Offset — %s (%s-level)', task_title, level), ...
            'FontSize', 16, 'FontWeight', 'bold');
        for il = 1:length(lat_labels)
            subplot(2, 2, il); hold on;
            mask = strcmp(T.latency_ms, lat_labels{il});
            for ie = 1:length(elec_labels)
                sub_mask = mask & strcmp(T.electrodes, elec_labels{ie});
                if ~any(sub_mask), continue, end
                scatter(T.aperiodic_offset(sub_mask), T.r_squared(sub_mask), 10, ...
                    colors_elec(ie, :), 'filled', 'MarkerFaceAlpha', 0.35, ...
                    'DisplayName', elec_labels{ie});
            end
            yline(r2_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
            xlabel('Aperiodic Offset'); ylabel('R^2');
            title(strrep(lat_labels{il}, '_', ' '), 'FontWeight', 'bold');
            if il == 1, legend('Location', 'southwest', 'FontSize', 9); end
            ylim([min(0, min(T.r_squared(mask)) - 0.05) 1.05]);
            hold off;
        end
        fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_vs_offset_%s%s.png', task, suffix));
        exportgraphics(fig3, fname, 'Resolution', 300);
        fprintf('  Saved: %s\n', fname);
        close(fig3);

        %% --- Figure 4: R² distribution histogram ---
        fig4 = figure('Position', [0 0 1512 982], 'Color', 'w');
        valid_r2 = T.r_squared(isfinite(T.r_squared));
        histogram(valid_r2, 50, 'FaceColor', [0.3 0.5 0.7], 'EdgeColor', 'none', ...
            'FaceAlpha', 0.8);
        xline(0.60, '--', 'Color', [0.8 0.5 0], 'LineWidth', 1.5);
        xline(r2_thresh, '--r', 'LineWidth', 1.5);
        xlabel('FOOOF R^2', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold');
        title(sprintf('FOOOF R^2 Distribution — %s (%s-level)  [%d / %d below %.2f (%.1f%%)]', ...
            task_title, level, n_below, n_valid, r2_thresh, 100*n_below/max(n_valid,1)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'FontSize', 12);
        fname = fullfile(fig_dir, sprintf('AOC_fooof_r2_distribution_%s%s.png', task, suffix));
        exportgraphics(fig4, fname, 'Resolution', 300);
        fprintf('  Saved: %s\n', fname);
        close(fig4);
    end
end

fprintf('\n=== FOOOF R² QC DONE ===\n');
