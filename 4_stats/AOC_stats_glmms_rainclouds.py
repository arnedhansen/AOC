# %% AOC Stats Rainclouds — Combined (Sternberg + N-Back)
# Loads merged CSV for both tasks, runs ANOVAs and mixed models, generates raincloud plots. Saves figures and model/ANOVA tables.
#
# Key outputs:
#   Raincloud figures; ANOVA/mixed-model tables

# %% Imports
import sys, os
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from statsmodels.stats.anova import AnovaRM
import statsmodels.formula.api as smf

from functions.stats_helpers import (
    iqr_outlier_filter, mixedlm_pairwise_contrasts, p_to_signif
)

from functions.rainclouds_plotting_helpers import add_stat_brackets

from functions.export_model_table import export_model_table

from functions.mixedlm_helpers import (
    fit_mixedlm, drop1_lrt, pairwise_condition_contrasts_at_mean_gaze, mixedlm_fixed_effects_to_df, lr_effect_sizes
)

# %% Parameters

# Colours
pal = ["#93B8C4", "#82AD82", "#D998A2"]  # AOC pastels

# Plot appearance
mpl.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 300,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "ps.fonttype": 42,
    "font.size": 15,
    "axes.titlesize": 15,
    "axes.labelsize": 15,
    "legend.fontsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "mathtext.default": "regular",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.edgecolor": "white",
    "axes.edgecolor": "white",
})

sns.set_style("white")  # clean white background, no grey panel

# %% I/O Paths
base_dir   = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
output_dir = f"{base_dir}/figures/stats/rainclouds"
output_dir_stats = f"{base_dir}/data/stats"
anova_dir = f"{base_dir}/data/stats/anova"

# %% Variables and labelling

variables  = ["Accuracy", "ReactionTime", "GazeDeviation", "MSRate", "Fixations", "Saccades", "PupilSize", "ScanPathLength", "AlphaPower", "IAF"]
titles     = ["Accuracy", "Reaction Time", "Gaze Deviation", "Microsaccade Rate", "Fixations", "Saccades", "Pupil Size", "Scan Path Length", "Alpha Power", "IAF"]
y_labels   = ["Accuracy [%]", "Reaction Time [s]", "Gaze Deviation [px]", "Microsaccade Rate [MS/s]", "Fixations", "Saccades", "Pupil Size [a.u.]", "Scan Path Length [px]", "Alpha Power [\u03BCV²/Hz]", "IAF [Hz]"]
save_names = ["acc", "rt", "gazedev", "ms", "fix", "sacc", "pupil", "spl", "pow", "iaf"]

# Manual y ticks and ylims per variable
yticks_map = {
    "Accuracy"      : np.arange(60, 101, 5),
    "ReactionTime"  : np.arange(0.3, 1.35, 0.1),
    "GazeDeviation" : np.arange(0, 65, 10),
    "MSRate"        : np.arange(0, 3.8, 0.5),
    "Fixations"     : np.arange(0, 6, 1),
    "Saccades"      : np.arange(0, 4.25, 1),
    "PupilSize"     : np.arange(0, 5.5, 1),
    "ScanPathLength": np.arange(0, 410, 50),
    "AlphaPower"    : np.arange(0, 1.52, 0.25),
    "IAF"           : np.arange(8, 14, 1),
}
ylims_map = {
    "Accuracy"      : (60, 102),
    "ReactionTime"  : (0.3, 1.35),
    "GazeDeviation" : (0, 65),
    "MSRate"        : (0, 3.8),
    "Fixations"     : (0, 6),
    "Saccades"      : (0, 4.25),
    "PupilSize"     : (0, 5.5),
    "ScanPathLength": (0, 410),
    "AlphaPower"    : (0, 1.6),
    "IAF"           : (8, 14.1),
}

# %% Task configurations

tasks = [
    {
        "name"       : "sternberg",
        "input_csv"  : f"{base_dir}/data/features/merged_data_sternberg.csv",
        # Accept numeric encodings {1,2,3} or {2,4,6}; otherwise normalise existing strings
        "cond_to_label_numeric": [{1: "WM load 2", 2: "WM load 4", 3: "WM load 6"},
                                  {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"}],
        "categories" : ["WM load 2", "WM load 4", "WM load 6"],
        "comparisons": [("WM load 2", "WM load 4"), ("WM load 2", "WM load 6"), ("WM load 4", "WM load 6")],
        "xlabel"     : "Condition"
    },
    {
        "name"       : "nback",
        "input_csv"  : f"{base_dir}/data/features/merged_data_nback.csv",
        "cond_to_label_numeric": [{1: "1-back", 2: "2-back", 3: "3-back"}],
        "categories" : ["1-back", "2-back", "3-back"],
        "comparisons": [("1-back", "2-back"), ("1-back", "3-back"), ("2-back", "3-back")],
        "xlabel"     : "Condition"
    }
]

# %% Pre-scan for global upper bounds (shared bracket baseline per variable)
global_upper = {var: np.nan for var in variables}

for _task in tasks:
    _dat = pd.read_csv(_task["input_csv"])
    _dat.loc[_dat["Accuracy"] > 100, "Accuracy"] = np.nan  # same impossible-value rule

    # normalise/label Condition exactly like in the main loop
    _cond = _dat["Condition"]
    if np.issubdtype(_cond.dtype, np.number):
        _uniq = sorted(pd.unique(_cond.dropna()).tolist())
        _applied_map = None
        for _cand in _task["cond_to_label_numeric"]:
            if set(_uniq).issubset(set(_cand.keys())):
                _applied_map = _cand
                break
        if _applied_map is None:
            _applied_map = {val: _task["categories"][i] for i, val in enumerate(_uniq[:len(_task["categories"])])}
        _dat["Condition"] = _dat["Condition"].map(_applied_map)
    else:
        _dat["Condition"] = _dat["Condition"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()

    _dat["Condition"] = pd.Categorical(_dat["Condition"], categories=_task["categories"], ordered=True)
    if _dat["ID"].dtype != "O":
        _dat["ID"] = _dat["ID"].astype(str)

    # apply the same outlier filter as the main pipeline
    _dat = iqr_outlier_filter(_dat, variables, by="Condition")

    # update global upper bounds
    for var in variables:
        vmax = pd.to_numeric(_dat[var], errors="coerce").max()
        if np.isfinite(vmax):
            if not np.isfinite(global_upper[var]) or (vmax > global_upper[var]):
                global_upper[var] = float(vmax)

# %% Processing + plotting
outputs = {}

for task in tasks:

    # %% Load
    dat = pd.read_csv(task["input_csv"])

    # Remove impossible values
    dat.loc[dat["Accuracy"] > 100, "Accuracy"] = np.nan

    # %% Condition labelling (robust to numeric or already-labelled strings)
    cond = dat["Condition"]
    if np.issubdtype(cond.dtype, np.number):
        uniq = sorted(pd.unique(cond.dropna()).tolist())
        applied_map = None
        for candidate_map in task["cond_to_label_numeric"]:
            if set(uniq).issubset(set(candidate_map.keys())):
                applied_map = candidate_map
                break
        if applied_map is None:
            # Fallback by rank order
            applied_map = {val: task["categories"][i] for i, val in enumerate(uniq[:len(task["categories"])])}
        dat["Condition"] = dat["Condition"].map(applied_map)
    else:
        dat["Condition"] = dat["Condition"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()

    dat["Condition"] = pd.Categorical(dat["Condition"], categories=task["categories"], ordered=True)

    # Ensure ID is string for grouping
    if dat["ID"].dtype != "O":
        dat["ID"] = dat["ID"].astype(str)

    # Outlier removal per condition & variable (1.5×IQR)
    dat = iqr_outlier_filter(dat, variables, by="Condition")

    # %% Descriptive statistics per condition (save as CSV)
    desc_rows = []

    for var in variables:
        if var not in dat.columns:
            continue

        subdat = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
        if subdat.empty:
            continue

        # Group by Condition (subjects already collapsed to one row per condition)
        for cond_lab, g in subdat.groupby("Condition", observed=True):
            if g.empty:
                continue

            vals = g[var].to_numpy()

            n     = vals.size
            mean  = np.nanmean(vals)
            sd    = np.nanstd(vals, ddof=1)
            sem   = sd / np.sqrt(n) if n > 0 else np.nan
            med   = np.nanmedian(vals)
            q1    = np.nanpercentile(vals, 25)
            q3    = np.nanpercentile(vals, 75)
            iqr   = q3 - q1

            desc_rows.append([
                task["name"], var, cond_lab, n,
                round(mean, 2), round(sd, 2), round(sem, 2),
                round(med, 2), round(q1, 2), round(q3, 2), round(iqr, 2)
            ])

    desc_table = pd.DataFrame(
        desc_rows,
        columns=[
            "Task", "Variable", "Condition", "N",
            "Mean", "SD", "SEM",
            "Median", "Q1", "Q3", "IQR"
        ]
    )

    # Save descriptive summary for the task
    out_csv = os.path.join(output_dir_stats, f"AOC_descriptives_{task['name']}.csv")
    desc_table.to_csv(out_csv, index=False)
    outputs[f"descriptives_{task['name']}"] = desc_table.copy()
    print(f"Saved descriptives → {out_csv}")

    # Category order / palette mapping
    condition_order = list(dat["Condition"].dropna().unique())
    pal_dict = dict(zip(condition_order, pal))

    anova_rows = []               # collects ANOVA rows for this task
    pairwise_effsize_rows = []    # collects pairwise effect sizes for this task

    # %% Loop variables
    for var, ttl, ylab, sname in zip(variables, titles, y_labels, save_names):

        dvar = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            continue

        # Ensure categorical ordering
        dvar["Condition"] = pd.Categorical(dvar["Condition"], categories=condition_order, ordered=True)

        # %% Repeated-measures ANOVA (within-subject: Condition)
        # balance subjects: keep only subjects present in all conditions for this var
        present = dvar.groupby('ID')['Condition'].nunique()
        keep_ids = present[present == len(condition_order)].index
        dvar_bal = dvar[dvar['ID'].isin(keep_ids)].copy()

        # run ANOVA on the balanced panel
        if dvar_bal['ID'].nunique() >= 2:
            aov = AnovaRM(data=dvar_bal, depvar=var, subject='ID', within=['Condition']).fit()
            aov_tab = aov.anova_table.reset_index().rename(columns={'index': 'Effect'})
            for _, r in aov_tab.iterrows():
                F   = float(r['F Value'])
                df1 = float(r['Num DF'])
                df2 = float(r['Den DF'])
                p   = float(r['Pr > F'])
                etap = (F * df1) / (F * df1 + df2) if np.isfinite(F) else np.nan
                anova_rows.append([task['name'], var, r['Effect'], df1, df2, F, p, etap])
        else:
            anova_rows.append([task['name'], var, 'Condition', np.nan, np.nan, np.nan, np.nan, np.nan])

        # Data bounds (used to size violins and brackets)
        lower_bound = float(dvar[var].min())
        upper_bound = float(dvar[var].max())

        # Mixed model + Bonferroni pairwise
        pw = mixedlm_pairwise_contrasts(
            dvar.rename(columns={var: "value"}),
            value_col="value",
            group_col="Condition",
            id_col="ID",
            p_adjust="bonferroni"
        )

        # %% Pairwise within-subject effect sizes (Cohen's dz) + 95% CI of mean difference
        # build wide table to get paired diffs
        wide = dvar.pivot(index='ID', columns='Condition', values=var)
        for (g1, g2) in task["comparisons"]:
            if (g1 in wide.columns) and (g2 in wide.columns):
                diffs = (wide[g2] - wide[g1]).dropna()
                n = diffs.shape[0]
                if n >= 2 and np.nanstd(diffs, ddof=1) > 0:
                    md   = float(np.nanmean(diffs))
                    sd_d = float(np.nanstd(diffs, ddof=1))
                    dz   = md / sd_d
                    # 95% CI for mean difference (paired t): md ± t*sd_d/sqrt(n)
                    from scipy import stats
                    tcrit = stats.t.ppf(0.975, df=n-1)
                    se_md = sd_d / np.sqrt(n)
                    ci_lo = md - tcrit * se_md
                    ci_hi = md + tcrit * se_md
                else:
                    md = dz = ci_lo = ci_hi = np.nan
                    n = int(n)
                # attach adjusted p if available
                padj = np.nan
                row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
                if not row.empty and "p_adj" in row:
                    try:
                        padj = float(row["p_adj"].iloc[0])
                    except Exception:
                        padj = np.nan
                pairwise_effsize_rows.append([
                    task["name"], var, g1, g2, n, md, dz, ci_lo, ci_hi, padj
                ])

        # %% Fit mixed model to export as a Word table
        # Random-intercept model (subjects), Condition as fixed effect.
        dvar_m = dvar.rename(columns={var: "value"}).copy()
        # Ensure Condition is categorical with your order (already set above)
        dvar_m["Condition"] = pd.Categorical(dvar_m["Condition"],
                                            categories=condition_order, ordered=True)

        # Use treatment coding with the FIRST level as reference (your ordered categories)
        # reml=False -> so the reported log-likelihood is comparable across models.
        # If it struggles to converge, it will fall back to random-intercepts only and try a different optimizer.
        model_result = None
        try:
            # Random intercepts only
            m = smf.mixedlm("value ~ C(Condition, Treatment(reference=condition_order[0]))",
                            data=dvar_m, groups=dvar_m["ID"])
            model_result = m.fit(method="lbfgs", reml=False, maxiter=200, disp=False)
        except Exception as e1:
            try:
                # Try Nelder-Mead as a fallback
                model_result = m.fit(method="nm", reml=False, maxiter=400, disp=False)
            except Exception as e2:
                # Final fallback: simple OLS (no random effects), so you still get a table
                m_ols = smf.ols("value ~ C(Condition, Treatment(reference=condition_order[0]))",
                                data=dvar_m).fit()
                model_result = m_ols

        # Export to Word
        doc_name = f"AOC_modeltable_{sname}_{task['name']}.docx"
        doc_path = os.path.join(output_dir_stats, doc_name)
        export_model_table(model_result, doc_path)
        print(f"Saved model table    → {doc_name}")
        
        # Also save fixed effects (β, SE, z/t, p, CI) to CSV and outputs
        from functions.mixedlm_helpers import mixedlm_fixed_effects_to_df

        fe_df = mixedlm_fixed_effects_to_df(
            model_result,
            task=task["name"],
            variable=var,
            model_label="DV ~ Condition + (1|ID)"
        )
        fe_csv = os.path.join(output_dir_stats, f"AOC_mixedlm_fixed_{sname}_{task['name']}.csv")
        fe_df.to_csv(fe_csv, index=False)
        outputs[f"mixedlm_fixed_{sname}_{task['name']}"] = fe_df.copy()
        print(f"Saved fixed effects  → {os.path.basename(fe_csv)}")

        # %% Figure
        fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
        fig.patch.set_alpha(1.0)
        ax.patch.set_alpha(1.0)
        ax.set_facecolor("white")

        # %% Manual raincloud parameters
        viol_alpha  = 0.60
        dot_alpha   = 0.50
        dot_size    = 30
        box_width   = 0.20
        cloud_offset = -0.20
        max_violsw  = 0.40
        bw_method   = 0.15
        if var == "Accuracy":
            cloud_offset = -0.20
            max_violsw   = 0.50
            bw_method    = 0.25

        # x positions for categories
        xpos = {c: i for i, c in enumerate(condition_order)}

        # Deterministic jitter
        rng = np.random.default_rng(12345)

        # %% Draw per condition
        for cond_lab in condition_order:
            yvals = dvar.loc[dvar["Condition"] == cond_lab, var].dropna().to_numpy()
            if yvals.size == 0:
                continue

            # VIOLIN (left half)
            # Determine hard cap for this variable
            ymax_cap = ylims_map[var][1] if var in ylims_map else float(dvar[var].max())
            ymin_cap = ylims_map[var][0] if var in ylims_map else float(dvar[var].min())
            yr = ymax_cap - ymin_cap

            pad_top = min(0.02 * yr, 0.3)
            y_grid_top = ymin_cap + yr + pad_top

            # build KDE
            kde = gaussian_kde(yvals, bw_method=bw_method)  # <-- REINSERT THIS LINE

            # grid
            y_grid = np.linspace(ymin_cap, y_grid_top, 400)

            dens = kde(y_grid)
            scale = (max_violsw / np.nanmax(dens)) if np.nanmax(dens) > 0 else 0.0

            x_left  = xpos[cond_lab] + cloud_offset - dens * scale
            x_right = np.full_like(y_grid, xpos[cond_lab] + cloud_offset)

            poly_x = np.concatenate([x_right, x_left[::-1]])
            poly_y = np.concatenate([y_grid, y_grid[::-1]])
            ax.fill(poly_x, poly_y,
                    facecolor=pal_dict[cond_lab], edgecolor="none",
                    alpha=viol_alpha, clip_on=True)

            # DOTS
            x_jit = xpos[cond_lab] + rng.uniform(-box_width / 2, box_width / 2, size=yvals.size)
            ax.scatter(x_jit, yvals, s=dot_size, alpha=dot_alpha, color=pal_dict[cond_lab], linewidths=0, zorder=3)

            # BOXPLOT
            bp = ax.boxplot(
                [yvals], positions=[xpos[cond_lab]], widths=box_width, vert=True,
                patch_artist=True, showfliers=False, whis=(5, 95),
                medianprops=dict(color="black", linewidth=1.5),
                boxprops=dict(linewidth=1.0, edgecolor="black"),
                whiskerprops=dict(linewidth=1.0, color="black"),
                capprops=dict(linewidth=1.0, color="black"),
                meanline=False, showmeans=False
            )
            for patch in bp["boxes"]:
                patch.set_facecolor(mpl.colors.to_rgba(pal_dict[cond_lab], 0.05))
                patch.set_edgecolor("black")
            for elem in ["whiskers", "caps", "medians"]:
                for artist in bp[elem]:
                    artist.set_color("black")

        # No legend
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

        # Clean spines and grid
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.yaxis.grid(True, linewidth=1, alpha=0.35)
        ax.xaxis.grid(False)

        # Labels and x-ticks
        ax.set_title("")
        ax.set_xticks(range(len(condition_order)))
        ax.set_xticklabels(condition_order)
        ax.set_xlabel("")
        if len(condition_order) >= 2:
            ax.annotate(
                task["xlabel"],
                xy=(xpos[condition_order[1]], 0),
                xycoords=("data", "axes fraction"),
                xytext=(0, -28),
                textcoords="offset points",
                ha="center", va="top"
            )

        # %% Bracket layout with shared (global) ymax per variable
        # If manual limits exist, use them for bracket layout;
        # otherwise fall back to the data/global upper.
        if var in ylims_map:
            ymin, ymax_cap = ylims_map[var]
        else:
            ymin = float(dvar[var].min()) if np.isfinite(dvar[var].min()) else 0.0
            ymax_data_local = float(dvar[var].max()) if np.isfinite(dvar[var].max()) else ymin
            ymax_cap = global_upper.get(var, np.nan)
            if not np.isfinite(ymax_cap):
                ymax_cap = ymax_data_local

        range_y = max(ymax_cap - ymin, 1.0)

        head = 0.02 * range_y
        step = 0.10 * range_y

        y_positions = []
        start = ymax_cap + 0.075 * range_y
        for i in range(len(task["comparisons"])):
            y_positions.append(start + i * step)

        # y-label at data midpoint
        ymin_cur, ymax_cur = ylims_map[var]
        ymid = 0.5 * (ymin_cur + ymax_cur)
        ax.set_ylabel("")
        ax.yaxis.get_label().set_visible(False)
        ax.text(
            -0.12,
            ymid,
            ylab,
            transform=ax.get_yaxis_transform(which='grid'),
            rotation=90,
            ha='center',
            va='center'
        )

        # Signif labels
        labels = []
        for (g1, g2) in task["comparisons"]:
            row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
            labels.append("n.s." if row.empty else p_to_signif(float(row["p_adj"].iloc[0])))

        # %% slightly increase bracket spacing only for Accuracy in N-back
        if (task["name"] == "nback") and (var == "Accuracy"):
            yr = ax.get_ylim()[1] - ax.get_ylim()[0]
            # assume you already computed y_positions; just spread them a bit more
            step_bump = 0.025 * yr   # +2.5% of axis range on each successive bracket
            y_positions = [y + i * step_bump for i, y in enumerate(y_positions)]

        add_stat_brackets(
            ax=ax,
            xcats=condition_order,
            comparisons=task["comparisons"],
            y_positions=y_positions,
            labels=labels,
            xmap=xpos
        )

        # %% Manual y-ticks (identical for both tasks)
        if var in yticks_map:
            ax.set_yticks(yticks_map[var])

        if var in ylims_map:
            ymin_set, ymax_set = ylims_map[var]
            # extend a bit to make sure brackets are inside
            extra = 0.12 * (ymax_set - ymin_set)
            ax.set_ylim(ymin_set, ymax_set + extra)
        else:
            ax.set_ylim(ymin, ymax_cap + 0.15 * (ymax_cap - ymin))


        # %% Save raincloud figure for each variable
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_dir, f"AOC_stats_rainclouds_{sname}_{task['name']}.png"),
            dpi=300,
            transparent=False,
            facecolor=fig.get_facecolor(),
            edgecolor=fig.get_edgecolor() if hasattr(fig, "get_edgecolor") else "white"
        )
        plt.close(fig)
        saveNameFig = os.path.join(f"AOC_stats_rainclouds_{sname}_{task['name']}.png")
        print(f"Saved raincloud fig. → {saveNameFig}")

    # Save ANOVA for this task
    anova_df = pd.DataFrame(
        anova_rows,
        columns=["Task", "Variable", "Effect", "DF_num", "DF_den", "F", "p", "eta_p2"]
    )
    anova_csv = os.path.join(anova_dir, f"AOC_anova_{task['name']}.csv")
    anova_df.to_csv(anova_csv, index=False)
    outputs[f"anova_{task['name']}"] = anova_df.copy()
    print(f"Saved ANOVA → {anova_csv}")

    # Save pairwise effect sizes for this task
    pw_eff_df = pd.DataFrame(
        pairwise_effsize_rows,
        columns=["Task", "Variable", "Group1", "Group2", "N", "MeanDiff", "Cohens_dz", "CI95_low", "CI95_high", "p_adj"]
    )
    pw_csv = os.path.join(output_dir_stats, f"AOC_pairwise_effectsizes_{task['name']}.csv")
    pw_eff_df.to_csv(pw_csv, index=False)
    outputs[f"pairwise_effect_sizes_{task['name']}"] = pw_eff_df.copy()
    print(f"Saved pairwise effect sizes → {pw_csv}")

# %% Alpha ~ Gaze * Condition + (1|ID) — run twice (GazeDeviation, MSRate)
print("\n=== Alpha ~ Gaze * Condition mixed models (per task) ===")

for task in tasks:
    print(f"\nTask: {task['name']}")
    dat = pd.read_csv(task["input_csv"])
    dat.loc[dat["Accuracy"] > 100, "Accuracy"] = np.nan

    # Harmonise Condition labels and order (same logic as above)
    cond = dat["Condition"]
    ref_level = task["categories"][0]
    if np.issubdtype(cond.dtype, np.number):
        uniq = sorted(pd.unique(cond.dropna()).tolist())
        applied_map = None
        for cand in task["cond_to_label_numeric"]:
            if set(uniq).issubset(set(cand.keys())):
                applied_map = cand
                break
        if applied_map is None:
            applied_map = {val: task["categories"][i] for i, val in enumerate(uniq[:len(task["categories"])])}
        dat["Condition"] = dat["Condition"].map(applied_map)
    else:
        dat["Condition"] = dat["Condition"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()

    dat["Condition"] = pd.Categorical(dat["Condition"], categories=task["categories"], ordered=True)
    if dat["ID"].dtype != "O":
        dat["ID"] = dat["ID"].astype(str)

    # Apply same outlier handling
    dat = iqr_outlier_filter(dat, variables, by="Condition")

    # Keep essential columns
    cols_needed = ["ID", "Condition", "AlphaPower", "GazeDeviation", "MSRate"]
    dat = dat[[c for c in cols_needed if c in dat.columns]].copy()

    # Drop rows with missing DV or predictors
    # We'll do this separately for each gaze predictor
    for gaze_var in ["GazeDeviation", "MSRate"]:
        if gaze_var not in dat.columns:
            print(f" - {gaze_var} not present. Skipping.")
            continue

        sub = dat[["ID", "Condition", "AlphaPower", gaze_var]].dropna().copy()

        # Ensure we have all condition levels for at least some subjects
        if sub.empty or sub["ID"].nunique() < 2:
            print(f" - Not enough data for {gaze_var}. Skipping.")
            continue

        # Mean-centre gaze within task (keeps units; recommended to reduce collinearity)
        sub[gaze_var + "_c"] = sub[gaze_var] - sub[gaze_var].mean()

        # MixedLM (REML=False for LRT comparability)
        # Full model with interaction
        formula_full = (
            f'AlphaPower ~ {gaze_var}_c * C(Condition, Treatment(reference="{ref_level}"))'
        )
        try:
            full_res = fit_mixedlm(formula_full, data=sub, group="ID", reml=False)
        except Exception as e:
            # Fallback: try different optimiser
            full_res = fit_mixedlm(formula_full, data=sub, group="ID", reml=False, method="nm", maxiter=800)

        # Reduced (drop interaction)
        formula_red = (
            f'AlphaPower ~ {gaze_var}_c + C(Condition, Treatment(reference="{ref_level}"))'
        )        
        try:
            red_res = fit_mixedlm(formula_red, data=sub, group="ID", reml=False)
        except Exception as e:
            red_res = fit_mixedlm(formula_red, data=sub, group="ID", reml=False, method="nm", maxiter=800)

        # drop1-style LRT for interaction
        lrt = drop1_lrt(full_res, red_res)
        lrt_df = pd.DataFrame([{
            "Task": task["name"],
            "GazePredictor": gaze_var,
            "LL_full": lrt["LL_full"], "LL_reduced": lrt["LL_reduced"],
            "df_full": lrt["df_full"], "df_reduced": lrt["df_reduced"],
            "df_diff": lrt["df_diff"], "LR": lrt["LR"], "p": lrt["p"]
        }])

        # Select final model based on LRT
        if np.isfinite(lrt["p"]) and (lrt["p"] < 0.05):
            final_res = full_res
            final_formula = formula_full
            interaction_kept = True
        else:
            final_res = red_res
            final_formula = formula_red
            interaction_kept = False

        # LRT "ANOVA" style tests (all based on the no-interaction model red_res)
        # - Interaction: already tested via lrt (full_res vs red_res)
        # - Condition main effect: red_res vs model without Condition
        # - Gaze main effect: red_res vs model without gaze
        #
        # red_res is the no-interaction model that is already fitted.
        # lrt is still the LRT for the interaction (full_res vs red_res).
        # lrt_cond and lrt_gaze are LRTs on the reduced model:
        # Condition: compare AlphaPower ~ gaze_c + Condition vs AlphaPower ~ gaze_c
        # Gaze: compare AlphaPower ~ gaze_c + Condition vs AlphaPower ~ Condition

        # Model without Condition (keeps gaze only)
        formula_nocond = (
            f'AlphaPower ~ {gaze_var}_c'
        )
        try:
            res_nocond = fit_mixedlm(formula_nocond, data=sub, group="ID", reml=False)
        except Exception:
            res_nocond = fit_mixedlm(formula_nocond, data=sub, group="ID",
                                     reml=False, method="nm", maxiter=800)

        lrt_cond = drop1_lrt(red_res, res_nocond)

        # Model without gaze (keeps Condition only)
        formula_nogaze = (
            f'AlphaPower ~ C(Condition, Treatment(reference="{ref_level}"))'
        )
        try:
            res_nogaze = fit_mixedlm(formula_nogaze, data=sub, group="ID", reml=False)
        except Exception:
            res_nogaze = fit_mixedlm(formula_nogaze, data=sub, group="ID",
                                     reml=False, method="nm", maxiter=800)

        lrt_gaze = drop1_lrt(red_res, res_nogaze)

        # Build LRT table with LR-based effect sizes
        rows_lrt = []
        n_obs = sub.shape[0]

        R2_cond, f2_cond = lr_effect_sizes(lrt_cond["LR"], lrt_cond["df_diff"], n_obs)
        rows_lrt.append({
            "Term": "Condition",
            "df":   lrt_cond["df_diff"],
            "LR_Chi2": lrt_cond["LR"],
            "p":    lrt_cond["p"],
            "R2_LR": R2_cond,
            "f2_LR": f2_cond
        })

        R2_gaze, f2_gaze = lr_effect_sizes(lrt_gaze["LR"], lrt_gaze["df_diff"], n_obs)
        rows_lrt.append({
            "Term": f"{gaze_var}_c",
            "df":   lrt_gaze["df_diff"],
            "LR_Chi2": lrt_gaze["LR"],
            "p":    lrt_gaze["p"],
            "R2_LR": R2_gaze,
            "f2_LR": f2_gaze
        })

        # Interaction row
        if interaction_kept and np.isfinite(lrt["p"]):
            R2_int, f2_int = lr_effect_sizes(lrt["LR"], lrt["df_diff"], n_obs)
            rows_lrt.append({
                "Term": "Interaction",
                "df":   lrt["df_diff"],
                "LR_Chi2": lrt["LR"],
                "p":    lrt["p"],
                "R2_LR": R2_int,
                "f2_LR": f2_int
            })

        for row in rows_lrt:
            row["N_obs"] = int(n_obs)

        lrtTbl_df = pd.DataFrame(rows_lrt, columns=["Term", "df", "LR_Chi2", "p", "R2_LR", "f2_LR"])
        lrtTbl_df.insert(0, "Task", task["name"])
        lrtTbl_df.insert(1, "GazePredictor", gaze_var)

        # Pairwise contrasts across Condition at mean gaze (centred = 0)
        cond_levels = list(sub["Condition"].cat.categories)
        design_prefix = f'C(Condition, Treatment(reference="{ref_level}"))'
        contrast_df = pairwise_condition_contrasts_at_mean_gaze(
            final_res, cond_levels, design_prefix=design_prefix
        )        
        contrast_df.insert(0, "Task", task["name"])
        contrast_df.insert(1, "GazePredictor", gaze_var)

        # Export model table (DOCX)
        doc_name = f"AOC_alpha_mixedlm_{gaze_var.lower()}_{task['name']}.docx"
        doc_path = os.path.join(output_dir_stats, doc_name)
        export_model_table(final_res, doc_path)
        print(f"Saved model table → {doc_name}")

        # Also save fixed effects for the alpha ~ gaze * Condition model
        fe_alpha_df = mixedlm_fixed_effects_to_df(
            final_res,
            task=task["name"],
            variable="AlphaPower",
            model_label=f"AlphaPower ~ {gaze_var}_c * Condition + (1|ID)" if interaction_kept
                        else f"AlphaPower ~ {gaze_var}_c + Condition + (1|ID)"
        )
        fe_alpha_csv = os.path.join(output_dir_stats, f"AOC_alpha_fixed_{gaze_var.lower()}_{task['name']}.csv")
        fe_alpha_df.to_csv(fe_alpha_csv, index=False)
        outputs[f"alpha_fixed_{gaze_var.lower()}_{task['name']}"] = fe_alpha_df.copy()
        print(f"Saved alpha fixed-effects → {os.path.basename(fe_alpha_csv)}")

        # Save LRT, lrtTbl table, and contrasts
        lrt_csv = os.path.join(output_dir_stats, f"AOC_alpha_drop1_{gaze_var.lower()}_{task['name']}.csv")
        lrtTbl_csv = os.path.join(output_dir_stats, f"AOC_alpha_lrtTbl_{gaze_var.lower()}_{task['name']}.csv")
        con_csv  = os.path.join(output_dir_stats, f"AOC_alpha_contrasts_{gaze_var.lower()}_{task['name']}.csv")
        lrt_df.to_csv(lrt_csv, index=False)
        lrtTbl_df.to_csv(lrtTbl_csv, index=False)
        contrast_df.to_csv(con_csv, index=False)
        outputs[f"alpha_drop1_{gaze_var.lower()}_{task['name']}"] = lrt_df.copy()
        outputs[f"alpha_lrtTbl_{gaze_var.lower()}_{task['name']}"] = lrtTbl_df.copy()
        outputs[f"alpha_contrasts_{gaze_var.lower()}_{task['name']}"] = contrast_df.copy()
        print(f"Saved LRT/LRT Table/Contrasts → {os.path.basename(lrt_csv)}, {os.path.basename(lrtTbl_csv)}, {os.path.basename(con_csv)}")

# %% Display all stored dataframes
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 0)
print("\n=== Output DataFrames ===\n")
for name, df in outputs.items():
    print(f"{name}: {df.shape[0]} rows × {df.shape[1]} columns")
    print(df.to_string(index=False))
    print("\n")
