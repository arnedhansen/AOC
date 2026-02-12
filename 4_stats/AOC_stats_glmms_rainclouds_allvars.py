# %% AOC Stats Rainclouds LONG — All Variables (Sternberg + N-Back)
# Loads merged_data_*.csv for both tasks, runs rm-ANOVAs and pairwise
# mixed-model contrasts, generates raincloud plots for every numeric DV.
# Saves figures to figures/stats/rainclouds/allvars.
#
# Key outputs:
#   One raincloud figure per variable × task
#   ANOVA summary CSV per task
#   Pairwise effect-size CSV per task

# %% Imports
import sys, os
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy import stats as sp_stats
from statsmodels.stats.anova import AnovaRM

from functions.stats_helpers import (
    iqr_outlier_filter, mixedlm_pairwise_contrasts, p_to_signif
)
from functions.rainclouds_plotting_helpers import add_stat_brackets

# %% Parameters

# Colours (AOC pastels)
pal = ["#93B8C4", "#82AD82", "#D998A2"]

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

import seaborn as sns
sns.set_style("white")

# %% I/O Paths
base_dir   = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
output_dir = f"{base_dir}/figures/stats/rainclouds/allvars"
os.makedirs(output_dir, exist_ok=True)

# %% Variables — every numeric DV in the LONG files
# (ID, Gender, Age, Handedness, OcularDominance, Condition are excluded)

variables = [
    # --- Original ---
    "Accuracy", "ReactionTime",
    "GazeDeviation", "GazeStdX", "GazeStdY",
    "PupilSize", "MSRate", "Blinks", "Fixations", "Saccades", "ScanPathLength",
    "BCEA", "BCEALateralization",
    "AlphaPower", "IAF", "Lateralization",
    # --- FOOOF alpha ---
    "AlphaPower_FOOOF",
    "AlphaPower_FOOOF_bl",
    "AlphaPower_FOOOF_bl_early",
    "AlphaPower_FOOOF_bl_late",
    # --- Baselined raw alpha ---
    "AlphaPower_bl",
    "AlphaPower_bl_early",
    "AlphaPower_bl_late",
    # --- Baselined gaze ---
    "GazeDeviationFullBL", "GazeDeviationEarlyBL", "GazeDeviationLateBL",
    "ScanPathLengthFullBL", "ScanPathLengthEarlyBL", "ScanPathLengthLateBL",
    "PupilSizeFullBL", "PupilSizeEarlyBL", "PupilSizeLateBL",
    "MSRateFullBL", "MSRateEarlyBL", "MSRateLateBL",
    # --- Baselined BCEA ---
    "BCEAFullBL", "BCEAEarlyBL", "BCEALateBL",
    "BCEALatFullBL", "BCEALatEarlyBL", "BCEALatLateBL",
]

# Pretty labels (auto-generated; override specific ones below)
_unit_map = {
    "Accuracy": "[%]",
    "ReactionTime": "[s]",
    "GazeDeviation": "[px]", "GazeStdX": "[px]", "GazeStdY": "[px]",
    "PupilSize": "[a.u.]",
    "MSRate": "[MS/s]",
    "Blinks": "", "Fixations": "", "Saccades": "",
    "ScanPathLength": "[px]",
    "BCEA": "[px²]",
    "BCEALateralization": "[L−R]",
    "AlphaPower": "[\u03BCV\u00B2/Hz]",
    "IAF": "[Hz]",
    "Lateralization": "[R\u2212L / R+L]",
    "AlphaPower_FOOOF": "[FOOOF log]",
    "AlphaPower_FOOOF_bl": "[FOOOF log, BL]",
    "AlphaPower_FOOOF_bl_early": "[FOOOF log, BL early]",
    "AlphaPower_FOOOF_bl_late": "[FOOOF log, BL late]",
    "AlphaPower_bl": "[dB]",
    "AlphaPower_bl_early": "[dB, early]",
    "AlphaPower_bl_late": "[dB, late]",
    "GazeDeviationFullBL": "[dB]", "GazeDeviationEarlyBL": "[dB]", "GazeDeviationLateBL": "[dB]",
    "ScanPathLengthFullBL": "[dB]", "ScanPathLengthEarlyBL": "[dB]", "ScanPathLengthLateBL": "[dB]",
    "PupilSizeFullBL": "[%\u0394]", "PupilSizeEarlyBL": "[%\u0394]", "PupilSizeLateBL": "[%\u0394]",
    "MSRateFullBL": "[dB]", "MSRateEarlyBL": "[dB]", "MSRateLateBL": "[dB]",
    "BCEA": "[px\u00B2]",
    "BCEALateralization": "[L\u2212R]",
    "BCEAFullBL": "[dB]", "BCEAEarlyBL": "[dB]", "BCEALateBL": "[dB]",
    "BCEALatFullBL": "[\u0394]", "BCEALatEarlyBL": "[\u0394]", "BCEALatLateBL": "[\u0394]",
}

def _pretty_label(var):
    """Human-readable axis label from variable name."""
    # Insert spaces before capitals and underscores
    nice = var.replace("_", " ")
    unit = _unit_map.get(var, "")
    return f"{nice} {unit}".strip()

def _save_name(var):
    """File-safe short name."""
    return var.lower()

# %% Task configurations
tasks = [
    {
        "name"       : "sternberg",
        "input_csv"  : f"{base_dir}/data/features/merged_data_sternberg.csv",
        "cond_to_label_numeric": [
            {1: "WM load 2", 2: "WM load 4", 3: "WM load 6"},
            {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"},
        ],
        "categories" : ["WM load 2", "WM load 4", "WM load 6"],
        "comparisons": [("WM load 2", "WM load 4"),
                        ("WM load 2", "WM load 6"),
                        ("WM load 4", "WM load 6")],
        "xlabel"     : "Condition",
    },
    {
        "name"       : "nback",
        "input_csv"  : f"{base_dir}/data/features/merged_data_nback.csv",
        "cond_to_label_numeric": [{1: "1-back", 2: "2-back", 3: "3-back"}],
        "categories" : ["1-back", "2-back", "3-back"],
        "comparisons": [("1-back", "2-back"),
                        ("1-back", "3-back"),
                        ("2-back", "3-back")],
        "xlabel"     : "Condition",
    },
]

# %% Helper: condition-label harmonisation
def harmonise_conditions(dat, task):
    """Map numeric or string Condition column to ordered categorical."""
    cond = dat["Condition"]
    if np.issubdtype(cond.dtype, np.number):
        uniq = sorted(pd.unique(cond.dropna()).tolist())
        applied_map = None
        for cand in task["cond_to_label_numeric"]:
            if set(uniq).issubset(set(cand.keys())):
                applied_map = cand
                break
        if applied_map is None:
            applied_map = {val: task["categories"][i]
                           for i, val in enumerate(uniq[:len(task["categories"])])}
        dat["Condition"] = dat["Condition"].map(applied_map)
    else:
        dat["Condition"] = (dat["Condition"].astype(str)
                            .str.replace(r"\s+", " ", regex=True).str.strip())
    dat["Condition"] = pd.Categorical(dat["Condition"],
                                      categories=task["categories"], ordered=True)
    if dat["ID"].dtype != "O":
        dat["ID"] = dat["ID"].astype(str)
    return dat


# %% Helper: draw one raincloud
def draw_raincloud(ax, dvar, var, condition_order, pal_dict, *, bw_method=0.15):
    """Draw half-violin + jittered dots + boxplot for each condition."""
    viol_alpha  = 0.60
    dot_alpha   = 0.50
    dot_size    = 30
    box_width   = 0.20
    cloud_offset = -0.20
    max_violsw  = 0.40

    if var == "Accuracy":
        cloud_offset = -0.20
        max_violsw   = 0.50
        bw_method    = 0.25

    xpos = {c: i for i, c in enumerate(condition_order)}
    rng = np.random.default_rng(12345)

    for cond_lab in condition_order:
        yvals = dvar.loc[dvar["Condition"] == cond_lab, var].dropna().to_numpy()
        yvals = yvals[np.isfinite(yvals)]
        if yvals.size < 3 or np.nanstd(yvals) == 0:
            continue

        # KDE (half-violin on the left)
        try:
            kde = gaussian_kde(yvals, bw_method=bw_method)
        except (np.linalg.LinAlgError, ValueError):
            continue
        y_lo, y_hi = yvals.min(), yvals.max()
        yr = y_hi - y_lo if y_hi > y_lo else 1.0
        y_grid = np.linspace(y_lo - 0.05 * yr, y_hi + 0.05 * yr, 400)
        dens = kde(y_grid)
        scale = (max_violsw / np.nanmax(dens)) if np.nanmax(dens) > 0 else 0.0

        x_left  = xpos[cond_lab] + cloud_offset - dens * scale
        x_right = np.full_like(y_grid, xpos[cond_lab] + cloud_offset)
        poly_x = np.concatenate([x_right, x_left[::-1]])
        poly_y = np.concatenate([y_grid, y_grid[::-1]])
        ax.fill(poly_x, poly_y, facecolor=pal_dict[cond_lab],
                edgecolor="none", alpha=viol_alpha, clip_on=True)

        # Jittered dots
        x_jit = xpos[cond_lab] + rng.uniform(-box_width / 2, box_width / 2,
                                               size=yvals.size)
        ax.scatter(x_jit, yvals, s=dot_size, alpha=dot_alpha,
                   color=pal_dict[cond_lab], linewidths=0, zorder=3)

        # Boxplot
        bp = ax.boxplot(
            [yvals], positions=[xpos[cond_lab]], widths=box_width, vert=True,
            patch_artist=True, showfliers=False, whis=(5, 95),
            medianprops=dict(color="black", linewidth=1.5),
            boxprops=dict(linewidth=1.0, edgecolor="black"),
            whiskerprops=dict(linewidth=1.0, color="black"),
            capprops=dict(linewidth=1.0, color="black"),
            meanline=False, showmeans=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(mpl.colors.to_rgba(pal_dict[cond_lab], 0.05))
            patch.set_edgecolor("black")
        for elem in ["whiskers", "caps", "medians"]:
            for artist in bp[elem]:
                artist.set_color("black")

    return xpos


# ============================================================
# %% Main loop
# ============================================================
all_anova_rows = []
all_effsize_rows = []

for task in tasks:
    print(f"\n{'='*60}")
    print(f"  Task: {task['name'].upper()}")
    print(f"{'='*60}")

    # Load & prep
    dat = pd.read_csv(task["input_csv"])
    dat.loc[dat.get("Accuracy", pd.Series(dtype=float)) > 100, "Accuracy"] = np.nan
    dat = harmonise_conditions(dat, task)

    # Determine which LONG variables are actually present in this CSV
    task_vars = [v for v in variables if v in dat.columns]

    # Replace inf / -inf with NaN (e.g. dB baselines where denominator was 0)
    dat[task_vars] = dat[task_vars].replace([np.inf, -np.inf], np.nan)

    # Outlier removal (1.5×IQR per condition)
    dat = iqr_outlier_filter(dat, task_vars, by="Condition")

    condition_order = list(dat["Condition"].cat.categories)
    pal_dict = dict(zip(condition_order, pal))
    xpos_map = {c: i for i, c in enumerate(condition_order)}

    # --------------------------------------------------------
    # Loop over variables
    # --------------------------------------------------------
    for var in task_vars:
        dvar = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
        if dvar.empty or dvar[var].nunique() < 2:
            continue

        dvar["Condition"] = pd.Categorical(dvar["Condition"],
                                           categories=condition_order, ordered=True)

        # ------ rm-ANOVA ------
        present = dvar.groupby("ID")["Condition"].nunique()
        keep_ids = present[present == len(condition_order)].index
        dvar_bal = dvar[dvar["ID"].isin(keep_ids)].copy()

        if dvar_bal["ID"].nunique() >= 2:
            try:
                aov = AnovaRM(data=dvar_bal, depvar=var,
                              subject="ID", within=["Condition"]).fit()
                aov_tab = aov.anova_table.reset_index().rename(columns={"index": "Effect"})
                for _, r in aov_tab.iterrows():
                    F   = float(r["F Value"])
                    df1 = float(r["Num DF"])
                    df2 = float(r["Den DF"])
                    p   = float(r["Pr > F"])
                    etap = (F * df1) / (F * df1 + df2) if np.isfinite(F) else np.nan
                    all_anova_rows.append([task["name"], var, r["Effect"],
                                           df1, df2, F, p, etap])
                    print(f"  rm-ANOVA: {var} [{task['name']}]  "
                          f"F({df1:.0f},{df2:.0f}) = {F:.3f}, p = {p:.4f}, η²p = {etap:.3f}")
            except Exception:
                all_anova_rows.append([task["name"], var, "Condition",
                                       np.nan, np.nan, np.nan, np.nan, np.nan])
        else:
            all_anova_rows.append([task["name"], var, "Condition",
                                   np.nan, np.nan, np.nan, np.nan, np.nan])

        # ------ Pairwise mixed-model contrasts ------
        try:
            pw = mixedlm_pairwise_contrasts(
                dvar.rename(columns={var: "value"}),
                value_col="value", group_col="Condition",
                id_col="ID", p_adjust="bonferroni",
            )
        except Exception:
            pw = pd.DataFrame(columns=["group1", "group2", "p_adj"])

        # Print pairwise to console
        if not pw.empty:
            for _, r in pw.iterrows():
                sig = p_to_signif(float(r['p_adj'])) if 'p_adj' in r else ''
                print(f"    {r['group1']} vs {r['group2']}: p_adj = {float(r['p_adj']):.4f} {sig}")

        # ------ Pairwise Cohen's dz ------
        wide = dvar.pivot(index="ID", columns="Condition", values=var)
        for g1, g2 in task["comparisons"]:
            if g1 in wide.columns and g2 in wide.columns:
                diffs = (wide[g2] - wide[g1]).dropna()
                n = diffs.shape[0]
                if n >= 2 and np.nanstd(diffs, ddof=1) > 0:
                    md   = float(np.nanmean(diffs))
                    sd_d = float(np.nanstd(diffs, ddof=1))
                    dz   = md / sd_d
                    tcrit = sp_stats.t.ppf(0.975, df=n - 1)
                    se_md = sd_d / np.sqrt(n)
                    ci_lo = md - tcrit * se_md
                    ci_hi = md + tcrit * se_md
                else:
                    md = dz = ci_lo = ci_hi = np.nan
                padj = np.nan
                row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
                if not row.empty and "p_adj" in row:
                    try:
                        padj = float(row["p_adj"].iloc[0])
                    except Exception:
                        pass
                all_effsize_rows.append([task["name"], var, g1, g2,
                                         int(n), md, dz, ci_lo, ci_hi, padj])

        # ======================================================
        # FIGURE
        # ======================================================
        fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
        fig.patch.set_alpha(1.0)
        ax.set_facecolor("white")

        xpos = draw_raincloud(ax, dvar, var, condition_order, pal_dict)

        # Remove legend if any
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

        # Spines & grid
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.yaxis.grid(True, linewidth=1, alpha=0.35)
        ax.xaxis.grid(False)

        # X axis
        ax.set_xticks(range(len(condition_order)))
        ax.set_xticklabels(condition_order)
        ax.set_xlabel("")
        ax.set_title("")
        if len(condition_order) >= 2:
            ax.annotate(
                task["xlabel"],
                xy=(xpos[condition_order[1]], 0),
                xycoords=("data", "axes fraction"),
                xytext=(0, -28), textcoords="offset points",
                ha="center", va="top",
            )

        # Y label
        ylab = _pretty_label(var)
        ax.set_ylabel(ylab)

        # ------ Auto y-limits ------
        vals_all = dvar[var].dropna()
        ymin_data = float(vals_all.min())
        ymax_data = float(vals_all.max())
        yr = ymax_data - ymin_data if ymax_data > ymin_data else 1.0
        pad_lo = 0.08 * yr
        pad_hi = 0.35 * yr  # extra room for brackets
        ax.set_ylim(ymin_data - pad_lo, ymax_data + pad_hi)

        # ------ Significance brackets ------
        labels = []
        for g1, g2 in task["comparisons"]:
            row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
            labels.append("n.s." if row.empty
                          else p_to_signif(float(row["p_adj"].iloc[0])))

        head = 0.02 * yr
        step = 0.10 * yr
        start = ymax_data + 0.075 * yr
        y_positions = [start + i * step for i in range(len(task["comparisons"]))]

        add_stat_brackets(
            ax=ax,
            xcats=condition_order,
            comparisons=task["comparisons"],
            y_positions=y_positions,
            labels=labels,
            xmap=xpos,
        )

        # ------ Save ------
        sname = _save_name(var)
        fname = f"AOC_rainclouds_{sname}_{task['name']}.png"
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, fname), dpi=300,
                    transparent=False, facecolor="white")
        plt.close(fig)
        print(f"  {fname}")

# ============================================================
# %% Save summary tables
# ============================================================
anova_df = pd.DataFrame(
    all_anova_rows,
    columns=["Task", "Variable", "Effect", "DF_num", "DF_den",
             "F", "p", "eta_p2"],
)
anova_csv = os.path.join(output_dir, "AOC_anova_allvars.csv")
anova_df.to_csv(anova_csv, index=False)
print(f"\nSaved ANOVA table → {anova_csv}")

effsize_df = pd.DataFrame(
    all_effsize_rows,
    columns=["Task", "Variable", "Group1", "Group2", "N",
             "MeanDiff", "Cohens_dz", "CI95_low", "CI95_high", "p_adj"],
)
effsize_csv = os.path.join(output_dir, "AOC_pairwise_effectsizes_allvars.csv")
effsize_df.to_csv(effsize_csv, index=False)
print(f"Saved pairwise effect sizes → {effsize_csv}")

# %% Done
print(f"\n{'='*60}")
print(f"  All done — {len(all_anova_rows)} ANOVA rows, "
      f"{len(all_effsize_rows)} pairwise rows")
n_figs = sum(1 for t in tasks
             for v in variables if v in pd.read_csv(t["input_csv"], nrows=0).columns)
print(f"  {n_figs} raincloud figures saved to {output_dir}")
print(f"{'='*60}")
