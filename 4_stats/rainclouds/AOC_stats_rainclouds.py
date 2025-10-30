# %% AOC Stats Rainclouds — Combined (Sternberg + N-back)

# %% Imports
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

from rainclouds_stats_helpers import (
    iqr_outlier_filter,
    mixedlm_pairwise_contrasts,
    p_to_signif
)

from rainclouds_plotting_helpers import (
    add_stat_brackets,
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

# %% Variables and labelling

variables  = ["Accuracy", "ReactionTime", "GazeDeviation", "MSRate", "Fixations", "Saccades", "AlphaPower", "IAF"]
titles     = ["Accuracy", "Reaction Time", "GazeDeviation", "Microsaccade Rate", "Fixations", "Saccades", "Alpha Power", "IAF"]
y_labels   = ["Accuracy [%]", "Reaction Time [s]", "Gaze Deviation [px]", "Microsaccade Rate [MS/s]", "Fixations", "Saccades", "Alpha Power [\u03BCV²/Hz]", "IAF [Hz]"]
save_names = ["acc", "rt", "gazedev", "ms", "fix", "sacc", "pow", "iaf"]

# Manual y ticks and ylims per variable
yticks_map = {
    "Accuracy"     : np.arange(60, 101, 5),
    "ReactionTime" : np.arange(0.3, 1.35, 0.1),
    "GazeDeviation": np.arange(0, 125, 10),
    "MSRate"       : np.arange(0, 4.1, 0.5),
    "Fixations"    : np.arange(0, 8.5, 1),
    "Saccades"     : np.arange(0, 4.25, 1),
    "AlphaPower"   : np.arange(0, 1.52, 0.25),
    "IAF"          : np.arange(8, 14, 1),
}
ylims_map = {
    "Accuracy"     : (60, 102),
    "ReactionTime" : (0.3, 1.35),
    "GazeDeviation": (0, 125),
    "MSRate"       : (0, 4.1),
    "Fixations"    : (0, 8.5),
    "Saccades"     : (0, 4.25),
    "AlphaPower"   : (0, 1.6),
    "IAF"          : (8, 14),
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

for task in tasks:

    # --- Load
    dat = pd.read_csv(task["input_csv"])

    # Remove impossible values
    dat.loc[dat["Accuracy"] > 100, "Accuracy"] = np.nan

    # --- Condition labelling (robust to numeric or already-labelled strings)
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

    # Category order / palette mapping
    condition_order = list(dat["Condition"].dropna().unique())
    pal_dict = dict(zip(condition_order, pal))

    # --- Loop variables
    for var, ttl, ylab, sname in zip(variables, titles, y_labels, save_names):

        dvar = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            continue

        dvar["Condition"] = pd.Categorical(dvar["Condition"], categories=condition_order, ordered=True)

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

        # === Figure
        fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
        fig.patch.set_alpha(1.0)
        ax.patch.set_alpha(1.0)
        ax.set_facecolor("white")

        # --- Manual raincloud parameters
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

        # --- Draw per condition
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

        # --- Bracket layout with shared (global) ymax per variable
        ymin = float(dvar[var].min()) if np.isfinite(dvar[var].min()) else 0.0
        ymax_data_local = float(dvar[var].max()) if np.isfinite(dvar[var].max()) else ymin
        ymax_cap = global_upper.get(var, np.nan)
        if not np.isfinite(ymax_cap):
            ymax_cap = ymax_data_local  # fallback if pre-scan found nothing

        # use the global cap for bracket baseline and ylim ceiling
        range_y = max(ymax_cap - ymin, 1.0)
        head = 0.02 * range_y
        step = 0.10 * range_y

        y_positions = []
        start = ymax_cap + 0.075 * range_y
        for i in range(len(task["comparisons"])):
            y_positions.append(start + i * step)

        # y-label at data midpoint
        ymin_cur, ymax_cur = ylims_map[var]
        ymid = (ymin_cur + ymax_cap) / 2.0
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

        # --- FIX: slightly increase bracket spacing only for Accuracy in N-back
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

        # --- Manual y-ticks (identical for both tasks)
        if var in yticks_map:
            ax.set_yticks(yticks_map[var])

        if var in ylims_map:
            ymin_set, ymax_set = ylims_map[var]
            ax.set_ylim(ymin_set, ymax_set)

        # --- Save
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_dir, f"AOC_stats_rainclouds_{sname}_{task['name']}_stats.png"),
            dpi=300,
            transparent=False,
            facecolor=fig.get_facecolor(),
            edgecolor=fig.get_edgecolor() if hasattr(fig, "get_edgecolor") else "white"
        )
        plt.close(fig)
