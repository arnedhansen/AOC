# AOC Stats Rainclouds Sternberg

#python -m pip install <package>
#conda install -n aoc-py311 -c conda-forge <package>
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt

from stats_helpers import (
    iqr_outlier_filter,
    mixedlm_pairwise_contrasts,
    p_to_signif
)
from plotting_helpers import (
    add_stat_brackets,
    ensure_dir
)


# I/O
input_csv  = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/merged_data_nback.csv"
output_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/rainclouds"
ensure_dir(output_dir)


# Appearance
pal = ["#93B8C4", "#82AD82", "#D998A2"]  # AOC pastels
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "axes.spines.right": False,
    "axes.spines.top": False
})

# Load & derive
dat = pd.read_csv(input_csv)

# harmonise columns like in R
dat["GazeStd"] = (dat["GazeStdX"] + dat["GazeStdY"]) / 2
dat["Condition"] = dat["Condition"].map({1: "1-back", 2: "2-back", 3: "3-back"})
dat["Condition"] = pd.Categorical(dat["Condition"], categories=["1-back", "2-back", "3-back"], ordered=True)
if dat["ID"].dtype != "O":
    dat["ID"] = dat["ID"].astype(str)

variables  = ["Accuracy", "ReactionTime"]
titles     = ["Accuracy", "Reaction Time"]
y_labels   = ["Accuracy [%]", "Reaction Time [s]"]
save_names = ["acc", "rt"]

# outlier removal per condition & variable (1.5×IQR), like your R code
dat = iqr_outlier_filter(dat, variables, by="Condition")

# comparisons (must match Condition levels)
comparisons = [("1-back", "2-back"), ("1-back", "3-back"), ("2-back", "3-back")]

# Loop variables
for var, ttl, ylab, sname in zip(variables, titles, y_labels, save_names):
    dvar = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
    if dvar.empty:
        continue

    # y-limits & tick scaffolding to mirror your R choices
    if var == "Accuracy":
        lower_bound = 64
        nominal_upper = 100
    elif var == "GazeDeviation":
        lower_bound = 5
        nominal_upper = None
    else:
        # data-driven
        vmin = dvar[var].min()
        vmax = dvar[var].max()
        lower_bound = vmin if np.isfinite(vmin) else None
        nominal_upper = vmax if np.isfinite(vmax) else None

    
    # Base raincloud
    fig, ax = plt.subplots(figsize=(8, 6))

    # RainCloud: one draw per condition
    # pt.RainCloud expects a long-form df; we already have it
    pt.RainCloud(
        x="Condition",
        y=var,
        data=dvar,
        palette=pal,
        bw=.2,
        width_viol=.6,
        width_box=.25,
        move=.2,
        pointplot=True,
        box_showfliers=False,
        ax=ax
    )

    ax.set_xlabel("Condition")
    ax.set_ylabel(ylab)
    ax.set_title(ttl)
    ax.set_ylim(bottom=lower_bound)  # top will be refined below

    # nice ticks for Accuracy / RT similar to your R code
    if var == "Accuracy":
        ax.set_yticks(np.arange(65, 101, 5))
        ymax = 100
    elif var == "ReactionTime":
        ax.set_yticks(np.round(np.arange(0.3, 1.31, 0.1), 2))
        ymax = max( dvar[var].max(), ax.get_ylim()[1] )
    else:
        ymax = max( dvar[var].max(), ax.get_ylim()[1] )

    ax.set_ylim(lower_bound, ymax)

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"AOC_stats_rainclouds_{sname}_nback.png"))
    plt.close(fig)

    
    # Mixed model + pairwise (emmeans-like) + brackets
    
    # Get pairwise contrasts from MixedLM with random intercept for ID
    # Bonferroni-adjust to mirror your R call
    pw = mixedlm_pairwise_contrasts(
        dvar.rename(columns={var: "value"}),
        value_col="value",
        group_col="Condition",
        id_col="ID",
        p_adjust="bonferroni"
    )

    # Build stats plot: start from a fresh raincloud to avoid artefacts
    fig, ax = plt.subplots(figsize=(8, 6))
    pt.RainCloud(
        x="Condition",
        y=var,
        data=dvar,
        palette=pal,
        bw=.2,
        width_viol=.6,
        width_box=.25,
        move=.2,
        pointplot=True,
        box_showfliers=False,
        ax=ax
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel(ylab)
    ax.set_title("")

    # Base ylim (room for brackets)
    ymin = lower_bound
    ymax = max( dvar[var].max(), nominal_upper if nominal_upper is not None else dvar[var].max() )
    rng  = ymax - ymin
    step = 0.10 * rng if rng > 0 else 1.0  # vertical step between brackets
    head = 0.02 * rng if rng > 0 else 0.2  # headroom above highest bracket

    # compute bracket positions in ascending order
    y_positions = []
    start = ymax + 0.00 * rng  # start right above the data top
    for i in range(len(comparisons)):
        y_positions.append(start + i * step)
    ax.set_ylim(ymin, y_positions[-1] + head if len(y_positions) else ymax)

    # labels as significance stars
    labels = []
    for (g1, g2) in comparisons:
        row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
        if row.empty:
            labels.append("n.s.")
        else:
            labels.append(p_to_signif(float(row["p_adj"].iloc[0])))

    # draw brackets
    add_stat_brackets(
        ax=ax,
        xcats=list(dvar["Condition"].cat.categories),
        comparisons=comparisons,
        y_positions=y_positions,
        labels=labels
    )

    # tidy y-ticks to match your R aesthetics
    if var == "Accuracy":
        ax.set_yticks(np.arange(65, 101, 5))
    elif var == "ReactionTime":
        ax.set_yticks(np.round(np.arange(0.3, 1.31, 0.1), 2))

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"AOC_stats_rainclouds_{sname}_nback_stats.png"))
    plt.close(fig)