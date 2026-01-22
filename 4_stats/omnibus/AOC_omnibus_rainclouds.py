# %% AOC Omnibus Rainclouds — Alpha Power
# Loads CSV from omnibus analysis, generates raincloud plots for Sternberg and N-back alpha power.
# Saves figures to figures/stats/omnibus.
#
# Key outputs:
#   Raincloud figures for alpha power by condition

# %% Imports
import sys, os
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy import stats
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

from functions.stats_helpers import (
    iqr_outlier_filter, p_to_signif
)

from functions.rainclouds_plotting_helpers import add_stat_brackets

# %% Parameters

# Colours (AOC pastels - matching AOC_stats_glmms_rainclouds.py)
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
base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
input_csv = f"{base_dir}/data/features/AOC_omnibus_raincloud_data.csv"
output_dir = f"{base_dir}/figures/stats/omnibus"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# %% Load data
print("Loading raincloud data...")
dat = pd.read_csv(input_csv)

# Ensure ID is string
if dat["ID"].dtype != "O":
    dat["ID"] = dat["ID"].astype(str)

# Ensure Condition is categorical with proper ordering
dat["Condition"] = pd.Categorical(dat["Condition"], ordered=True)

print(f"Loaded {len(dat)} rows")
print(f"Tasks: {dat['Task'].unique()}")
print(f"Conditions: {dat['Condition'].unique()}")

# %% Task configurations
tasks = [
    {
        "name": "sternberg",
        "categories": ["WM load 2", "WM load 4", "WM load 6"],
        "comparisons": [("WM load 2", "WM load 4"), ("WM load 2", "WM load 6"), ("WM load 4", "WM load 6")],
        "xlabel": "Condition",
        "palette": pal,
        "title": "Sternberg"
    },
    {
        "name": "nback",
        "categories": ["1-back", "2-back", "3-back"],
        "comparisons": [("1-back", "2-back"), ("1-back", "3-back"), ("2-back", "3-back")],
        "xlabel": "Condition",
        "palette": pal,
        "title": "N-back"
    }
]

# Variable configuration
var = "AlphaPower"
ttl = "Alpha Power"
ylab = "Alpha power [change from bl]"
sname = "alpha"

# Manual y ticks and ylims
yticks_alpha = np.arange(-0.25, 0.3, 0.05)  # Include 0.25 (arange is exclusive of end, so use 0.3)
ylims_alpha = (-0.3, 0.3)

# %% Calculate global data maximum across all tasks (for consistent bracket positioning)
global_data_max = float(dat[var].max()) if not dat[var].isna().all() else ylims_alpha[1]
print(f"Global data maximum: {global_data_max:.4f}")

# %% Processing + plotting
for task in tasks:
    print(f"\n=== Processing {task['name']} ===")
    
    # Filter data for this task
    dvar = dat.loc[(dat["Task"] == task["name"]) & (~dat[var].isna()), ["ID", "Condition", var]].copy()
    if dvar.empty:
        print(f"  No data for {task['name']}. Skipping.")
        continue
    
    # Ensure categorical ordering
    dvar["Condition"] = pd.Categorical(dvar["Condition"], categories=task["categories"], ordered=True)
    condition_order = list(dvar["Condition"].dropna().unique())
    
    # Category order / palette mapping
    pal_dict = dict(zip(condition_order, task["palette"]))
    
    # %% Repeated-measures ANOVA (within-subject: Condition)
    # balance subjects: keep only subjects present in all conditions for this var
    present = dvar.groupby('ID')['Condition'].nunique()
    keep_ids = present[present == len(condition_order)].index
    dvar_bal = dvar[dvar['ID'].isin(keep_ids)].copy()
    
    # run ANOVA on the balanced panel
    if dvar_bal['ID'].nunique() >= 2:
        aov = AnovaRM(data=dvar_bal, depvar=var, subject='ID', within=['Condition']).fit()
        aov_tab = aov.anova_table.reset_index().rename(columns={'index': 'Effect'})
        print(f"  ANOVA: F = {aov_tab.loc[0, 'F Value']:.3f}, p = {aov_tab.loc[0, 'Pr > F']:.4f}")
    else:
        print(f"  Warning: Not enough balanced subjects for ANOVA")
    
    # Paired t-tests with FDR correction
    # Build wide table for paired comparisons
    wide = dvar.pivot(index='ID', columns='Condition', values=var)
    
    # Perform paired t-tests for each comparison
    pw_rows = []
    p_values = []
    for (g1, g2) in task["comparisons"]:
        if (g1 in wide.columns) and (g2 in wide.columns):
            # Get paired data (drop NaNs in either condition)
            paired_data = wide[[g1, g2]].dropna()
            if len(paired_data) >= 2:
                # Paired t-test (like MATLAB's ttest)
                t_stat, p_val = stats.ttest_rel(paired_data[g2], paired_data[g1])
                pw_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "p": p_val
                })
                p_values.append(p_val)
            else:
                pw_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "p": np.nan
                })
                p_values.append(np.nan)
    
    # Apply FDR correction (Benjamini-Hochberg) - less conservative than Bonferroni
    if len(p_values) > 0 and not all(np.isnan(p_values)):
        # Filter out NaNs for correction
        valid_mask = ~np.isnan(p_values)
        if np.sum(valid_mask) > 0:
            p_valid = np.array(p_values)[valid_mask]
            _, p_adj_valid, _, _ = multipletests(p_valid, method='fdr_bh')
            # Map back to original positions
            p_adj_all = np.full(len(p_values), np.nan)
            p_adj_all[valid_mask] = p_adj_valid
        else:
            p_adj_all = np.array(p_values)
    else:
        p_adj_all = np.array(p_values)
    
    # Add adjusted p-values to dataframe
    for i, row in enumerate(pw_rows):
        row["p_adj"] = p_adj_all[i] if i < len(p_adj_all) else np.nan
    
    pw = pd.DataFrame(pw_rows)
    
    # Print pairwise results
    print("  Pairwise comparisons (paired t-tests with FDR correction):")
    for _, row in pw.iterrows():
        print(f"    {row['group1']} vs {row['group2']}: p = {row['p']:.4f}, p_adj = {row['p_adj']:.4f}")
    
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
        ymin_cap, ymax_cap = ylims_alpha
        yr = ymax_cap - ymin_cap
        
        pad_top = min(0.02 * yr, 0.3)
        y_grid_top = ymin_cap + yr + pad_top
        
        # build KDE
        kde = gaussian_kde(yvals, bw_method=bw_method)
        
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
    
    # %% Bracket layout (closer to data points with moderate spacing)
    ymin, ymax_cap = ylims_alpha
    range_y = max(ymax_cap - ymin, 1.0)
    
    # Use global data max for consistent bracket positioning across tasks
    head = 0.02 * range_y
    step = 0.07 * range_y  # Moderate spacing between brackets
    
    y_positions = []
    start = 0.26  # Use global max for consistent positioning
    for i in range(len(task["comparisons"])):
        y_positions.append(start + i * step)
    
    # y-label at data midpoint
    ymid = 0.5 * (ymin + ymax_cap)
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
    
    add_stat_brackets(
        ax=ax,
        xcats=condition_order,
        comparisons=task["comparisons"],
        y_positions=y_positions,
        labels=labels,
        xmap=xpos
    )
    
    # %% Manual y-ticks and ylims
    # Ensure y=0 is always included in ticks
    yticks_final = list(yticks_alpha)
    if 0.0 not in yticks_final:
        yticks_final.append(0.0)
        yticks_final = sorted(yticks_final)
    ax.set_yticks(yticks_final)
    ymin_set, ymax_set = ylims_alpha
    # extend a bit to make sure brackets are inside
    extra = 0.08 * (ymax_set - ymin_set)
    ax.set_ylim(ymin_set, ymax_set + extra)
    
    # No title (removed as requested)
    
    # %% Save raincloud figure
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, f"AOC_omnibus_raincloud_{sname}_{task['name']}.png"),
        dpi=300,
        transparent=False,
        facecolor=fig.get_facecolor(),
        edgecolor=fig.get_edgecolor() if hasattr(fig, "get_edgecolor") else "white"
    )
    plt.close(fig)
    saveNameFig = f"AOC_omnibus_raincloud_{sname}_{task['name']}.png"
    print(f"  Saved raincloud fig. → {saveNameFig}")

print("\n=== Done ===")
