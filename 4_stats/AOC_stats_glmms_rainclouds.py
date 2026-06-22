# %% AOC Stats Rainclouds — Combined (Sternberg + N-Back)
# Loads merged CSVs, maps window-specific columns, applies IQR filtering, and saves raincloud figures.
# Significance brackets on all rainclouds: R pairwise CSV when available (confirmatory DVs),
# otherwise figure-only MixedLM contrasts in Python.
#
# Run order:
#   Rscript AOC_stats_glmm_nback.R
#   Rscript AOC_stats_glmm_sternberg.R
#   Rscript AOC_stats_glmm_export.R
#   python AOC_stats_glmms_rainclouds.py
#
# Key outputs:
#   Raincloud figures in figures/stats/rainclouds/

# %% Imports
import sys, os
import warnings
from typing import Optional
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from functions.stats_helpers import iqr_outlier_filter, p_to_signif, mixedlm_pairwise_contrasts
from functions.rainclouds_plotting_helpers import add_stat_brackets

warnings.filterwarnings("ignore", category=ConvergenceWarning, module="statsmodels")

# %% Task window sources (merged tables from AOC_master_matrix_*.m)
# n-back: full epoch [0 2]s for EEG alpha, ERSD, and baselined gaze (% change from FullBL).
# Sternberg: late epoch [1 2]s for EEG alpha, ERSD, and baselined gaze (LateBL).

RAW_ALPHA_BY_TASK = {
    "nback": "AlphaPower_raw_full",
    "sternberg": "AlphaPower_raw_late",
}

ALPHA_BL_BY_TASK = {
    "nback": "AlphaPower_bl_full",
    "sternberg": "AlphaPower_bl_late",
}

ERSD_BY_TASK = {
    "nback": "ERSD_full",
    "sternberg": "ERSD_late",
}

GAZE_BL_SOURCE_BY_TASK = {
    "nback": {
        "GazeDeviationBL": "GazeDeviationFullBL",
        "MSRateBL": "MSRateFullBL",
    },
    "sternberg": {
        "GazeDeviationBL": "GazeDeviationLateBL",
        "MSRateBL": "MSRateLateBL",
    },
}

# DVs with confirmatory R GLMMs; use R pairwise p-values for brackets when CSV is present
R_BRACKET_VARS = {"Accuracy", "ReactionTime", "GazeDeviation", "MSRate", "ERSD"}


def enforce_task_raw_alpha(df: pd.DataFrame, task_name: str) -> pd.DataFrame:
    """Ensure canonical `AlphaPower` follows the task-specific raw-window definition."""
    source_col = RAW_ALPHA_BY_TASK.get(task_name)
    if source_col is None:
        print(f"[{task_name}] WARNING: no raw-alpha mapping configured; using existing AlphaPower.")
        return df

    if source_col in df.columns:
        df["AlphaPower"] = pd.to_numeric(df[source_col], errors="coerce")
    else:
        print(f"[{task_name}] WARNING: {source_col} missing; using existing AlphaPower.")
    return df


def enforce_task_window_columns(df: pd.DataFrame, task_name: str) -> pd.DataFrame:
    """Map window-specific source columns into canonical names used in `variables` (in-place)."""
    df = enforce_task_raw_alpha(df, task_name)

    bl_src = ALPHA_BL_BY_TASK.get(task_name)
    if bl_src and bl_src in df.columns:
        df["AlphaPower_bl"] = pd.to_numeric(df[bl_src], errors="coerce")
    elif bl_src:
        print(f"[{task_name}] WARNING: {bl_src} missing; AlphaPower_bl not set.")

    ersd_src = ERSD_BY_TASK.get(task_name)
    if ersd_src and ersd_src in df.columns:
        df["ERSD"] = pd.to_numeric(df[ersd_src], errors="coerce")
    elif ersd_src:
        print(f"[{task_name}] WARNING: {ersd_src} missing; ERSD not set.")

    gaze_map = GAZE_BL_SOURCE_BY_TASK.get(task_name, {})
    for canonical, src in gaze_map.items():
        if src in df.columns:
            df[canonical] = pd.to_numeric(df[src], errors="coerce")
        else:
            print(f"[{task_name}] WARNING: {src} missing; {canonical} not set.")

    return df


def label_condition(dat: pd.DataFrame, task: dict) -> pd.DataFrame:
    """Map numeric or string Condition values to ordered categorical labels."""
    out = dat.copy()
    cond = out["Condition"]
    if np.issubdtype(cond.dtype, np.number):
        uniq = sorted(pd.unique(cond.dropna()).tolist())
        applied_map = None
        for candidate_map in task["cond_to_label_numeric"]:
            if set(uniq).issubset(set(candidate_map.keys())):
                applied_map = candidate_map
                break
        if applied_map is None:
            applied_map = {
                val: task["categories"][i]
                for i, val in enumerate(uniq[: len(task["categories"])])
            }
        out["Condition"] = out["Condition"].map(applied_map)
    else:
        out["Condition"] = (
            out["Condition"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
        )

    out["Condition"] = pd.Categorical(out["Condition"], categories=task["categories"], ordered=True)
    if out["ID"].dtype != "O":
        out["ID"] = out["ID"].astype(str)
    return out


def load_task_data(task: dict, variables: list) -> pd.DataFrame:
    """Load, map columns, label conditions, and apply IQR outlier filtering."""
    dat = pd.read_csv(task["input_csv"])
    dat = enforce_task_window_columns(dat, task["name"])
    dat.loc[dat["Accuracy"] > 100, "Accuracy"] = np.nan
    dat = label_condition(dat, task)
    vars_present = [v for v in variables if v in dat.columns]
    dat = iqr_outlier_filter(dat, vars_present, by="Condition")
    return dat


def load_pairwise_contrasts(stats_dir: str, task_name: str) -> Optional[pd.DataFrame]:
    """Load R-exported pairwise contrasts for bracket annotations."""
    path = os.path.join(stats_dir, f"AOC_pairwise_mixedlm_{task_name}.csv")
    if not os.path.isfile(path):
        print(f"WARNING: pairwise contrasts not found → {path}")
        print("         Figure brackets will use Python MixedLM contrasts instead.")
        return None
    return pd.read_csv(path)


def bracket_labels_for_var(
    pw: Optional[pd.DataFrame],
    var: str,
    comparisons: list,
) -> list[str]:
    """Map model-based pairwise p_adj values to significance bracket labels."""
    labels = []
    for g1, g2 in comparisons:
        if pw is None or pw.empty:
            labels.append("n.s.")
            continue
        sub = pw.loc[pw["Variable"] == var]
        row = sub.loc[
            ((sub["Group1"] == g1) & (sub["Group2"] == g2))
            | ((sub["Group1"] == g2) & (sub["Group2"] == g1))
        ]
        if row.empty:
            labels.append("n.s.")
        else:
            labels.append(p_to_signif(float(row["p_adj"].iloc[0])))
    return labels


def bracket_labels_from_python_mixedlm(
    dvar: pd.DataFrame,
    var: str,
    comparisons: list,
    task_name: str,
) -> list[str]:
    """Figure-only Load contrasts for raincloud brackets (not written to stats CSVs)."""
    df = dvar.rename(columns={var: "value"}).copy()
    df["Condition"] = pd.Categorical(df["Condition"], categories=df["Condition"].cat.categories, ordered=True)
    try:
        pw = mixedlm_pairwise_contrasts(
            df, value_col="value", group_col="Condition", id_col="ID", p_adjust="fdr_bh"
        )
    except Exception as exc:
        print(f"WARNING: figure-only MixedLM brackets failed for {var} [{task_name}]: {exc}")
        return ["n.s."] * len(comparisons)

    labels = []
    for g1, g2 in comparisons:
        row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
        labels.append("n.s." if row.empty else p_to_signif(float(row["p_adj"].iloc[0])))
    return labels


def bracket_labels_for_raincloud(
    dvar: pd.DataFrame,
    var: str,
    comparisons: list,
    pw: Optional[pd.DataFrame],
    task_name: str,
) -> list[str]:
    """R pairwise contrasts for confirmatory DVs when available; else Python MixedLM."""
    if pw is not None and not pw.empty and var in R_BRACKET_VARS:
        sub = pw.loc[pw["Variable"] == var]
        if not sub.empty:
            return bracket_labels_for_var(pw, var, comparisons)
    return bracket_labels_from_python_mixedlm(dvar, var, comparisons, task_name)


def plot_raincloud(
    dvar: pd.DataFrame,
    var: str,
    ylab: str,
    sname: str,
    task: dict,
    condition_order: list,
    pal_dict: dict,
    global_upper: dict,
    yticks_map: dict,
    ylims_map: dict,
    output_dir: str,
    pw: Optional[pd.DataFrame],
    figure_save_dpi: int,
) -> None:
    """Draw and save one raincloud figure."""
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    fig.patch.set_alpha(1.0)
    ax.patch.set_alpha(1.0)
    ax.set_facecolor("white")

    viol_alpha = 0.60
    dot_alpha = 0.50
    dot_size = 30
    box_width = 0.20
    cloud_offset = -0.20
    max_violsw = 0.40
    bw_method = 0.15
    if var == "Accuracy":
        cloud_offset = -0.20
        max_violsw = 0.50
        bw_method = 0.25

    xpos = {c: i for i, c in enumerate(condition_order)}
    rng = np.random.default_rng(12345)

    for cond_lab in condition_order:
        yvals = dvar.loc[dvar["Condition"] == cond_lab, var].dropna().to_numpy()
        if yvals.size == 0:
            continue

        ymax_cap = ylims_map[var][1] if var in ylims_map else float(dvar[var].max())
        ymin_cap = ylims_map[var][0] if var in ylims_map else float(dvar[var].min())
        yr = ymax_cap - ymin_cap

        pad_top = min(0.02 * yr, 0.3)
        y_grid_top = ymin_cap + yr + pad_top

        kde = gaussian_kde(yvals, bw_method=bw_method)
        y_grid = np.linspace(ymin_cap, y_grid_top, 400)

        dens = kde(y_grid)
        scale = (max_violsw / np.nanmax(dens)) if np.nanmax(dens) > 0 else 0.0

        x_left = xpos[cond_lab] + cloud_offset - dens * scale
        x_right = np.full_like(y_grid, xpos[cond_lab] + cloud_offset)

        poly_x = np.concatenate([x_right, x_left[::-1]])
        poly_y = np.concatenate([y_grid, y_grid[::-1]])
        ax.fill(
            poly_x, poly_y,
            facecolor=pal_dict[cond_lab], edgecolor="none",
            alpha=viol_alpha, clip_on=True,
        )

        x_jit = xpos[cond_lab] + rng.uniform(-box_width / 2, box_width / 2, size=yvals.size)
        ax.scatter(
            x_jit, yvals, s=dot_size, alpha=dot_alpha,
            color=pal_dict[cond_lab], linewidths=0, zorder=3,
        )

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

    leg = ax.get_legend()
    if leg is not None:
        leg.remove()

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.yaxis.grid(True, linewidth=1, alpha=0.35)
    ax.xaxis.grid(False)

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
            ha="center", va="top",
        )

    if var in ylims_map:
        ymin, ymax_cap = ylims_map[var]
    else:
        ymin = float(dvar[var].min()) if np.isfinite(dvar[var].min()) else 0.0
        ymax_data_local = float(dvar[var].max()) if np.isfinite(dvar[var].max()) else ymin
        ymax_cap = global_upper.get(var, np.nan)
        if not np.isfinite(ymax_cap):
            ymax_cap = ymax_data_local

    range_y = max(ymax_cap - ymin, 1.0)
    step = 0.10 * range_y
    y_positions = []
    start = ymax_cap + 0.075 * range_y
    for i in range(len(task["comparisons"])):
        y_positions.append(start + i * step)

    if var in ylims_map:
        ymin_cur, ymax_cur = ylims_map[var]
    else:
        ymin_cur = float(dvar[var].min()) if np.isfinite(dvar[var].min()) else 0.0
        ymax_cur = float(dvar[var].max()) if np.isfinite(dvar[var].max()) else 1.0
    ymid = 0.5 * (ymin_cur + ymax_cur)
    ax.set_ylabel("")
    ax.yaxis.get_label().set_visible(False)
    ax.text(
        YLABEL_GRID_X,
        ymid,
        ylab,
        transform=ax.get_yaxis_transform(which="grid"),
        rotation=90,
        ha="center",
        va="center",
    )

    labels = bracket_labels_for_raincloud(
        dvar, var, task["comparisons"], pw, task["name"]
    )
    if (task["name"] == "nback") and (var == "Accuracy"):
        yr_ax = ax.get_ylim()[1] - ax.get_ylim()[0]
        step_bump = 0.025 * yr_ax
        y_positions = [y + i * step_bump for i, y in enumerate(y_positions)]

    add_stat_brackets(
        ax=ax,
        xcats=condition_order,
        comparisons=task["comparisons"],
        y_positions=y_positions,
        labels=labels,
        xmap=xpos,
    )

    if var in yticks_map:
        ax.set_yticks(yticks_map[var])

    if var in ylims_map:
        ymin_set, ymax_set = ylims_map[var]
        extra = 0.12 * (ymax_set - ymin_set)
        ax.set_ylim(ymin_set, ymax_set + extra)
    else:
        ax.set_ylim(ymin, ymax_cap + 0.15 * (ymax_cap - ymin))

    fig.tight_layout()
    fig.subplots_adjust(left=0.17)
    out_path = os.path.join(output_dir, f"AOC_stats_rainclouds_{sname}_{task['name']}.png")
    fig.savefig(
        out_path,
        dpi=figure_save_dpi,
        transparent=False,
        facecolor=fig.get_facecolor(),
        edgecolor=fig.get_edgecolor() if hasattr(fig, "get_edgecolor") else "white",
    )
    plt.close(fig)
    print(f"Saved raincloud fig. → {os.path.basename(out_path)}")


# %% Parameters

pal = ["#93B8C4", "#82AD82", "#D998A2"]  # AOC pastels

FIGURE_SAVE_DPI = 600
mpl.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": FIGURE_SAVE_DPI,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "ps.fonttype": 42,
    "font.size": 25,
    "axes.titlesize": 25,
    "axes.labelsize": 20,
    "legend.fontsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "mathtext.default": "regular",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.edgecolor": "white",
    "axes.edgecolor": "white",
})

YLABEL_GRID_X = -0.125

sns.set_style("white")

# %% I/O Paths
base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
output_dir = f"{base_dir}/figures/stats/rainclouds"
stats_dir = f"{base_dir}/data/stats"
os.makedirs(output_dir, exist_ok=True)

MERGED_CSV_STERNBERG = f"{base_dir}/data/features/AOC_merged_data_sternberg.csv"
MERGED_CSV_NBACK = f"{base_dir}/data/features/AOC_merged_data_nback.csv"

# %% Variables and labelling
variables = [
    "Accuracy", "ReactionTime",
    "GazeDeviation", "MSRate",
    "GazeDeviationBL", "MSRateBL",
    "AlphaPower", "AlphaPower_bl",
    "IAF", "ERSD",
]
y_labels = [
    "Accuracy [%]", "Reaction Time [s]",
    "Gaze Deviation [px]", "Microsaccade Rate [MS/s]",
    "Gaze Deviation [%]", "Microsaccade Rate [%]",
    "Alpha Power [\u03BCV²/Hz]",
    "Alpha Power [dB]",
    "IAF [Hz]",
    "ERS/ERD [dB]",
]
save_names = [
    "acc", "rt",
    "gazedev", "ms",
    "gazedev_bl", "ms_bl",
    "pow_raw", "pow_bl",
    "iaf", "ersd",
]

yticks_map = {
    "Accuracy": np.arange(60, 110, 10),
    "ReactionTime": np.arange(0.25, 1.5, 0.25),
    "GazeDeviation": np.arange(0, 80, 10),
    "MSRate": np.arange(0, 4, 1),
    "AlphaPower": np.arange(0, 12, 1),
    "GazeDeviationBL": np.arange(-50, 350, 50),
    "MSRateBL": np.arange(-100, 125, 25),
    "AlphaPower_bl": np.arange(-4, 4.1, 1),
    "IAF": np.arange(8, 15, 1),
    "ERSD": np.arange(-3, 4, 1),
}
ylims_map = {
    "Accuracy": (58, 101),
    "ReactionTime": (0.25, 1.35),
    "GazeDeviation": (0, 70),
    "MSRate": (0, 3.75),
    "AlphaPower": (0, 11),
    "GazeDeviationBL": (-50, 300),
    "MSRateBL": (-110, 110),
    "AlphaPower_bl": (-4, 4),
    "IAF": (8, 14),
    "ERSD": (-3.75, 3.75),
}

# %% Task configurations
tasks = [
    {
        "name": "sternberg",
        "input_csv": MERGED_CSV_STERNBERG,
        "cond_to_label_numeric": [
            {1: "WM load 2", 2: "WM load 4", 3: "WM load 6"},
            {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"},
        ],
        "categories": ["WM load 2", "WM load 4", "WM load 6"],
        "comparisons": [
            ("WM load 2", "WM load 4"),
            ("WM load 2", "WM load 6"),
            ("WM load 4", "WM load 6"),
        ],
        "xlabel": "Condition",
    },
    {
        "name": "nback",
        "input_csv": MERGED_CSV_NBACK,
        "cond_to_label_numeric": [{1: "1-back", 2: "2-back", 3: "3-back"}],
        "categories": ["1-back", "2-back", "3-back"],
        "comparisons": [
            ("1-back", "2-back"),
            ("1-back", "3-back"),
            ("2-back", "3-back"),
        ],
        "xlabel": "Condition",
    },
]

# %% Pre-scan for global upper bounds (shared bracket baseline per variable)
global_upper = {var: np.nan for var in variables}

for _task in tasks:
    _dat = load_task_data(_task, variables)
    for var in variables:
        if var not in _dat.columns:
            continue
        vmax = pd.to_numeric(_dat[var], errors="coerce").max()
        if np.isfinite(vmax):
            if not np.isfinite(global_upper[var]) or (vmax > global_upper[var]):
                global_upper[var] = float(vmax)

# %% Raincloud plotting
for task in tasks:
    dat = load_task_data(task, variables)
    pw = load_pairwise_contrasts(stats_dir, task["name"])

    condition_order = list(dat["Condition"].dropna().unique())
    pal_dict = dict(zip(condition_order, pal))

    for var, ylab, sname in zip(variables, y_labels, save_names):
        if var not in dat.columns:
            continue

        dvar = dat.loc[~dat[var].isna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            continue

        dvar["Condition"] = pd.Categorical(dvar["Condition"], categories=condition_order, ordered=True)

        plot_raincloud(
            dvar=dvar,
            var=var,
            ylab=ylab,
            sname=sname,
            task=task,
            condition_order=condition_order,
            pal_dict=pal_dict,
            global_upper=global_upper,
            yticks_map=yticks_map,
            ylims_map=ylims_map,
            output_dir=output_dir,
            pw=pw,
            figure_save_dpi=FIGURE_SAVE_DPI,
        )

print("\nRaincloud figures saved to:", output_dir)
