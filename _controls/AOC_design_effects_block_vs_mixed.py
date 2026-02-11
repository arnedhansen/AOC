# %% AOC Design Effects — Block vs Mixed Design Influence on Alpha Power
# Analyses whether the different experimental designs (block design in N-back
# vs. random/mixed design in Sternberg) introduce systematic order or carry-over
# effects on alpha power.
#
# N-back: 6 blocks × 75 trials, each block = one condition (1-/2-/3-back).
#         Block order randomised per subject → "starting condition" varies.
# Sternberg: 6 blocks × 50 trials, WM load (2/4/6) randomised per trial
#            within each block → no block-level condition confound.
#
# MAIN QUESTION: Does it matter which N-back condition comes first?
#   → Do participants who start with 1-back, 2-back, or 3-back show
#     different alpha power patterns across conditions?
#
# Additional questions:
#   1. N-back: Does alpha power drift across block positions (1–6)?
#   2. N-back: Does the 1st vs 2nd occurrence of the same condition differ?
#   3. Sternberg: Does the preceding trial's condition influence alpha power?
#   4. Sternberg: Do condition switches vs repeats affect alpha?
#   5. Cross-task: How stable is alpha across the experiment in each design?
#
# Key outputs:
#   Figures (PNG); statistical summary table (CSV); console output

# %% Imports
import sys, os
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sp_stats
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# %% Plot settings (match project style)
pal = ["#93B8C4", "#82AD82", "#D998A2"]  # AOC pastels: blue, green, pink
pal_start = {"1-back": "#93B8C4", "2-back": "#82AD82", "3-back": "#D998A2"}

mpl.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 300,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "ps.fonttype": 42,
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "mathtext.default": "regular",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})
sns.set_style("white")

# %% I/O paths
base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
feature_dir = f"{base_dir}/data/features"
fig_dir = f"{base_dir}/data/controls/block-design-effect"
stats_dir = f"{base_dir}/data/controls/block-design-effect"
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(stats_dir, exist_ok=True)

# %% Configuration
NBACK_TRIALS_PER_BLOCK = 75
NBACK_N_BLOCKS = 6
STERN_TRIALS_PER_BLOCK = 50
STERN_N_BLOCKS = 6

# Alpha power columns (determined by data availability)
# N-back trials:    only AlphaPowerEarly has data (Late & Full are NaN)
# N-back subjects:  AlphaPower (aggregated)
# Sternberg trials: AlphaPowerLate has data
ALPHA_COL_NB_TRIALS = "AlphaPowerEarly"  # 0–1 s post-stimulus
ALPHA_COL_ST_TRIALS = "AlphaPowerLate"   # 1–2 s retention interval

# Condition labels
NBACK_COND_MAP = {1: "1-back", 2: "2-back", 3: "3-back"}
STERN_COND_MAP = {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"}

# %% ──────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ──────────────────────────────────────────────────────────────────────────────

# %% Load N-back subject-level data (has AlphaPower)
nb_subj_path = os.path.join(feature_dir, "merged_data_nback.csv")
nb_subj_raw = pd.read_csv(nb_subj_path)
nb_subj_raw["ID"] = nb_subj_raw["ID"].astype(str)
nb_subj_raw["Condition"] = nb_subj_raw["Condition"].astype(int)

# %% Load N-back trial-level data
nback_path = os.path.join(feature_dir, "merged_data_nback_trials.csv")
nb = pd.read_csv(nback_path)
nb["ID"] = nb["ID"].astype(str)
nb = nb.dropna(subset=[ALPHA_COL_NB_TRIALS, "Trial", "Condition"])
nb["Condition"] = nb["Condition"].astype(int)
nb["Trial"] = nb["Trial"].astype(int)

# Derive block position (1–6) from trial number
nb["Block"] = np.ceil(nb["Trial"] / NBACK_TRIALS_PER_BLOCK).astype(int)
nb["TrialInBlock"] = ((nb["Trial"] - 1) % NBACK_TRIALS_PER_BLOCK) + 1

# %% Load Sternberg trial-level data
stern_path = os.path.join(feature_dir, "merged_data_sternberg_trials.csv")
st = pd.read_csv(stern_path)
st["ID"] = st["ID"].astype(str)
st = st.dropna(subset=[ALPHA_COL_ST_TRIALS, "Trial", "Condition"])
st["Condition"] = st["Condition"].astype(int)
st["Trial"] = st["Trial"].astype(int)

# Derive block position
st["Block"] = np.ceil(st["Trial"] / STERN_TRIALS_PER_BLOCK).astype(int)
st["TrialInBlock"] = ((st["Trial"] - 1) % STERN_TRIALS_PER_BLOCK) + 1

print(f"N-back (trials):    {nb['ID'].nunique()} subjects, {len(nb)} trials")
print(f"N-back (subjects):  {nb_subj_raw['ID'].nunique()} subjects")
print(f"Sternberg (trials): {st['ID'].nunique()} subjects, {len(st)} trials")

# %% ──────────────────────────────────────────────────────────────────────────
# N-BACK: RECONSTRUCT BLOCK ORDER AND STARTING CONDITION
# ──────────────────────────────────────────────────────────────────────────────

# %% From trial data: all trials in a block share the same condition
nback_block_order = (
    nb.groupby(["ID", "Block"])["Condition"]
    .first()
    .reset_index()
    .rename(columns={"Condition": "BlockCondition"})
    .sort_values(["ID", "Block"])
)

# Starting condition = condition of block 1
start_cond = (
    nback_block_order[nback_block_order["Block"] == 1][["ID", "BlockCondition"]]
    .rename(columns={"BlockCondition": "StartCond"})
)

# Merge into trial data
nb = nb.merge(start_cond, on="ID", how="left")
nb["StartCondLabel"] = nb["StartCond"].map(NBACK_COND_MAP)

# Merge into subject-level data
nb_subj_raw = nb_subj_raw.merge(start_cond, on="ID", how="left")
nb_subj_raw["StartCondLabel"] = nb_subj_raw["StartCond"].map(NBACK_COND_MAP)
nb_subj_raw["CondLabel"] = nb_subj_raw["Condition"].map(NBACK_COND_MAP)
nb_subj_raw["CondLabel"] = pd.Categorical(
    nb_subj_raw["CondLabel"], categories=["1-back", "2-back", "3-back"], ordered=True
)

# Condition repetition number (1st or 2nd occurrence)
nback_block_order["Repetition"] = (
    nback_block_order.groupby(["ID", "BlockCondition"]).cumcount() + 1
)
nb = nb.merge(
    nback_block_order[["ID", "Block", "Repetition"]],
    on=["ID", "Block"], how="left"
)

# Condition labels for trial data
nb["CondLabel"] = nb["Condition"].map(NBACK_COND_MAP)
nb["CondLabel"] = pd.Categorical(
    nb["CondLabel"], categories=["1-back", "2-back", "3-back"], ordered=True
)

print(f"\nN-back starting condition distribution:")
print(start_cond["StartCond"].value_counts().sort_index().to_string())

# %% ──────────────────────────────────────────────────────────────────────────
# STERNBERG: DERIVE SEQUENTIAL TRIAL FEATURES
# ──────────────────────────────────────────────────────────────────────────────

st = st.sort_values(["ID", "Trial"])
st["PrevCondition"] = st.groupby("ID")["Condition"].shift(1)
st["ConditionSwitch"] = (st["Condition"] != st["PrevCondition"]).astype(float)
st.loc[st["PrevCondition"].isna(), "ConditionSwitch"] = np.nan
st["PrevCondition"] = st["PrevCondition"].astype("Int64")

st["CondLabel"] = st["Condition"].map(STERN_COND_MAP)
st["CondLabel"] = pd.Categorical(
    st["CondLabel"], categories=["WM load 2", "WM load 4", "WM load 6"], ordered=True
)
st["PrevCondLabel"] = st["PrevCondition"].map(STERN_COND_MAP)

# %% Helper
def subj_mean_sem(df, alpha_col, groupby_cols):
    """Aggregate to subject means, then compute group mean ± SEM."""
    subj = df.groupby(["ID"] + groupby_cols)[alpha_col].mean().reset_index()
    out = subj.groupby(groupby_cols)[alpha_col].agg(["mean", "sem", "count"]).reset_index()
    out.columns = groupby_cols + ["mean", "sem", "n"]
    return out

# %% ══════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════
#                    MAIN QUESTION: STARTING CONDITION EFFECT
# ══════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 1 — N-BACK: ALPHA BY CONDITION × STARTING CONDITION (interaction plot)
# ──────────────────────────────────────────────────────────────────────────────

# %% Use subject-level data (AlphaPower) for cleanest comparison
nb_subj_clean = nb_subj_raw.dropna(subset=["AlphaPower", "StartCond"]).copy()

fig1, ax1 = plt.subplots(figsize=(9, 6), facecolor="white")

start_labels = ["1-back", "2-back", "3-back"]
cond_positions = np.arange(3)
offsets = [-0.2, 0.0, 0.2]
markers = ["o", "s", "D"]

for si, start_lab in enumerate(start_labels):
    sub = nb_subj_clean[nb_subj_clean["StartCondLabel"] == start_lab]
    if sub.empty:
        continue

    means = sub.groupby("CondLabel")["AlphaPower"].mean()
    sems = sub.groupby("CondLabel")["AlphaPower"].sem()
    n = sub["ID"].nunique()

    ax1.errorbar(
        cond_positions + offsets[si], means.values, yerr=sems.values,
        color=pal[si], marker=markers[si], markersize=9,
        capsize=4, linewidth=2, label=f"Start: {start_lab} (n={n})"
    )

ax1.set_xticks(cond_positions)
ax1.set_xticklabels(["1-back", "2-back", "3-back"])
ax1.set_xlabel("Condition")
ax1.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax1.set_title("N-back — Does the starting condition affect alpha?")
ax1.legend(frameon=False, loc="best")
ax1.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax1.spines.values():
    spine.set_visible(False)
fig1.tight_layout()
fig1.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_startcond_interaction.png"))
plt.close(fig1)
print("Saved: Fig 1 — N-back alpha × starting condition interaction")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 2 — N-BACK: ALPHA BY CONDITION, FACETED BY STARTING CONDITION
# ──────────────────────────────────────────────────────────────────────────────

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5), sharey=True, facecolor="white")

for idx, (start_val, start_lab) in enumerate(NBACK_COND_MAP.items()):
    ax = axes2[idx]
    sub = nb_subj_clean[nb_subj_clean["StartCond"] == start_val]

    if sub.empty:
        ax.set_title(f"Start: {start_lab}\n(n=0)")
        continue

    n_subj = sub["ID"].nunique()

    sns.boxplot(
        data=sub, x="CondLabel", y="AlphaPower",
        palette=pal, width=0.5, fliersize=0,
        boxprops=dict(alpha=0.4), ax=ax
    )
    sns.stripplot(
        data=sub, x="CondLabel", y="AlphaPower",
        palette=pal, size=5, alpha=0.6, jitter=0.15, ax=ax
    )
    ax.set_title(f"Start: {start_lab}  (n={n_subj})")
    ax.set_xlabel("Condition")
    if idx == 0:
        ax.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
    else:
        ax.set_ylabel("")
    ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(False)

fig2.suptitle("N-back — Alpha power by condition, split by starting condition",
              fontsize=16, y=1.02)
fig2.tight_layout()
fig2.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_by_startcond.png"))
plt.close(fig2)
print("Saved: Fig 2 — N-back alpha faceted by starting condition")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 3 — N-BACK: ALPHA CHANGE (3-back minus 1-back) BY STARTING CONDITION
# ──────────────────────────────────────────────────────────────────────────────

# %% Pivot to wide: one row per subject, columns = conditions
nb_wide = nb_subj_clean.pivot(index="ID", columns="CondLabel", values="AlphaPower")
nb_wide = nb_wide.merge(start_cond, on="ID", how="left")
nb_wide["StartCondLabel"] = nb_wide["StartCond"].map(NBACK_COND_MAP)
nb_wide["AlphaChange"] = nb_wide["3-back"] - nb_wide["1-back"]

fig3, ax3 = plt.subplots(figsize=(8, 5), facecolor="white")

for si, start_lab in enumerate(start_labels):
    sub = nb_wide[nb_wide["StartCondLabel"] == start_lab].dropna(subset=["AlphaChange"])
    if sub.empty:
        continue
    vals = sub["AlphaChange"].values
    n = len(vals)
    m = np.mean(vals)
    se = sp_stats.sem(vals)

    # Individual dots
    jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=n)
    ax3.scatter(np.full(n, si) + jitter, vals,
                color=pal[si], alpha=0.45, s=35, zorder=3)
    # Mean ± SEM
    ax3.errorbar(si, m, yerr=se, color=pal[si], marker="D", markersize=10,
                  capsize=5, linewidth=2.5, zorder=4, markeredgecolor="black",
                  markeredgewidth=0.5)

ax3.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
ax3.set_xticks(range(3))
ax3.set_xticklabels([f"Start: {s}" for s in start_labels])
ax3.set_xlabel("Starting condition")
ax3.set_ylabel(r"$\Delta$ Alpha power (3-back $-$ 1-back) [$\mu$V²/Hz]")
ax3.set_title("N-back — Alpha change by starting condition")
ax3.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax3.spines.values():
    spine.set_visible(False)
fig3.tight_layout()
fig3.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_change_by_startcond.png"))
plt.close(fig3)
print("Saved: Fig 3 — N-back alpha change by starting condition")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 4 — N-BACK: ALPHA TIME COURSE ACROSS BLOCKS
# ──────────────────────────────────────────────────────────────────────────────

fig4, ax4 = plt.subplots(figsize=(14, 5), facecolor="white")

BIN_SIZE_NB = 5
nb["TrialBin"] = ((nb["Trial"] - 1) // BIN_SIZE_NB) * BIN_SIZE_NB + BIN_SIZE_NB / 2

gm_nb = nb.groupby(["TrialBin", "CondLabel"])[ALPHA_COL_NB_TRIALS].agg(["mean", "sem"]).reset_index()
gm_nb.columns = ["TrialBin", "CondLabel", "mean", "sem"]

for i, cond in enumerate(["1-back", "2-back", "3-back"]):
    sub = gm_nb[gm_nb["CondLabel"] == cond]
    ax4.plot(sub["TrialBin"], sub["mean"], color=pal[i], linewidth=1.5, label=cond)
    ax4.fill_between(sub["TrialBin"],
                      sub["mean"] - sub["sem"],
                      sub["mean"] + sub["sem"],
                      color=pal[i], alpha=0.2)

for b in range(1, NBACK_N_BLOCKS):
    ax4.axvline(b * NBACK_TRIALS_PER_BLOCK + 0.5, color="grey", linestyle="--",
                linewidth=0.8, alpha=0.6)

ax4.set_xlabel("Trial number")
ax4.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax4.set_title("N-back — Alpha power time course (block design)")
ax4.legend(frameon=False)
ax4.yaxis.grid(True, linewidth=0.5, alpha=0.3)
fig4.tight_layout()
fig4.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_timecourse.png"))
plt.close(fig4)
print("Saved: Fig 4 — N-back alpha time course")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 5 — N-BACK: ALPHA BY BLOCK POSITION
# ──────────────────────────────────────────────────────────────────────────────

fig5, ax5 = plt.subplots(figsize=(10, 5), facecolor="white")

nb_block = subj_mean_sem(nb, ALPHA_COL_NB_TRIALS, ["Block", "CondLabel"])

for i, cond in enumerate(["1-back", "2-back", "3-back"]):
    sub = nb_block[nb_block["CondLabel"] == cond]
    ax5.errorbar(sub["Block"], sub["mean"], yerr=sub["sem"],
                  color=pal[i], marker="o", markersize=6,
                  capsize=3, linewidth=1.5, label=cond)

ax5.set_xlabel("Block position (1–6)")
ax5.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax5.set_title("N-back — Alpha power by block position")
ax5.set_xticks(range(1, 7))
ax5.legend(frameon=False)
ax5.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax5.spines.values():
    spine.set_visible(False)
fig5.tight_layout()
fig5.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_blockposition.png"))
plt.close(fig5)
print("Saved: Fig 5 — N-back alpha by block position")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 6 — N-BACK: 1ST VS 2ND REPETITION
# ──────────────────────────────────────────────────────────────────────────────

nb_rep = nb.groupby(["ID", "CondLabel", "Repetition"])[ALPHA_COL_NB_TRIALS].mean().reset_index()

fig6, ax6 = plt.subplots(figsize=(9, 5), facecolor="white")

rep_colors = ["#6b9dba", "#c4785a"]

for i, cond in enumerate(["1-back", "2-back", "3-back"]):
    for rep in [1, 2]:
        sub = nb_rep[(nb_rep["CondLabel"] == cond) & (nb_rep["Repetition"] == rep)]
        jit = np.random.default_rng(42 + rep).uniform(-0.08, 0.08, size=len(sub))
        ax6.scatter(
            np.full(len(sub), i) + (rep - 1.5) * 0.12 + jit,
            sub[ALPHA_COL_NB_TRIALS],
            color=rep_colors[rep - 1], alpha=0.35, s=20, zorder=3
        )
        m = sub[ALPHA_COL_NB_TRIALS].mean()
        se = sub[ALPHA_COL_NB_TRIALS].sem()
        ax6.errorbar(i + (rep - 1.5) * 0.2, m, yerr=se,
                      color=rep_colors[rep - 1], marker="D", markersize=8,
                      capsize=4, linewidth=2, zorder=4,
                      label=f"Rep {rep}" if i == 0 else None)

ax6.set_xticks(range(3))
ax6.set_xticklabels(["1-back", "2-back", "3-back"])
ax6.set_xlabel("Condition")
ax6.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax6.set_title("N-back — 1st vs 2nd block occurrence per condition")
ax6.legend(frameon=False)
ax6.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax6.spines.values():
    spine.set_visible(False)
fig6.tight_layout()
fig6.savefig(os.path.join(fig_dir, "AOC_design_nback_alpha_repetition.png"))
plt.close(fig6)
print("Saved: Fig 6 — N-back alpha 1st vs 2nd repetition")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 7 — STERNBERG: ALPHA TIME COURSE
# ──────────────────────────────────────────────────────────────────────────────

fig7, ax7 = plt.subplots(figsize=(14, 5), facecolor="white")

BIN_SIZE_ST = 5
st["TrialBin"] = ((st["Trial"] - 1) // BIN_SIZE_ST) * BIN_SIZE_ST + BIN_SIZE_ST / 2

gm_st = st.groupby(["TrialBin", "CondLabel"])[ALPHA_COL_ST_TRIALS].agg(["mean", "sem"]).reset_index()
gm_st.columns = ["TrialBin", "CondLabel", "mean", "sem"]

for i, cond in enumerate(["WM load 2", "WM load 4", "WM load 6"]):
    sub = gm_st[gm_st["CondLabel"] == cond]
    ax7.plot(sub["TrialBin"], sub["mean"], color=pal[i], linewidth=1.5, label=cond)
    ax7.fill_between(sub["TrialBin"],
                      sub["mean"] - sub["sem"],
                      sub["mean"] + sub["sem"],
                      color=pal[i], alpha=0.2)

for b in range(1, STERN_N_BLOCKS):
    ax7.axvline(b * STERN_TRIALS_PER_BLOCK + 0.5, color="grey", linestyle="--",
                linewidth=0.8, alpha=0.6)

ax7.set_xlabel("Trial number")
ax7.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax7.set_title("Sternberg — Alpha power time course (mixed/random design)")
ax7.legend(frameon=False)
ax7.yaxis.grid(True, linewidth=0.5, alpha=0.3)
fig7.tight_layout()
fig7.savefig(os.path.join(fig_dir, "AOC_design_stern_alpha_timecourse.png"))
plt.close(fig7)
print("Saved: Fig 7 — Sternberg alpha time course")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 8 — STERNBERG: PREVIOUS TRIAL CONDITION EFFECT
# ──────────────────────────────────────────────────────────────────────────────

st_seq = st.dropna(subset=["PrevCondition"]).copy()
st_seq["PrevCondLabel"] = pd.Categorical(
    st_seq["PrevCondition"].map(STERN_COND_MAP),
    categories=["WM load 2", "WM load 4", "WM load 6"], ordered=True
)

st_seq_subj = (
    st_seq.groupby(["ID", "CondLabel", "PrevCondLabel"])[ALPHA_COL_ST_TRIALS]
    .mean().reset_index()
)

fig8, axes8 = plt.subplots(1, 3, figsize=(16, 5), sharey=True, facecolor="white")

for idx, cond in enumerate(["WM load 2", "WM load 4", "WM load 6"]):
    ax = axes8[idx]
    sub = st_seq_subj[st_seq_subj["CondLabel"] == cond]

    sns.boxplot(
        data=sub, x="PrevCondLabel", y=ALPHA_COL_ST_TRIALS,
        palette=pal, width=0.5, fliersize=0,
        boxprops=dict(alpha=0.4), ax=ax
    )
    sns.stripplot(
        data=sub, x="PrevCondLabel", y=ALPHA_COL_ST_TRIALS,
        palette=pal, size=4, alpha=0.5, jitter=0.15, ax=ax
    )
    ax.set_title(f"Current: {cond}")
    ax.set_xlabel("Previous trial condition")
    if idx == 0:
        ax.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
    else:
        ax.set_ylabel("")
    ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(False)

fig8.suptitle("Sternberg — Alpha power by previous trial condition", fontsize=16, y=1.02)
fig8.tight_layout()
fig8.savefig(os.path.join(fig_dir, "AOC_design_stern_alpha_prev_condition.png"))
plt.close(fig8)
print("Saved: Fig 8 — Sternberg alpha by previous condition")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 9 — STERNBERG: SWITCH VS REPEAT
# ──────────────────────────────────────────────────────────────────────────────

st_switch = st.dropna(subset=["ConditionSwitch"]).copy()
st_switch["SwitchLabel"] = st_switch["ConditionSwitch"].map({0.0: "Repeat", 1.0: "Switch"})

st_sw_subj = (
    st_switch.groupby(["ID", "CondLabel", "SwitchLabel"])[ALPHA_COL_ST_TRIALS]
    .mean().reset_index()
)

fig9, ax9 = plt.subplots(figsize=(9, 5), facecolor="white")

sns.boxplot(
    data=st_sw_subj, x="CondLabel", y=ALPHA_COL_ST_TRIALS, hue="SwitchLabel",
    palette={"Repeat": "#82AD82", "Switch": "#D998A2"},
    width=0.6, fliersize=0, boxprops=dict(alpha=0.4), ax=ax9
)
sns.stripplot(
    data=st_sw_subj, x="CondLabel", y=ALPHA_COL_ST_TRIALS, hue="SwitchLabel",
    palette={"Repeat": "#82AD82", "Switch": "#D998A2"},
    size=4, alpha=0.5, jitter=0.1, dodge=True, ax=ax9
)
handles, labels = ax9.get_legend_handles_labels()
n_unique = len(set(labels))
ax9.legend(handles[:n_unique], labels[:n_unique], frameon=False)

ax9.set_xlabel("Condition")
ax9.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax9.set_title("Sternberg — Alpha power: condition switch vs repeat")
ax9.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax9.spines.values():
    spine.set_visible(False)
fig9.tight_layout()
fig9.savefig(os.path.join(fig_dir, "AOC_design_stern_alpha_switch_repeat.png"))
plt.close(fig9)
print("Saved: Fig 9 — Sternberg alpha switch vs repeat")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 10 — CROSS-TASK: ALPHA STABILITY COMPARISON
# ──────────────────────────────────────────────────────────────────────────────

fig10, (ax10a, ax10b) = plt.subplots(1, 2, figsize=(16, 5), facecolor="white")

# N-back: grand average per trial bin (collapsed across conditions)
gm_nb_all = nb.groupby("TrialBin")[ALPHA_COL_NB_TRIALS].agg(["mean", "sem"]).reset_index()
ax10a.plot(gm_nb_all["TrialBin"], gm_nb_all["mean"], color="#555555", linewidth=1.5)
ax10a.fill_between(gm_nb_all["TrialBin"],
                    gm_nb_all["mean"] - gm_nb_all["sem"],
                    gm_nb_all["mean"] + gm_nb_all["sem"],
                    color="#555555", alpha=0.2)
for b in range(1, NBACK_N_BLOCKS):
    ax10a.axvline(b * NBACK_TRIALS_PER_BLOCK + 0.5, color="grey",
                  linestyle="--", linewidth=0.6, alpha=0.5)
ax10a.set_xlabel("Trial number")
ax10a.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax10a.set_title("N-back (block design)")
ax10a.yaxis.grid(True, linewidth=0.5, alpha=0.3)

# Sternberg: grand average per trial bin
gm_st_all = st.groupby("TrialBin")[ALPHA_COL_ST_TRIALS].agg(["mean", "sem"]).reset_index()
ax10b.plot(gm_st_all["TrialBin"], gm_st_all["mean"], color="#555555", linewidth=1.5)
ax10b.fill_between(gm_st_all["TrialBin"],
                    gm_st_all["mean"] - gm_st_all["sem"],
                    gm_st_all["mean"] + gm_st_all["sem"],
                    color="#555555", alpha=0.2)
for b in range(1, STERN_N_BLOCKS):
    ax10b.axvline(b * STERN_TRIALS_PER_BLOCK + 0.5, color="grey",
                  linestyle="--", linewidth=0.6, alpha=0.5)
ax10b.set_xlabel("Trial number")
ax10b.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax10b.set_title("Sternberg (mixed/random design)")
ax10b.yaxis.grid(True, linewidth=0.5, alpha=0.3)

fig10.suptitle("Cross-task — Alpha power stability across the experiment",
               fontsize=16, y=1.02)
fig10.tight_layout()
fig10.savefig(os.path.join(fig_dir, "AOC_design_crosstask_alpha_stability.png"))
plt.close(fig10)
print("Saved: Fig 10 — Cross-task alpha stability")

# %% ──────────────────────────────────────────────────────────────────────────
# FIGURE 11 — STERNBERG: ALPHA BY BLOCK POSITION
# ──────────────────────────────────────────────────────────────────────────────

fig11, ax11 = plt.subplots(figsize=(10, 5), facecolor="white")

st_block = subj_mean_sem(st, ALPHA_COL_ST_TRIALS, ["Block", "CondLabel"])

for i, cond in enumerate(["WM load 2", "WM load 4", "WM load 6"]):
    sub = st_block[st_block["CondLabel"] == cond]
    ax11.errorbar(sub["Block"], sub["mean"], yerr=sub["sem"],
                  color=pal[i], marker="o", markersize=6,
                  capsize=3, linewidth=1.5, label=cond)

ax11.set_xlabel("Block position (1–6)")
ax11.set_ylabel(r"Alpha power [$\mu$V²/Hz]")
ax11.set_title("Sternberg — Alpha power by block position")
ax11.set_xticks(range(1, 7))
ax11.legend(frameon=False)
ax11.yaxis.grid(True, linewidth=0.5, alpha=0.3)
for spine in ax11.spines.values():
    spine.set_visible(False)
fig11.tight_layout()
fig11.savefig(os.path.join(fig_dir, "AOC_design_stern_alpha_blockposition.png"))
plt.close(fig11)
print("Saved: Fig 11 — Sternberg alpha by block position")

# %% ══════════════════════════════════════════════════════════════════════════
#                         STATISTICAL ANALYSES
# ══════════════════════════════════════════════════════════════════════════

all_stats = []

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 1 — MAIN: N-BACK: AlphaPower ~ Condition * StartCond + (1|ID)
#   Does starting condition moderate the WM-load effect on alpha?
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 1 [MAIN]: N-back — AlphaPower ~ Condition * StartCond + (1|ID)")
print("  → Does it matter which condition comes first?")
print("=" * 70)

nb_m1 = nb_subj_clean[["ID", "AlphaPower", "Condition", "StartCond"]].dropna().copy()
nb_m1["Condition"] = nb_m1["Condition"].astype(str)
nb_m1["StartCond"] = nb_m1["StartCond"].astype(str)

try:
    m1_full = smf.mixedlm(
        "AlphaPower ~ C(Condition) * C(StartCond)",
        data=nb_m1, groups=nb_m1["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m1_red = smf.mixedlm(
        "AlphaPower ~ C(Condition) + C(StartCond)",
        data=nb_m1, groups=nb_m1["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m1_nocond = smf.mixedlm(
        "AlphaPower ~ C(StartCond)",
        data=nb_m1, groups=nb_m1["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m1_nostart = smf.mixedlm(
        "AlphaPower ~ C(Condition)",
        data=nb_m1, groups=nb_m1["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    # LRT: Interaction
    lr1_int = -2 * (m1_red.llf - m1_full.llf)
    df1_int = m1_full.df_modelwc - m1_red.df_modelwc
    p1_int = sp_stats.chi2.sf(lr1_int, df1_int) if df1_int > 0 else np.nan

    # LRT: Main effect of StartCond
    lr1_start = -2 * (m1_nostart.llf - m1_red.llf)
    df1_start = m1_red.df_modelwc - m1_nostart.df_modelwc
    p1_start = sp_stats.chi2.sf(lr1_start, df1_start) if df1_start > 0 else np.nan

    # LRT: Main effect of Condition
    lr1_cond = -2 * (m1_nocond.llf - m1_red.llf)
    df1_cond = m1_red.df_modelwc - m1_nocond.df_modelwc
    p1_cond = sp_stats.chi2.sf(lr1_cond, df1_cond) if df1_cond > 0 else np.nan

    print(f"\n  Condition main effect:    χ²({df1_cond}) = {lr1_cond:.3f}, p = {p1_cond:.4f}")
    print(f"  StartCond main effect:   χ²({df1_start}) = {lr1_start:.3f}, p = {p1_start:.4f}")
    print(f"  Condition × StartCond:   χ²({df1_int}) = {lr1_int:.3f}, p = {p1_int:.4f}")
    print("\n  Full model summary:")
    print(m1_full.summary())

    all_stats.append({
        "Analysis": "N-back [MAIN]: Condition main effect",
        "Model": "AlphaPower ~ Condition + StartCond + (1|ID) vs ~ StartCond + (1|ID)",
        "LR_Chi2": round(lr1_cond, 3), "df": df1_cond, "p": round(p1_cond, 4),
        "Note": "WM load effect on alpha"
    })
    all_stats.append({
        "Analysis": "N-back [MAIN]: StartCond main effect",
        "Model": "AlphaPower ~ Condition + StartCond + (1|ID) vs ~ Condition + (1|ID)",
        "LR_Chi2": round(lr1_start, 3), "df": df1_start, "p": round(p1_start, 4),
        "Note": "Does starting condition shift alpha overall?"
    })
    all_stats.append({
        "Analysis": "N-back [MAIN]: Condition × StartCond interaction",
        "Model": "AlphaPower ~ Condition * StartCond + (1|ID)",
        "LR_Chi2": round(lr1_int, 3), "df": df1_int, "p": round(p1_int, 4),
        "Note": "Does starting condition moderate the WM-load effect?"
    })
except Exception as e:
    print(f"  Model 1 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 2 — N-BACK: One-way test of alpha change (3-back − 1-back) by StartCond
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 2: Kruskal–Wallis on alpha change (3-back − 1-back) by StartCond")
print("=" * 70)

groups_change = []
for s_val in [1, 2, 3]:
    g = nb_wide.loc[nb_wide["StartCond"] == s_val, "AlphaChange"].dropna().values
    groups_change.append(g)
    print(f"  Start {s_val}-back: n={len(g)}, mean={np.mean(g):.4f}, SD={np.std(g, ddof=1):.4f}")

if all(len(g) >= 2 for g in groups_change):
    H, p_kw = sp_stats.kruskal(*groups_change)
    print(f"\n  Kruskal–Wallis: H = {H:.3f}, p = {p_kw:.4f}")
    all_stats.append({
        "Analysis": "N-back: Alpha change by StartCond (Kruskal-Wallis)",
        "Model": "Kruskal-Wallis on (3-back − 1-back) grouped by StartCond",
        "LR_Chi2": round(H, 3), "df": 2, "p": round(p_kw, 4),
        "Note": "Non-parametric group comparison"
    })

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 3 — N-BACK: AlphaPower ~ Condition * Block_c + (1|ID)
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 3: N-back — AlphaPower ~ Condition * Block + (1|ID)")
print("=" * 70)

nb_m3 = nb[["ID", ALPHA_COL_NB_TRIALS, "Condition", "Block"]].dropna().copy()
nb_m3["Condition"] = nb_m3["Condition"].astype(str)
nb_m3["Block_c"] = nb_m3["Block"] - nb_m3["Block"].mean()

try:
    m3_full = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition) * Block_c",
        data=nb_m3, groups=nb_m3["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m3_red = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition) + Block_c",
        data=nb_m3, groups=nb_m3["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr3 = -2 * (m3_red.llf - m3_full.llf)
    df3 = m3_full.df_modelwc - m3_red.df_modelwc
    p3 = sp_stats.chi2.sf(lr3, df3) if df3 > 0 else np.nan

    print(f"\n  Interaction LRT: χ²({df3}) = {lr3:.3f}, p = {p3:.4f}")

    # Main effect of Block
    m3_noblk = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition)",
        data=nb_m3, groups=nb_m3["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr3b = -2 * (m3_noblk.llf - m3_red.llf)
    df3b = m3_red.df_modelwc - m3_noblk.df_modelwc
    p3b = sp_stats.chi2.sf(lr3b, df3b) if df3b > 0 else np.nan
    print(f"  Main effect of Block: χ²({df3b}) = {lr3b:.3f}, p = {p3b:.4f}")
    print("\n  Full model summary:")
    print(m3_full.summary())

    all_stats.append({
        "Analysis": "N-back: Condition × Block interaction",
        "Model": f"{ALPHA_COL_NB_TRIALS} ~ Condition * Block_c + (1|ID)",
        "LR_Chi2": round(lr3, 3), "df": df3, "p": round(p3, 4),
        "Note": "Does the condition effect change across blocks?"
    })
    all_stats.append({
        "Analysis": "N-back: Block position main effect",
        "Model": f"{ALPHA_COL_NB_TRIALS} ~ Condition + Block_c vs ~ Condition",
        "LR_Chi2": round(lr3b, 3), "df": df3b, "p": round(p3b, 4),
        "Note": "Temporal drift across blocks"
    })
except Exception as e:
    print(f"  Model 3 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 4 — N-BACK: AlphaPower ~ Condition * Repetition + (1|ID)
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 4: N-back — AlphaPower ~ Condition * Repetition + (1|ID)")
print("=" * 70)

nb_m4 = nb[["ID", ALPHA_COL_NB_TRIALS, "Condition", "Repetition"]].dropna().copy()
nb_m4["Condition"] = nb_m4["Condition"].astype(str)
nb_m4["Repetition"] = nb_m4["Repetition"].astype(str)

try:
    m4_full = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition) * C(Repetition)",
        data=nb_m4, groups=nb_m4["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m4_red = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition) + C(Repetition)",
        data=nb_m4, groups=nb_m4["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr4 = -2 * (m4_red.llf - m4_full.llf)
    df4 = m4_full.df_modelwc - m4_red.df_modelwc
    p4 = sp_stats.chi2.sf(lr4, df4) if df4 > 0 else np.nan

    m4_norep = smf.mixedlm(
        f"{ALPHA_COL_NB_TRIALS} ~ C(Condition)",
        data=nb_m4, groups=nb_m4["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr4b = -2 * (m4_norep.llf - m4_red.llf)
    df4b = m4_red.df_modelwc - m4_norep.df_modelwc
    p4b = sp_stats.chi2.sf(lr4b, df4b) if df4b > 0 else np.nan

    print(f"\n  Interaction LRT: χ²({df4}) = {lr4:.3f}, p = {p4:.4f}")
    print(f"  Main effect of Repetition: χ²({df4b}) = {lr4b:.3f}, p = {p4b:.4f}")
    print("\n  Full model summary:")
    print(m4_full.summary())

    all_stats.append({
        "Analysis": "N-back: Condition × Repetition interaction",
        "Model": f"{ALPHA_COL_NB_TRIALS} ~ Condition * Repetition + (1|ID)",
        "LR_Chi2": round(lr4, 3), "df": df4, "p": round(p4, 4),
        "Note": "1st vs 2nd block of same condition"
    })
    all_stats.append({
        "Analysis": "N-back: Repetition main effect",
        "Model": f"{ALPHA_COL_NB_TRIALS} ~ Condition + Repetition vs ~ Condition",
        "LR_Chi2": round(lr4b, 3), "df": df4b, "p": round(p4b, 4),
        "Note": "1st vs 2nd occurrence overall"
    })
except Exception as e:
    print(f"  Model 4 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 5 — STERNBERG: AlphaPower ~ Condition * PrevCondition + (1|ID)
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 5: Sternberg — AlphaPower ~ Condition * PrevCondition + (1|ID)")
print("=" * 70)

st_m5 = st_seq[["ID", ALPHA_COL_ST_TRIALS, "Condition", "PrevCondition"]].dropna().copy()
st_m5["Condition"] = st_m5["Condition"].astype(str)
st_m5["PrevCondition"] = st_m5["PrevCondition"].astype(str)

try:
    m5_full = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) * C(PrevCondition)",
        data=st_m5, groups=st_m5["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m5_red = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) + C(PrevCondition)",
        data=st_m5, groups=st_m5["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr5 = -2 * (m5_red.llf - m5_full.llf)
    df5 = m5_full.df_modelwc - m5_red.df_modelwc
    p5 = sp_stats.chi2.sf(lr5, df5) if df5 > 0 else np.nan

    m5_noprev = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition)",
        data=st_m5, groups=st_m5["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr5b = -2 * (m5_noprev.llf - m5_red.llf)
    df5b = m5_red.df_modelwc - m5_noprev.df_modelwc
    p5b = sp_stats.chi2.sf(lr5b, df5b) if df5b > 0 else np.nan

    print(f"\n  Interaction LRT: χ²({df5}) = {lr5:.3f}, p = {p5:.4f}")
    print(f"  Main effect of PrevCondition: χ²({df5b}) = {lr5b:.3f}, p = {p5b:.4f}")
    print("\n  Full model summary:")
    print(m5_full.summary())

    all_stats.append({
        "Analysis": "Sternberg: Condition × PrevCondition interaction",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition * PrevCondition + (1|ID)",
        "LR_Chi2": round(lr5, 3), "df": df5, "p": round(p5, 4),
        "Note": "Sequential dependency"
    })
    all_stats.append({
        "Analysis": "Sternberg: PrevCondition main effect",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition + PrevCond vs ~ Condition",
        "LR_Chi2": round(lr5b, 3), "df": df5b, "p": round(p5b, 4),
        "Note": "Previous trial carry-over"
    })
except Exception as e:
    print(f"  Model 5 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 6 — STERNBERG: AlphaPower ~ Condition * Switch + (1|ID)
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 6: Sternberg — AlphaPower ~ Condition * Switch + (1|ID)")
print("=" * 70)

st_m6 = st_switch[["ID", ALPHA_COL_ST_TRIALS, "Condition", "ConditionSwitch"]].dropna().copy()
st_m6["Condition"] = st_m6["Condition"].astype(str)
st_m6["ConditionSwitch"] = st_m6["ConditionSwitch"].astype(str)

try:
    m6_full = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) * C(ConditionSwitch)",
        data=st_m6, groups=st_m6["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m6_red = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) + C(ConditionSwitch)",
        data=st_m6, groups=st_m6["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr6 = -2 * (m6_red.llf - m6_full.llf)
    df6 = m6_full.df_modelwc - m6_red.df_modelwc
    p6 = sp_stats.chi2.sf(lr6, df6) if df6 > 0 else np.nan

    m6_nosw = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition)",
        data=st_m6, groups=st_m6["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr6b = -2 * (m6_nosw.llf - m6_red.llf)
    df6b = m6_red.df_modelwc - m6_nosw.df_modelwc
    p6b = sp_stats.chi2.sf(lr6b, df6b) if df6b > 0 else np.nan

    print(f"\n  Interaction LRT: χ²({df6}) = {lr6:.3f}, p = {p6:.4f}")
    print(f"  Main effect of Switch: χ²({df6b}) = {lr6b:.3f}, p = {p6b:.4f}")
    print("\n  Full model summary:")
    print(m6_full.summary())

    all_stats.append({
        "Analysis": "Sternberg: Condition × Switch interaction",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition * Switch + (1|ID)",
        "LR_Chi2": round(lr6, 3), "df": df6, "p": round(p6, 4),
        "Note": "Switch vs repeat × WM load"
    })
    all_stats.append({
        "Analysis": "Sternberg: Switch main effect",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition + Switch vs ~ Condition",
        "LR_Chi2": round(lr6b, 3), "df": df6b, "p": round(p6b, 4),
        "Note": "Condition switch cost on alpha"
    })
except Exception as e:
    print(f"  Model 6 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# STAT 7 — STERNBERG: AlphaPower ~ Condition * Block + (1|ID)
# ──────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STAT 7: Sternberg — AlphaPower ~ Condition * Block + (1|ID)")
print("=" * 70)

st_m7 = st[["ID", ALPHA_COL_ST_TRIALS, "Condition", "Block"]].dropna().copy()
st_m7["Condition"] = st_m7["Condition"].astype(str)
st_m7["Block_c"] = st_m7["Block"] - st_m7["Block"].mean()

try:
    m7_full = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) * Block_c",
        data=st_m7, groups=st_m7["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    m7_red = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition) + Block_c",
        data=st_m7, groups=st_m7["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr7 = -2 * (m7_red.llf - m7_full.llf)
    df7 = m7_full.df_modelwc - m7_red.df_modelwc
    p7 = sp_stats.chi2.sf(lr7, df7) if df7 > 0 else np.nan

    m7_noblk = smf.mixedlm(
        f"{ALPHA_COL_ST_TRIALS} ~ C(Condition)",
        data=st_m7, groups=st_m7["ID"]
    ).fit(reml=False, method="lbfgs", maxiter=300, disp=False)

    lr7b = -2 * (m7_noblk.llf - m7_red.llf)
    df7b = m7_red.df_modelwc - m7_noblk.df_modelwc
    p7b = sp_stats.chi2.sf(lr7b, df7b) if df7b > 0 else np.nan

    print(f"\n  Interaction LRT: χ²({df7}) = {lr7:.3f}, p = {p7:.4f}")
    print(f"  Main effect of Block: χ²({df7b}) = {lr7b:.3f}, p = {p7b:.4f}")
    print("\n  Full model summary:")
    print(m7_full.summary())

    all_stats.append({
        "Analysis": "Sternberg: Condition × Block interaction",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition * Block_c + (1|ID)",
        "LR_Chi2": round(lr7, 3), "df": df7, "p": round(p7, 4),
        "Note": "Temporal drift × condition"
    })
    all_stats.append({
        "Analysis": "Sternberg: Block main effect",
        "Model": f"{ALPHA_COL_ST_TRIALS} ~ Condition + Block_c vs ~ Condition",
        "LR_Chi2": round(lr7b, 3), "df": df7b, "p": round(p7b, 4),
        "Note": "Temporal drift in mixed design"
    })
except Exception as e:
    print(f"  Model 7 failed: {e}")

# %% ──────────────────────────────────────────────────────────────────────────
# SAVE SUMMARY TABLE
# ──────────────────────────────────────────────────────────────────────────────

if all_stats:
    stats_df = pd.DataFrame(all_stats)
    stats_csv = os.path.join(stats_dir, "AOC_design_effects_summary.csv")
    stats_df.to_csv(stats_csv, index=False)
    print(f"\nSaved statistical summary → {stats_csv}")

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 0)
    print("\n" + "=" * 70)
    print("SUMMARY OF ALL DESIGN-EFFECT ANALYSES")
    print("=" * 70)
    print(stats_df.to_string(index=False))
else:
    print("No statistical results to save.")

print("\n✓ All design-effect analyses complete.")
