# %% AOC Keyboard Switch Effects — Control Analysis
# Checks whether accuracy and reaction time depend on subject ID (odd vs even),
# because answer keys on the keyboard are switched between subjects in the paradigm.
#
# Odd-ID subjects and even-ID subjects have reversed key mappings for
# "match" / "no match" responses. This script tests whether that counterbalancing
# introduces any systematic bias in behavioural performance, and whether
# participants' handedness interacts with the key mapping.
#
# Analyses (per task: Sternberg + N-back):
#   1. Group comparison: Accuracy and RT for odd vs even IDs (t-test / Mann-Whitney)
#   2. Interaction: DV ~ Condition * KeyGroup + (1|ID)  (mixed model + LRT)
#   3. Per-condition group comparisons
#   4. Handedness link: cross-tab KeyGroup × Handedness; DV ~ Condition * Handedness;
#      DV ~ Condition * KeyGroup * Handedness (3-way interaction)
#   5. Figures: Accuracy and RT by key group, by handedness, and their combination
#
# Key outputs:
#   Figures (PNG); statistical summary table (CSV); console output

# %% Imports
import sys, os, warnings
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sp_stats
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*onvergence.*")

# %% Plot settings (match project style)
pal_group = {"odd": "#93B8C4", "even": "#D998A2"}  # blue = odd, pink = even
pal_cond  = ["#93B8C4", "#82AD82", "#D998A2"]       # AOC pastels

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
base_dir   = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
fig_dir    = f"{base_dir}/data/controls/keyboard-switch"
stats_dir  = fig_dir
os.makedirs(fig_dir, exist_ok=True)

# %% Task configuration
tasks = {
    "sternberg": {
        "csv":        f"{base_dir}/data/features/merged_data_sternberg_LONG.csv",
        "cond_map":   {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"},
        "cond_order": ["WM load 2", "WM load 4", "WM load 6"],
    },
    "nback": {
        "csv":        f"{base_dir}/data/features/merged_data_nback_LONG.csv",
        "cond_map":   {1: "1-back", 2: "2-back", 3: "3-back"},
        "cond_order": ["1-back", "2-back", "3-back"],
    },
}

DVs = ["Accuracy", "ReactionTime"]
DV_labels = {"Accuracy": "Accuracy [%]", "ReactionTime": "Reaction Time [s]"}
pal_hand = {"R": "#82AD82", "L": "#D998A2", "A": "#C4A882"}  # green/pink/amber

all_stats = []


# %% Helper: format p value
def fmt_p(p):
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


# %% Helper: load and prep data
def load_and_prep(task_name):
    cfg = tasks[task_name]
    dat = pd.read_csv(cfg["csv"])
    dat["ID"] = dat["ID"].astype(int)

    # Harmonise condition labels
    if np.issubdtype(dat["Condition"].dtype, np.number):
        dat["Condition"] = dat["Condition"].map(cfg["cond_map"])
    dat["Condition"] = pd.Categorical(
        dat["Condition"], categories=cfg["cond_order"], ordered=True
    )

    # Key group: odd vs even subject ID
    dat["KeyGroup"] = np.where(dat["ID"] % 2 == 1, "odd", "even")

    # Handedness — harmonise to R / L / A
    if "Handedness" in dat.columns:
        h = dat["Handedness"].astype(str).str.strip().str.upper()
        h = h.replace({
            "RECHTS": "R", "RIGHT": "R", "R": "R",
            "LINKS": "L",  "LEFT": "L",  "L": "L",
            "AMBIDEXTROUS": "A", "BEIDH": "A", "A": "A",
        })
        dat["Hand"] = h
    else:
        dat["Hand"] = np.nan

    # Clean
    dat.loc[dat.get("Accuracy", pd.Series(dtype=float)) > 100, "Accuracy"] = np.nan
    dat.replace([np.inf, -np.inf], np.nan, inplace=True)

    return dat


# ============================================================
# %% Main loop over tasks
# ============================================================
for task_name in ["sternberg", "nback"]:
    cfg = tasks[task_name]

    print(f"\n{'=' * 70}")
    print(f"  Task: {task_name.upper()} — Keyboard switch control")
    print(f"{'=' * 70}")

    dat = load_and_prep(task_name)

    # Report group sizes
    subj_groups = dat.groupby("ID")["KeyGroup"].first()
    n_odd  = (subj_groups == "odd").sum()
    n_even = (subj_groups == "even").sum()
    print(f"  Subjects: {n_odd} odd-ID, {n_even} even-ID  "
          f"(total {subj_groups.shape[0]})")

    # ──────────────────────────────────────────────────────────────
    # 1. Overall group comparison per DV (collapsed across conditions)
    # ──────────────────────────────────────────────────────────────
    print(f"\n  {'─' * 60}")
    print(f"  1. Overall group comparison (odd vs even, collapsed)")
    print(f"  {'─' * 60}")

    for dv in DVs:
        if dv not in dat.columns:
            continue

        # Subject-level means (average across conditions)
        subj_means = dat.groupby(["ID", "KeyGroup"])[dv].mean().reset_index()
        odd_vals  = subj_means.loc[subj_means["KeyGroup"] == "odd",  dv].dropna().values
        even_vals = subj_means.loc[subj_means["KeyGroup"] == "even", dv].dropna().values

        if len(odd_vals) < 2 or len(even_vals) < 2:
            continue

        # Independent t-test
        t_val, p_t = sp_stats.ttest_ind(odd_vals, even_vals)
        # Mann-Whitney U (non-parametric)
        u_val, p_u = sp_stats.mannwhitneyu(odd_vals, even_vals, alternative="two-sided")
        # Cohen's d
        pooled_sd = np.sqrt(
            ((len(odd_vals) - 1) * np.std(odd_vals, ddof=1)**2
             + (len(even_vals) - 1) * np.std(even_vals, ddof=1)**2)
            / (len(odd_vals) + len(even_vals) - 2)
        )
        d = (np.mean(odd_vals) - np.mean(even_vals)) / pooled_sd if pooled_sd > 0 else np.nan

        print(f"\n    {dv}:")
        print(f"      Odd  (n={len(odd_vals)}):  M = {np.mean(odd_vals):.4f}, "
              f"SD = {np.std(odd_vals, ddof=1):.4f}")
        print(f"      Even (n={len(even_vals)}): M = {np.mean(even_vals):.4f}, "
              f"SD = {np.std(even_vals, ddof=1):.4f}")
        print(f"      t-test:       t({len(odd_vals)+len(even_vals)-2}) = {t_val:.3f}, "
              f"p = {fmt_p(p_t)}")
        print(f"      Mann-Whitney: U = {u_val:.1f}, p = {fmt_p(p_u)}")
        print(f"      Cohen's d = {d:.3f}")

        all_stats.append({
            "Task": task_name, "Analysis": f"Overall {dv}: odd vs even",
            "Test": "t-test", "Statistic": round(t_val, 3),
            "df": len(odd_vals) + len(even_vals) - 2,
            "p": round(p_t, 4), "Effect_size_d": round(d, 3),
            "M_odd": round(np.mean(odd_vals), 4),
            "M_even": round(np.mean(even_vals), 4),
            "Note": "Collapsed across conditions"
        })

    # ──────────────────────────────────────────────────────────────
    # 2. Interaction: DV ~ Condition * KeyGroup + (1|ID)
    # ──────────────────────────────────────────────────────────────
    print(f"\n  {'─' * 60}")
    print(f"  2. Mixed model: DV ~ Condition * KeyGroup + (1|ID)")
    print(f"  {'─' * 60}")

    ref_cond = cfg["cond_order"][0]

    for dv in DVs:
        if dv not in dat.columns:
            continue

        sub = dat[["ID", "Condition", "KeyGroup", dv]].dropna().copy()
        sub["ID"] = sub["ID"].astype(str)
        if sub.shape[0] < 10:
            continue

        # Full model: interaction
        f_full = f'{dv} ~ C(Condition, Treatment("{ref_cond}")) * C(KeyGroup, Treatment("odd"))'
        # Additive: no interaction
        f_add  = f'{dv} ~ C(Condition, Treatment("{ref_cond}")) + C(KeyGroup, Treatment("odd"))'
        # No KeyGroup
        f_nocg = f'{dv} ~ C(Condition, Treatment("{ref_cond}"))'

        res_full = res_add = res_nocg = None
        for formula, label in [(f_full, "full"), (f_add, "add"), (f_nocg, "nocg")]:
            for method in ["lbfgs", "powell", "nm"]:
                try:
                    m = smf.mixedlm(formula, data=sub, groups=sub["ID"])
                    r = m.fit(reml=False, method=method, maxiter=600, disp=False)
                    if label == "full":
                        res_full = r
                    elif label == "add":
                        res_add = r
                    else:
                        res_nocg = r
                    break
                except Exception:
                    continue

        if res_full is None or res_add is None or res_nocg is None:
            print(f"\n    {dv}: model convergence issue — skipping.")
            continue

        # LRT: interaction
        lr_int = -2 * (res_add.llf - res_full.llf)
        df_int = res_full.df_modelwc - res_add.df_modelwc
        p_int  = sp_stats.chi2.sf(max(lr_int, 0), max(df_int, 1))

        # LRT: main effect of KeyGroup
        lr_kg = -2 * (res_nocg.llf - res_add.llf)
        df_kg = res_add.df_modelwc - res_nocg.df_modelwc
        p_kg  = sp_stats.chi2.sf(max(lr_kg, 0), max(df_kg, 1))

        print(f"\n    {dv}:")
        print(f"      KeyGroup main effect:   χ²({df_kg}) = {lr_kg:.3f}, p = {fmt_p(p_kg)}")
        print(f"      Condition × KeyGroup:   χ²({df_int}) = {lr_int:.3f}, p = {fmt_p(p_int)}")
        print(f"      Full model fixed effects:")
        for term in res_full.params.index:
            b  = res_full.params[term]
            se = res_full.bse.get(term, np.nan)
            p  = res_full.pvalues.get(term, np.nan)
            print(f"        {term:55s}  β={b:+.5f}  SE={se:.5f}  p={fmt_p(p)}")

        all_stats.append({
            "Task": task_name, "Analysis": f"{dv}: KeyGroup main effect (LRT)",
            "Test": "LRT", "Statistic": round(lr_kg, 3),
            "df": df_kg, "p": round(p_kg, 4), "Effect_size_d": np.nan,
            "M_odd": np.nan, "M_even": np.nan,
            "Note": f"{dv} ~ Condition + KeyGroup vs ~ Condition"
        })
        all_stats.append({
            "Task": task_name, "Analysis": f"{dv}: Condition × KeyGroup interaction (LRT)",
            "Test": "LRT", "Statistic": round(lr_int, 3),
            "df": df_int, "p": round(p_int, 4), "Effect_size_d": np.nan,
            "M_odd": np.nan, "M_even": np.nan,
            "Note": f"{dv} ~ Condition * KeyGroup vs ~ Condition + KeyGroup"
        })

    # ──────────────────────────────────────────────────────────────
    # 3. Per-condition group comparisons
    # ──────────────────────────────────────────────────────────────
    print(f"\n  {'─' * 60}")
    print(f"  3. Per-condition comparisons (odd vs even)")
    print(f"  {'─' * 60}")

    for dv in DVs:
        if dv not in dat.columns:
            continue
        print(f"\n    {dv}:")
        for cond in cfg["cond_order"]:
            csub = dat[dat["Condition"] == cond]
            odd_c  = csub.loc[csub["KeyGroup"] == "odd",  dv].dropna().values
            even_c = csub.loc[csub["KeyGroup"] == "even", dv].dropna().values
            if len(odd_c) < 2 or len(even_c) < 2:
                continue
            t_c, p_c = sp_stats.ttest_ind(odd_c, even_c)
            print(f"      {cond:15s}  odd M={np.mean(odd_c):.4f}  even M={np.mean(even_c):.4f}  "
                  f"t={t_c:.3f}  p={fmt_p(p_c)}")

    # ──────────────────────────────────────────────────────────────
    # 4. Handedness analyses
    # ──────────────────────────────────────────────────────────────
    has_hand = dat["Hand"].notna().any()

    if has_hand:
        print(f"\n  {'─' * 60}")
        print(f"  4. Handedness analyses")
        print(f"  {'─' * 60}")

        # 4a. Cross-tabulation: KeyGroup × Handedness
        xtab = pd.crosstab(
            dat.groupby("ID")["Hand"].first(),
            dat.groupby("ID")["KeyGroup"].first(),
            margins=True,
        )
        print(f"\n    Cross-tabulation (KeyGroup × Handedness):")
        print(f"    {xtab.to_string()}")

        # 4b. DV ~ Condition * Handedness + (1|ID)
        print(f"\n    Mixed models: DV ~ Condition * Handedness + (1|ID)")

        for dv in DVs:
            if dv not in dat.columns:
                continue
            sub_h = dat[["ID", "Condition", "Hand", dv]].dropna().copy()
            sub_h["ID"] = sub_h["ID"].astype(str)
            hand_vals = sub_h["Hand"].unique()
            if len(hand_vals) < 2 or sub_h.shape[0] < 10:
                print(f"\n      {dv}: insufficient handedness variability — skipping.")
                continue

            ref_hand = "R"  # right-handed as reference
            f_h_full = f'{dv} ~ C(Condition, Treatment("{ref_cond}")) * C(Hand, Treatment("{ref_hand}"))'
            f_h_add  = f'{dv} ~ C(Condition, Treatment("{ref_cond}")) + C(Hand, Treatment("{ref_hand}"))'
            f_h_noh  = f'{dv} ~ C(Condition, Treatment("{ref_cond}"))'

            res_h_full = res_h_add = res_h_noh = None
            for formula, label in [(f_h_full, "full"), (f_h_add, "add"), (f_h_noh, "noh")]:
                for method in ["lbfgs", "powell", "nm"]:
                    try:
                        m = smf.mixedlm(formula, data=sub_h, groups=sub_h["ID"])
                        r = m.fit(reml=False, method=method, maxiter=600, disp=False)
                        if label == "full":
                            res_h_full = r
                        elif label == "add":
                            res_h_add = r
                        else:
                            res_h_noh = r
                        break
                    except Exception:
                        continue

            if res_h_full is not None and res_h_add is not None and res_h_noh is not None:
                lr_h_int = -2 * (res_h_add.llf - res_h_full.llf)
                df_h_int = res_h_full.df_modelwc - res_h_add.df_modelwc
                p_h_int  = sp_stats.chi2.sf(max(lr_h_int, 0), max(df_h_int, 1))

                lr_h_me = -2 * (res_h_noh.llf - res_h_add.llf)
                df_h_me = res_h_add.df_modelwc - res_h_noh.df_modelwc
                p_h_me  = sp_stats.chi2.sf(max(lr_h_me, 0), max(df_h_me, 1))

                print(f"\n      {dv}:")
                print(f"        Handedness main effect:     χ²({df_h_me}) = {lr_h_me:.3f}, "
                      f"p = {fmt_p(p_h_me)}")
                print(f"        Condition × Handedness:     χ²({df_h_int}) = {lr_h_int:.3f}, "
                      f"p = {fmt_p(p_h_int)}")

                all_stats.append({
                    "Task": task_name, "Analysis": f"{dv}: Handedness main effect (LRT)",
                    "Test": "LRT", "Statistic": round(lr_h_me, 3),
                    "df": df_h_me, "p": round(p_h_me, 4), "Effect_size_d": np.nan,
                    "M_odd": np.nan, "M_even": np.nan,
                    "Note": f"{dv} ~ Condition + Handedness vs ~ Condition"
                })
                all_stats.append({
                    "Task": task_name,
                    "Analysis": f"{dv}: Condition × Handedness interaction (LRT)",
                    "Test": "LRT", "Statistic": round(lr_h_int, 3),
                    "df": df_h_int, "p": round(p_h_int, 4), "Effect_size_d": np.nan,
                    "M_odd": np.nan, "M_even": np.nan,
                    "Note": f"{dv} ~ Condition * Handedness vs ~ Condition + Handedness"
                })
            else:
                print(f"\n      {dv}: Handedness model convergence issue — skipping.")

        # 4c. DV ~ Condition * KeyGroup * Handedness + (1|ID)  (3-way)
        print(f"\n    3-way model: DV ~ Condition * KeyGroup * Handedness + (1|ID)")

        for dv in DVs:
            if dv not in dat.columns:
                continue
            sub_3 = dat[["ID", "Condition", "KeyGroup", "Hand", dv]].dropna().copy()
            sub_3["ID"] = sub_3["ID"].astype(str)
            hand_vals = sub_3["Hand"].unique()
            if len(hand_vals) < 2 or sub_3.shape[0] < 10:
                continue

            # Full 3-way
            f_3way = (f'{dv} ~ C(Condition, Treatment("{ref_cond}")) * '
                      f'C(KeyGroup, Treatment("odd")) * '
                      f'C(Hand, Treatment("R"))')
            # All 2-way (drop 3-way interaction)
            f_2way = (f'{dv} ~ (C(Condition, Treatment("{ref_cond}")) + '
                      f'C(KeyGroup, Treatment("odd")) + '
                      f'C(Hand, Treatment("R"))) ** 2')
            # KeyGroup × Handedness only (no condition interaction)
            f_kg_h = (f'{dv} ~ C(Condition, Treatment("{ref_cond}")) + '
                      f'C(KeyGroup, Treatment("odd")) * '
                      f'C(Hand, Treatment("R"))')
            # No KeyGroup × Handedness
            f_no_kgh = (f'{dv} ~ C(Condition, Treatment("{ref_cond}")) + '
                        f'C(KeyGroup, Treatment("odd")) + '
                        f'C(Hand, Treatment("R"))')

            results_3w = {}
            for formula, label in [(f_3way, "3way"), (f_2way, "2way"),
                                   (f_kg_h, "kg_h"), (f_no_kgh, "no_kgh")]:
                for method in ["lbfgs", "powell", "nm"]:
                    try:
                        m = smf.mixedlm(formula, data=sub_3, groups=sub_3["ID"])
                        results_3w[label] = m.fit(reml=False, method=method,
                                                  maxiter=600, disp=False)
                        break
                    except Exception:
                        continue

            if "kg_h" in results_3w and "no_kgh" in results_3w:
                lr_kgh = -2 * (results_3w["no_kgh"].llf - results_3w["kg_h"].llf)
                df_kgh = results_3w["kg_h"].df_modelwc - results_3w["no_kgh"].df_modelwc
                p_kgh  = sp_stats.chi2.sf(max(lr_kgh, 0), max(df_kgh, 1))

                print(f"\n      {dv}:")
                print(f"        KeyGroup × Handedness:         χ²({df_kgh}) = {lr_kgh:.3f}, "
                      f"p = {fmt_p(p_kgh)}")

                all_stats.append({
                    "Task": task_name,
                    "Analysis": f"{dv}: KeyGroup × Handedness interaction (LRT)",
                    "Test": "LRT", "Statistic": round(lr_kgh, 3),
                    "df": df_kgh, "p": round(p_kgh, 4), "Effect_size_d": np.nan,
                    "M_odd": np.nan, "M_even": np.nan,
                    "Note": "Does handedness moderate the key-mapping effect?"
                })

            if "3way" in results_3w and "2way" in results_3w:
                lr_3w = -2 * (results_3w["2way"].llf - results_3w["3way"].llf)
                df_3w = results_3w["3way"].df_modelwc - results_3w["2way"].df_modelwc
                p_3w  = sp_stats.chi2.sf(max(lr_3w, 0), max(df_3w, 1))

                print(f"        Cond × KeyGroup × Handedness: χ²({df_3w}) = {lr_3w:.3f}, "
                      f"p = {fmt_p(p_3w)}")

                all_stats.append({
                    "Task": task_name,
                    "Analysis": f"{dv}: Condition × KeyGroup × Handedness 3-way (LRT)",
                    "Test": "LRT", "Statistic": round(lr_3w, 3),
                    "df": df_3w, "p": round(p_3w, 4), "Effect_size_d": np.nan,
                    "M_odd": np.nan, "M_even": np.nan,
                    "Note": "3-way interaction (full model)"
                })

    else:
        print(f"\n  4. Handedness: column not available — skipping.")

    # ──────────────────────────────────────────────────────────────
    # FIGURES
    # ──────────────────────────────────────────────────────────────
    for dv in DVs:
        if dv not in dat.columns:
            continue

        # ---- Figure A: Overall box + strip by KeyGroup ----
        subj_means = dat.groupby(["ID", "KeyGroup"])[dv].mean().reset_index()

        fig, ax = plt.subplots(figsize=(6, 5), facecolor="white")
        sns.boxplot(
            data=subj_means, x="KeyGroup", y=dv, order=["odd", "even"],
            palette=pal_group, width=0.5, fliersize=0,
            boxprops=dict(alpha=0.4), ax=ax
        )
        sns.stripplot(
            data=subj_means, x="KeyGroup", y=dv, order=["odd", "even"],
            palette=pal_group, size=6, alpha=0.6, jitter=0.15, ax=ax
        )
        ax.set_xlabel("Key group (subject ID parity)")
        ax.set_ylabel(DV_labels.get(dv, dv))
        ax.set_title(f"{task_name.upper()} — {dv} by key group")
        ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout()
        fname = f"AOC_keyswitch_{dv.lower()}_{task_name}_overall.png"
        fig.savefig(os.path.join(fig_dir, fname))
        plt.close(fig)
        print(f"\n  Saved: {fname}")

        # ---- Figure B: Condition × KeyGroup interaction plot ----
        fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")

        for gi, (grp, col) in enumerate(pal_group.items()):
            gsub = dat[dat["KeyGroup"] == grp]
            # Subject means per condition
            smeans = gsub.groupby(["ID", "Condition"])[dv].mean().reset_index()
            gmeans = smeans.groupby("Condition")[dv].agg(["mean", "sem"]).reset_index()

            n_grp = gsub["ID"].nunique()
            offsets = [-0.12, 0.12]
            ax.errorbar(
                np.arange(len(cfg["cond_order"])) + offsets[gi],
                gmeans["mean"].values, yerr=gmeans["sem"].values,
                color=col, marker="os"[gi], markersize=9, capsize=4,
                linewidth=2, label=f"{grp} IDs (n={n_grp})"
            )

        ax.set_xticks(range(len(cfg["cond_order"])))
        ax.set_xticklabels(cfg["cond_order"])
        ax.set_xlabel("Condition")
        ax.set_ylabel(DV_labels.get(dv, dv))
        ax.set_title(f"{task_name.upper()} — {dv} by condition × key group")
        ax.legend(frameon=False, loc="best")
        ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout()
        fname = f"AOC_keyswitch_{dv.lower()}_{task_name}_interaction.png"
        fig.savefig(os.path.join(fig_dir, fname))
        plt.close(fig)
        print(f"  Saved: {fname}")

        # ---- Figure C: Faceted box+strip per condition ----
        fig, axes = plt.subplots(
            1, len(cfg["cond_order"]),
            figsize=(5 * len(cfg["cond_order"]), 5),
            sharey=True, facecolor="white"
        )
        if len(cfg["cond_order"]) == 1:
            axes = [axes]

        for idx, cond in enumerate(cfg["cond_order"]):
            ax = axes[idx]
            csub = dat[dat["Condition"] == cond]

            sns.boxplot(
                data=csub, x="KeyGroup", y=dv, order=["odd", "even"],
                palette=pal_group, width=0.5, fliersize=0,
                boxprops=dict(alpha=0.4), ax=ax
            )
            sns.stripplot(
                data=csub, x="KeyGroup", y=dv, order=["odd", "even"],
                palette=pal_group, size=5, alpha=0.5, jitter=0.15, ax=ax
            )
            ax.set_title(cond)
            ax.set_xlabel("Key group")
            if idx == 0:
                ax.set_ylabel(DV_labels.get(dv, dv))
            else:
                ax.set_ylabel("")
            ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
            for spine in ax.spines.values():
                spine.set_visible(False)

        fig.suptitle(
            f"{task_name.upper()} — {dv} by key group per condition",
            fontsize=16, y=1.02
        )
        fig.tight_layout()
        fname = f"AOC_keyswitch_{dv.lower()}_{task_name}_bycondition.png"
        fig.savefig(os.path.join(fig_dir, fname))
        plt.close(fig)
        print(f"  Saved: {fname}")


    # ──────────────────────────────────────────────────────────────
    # HANDEDNESS FIGURES
    # ──────────────────────────────────────────────────────────────
    if has_hand:
        hand_order = [h for h in ["R", "L", "A"] if h in dat["Hand"].unique()]
        pal_h_cur = {h: pal_hand[h] for h in hand_order}

        for dv in DVs:
            if dv not in dat.columns:
                continue

            # ---- Figure D: DV by Handedness (overall) ----
            subj_hand = (dat.groupby(["ID", "Hand"])[dv].mean()
                         .reset_index().dropna(subset=["Hand"]))

            fig, ax = plt.subplots(figsize=(6, 5), facecolor="white")
            sns.boxplot(
                data=subj_hand, x="Hand", y=dv, order=hand_order,
                palette=pal_h_cur, width=0.5, fliersize=0,
                boxprops=dict(alpha=0.4), ax=ax
            )
            sns.stripplot(
                data=subj_hand, x="Hand", y=dv, order=hand_order,
                palette=pal_h_cur, size=6, alpha=0.6, jitter=0.15, ax=ax
            )
            ax.set_xlabel("Handedness")
            ax.set_ylabel(DV_labels.get(dv, dv))
            ax.set_title(f"{task_name.upper()} — {dv} by handedness")
            ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.tight_layout()
            fname = f"AOC_keyswitch_{dv.lower()}_{task_name}_handedness.png"
            fig.savefig(os.path.join(fig_dir, fname))
            plt.close(fig)
            print(f"\n  Saved: {fname}")

            # ---- Figure E: Condition × Handedness interaction plot ----
            fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")
            markers_h = {"R": "o", "L": "s", "A": "D"}

            for hi, hval in enumerate(hand_order):
                hsub = dat[dat["Hand"] == hval]
                smeans = hsub.groupby(["ID", "Condition"])[dv].mean().reset_index()
                gmeans = (smeans.groupby("Condition")[dv]
                          .agg(["mean", "sem"]).reset_index())
                n_h = hsub["ID"].nunique()
                offsets_h = np.linspace(-0.15, 0.15, len(hand_order))
                ax.errorbar(
                    np.arange(len(cfg["cond_order"])) + offsets_h[hi],
                    gmeans["mean"].values, yerr=gmeans["sem"].values,
                    color=pal_hand[hval], marker=markers_h[hval],
                    markersize=9, capsize=4, linewidth=2,
                    label=f"{hval}-handed (n={n_h})"
                )

            ax.set_xticks(range(len(cfg["cond_order"])))
            ax.set_xticklabels(cfg["cond_order"])
            ax.set_xlabel("Condition")
            ax.set_ylabel(DV_labels.get(dv, dv))
            ax.set_title(f"{task_name.upper()} — {dv} by condition × handedness")
            ax.legend(frameon=False, loc="best")
            ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.tight_layout()
            fname = f"AOC_keyswitch_{dv.lower()}_{task_name}_handedness_interaction.png"
            fig.savefig(os.path.join(fig_dir, fname))
            plt.close(fig)
            print(f"  Saved: {fname}")

            # ---- Figure F: KeyGroup × Handedness combined ----
            sub_kh = (dat[dat["Hand"].notna()]
                      .groupby(["ID", "KeyGroup", "Hand"])[dv]
                      .mean().reset_index())
            sub_kh["KG_Hand"] = sub_kh["KeyGroup"] + " / " + sub_kh["Hand"]

            kg_hand_order = [f"{kg} / {h}"
                             for kg in ["odd", "even"] for h in hand_order]
            kg_hand_order = [x for x in kg_hand_order if x in sub_kh["KG_Hand"].values]

            if len(kg_hand_order) >= 2:
                fig, ax = plt.subplots(figsize=(max(8, 2.5 * len(kg_hand_order)), 5),
                                       facecolor="white")
                sns.boxplot(
                    data=sub_kh, x="KG_Hand", y=dv, order=kg_hand_order,
                    palette=[pal_group.get(x.split(" / ")[0], "#999")
                             for x in kg_hand_order],
                    width=0.55, fliersize=0, boxprops=dict(alpha=0.4), ax=ax
                )
                sns.stripplot(
                    data=sub_kh, x="KG_Hand", y=dv, order=kg_hand_order,
                    palette=[pal_group.get(x.split(" / ")[0], "#999")
                             for x in kg_hand_order],
                    size=5, alpha=0.55, jitter=0.15, ax=ax
                )
                ax.set_xlabel("Key group / Handedness")
                ax.set_ylabel(DV_labels.get(dv, dv))
                ax.set_title(f"{task_name.upper()} — {dv} by key group × handedness")
                ax.yaxis.grid(True, linewidth=0.5, alpha=0.3)
                for spine in ax.spines.values():
                    spine.set_visible(False)
                plt.xticks(rotation=25, ha="right")
                fig.tight_layout()
                fname = (f"AOC_keyswitch_{dv.lower()}_{task_name}"
                         f"_keygroup_x_handedness.png")
                fig.savefig(os.path.join(fig_dir, fname))
                plt.close(fig)
                print(f"  Saved: {fname}")


# ============================================================
# %% Save summary table
# ============================================================
if all_stats:
    stats_df = pd.DataFrame(all_stats)
    stats_csv = os.path.join(stats_dir, "AOC_keyboard_switch_effects_summary.csv")
    stats_df.to_csv(stats_csv, index=False)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 0)
    print(f"\n{'=' * 70}")
    print("SUMMARY OF ALL KEYBOARD-SWITCH ANALYSES")
    print(f"{'=' * 70}")
    print(stats_df.to_string(index=False))
    print(f"\nSaved → {stats_csv}")
else:
    print("\nNo statistical results to save.")

print(f"\n✓ Keyboard switch control analysis complete.")
