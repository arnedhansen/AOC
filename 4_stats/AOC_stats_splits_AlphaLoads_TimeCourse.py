#!/usr/bin/env python3
"""AOC stats for split AlphaLoads TimeCourse (Sternberg).

Input:
  data/stats/splits/AOC_splitAlphaLoads_TimeCourse_sternberg_stats_input.csv
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm


BASE_DIR = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
INPUT_CSV = f"{BASE_DIR}/data/stats/splits/AOC_splitAlphaLoads_TimeCourse_sternberg_stats_input.csv"

LOAD_ORDER = [2, 4, 6]
LOAD_LABELS = {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"}
GROUP_ORDER = ["Decrease", "Increase"]


def pooled_cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    nx, ny = len(x), len(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if not np.isfinite(pooled) or pooled == 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / pooled


def fit_mixedlm_with_fallback(formula_fixed: str, data: pd.DataFrame):
    formulas = [
        f"{formula_fixed} + (LoadValue|Subject)",
        f"{formula_fixed} + (1|Subject)",
    ]
    last_err = None
    for fml in formulas:
        try:
            model = smf.mixedlm(fml, data=data, groups=data["Subject"])
            res = model.fit(reml=False, method="lbfgs", maxiter=500, disp=False)
            return res, fml
        except Exception as exc:
            last_err = exc
            continue
    raise RuntimeError(f"LMM failed for all candidate formulas: {last_err}")


def run_group_comparison(df_incl: pd.DataFrame) -> None:
    print("\n=== Parametric group comparison (Decrease vs Increase) ===")
    subj = (
        df_incl.groupby(["ID", "Group"], as_index=False)["DevMetric"]
        .mean(numeric_only=True)
        .rename(columns={"DevMetric": "SubjMeanDev"})
    )

    dec = subj.loc[subj["Group"] == "Decrease", "SubjMeanDev"].dropna().to_numpy()
    inc = subj.loc[subj["Group"] == "Increase", "SubjMeanDev"].dropna().to_numpy()

    if len(dec) < 2 or len(inc) < 2:
        print(f"Dev: insufficient data (n_dec={len(dec)}, n_inc={len(inc)})")
        return

    t_two = stats.ttest_ind(dec, inc, equal_var=True, alternative="two-sided")
    d_val = pooled_cohens_d(dec, inc)
    dfree = len(dec) + len(inc) - 2

    print(
        "Dev: "
        f"Dec mean={np.mean(dec):.2f} (SD={np.std(dec, ddof=1):.2f}), "
        f"Inc mean={np.mean(inc):.2f} (SD={np.std(inc, ddof=1):.2f}); "
        f"diff={np.mean(dec) - np.mean(inc):.2f}; "
        f"two-tailed p={t_two.pvalue:.4f}; d={d_val:.3f}; "
        f"t({dfree})={t_two.statistic:.2f}; n_dec={len(dec)}, n_inc={len(inc)}"
    )


def run_group_load_anova_like(df_metric: pd.DataFrame, metric_name: str) -> None:
    print(f"\n=== Group x Load analysis ({metric_name}) ===")
    try:
        fit = smf.ols(
            "Value ~ C(Group) * C(LoadLabel) + C(Subject)",
            data=df_metric,
        ).fit()
        aov = anova_lm(fit, typ=2)

        p_group = aov.loc["C(Group)", "PR(>F)"] if "C(Group)" in aov.index else np.nan
        p_load = aov.loc["C(LoadLabel)", "PR(>F)"] if "C(LoadLabel)" in aov.index else np.nan
        int_key = "C(Group):C(LoadLabel)"
        p_int = aov.loc[int_key, "PR(>F)"] if int_key in aov.index else np.nan
        n_subj = df_metric["Subject"].nunique()

        print(
            f"{metric_name}: Load p={p_load:.4f}; Group p={p_group:.4f}; "
            f"LoadxGroup p={p_int:.4f}; n={n_subj}"
        )
        print("ANOVA table:")
        print(aov.to_string())
    except Exception as exc:
        print(f"{metric_name}: Group x Load analysis failed ({exc})")


def run_lmm(df_metric: pd.DataFrame, metric_name: str) -> None:
    print(f"\n=== LMM ({metric_name}) ===")
    if len(df_metric) < 15:
        print(f"{metric_name}: insufficient data for LMM (n_obs={len(df_metric)})")
        return
    try:
        res, formula_used = fit_mixedlm_with_fallback("Value ~ C(Group) * C(LoadLabel)", df_metric)
        print(f"{metric_name}: Formula: {formula_used}")
        print(res.summary())
        print(f"n_obs={len(df_metric)}, n_subj={df_metric['Subject'].nunique()}")
    except Exception as exc:
        print(f"{metric_name}: LMM failed ({exc})")


def print_descriptives(df_incl: pd.DataFrame, value_col: str, label: str) -> None:
    print(f"\n{label}")
    for grp in GROUP_ORDER:
        print(f"  {grp}:")
        dfg = df_incl[df_incl["Group"] == grp]
        for lv in LOAD_ORDER:
            vals = dfg.loc[dfg["LoadValue"] == lv, value_col].dropna().to_numpy()
            if len(vals) == 0:
                print(f"    {LOAD_LABELS[lv]}: n=0")
                continue
            sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else np.nan
            print(f"    {LOAD_LABELS[lv]}: n={len(vals)}, mean={np.mean(vals):.6g}, SEM={sem:.6g}")


def main() -> None:
    print("=== AOC stats split AlphaLoads TimeCourse (Python) ===")
    print(f"Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)

    req = ["ID", "LoadValue", "LoadLabel", "Group", "Included", "AlphaPower_FOOOF_bl"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if "GazeDeviationFullBL" in df.columns:
        df["DevMetric"] = pd.to_numeric(df["GazeDeviationFullBL"], errors="coerce")
    elif "Dev_dB" in df.columns:
        df["DevMetric"] = pd.to_numeric(df["Dev_dB"], errors="coerce")
    else:
        raise ValueError("Missing gaze deviation column: expected GazeDeviationFullBL or Dev_dB.")

    if "GazeDev_pct_0_2s" in df.columns:
        df["GazeFollowup"] = pd.to_numeric(df["GazeDev_pct_0_2s"], errors="coerce")
        gaze_followup_label = "GazeDev_pct_0_2s"
    elif "GazeDev_dB_0p5_1p5s" in df.columns:
        df["GazeFollowup"] = pd.to_numeric(df["GazeDev_dB_0p5_1p5s"], errors="coerce")
        gaze_followup_label = "GazeDev_dB_0p5_1p5s"
    else:
        raise ValueError("Missing follow-up gaze column: expected GazeDev_pct_0_2s or GazeDev_dB_0p5_1p5s.")

    df = df.copy()
    df = df[df["Included"] == True].copy()
    df = df[df["Group"].isin(GROUP_ORDER)].copy()
    df["LoadValue"] = pd.to_numeric(df["LoadValue"], errors="coerce")
    df = df[df["LoadValue"].isin(LOAD_ORDER)].copy()
    df["LoadLabel"] = pd.Categorical(df["LoadLabel"], categories=[LOAD_LABELS[x] for x in LOAD_ORDER], ordered=True)
    df["Group"] = pd.Categorical(df["Group"], categories=GROUP_ORDER, ordered=True)
    df["Subject"] = df["ID"].astype(str)

    print(f"Included observations: {len(df)}")
    print(f"Included subjects: {df['Subject'].nunique()}")
    print("Group sizes (subjects):")
    print(df.groupby("Group")["Subject"].nunique().to_string())

    run_group_comparison(df)

    dev_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "DevMetric"]].rename(columns={"DevMetric": "Value"}).dropna()
    run_group_load_anova_like(dev_df, "Dev")
    run_lmm(dev_df, "Dev")

    print("\n=== Follow-up: Alpha and gaze deviation vs load (Group x Load) ===")
    print("Alpha: mean AlphaPower_FOOOF_bl [dB] per load.")
    print(f"Gaze: mean {gaze_followup_label} per load.")
    print_descriptives(df, "AlphaPower_FOOOF_bl", "Alpha (AlphaPower_FOOOF_bl [dB])")
    print_descriptives(df, "GazeFollowup", gaze_followup_label)

    alpha_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "AlphaPower_FOOOF_bl"]].rename(
        columns={"AlphaPower_FOOOF_bl": "Value"}
    ).dropna()
    run_group_load_anova_like(alpha_df, "Alpha")
    run_lmm(alpha_df, "Alpha")

    gaze_follow_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "GazeFollowup"]].rename(columns={"GazeFollowup": "Value"}).dropna()
    run_group_load_anova_like(gaze_follow_df, gaze_followup_label)
    run_lmm(gaze_follow_df, gaze_followup_label)


if __name__ == "__main__":
    main()
