#!/usr/bin/env python3
"""AOC stats for split Alpha Amp/Red (Sternberg).

Inputs:
  1) Base split stats input (legacy):
     data/stats/splits/AOC_splitAmpRed_sternberg_stats_input.csv
  2) Follow-up model input (from MATLAB preprocessing):
     data/stats/splits/AOC_splitAmpRed_followup_model_input.csv

This script runs:
  - legacy split analyses (gaze group comparison, Group x Load, LMM),
  - follow-up models:
      AlphaPower_FOOOF_bl ~ GazeTOI_z * Load (+ random effects),
      ReactionTime ~ GazeTOI_z * Load + Group (+ random effects),
      Accuracy GLMM-like model with binomial family.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import statsmodels.api as sm


BASE_DIR = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
INPUT_CSV = f"{BASE_DIR}/data/stats/splits/AOC_splitAmpRed_sternberg_stats_input.csv"
FOLLOWUP_CSV = f"{BASE_DIR}/data/stats/splits/AOC_splitAmpRed_followup_model_input.csv"

LOAD_ORDER = [2, 4, 6]
LOAD_LABELS = {2: "WM load 2", 4: "WM load 4", 6: "WM load 6"}
GROUP_ORDER = ["Reduction", "Amplification"]


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
    print("\n=== Parametric group comparison (Reduction > Amplification) ===")
    subj = (
        df_incl.groupby(["ID", "Group"], as_index=False)["GazeDeviationFullBL"]
        .mean(numeric_only=True)
        .rename(columns={"GazeDeviationFullBL": "SubjMeanDev"})
    )

    red = subj.loc[subj["Group"] == "Reduction", "SubjMeanDev"].dropna().to_numpy()
    amp = subj.loc[subj["Group"] == "Amplification", "SubjMeanDev"].dropna().to_numpy()

    if len(red) < 2 or len(amp) < 2:
        print(f"Dev: insufficient data (n_red={len(red)}, n_amp={len(amp)})")
        return

    t_two = stats.ttest_ind(red, amp, equal_var=True, alternative="two-sided")
    t_one = stats.ttest_ind(red, amp, equal_var=True, alternative="greater")
    d_val = pooled_cohens_d(red, amp)
    dfree = len(red) + len(amp) - 2

    print(
        "Dev: "
        f"Red mean={np.mean(red):.2f} (SD={np.std(red, ddof=1):.2f}), "
        f"Amp mean={np.mean(amp):.2f} (SD={np.std(amp, ddof=1):.2f}); "
        f"diff={np.mean(red) - np.mean(amp):.2f}; "
        f"one-tailed p(Red>Amp)={t_one.pvalue:.4f}; "
        f"two-tailed p={t_two.pvalue:.4f}; d={d_val:.3f}; "
        f"t({dfree})={t_two.statistic:.2f}; n_red={len(red)}, n_amp={len(amp)}"
    )


def run_group_load_anova_like(df_metric: pd.DataFrame, metric_name: str) -> None:
    print(f"\n=== Group x Load analysis ({metric_name}) ===")
    # OLS with subject fixed effect reproduces repeated-measures structure.
    # Type-II ANOVA is used to report Group, Load, Group:Load terms.
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
    print("=== AOC stats splits Amp/Red (Python) ===")
    print(f"Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)

    req = [
        "ID",
        "LoadValue",
        "LoadLabel",
        "Group",
        "Included",
        "AlphaPower_FOOOF_bl",
        "GazeDeviationFullBL",
        "GazeDev_pct_0_2s",
    ]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

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

    # 1) Group comparison on subject means of gaze deviation.
    run_group_comparison(df)

    # 2) MATLAB block: Mixed ANOVA-style output for gaze deviation.
    dev_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "GazeDeviationFullBL"]].rename(
        columns={"GazeDeviationFullBL": "Value"}
    ).dropna()
    run_group_load_anova_like(dev_df, "Dev")

    # 3) MATLAB block: LMM for gaze deviation.
    run_lmm(dev_df, "Dev")

    # 4) Follow-up metrics (Alpha and GazeDev_pct_0_2s).
    print("\n=== Follow-up: Alpha and gaze deviation (%) vs load (Group x Load) ===")
    print("Alpha: mean AlphaPower_FOOOF_bl [dB] per load.")
    print("Gaze: mean GazeDev_pct_0_2s [%] per load.")

    print_descriptives(df, "AlphaPower_FOOOF_bl", "Alpha (AlphaPower_FOOOF_bl [dB])")
    print_descriptives(df, "GazeDev_pct_0_2s", "GazeDev_pct_0_2s [%]")

    alpha_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "AlphaPower_FOOOF_bl"]].rename(
        columns={"AlphaPower_FOOOF_bl": "Value"}
    ).dropna()
    run_group_load_anova_like(alpha_df, "Alpha")
    run_lmm(alpha_df, "Alpha")

    gaze_pct_df = df[["Subject", "Group", "LoadValue", "LoadLabel", "GazeDev_pct_0_2s"]].rename(
        columns={"GazeDev_pct_0_2s": "Value"}
    ).dropna()
    run_group_load_anova_like(gaze_pct_df, "GazeDev_pct_0_2s")
    run_lmm(gaze_pct_df, "GazeDev_pct_0_2s")

    # Follow-up models from dedicated CSV generated by MATLAB.
    run_followup_models()


def run_followup_models() -> None:
    print("\n=== Follow-up models from AOC_splitAmpRed_followup_model_input.csv ===")
    print(f"Loading: {FOLLOWUP_CSV}")
    dff = pd.read_csv(FOLLOWUP_CSV)

    req = [
        "Subject", "Group", "LoadValue", "GazeTOI_z",
        "AlphaPower_FOOOF_bl", "ReactionTime", "Accuracy", "NTrials",
    ]
    missing = [c for c in req if c not in dff.columns]
    if missing:
        raise ValueError(f"Missing required follow-up columns: {missing}")

    dff = dff.copy()
    dff["Subject"] = dff["Subject"].astype(str)
    dff["Group"] = pd.Categorical(dff["Group"], categories=GROUP_ORDER, ordered=True)
    dff["LoadValue"] = pd.to_numeric(dff["LoadValue"], errors="coerce")
    dff = dff[dff["LoadValue"].isin(LOAD_ORDER)].copy()
    dff["LoadLabel"] = pd.Categorical(
        dff["LoadValue"].map(LOAD_LABELS),
        categories=[LOAD_LABELS[x] for x in LOAD_ORDER],
        ordered=True,
    )

    # Model 1: AlphaPower_FOOOF_bl ~ GazeTOI_z * Load
    print("\n--- Follow-up Model 1: AlphaPower_FOOOF_bl ~ GazeTOI_z * Load ---")
    d1 = dff[["Subject", "LoadValue", "LoadLabel", "GazeTOI_z", "AlphaPower_FOOOF_bl"]].dropna().copy()
    if len(d1) < 15:
        print("Insufficient observations for Model 1.")
    else:
        try:
            m1, fml1 = fit_mixedlm_with_fallback("AlphaPower_FOOOF_bl ~ GazeTOI_z * C(LoadLabel)", d1)
            print(f"Formula: {fml1}")
            print(m1.summary())
        except Exception as exc:
            print(f"Model 1 failed: {exc}")

    # Model 2: ReactionTime ~ GazeTOI_z * Load + Group
    print("\n--- Follow-up Model 2: ReactionTime ~ GazeTOI_z * Load + Group ---")
    d2 = dff[["Subject", "Group", "LoadValue", "LoadLabel", "GazeTOI_z", "ReactionTime"]].dropna().copy()
    if len(d2) < 15:
        print("Insufficient observations for Model 2.")
    else:
        try:
            m2, fml2 = fit_mixedlm_with_fallback("ReactionTime ~ GazeTOI_z * C(LoadLabel) + C(Group)", d2)
            print(f"Formula: {fml2}")
            print(m2.summary())
        except Exception as exc:
            print(f"Model 2 failed: {exc}")

    # Model 3: Accuracy GLMM-like (binomial family, subject-cluster robust SE)
    print("\n--- Follow-up Model 3: Accuracy GLMM (binomial-logit) ---")
    d3 = dff[["Subject", "Group", "LoadLabel", "GazeTOI_z", "Accuracy", "NTrials"]].dropna().copy()
    if len(d3) < 15:
        print("Insufficient observations for Model 3.")
        return

    # Accuracy normalization
    acc = d3["Accuracy"].astype(float).to_numpy()
    if np.nanmax(acc) > 1.5:
        d3["AccProp"] = np.clip(d3["Accuracy"] / 100.0, 0, 1)
    else:
        d3["AccProp"] = np.clip(d3["Accuracy"], 0, 1)

    # Trial weights; fallback is 100 when unavailable
    ntr = pd.to_numeric(d3["NTrials"], errors="coerce")
    if ntr.notna().sum() == 0:
        d3["BinomialSize"] = 100.0
        print("No usable NTrials found; BinomialSize set to 100.")
    else:
        med_n = float(np.nanmedian(ntr))
        if not np.isfinite(med_n) or med_n <= 0:
            med_n = 100.0
        ntr = ntr.where((ntr > 0) & np.isfinite(ntr), med_n)
        d3["BinomialSize"] = ntr
        print(f"BinomialSize uses NTrials with median fallback={med_n:.1f}.")

    try:
        # statsmodels has no frequentist GLMM equivalent to lme4::glmer.
        # This uses binomial GLM with subject fixed effects and cluster-robust SE.
        # It is the closest stable alternative in this stack.
        glm = smf.glm(
            "AccProp ~ GazeTOI_z * C(LoadLabel) + C(Group) + C(Subject)",
            data=d3,
            family=sm.families.Binomial(),
            freq_weights=d3["BinomialSize"],
        ).fit(cov_type="cluster", cov_kwds={"groups": d3["Subject"]})
        print("Formula: AccProp ~ GazeTOI_z * C(LoadLabel) + C(Group) + C(Subject)")
        print(glm.summary())
        print("Note: GLM with subject fixed effects + clustered SE used as GLMM approximation.")
    except Exception as exc:
        print(f"Model 3 failed: {exc}")


if __name__ == "__main__":
    main()
