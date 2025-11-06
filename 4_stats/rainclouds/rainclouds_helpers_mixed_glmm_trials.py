"""
Rainclouds — Mixed-effects helpers for TRIAL-LEVEL stats.

This module fits *Gaussian* mixed-effects models using statsmodels.MixedLM
with:
  - random intercepts for subject (ID),
  - optional random slopes for Condition within subject,
  - a nested Trial random effect encoded as a variance component (ID:Trial).

It then produces pairwise contrasts between Condition levels via Wald tests,
with multiple-comparisons correction.

If you later need non-Gaussian GLMMs (e.g., binomial for Accuracy=0/1,
Poisson/neg-bin for counts) *together with* random slopes and nested trial
effects, consider switching the backend to R (lme4/glmmTMB via rpy2).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


# ---------------------------
# Utilities
# ---------------------------

def _ensure_categorical_in_order(df: pd.DataFrame, col: str, categories) -> pd.DataFrame:
    """
    Ensure `col` is an ordered categorical with the provided categories.
    """
    out = df.copy()
    out[col] = pd.Categorical(out[col], categories=categories, ordered=True)
    return out


def _build_id_trial_factor(df: pd.DataFrame, id_col: str, trial_col: str) -> pd.DataFrame:
    """
    Create a nested 'ID_Trial' factor suitable for vc_formula in MixedLM.
    """
    out = df.copy()
    out["ID"] = out[id_col].astype(str)
    out["Trial"] = out[trial_col].astype(str)
    out["ID_Trial"] = out["ID"] + "#" + out["Trial"]
    return out


# ---------------------------
# Model fitting
# ---------------------------

def _mixedlm_fit_with_vc(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    id_col: str,
    trial_col: str,
    add_random_slopes: bool = True
):
    """
    Fit Gaussian MixedLM with:
      fixed:           value ~ C(group_col)
      random (ID):     intercept (+ optional slopes for C(group_col))
      variance comps:  trial nested within ID via C(ID_Trial)

    Returns
    -------
    res : statsmodels MixedLMResults
    dfc : DataFrame (copy used internally with ensured types/columns)
    """
    dfc = _build_id_trial_factor(df, id_col=id_col, trial_col=trial_col)

    # Fixed-effects: treatment coding for Condition
    # Ensure categorical ordering is already set by caller
    if not isinstance(dfc[group_col].dtype, CategoricalDtype):
        cats = sorted(dfc[group_col].dropna().unique().tolist())
        dfc = _ensure_categorical_in_order(dfc, group_col, cats)

    fe_formula = f"{value_col} ~ C({group_col})"
    re_formula = f"1 + C({group_col})" if add_random_slopes else "1"

    # Encode nested trial as a variance component
    vc_formula = {"trial": "0 + C(ID_Trial)"}

    model = smf.mixedlm(
        fe_formula,
        data=dfc,
        groups=dfc["ID"],
        re_formula=re_formula,
        vc_formula=vc_formula
    )

    # LBFGS is typically stable here; fall back to Nelder-Mead if needed.
    res = model.fit(reml=True, method="lbfgs")
    return res, dfc


# ---------------------------
# Contrasts and pairwise tests
# ---------------------------

def _predicted_means_and_cov_from_res(res, dfc: pd.DataFrame, group_col: str):
    """
    With treatment coding:
      mean(L1) = Intercept
      mean(Li) = Intercept + C(group)[T.Li], i>=2
    """
    levels = list(dfc[group_col].cat.categories)
    fe = res.params
    cov = res.cov_params()

    means = {}
    L = {}
    names = list(fe.index)
    i_int = names.index("Intercept")

    for i, lev in enumerate(levels):
        if i == 0:
            means[lev] = fe["Intercept"]
            c = np.zeros(len(fe))
            c[i_int] = 1.0
            L[lev] = c
        else:
            pname = f"C({group_col})[T.{lev}]"
            beta = fe.get(pname, 0.0)
            means[lev] = fe["Intercept"] + beta
            c = np.zeros(len(fe))
            c[i_int] = 1.0
            if pname in fe.index:
                c[names.index(pname)] = 1.0
            L[lev] = c

    return levels, means, L, cov, fe


def _wald_contrast(c_vec: np.ndarray, fe: np.ndarray, cov: np.ndarray):
    """
    Wald test for a linear contrast c'β.
    Returns (estimate, se, z, p).
    """
    est = float(np.dot(c_vec, fe))
    se = float(np.sqrt(np.dot(c_vec, np.dot(cov, c_vec))))
    z = est / se if se > 0 else np.nan
    p = 2 * (1 - norm.cdf(abs(z)))
    return est, se, z, p


def mixedlm_pairwise_contrasts(
    df: pd.DataFrame,
    value_col: str = "value",
    group_col: str = "Condition",
    id_col: str = "ID",
    trial_col: str = "Trial",
    add_random_slopes: bool = True,
    p_adjust: str | None = "bonferroni"
) -> pd.DataFrame:
    """
    Fit MixedLM with nested Trial VC and return pairwise contrasts for Condition.

    Parameters
    ----------
    df : DataFrame
        Must contain value_col, group_col, id_col, trial_col.
        group_col should already be a Categorical with the desired order.
    value_col : str
    group_col : str
    id_col : str
    trial_col : str
    add_random_slopes : bool
        If True, random slopes for Condition within ID are estimated.
    p_adjust : {'bonferroni', 'fdr_bh', None}

    Returns
    -------
    DataFrame with columns:
        ['group1', 'group2', 'estimate', 'se', 'z', 'p', 'p_adj']
    """
    # Ensure categorical ordering if not already set
    if not isinstance(df[group_col].dtype, CategoricalDtype):
        cats = sorted(df[group_col].dropna().unique().tolist())
        df = _ensure_categorical_in_order(df, group_col, cats)

    res, dfc = _mixedlm_fit_with_vc(
        df=df,
        value_col=value_col,
        group_col=group_col,
        id_col=id_col,
        trial_col=trial_col,
        add_random_slopes=add_random_slopes
    )

    levels, means, L, cov, fe = _predicted_means_and_cov_from_res(res, dfc, group_col)

    pairs = []
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            g1, g2 = levels[i], levels[j]
            c_vec = L[g2] - L[g1]
            est, se, z, p = _wald_contrast(c_vec, fe, cov)
            pairs.append((g1, g2, est, se, z, p))

    out = pd.DataFrame(pairs, columns=["group1", "group2", "estimate", "se", "z", "p"])

    if p_adjust is not None and len(out) > 0:
        method = "bonferroni" if p_adjust == "bonferroni" else "fdr_bh"
        out["p_adj"] = multipletests(out["p"].to_numpy(), method=method)[1]
    else:
        out["p_adj"] = out["p"]

    return out
