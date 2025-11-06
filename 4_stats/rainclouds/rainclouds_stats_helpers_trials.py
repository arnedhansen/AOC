"""
Rainclouds — generic statistical helpers for TRIAL-LEVEL pipelines.

Contents
--------
- p_to_signif(p):      map p-values to significance stars.
- iqr_outlier_filter:  1.5×IQR filtering within groups (sets outliers to NaN).
"""

from __future__ import annotations
import numpy as np
import pandas as pd


def p_to_signif(p: float) -> str:
    """
    Convert a p-value to significance stars.

    Parameters
    ----------
    p : float

    Returns
    -------
    str
        '***' (p<.001), '**' (p<.01), '*' (p<.05), else 'n.s.'
    """
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "n.s."


def _mask_series(s: pd.Series) -> pd.Series:
    """
    Apply 1.5×IQR outlier masking within a vector.
    Outliers are set to NaN; original dtype/indices are preserved by caller.

    Notes
    -----
    - If fewer than 3 non-NaN values, returns unchanged.
    - If IQR==0 or non-finite, returns unchanged.
    """
    if s.size == 0:
        return s

    # Work on float array; keep NaNs intact
    x = s.to_numpy(dtype=float, copy=False)
    if np.sum(~np.isnan(x)) < 3:
        return s

    q1 = np.nanpercentile(x, 25)
    q3 = np.nanpercentile(x, 75)
    iqr = q3 - q1
    if not np.isfinite(iqr) or iqr == 0:
        return s

    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    keep = s.isna() | ((s >= lower) & (s <= upper))
    return s.where(keep)


def iqr_outlier_filter(df: pd.DataFrame, variables, by) -> pd.DataFrame:
    """
    Set outliers to NaN using 1.5×IQR, performed separately within groups.

    Parameters
    ----------
    df : DataFrame
        Input table (trial-level or otherwise).
    variables : list[str]
        Column names to filter.
    by : str | list[str]
        Column(s) defining groups (e.g., 'Condition').

    Returns
    -------
    DataFrame
        Copy of df with outliers set to NaN in the specified variables.
    """
    out = df.copy()
    by_cols = [by] if isinstance(by, str) else list(by)

    g = out.groupby(by_cols, observed=True, sort=False, dropna=False, group_keys=False)
    for v in variables:
        if v not in out.columns:
            continue
        masked = g[v].apply(_mask_series)      # keeps original row order (group_keys=False)
        out[v] = masked.reindex(out.index)     # align defensively

    return out
