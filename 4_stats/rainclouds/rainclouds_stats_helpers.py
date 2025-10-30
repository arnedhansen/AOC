import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def p_to_signif(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "n.s."
def _mask_series(s: pd.Series) -> pd.Series:
    """1.5×IQR mask within one group; outliers -> NaN."""
    if s.size == 0:
        return s
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

def iqr_outlier_filter(df: pd.DataFrame, variables, by):
    """
    Set outliers to NaN per group for each variable via 1.5×IQR.
    Uses SeriesGroupBy.apply with group_keys=False to keep the original row index.
    """
    out = df.copy()
    if isinstance(by, str):
        by = [by]

    g = out.groupby(by, observed=True, sort=False, dropna=False, group_keys=False)

    for v in variables:
        masked = g[v].apply(_mask_series)                 # no FutureWarning
        out[v] = masked.reindex(out.index)                # align defensively

    return out

def _mixedlm_fit(df, value_col, group_col, id_col):
    # Treatment coding with first category as baseline (like R’s default)
    df = df.copy()
    df[group_col] = pd.Categorical(df[group_col], ordered=True)
    formula = f"{value_col} ~ C({group_col})"
    model = smf.mixedlm(formula, data=df, groups=df[id_col], re_formula="1")
    res = model.fit(reml=True, method="lbfgs")
    return res, df

def _predicted_means_and_cov(res, df, group_col):
    # With treatment coding and k levels (L1 baseline),
    # mean(L1) = b0
    # mean(Li) = b0 + b_C[Li] for i>=2
    levels = list(df[group_col].cat.categories)
    fe = res.params  # fixed effects
    cov = res.cov_params()  # covariance of fixed effects
    # map param names
    # Expected names: Intercept, C(group)[T.level2], C(group)[T.level3], ...
    means = {}
    L = {}
    for i, lev in enumerate(levels):
        if i == 0:
            # intercept
            means[lev] = fe["Intercept"]
            c = np.zeros(len(fe))
            c[list(fe.index).index("Intercept")] = 1.0
            L[lev] = c
        else:
            name = f"C({group_col})[T.{lev}]"
            beta = fe.get(name, 0.0)
            means[lev] = fe["Intercept"] + beta
            c = np.zeros(len(fe))
            c[list(fe.index).index("Intercept")] = 1.0
            if name in fe.index:
                c[list(fe.index).index(name)] = 1.0
            L[lev] = c
    return levels, means, L, cov

def _contrast_test(c_vec, fe, cov):
    # c_vec: contrast vector (1D array-like)
    # fe:    fixed-effects coefficients (1D array-like)
    # cov:   covariance matrix of fixed effects (2D array)
    est = float(np.dot(c_vec, fe))
    se = float(np.sqrt(np.dot(c_vec, np.dot(cov, c_vec))))
    z = est / se if se > 0 else np.nan
    p = 2 * (1 - norm.cdf(abs(z)))  # two-sided p-value from standard normal
    return est, se, z, p

def mixedlm_pairwise_contrasts(df, value_col="value", group_col="Condition", id_col="ID", p_adjust="bonferroni"):
    res, dfc = _mixedlm_fit(df, value_col, group_col, id_col)
    levels, means, L, cov = _predicted_means_and_cov(res, dfc, group_col)
    fe = res.params

    # build pairwise list in canonical order
    pairs = []
    for i in range(len(levels)):
        for j in range(i+1, len(levels)):
            g1, g2 = levels[i], levels[j]
            # contrast mean(g2) - mean(g1) = (L[g2] - L[g1])' * beta
            c_vec = L[g2] - L[g1]
            est, se, z, p = _contrast_test(c_vec, fe, cov)
            pairs.append((g1, g2, est, se, z, p))

    out = pd.DataFrame(pairs, columns=["group1", "group2", "estimate", "se", "z", "p"])
    if p_adjust is not None and len(out) > 0:
        padj = multipletests(out["p"].to_numpy(), method=("bonferroni" if p_adjust == "bonferroni" else "fdr_bh"))[1]
        out["p_adj"] = padj
    else:
        out["p_adj"] = out["p"]
    return out
