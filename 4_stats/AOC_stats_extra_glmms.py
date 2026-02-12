# %% AOC Extra GLMMs — Recommended Supplementary Analyses
# Runs the additional mixed models recommended for the Stage-2 manuscript:
#
#   1. Trial-level GLMM for H5 (Alpha ~ Gaze * Condition, trial-level)
#   2. Random-slopes models (Condition slopes per subject)
#   3. Joint Task model  (Alpha ~ Gaze * Load * Task)
#   4. Pupil-size covariate (Alpha ~ Gaze * Condition + PupilSize)
#   5. Within-subject trial correlations (per-subject r, Fisher z, group t)
#   6. Formal mediation  (Load → Gaze → Alpha, Baron-Kenny + Sobel)
#
# Reads the same CSVs as the main stats scripts.
# Prints all results to console; saves a summary CSV.

# %% Imports
import sys, os, warnings
sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from scipy.stats import pearsonr, spearmanr, norm, chi2
import statsmodels.formula.api as smf

from functions.stats_helpers import iqr_outlier_filter
from functions.mixedlm_helpers import fit_mixedlm, drop1_lrt

warnings.filterwarnings("ignore", message=".*onvergence.*")
warnings.filterwarnings("ignore", message=".*The Hessian.*")
warnings.filterwarnings("ignore", message=".*delta_grad.*")

# %% Paths
base_dir   = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
output_dir = os.path.join(base_dir, "data", "stats", "extra_glmms")
os.makedirs(output_dir, exist_ok=True)

# %% Data sources
LONG_FILES = {
    "sternberg": os.path.join(base_dir, "data/features/merged_data_sternberg.csv"),
    "nback":     os.path.join(base_dir, "data/features/merged_data_nback.csv"),
}
TRIAL_FILES = {
    "sternberg": os.path.join(base_dir, "data/features/merged_data_sternberg_trials.csv"),
    "nback":     os.path.join(base_dir, "data/features/merged_data_nback_trials.csv"),
}

COND_MAP = {
    "sternberg": {2: "WM2", 4: "WM4", 6: "WM6"},
    "nback":     {1: "NB1", 2: "NB2", 3: "NB3"},
}
COND_ORDER = {
    "sternberg": ["WM2", "WM4", "WM6"],
    "nback":     ["NB1", "NB2", "NB3"],
}

# Gaze predictors to test in the alpha models
GAZE_VARS  = ["GazeDeviation", "MSRate", "ScanPathLength", "BCEA", "BCEALateralization"]
GAZE_VARS_TRIAL = {
    "sternberg": ["GazeDeviationFull", "ScanPathLengthFull", "MSRateFull", "BCEAFull", "BCEALatFull"],
    "nback":     ["GazeDeviationFull", "ScanPathLengthFull", "MSRateFull", "BCEAFull", "BCEALatFull"],
}
# Alpha DVs to iterate in trial-level analyses (sections 1 & 5).
# Sternberg: run twice — full epoch (0-2 s) and late/retention (1-2 s).
# N-back:    run once  — early epoch (0-1 s) only (full/late not yet extracted).
ALPHA_VARS_TRIAL = {
    "sternberg": [
        ("AlphaPowerFull", "alpha_full_0-2s"),
        ("AlphaPowerLate", "alpha_late_1-2s"),
    ],
    "nback": [
        ("AlphaPowerEarly", "alpha_early_0-1s"),
    ],
}

# Collector for all results
all_results = []

def log(section, task, model, term, stat_name, value, note=""):
    all_results.append({
        "Section": section, "Task": task, "Model": model,
        "Term": term, "Statistic": stat_name,
        "Value": value, "Note": note,
    })


# ============================================================
# %% Helper: load + prep LONG data
# ============================================================
def load_long(task):
    dat = pd.read_csv(LONG_FILES[task])
    dat.replace([np.inf, -np.inf], np.nan, inplace=True)
    if dat["Condition"].dtype in [np.float64, np.int64]:
        dat["Condition"] = dat["Condition"].map(COND_MAP[task])
    dat["Condition"] = pd.Categorical(
        dat["Condition"], categories=COND_ORDER[task], ordered=True
    )
    dat["ID"] = dat["ID"].astype(str)
    return dat


def load_trials(task):
    dat = pd.read_csv(TRIAL_FILES[task])
    dat.replace([np.inf, -np.inf], np.nan, inplace=True)
    if dat["Condition"].dtype in [np.float64, np.int64]:
        dat["Condition"] = dat["Condition"].map(COND_MAP[task])
    dat["Condition"] = pd.Categorical(
        dat["Condition"], categories=COND_ORDER[task], ordered=True
    )
    dat["ID"] = dat["ID"].astype(str)
    return dat


def safe_fit(formula, data, group="ID", reml=False, re_formula="1"):
    """Fit MixedLM with fallback optimisers."""
    for method in ["lbfgs", "powell", "nm"]:
        try:
            m = smf.mixedlm(formula, data=data, groups=data[group],
                            re_formula=re_formula)
            return m.fit(reml=reml, method=method, maxiter=600, disp=False)
        except Exception:
            continue
    return None


def fmt_p(p):
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


# ============================================================
# %% 1. TRIAL-LEVEL GLMM  (Alpha ~ Gaze * Condition + (1|ID))
# ============================================================
print("\n" + "=" * 70)
print("  1. TRIAL-LEVEL GLMMs: Alpha ~ Gaze_c * Condition + (1|ID)")
print("=" * 70)

for task in ["sternberg", "nback"]:
    trl = load_trials(task)
    ref = COND_ORDER[task][0]

    for alpha_col, alpha_label in ALPHA_VARS_TRIAL[task]:
        if alpha_col not in trl.columns:
            print(f"  [{task}] {alpha_col}: column not found. Skipping.")
            continue

        for gaze_col in GAZE_VARS_TRIAL[task]:
            sub = trl[["ID", "Condition", alpha_col, gaze_col]].dropna().copy()
            if sub.shape[0] < 50:
                print(f"  [{task}] {alpha_col} × {gaze_col}: "
                      f"too few rows ({sub.shape[0]}). Skipping.")
                continue

            # Mean-centre gaze
            sub["Gaze_c"] = sub[gaze_col] - sub[gaze_col].mean()

            # Full model (with interaction)
            f_full = f'{alpha_col} ~ Gaze_c * C(Condition, Treatment("{ref}"))'
            f_red  = f'{alpha_col} ~ Gaze_c + C(Condition, Treatment("{ref}"))'

            res_full = safe_fit(f_full, sub)
            res_red  = safe_fit(f_red, sub)

            if res_full is None or res_red is None:
                print(f"  [{task}] {alpha_col} × {gaze_col}: "
                      f"model did not converge.")
                continue

            lrt = drop1_lrt(res_full, res_red)

            # Choose final model
            if np.isfinite(lrt["p"]) and lrt["p"] < 0.05:
                final = res_full
                kept = "interaction"
            else:
                final = res_red
                kept = "no interaction"

            n_obs = sub.shape[0]
            n_subj = sub["ID"].nunique()

            print(f"\n  [{task.upper()}] {gaze_col} → {alpha_col}  ({alpha_label})")
            print(f"    N = {n_obs} trials, {n_subj} subjects")
            print(f"    Interaction LRT: χ²({lrt['df_diff']}) = {lrt['LR']:.2f}, "
                  f"p = {fmt_p(lrt['p'])}  → kept: {kept}")
            print(f"    Final model fixed effects:")
            for term in final.params.index:
                b  = final.params[term]
                se = final.bse[term] if term in final.bse.index else np.nan
                p  = final.pvalues[term] if term in final.pvalues.index else np.nan
                print(f"      {term:50s}  β={b:+.5f}  SE={se:.5f}  p={fmt_p(p)}")
                log("1_TrialGLMM", task,
                    f"{alpha_label}~{gaze_col}*Cond", term, "beta", b)
                log("1_TrialGLMM", task,
                    f"{alpha_label}~{gaze_col}*Cond", term, "p", p)

            log("1_TrialGLMM", task, f"{alpha_label}~{gaze_col}*Cond",
                "Interaction_LRT", "chi2", lrt["LR"], f"df={lrt['df_diff']}")
            log("1_TrialGLMM", task, f"{alpha_label}~{gaze_col}*Cond",
                "Interaction_LRT", "p", lrt["p"])


# ============================================================
# %% 2. RANDOM SLOPES  (DV ~ Condition + (1 + Condition | ID))
# ============================================================
print("\n" + "=" * 70)
print("  2. RANDOM SLOPES: DV ~ Condition + (1 + Condition | ID)")
print("=" * 70)

key_dvs = ["AlphaPower", "GazeDeviation", "ScanPathLength", "MSRate"]

for task in ["sternberg", "nback"]:
    dat = load_long(task)
    ref = COND_ORDER[task][0]

    for dv in key_dvs:
        if dv not in dat.columns:
            continue
        sub = dat[["ID", "Condition", dv]].dropna().copy()
        if sub["ID"].nunique() < 5:
            continue

        # Random-intercept only (baseline)
        f_ri = f'{dv} ~ C(Condition, Treatment("{ref}"))'
        res_ri = safe_fit(f_ri, sub, re_formula="1")

        # Random-intercept + random slope for Condition
        # statsmodels encodes random slopes via re_formula
        f_rs = f_ri  # same fixed effects
        res_rs = None
        try:
            m = smf.mixedlm(f_rs, data=sub, groups=sub["ID"],
                            re_formula=f'C(Condition, Treatment("{ref}"))')
            res_rs = m.fit(reml=False, method="lbfgs", maxiter=600, disp=False)
        except Exception:
            try:
                m = smf.mixedlm(f_rs, data=sub, groups=sub["ID"],
                                re_formula=f'C(Condition, Treatment("{ref}"))')
                res_rs = m.fit(reml=False, method="powell", maxiter=800, disp=False)
            except Exception:
                pass

        if res_rs is not None and res_ri is not None:
            # Compare via LRT
            lrt = drop1_lrt(res_rs, res_ri)
            status = "CONVERGED"
            print(f"\n  [{task.upper()}] {dv}: random slopes {status}")
            print(f"    LRT vs intercept-only: χ²({lrt['df_diff']}) = {lrt['LR']:.2f}, "
                  f"p = {fmt_p(lrt['p'])}")

            # Check if fixed effects change materially
            for term in res_ri.params.index:
                b_ri = res_ri.params[term]
                b_rs = res_rs.params.get(term, np.nan)
                print(f"    {term:50s}  β_RI={b_ri:+.5f}  β_RS={b_rs:+.5f}")

            log("2_RandSlopes", task, f"{dv}~Cond+(1+Cond|ID)",
                "LRT_vs_RI", "chi2", lrt["LR"], f"df={lrt['df_diff']}")
            log("2_RandSlopes", task, f"{dv}~Cond+(1+Cond|ID)",
                "LRT_vs_RI", "p", lrt["p"])
        else:
            print(f"\n  [{task.upper()}] {dv}: random slopes DID NOT CONVERGE")
            log("2_RandSlopes", task, f"{dv}~Cond+(1+Cond|ID)",
                "convergence", "status", 0, "did not converge")


# ============================================================
# %% 3. JOINT TASK MODEL  (Alpha ~ Gaze_c * Load_c * Task + (1|ID))
# ============================================================
print("\n" + "=" * 70)
print("  3. JOINT TASK MODEL: Alpha ~ Gaze_c * Load_c * Task + (1|ID)")
print("=" * 70)

frames = []
for task in ["sternberg", "nback"]:
    dat = load_long(task)
    # Numeric load for a common scale (1, 2, 3)
    cond_to_num = {c: i + 1 for i, c in enumerate(COND_ORDER[task])}
    dat["Load_num"] = dat["Condition"].map(cond_to_num).astype(float)
    dat["Task"] = task
    frames.append(dat)

joint = pd.concat(frames, ignore_index=True)

for gaze_col in ["GazeDeviation", "ScanPathLength", "MSRate"]:
    if gaze_col not in joint.columns:
        continue
    sub = joint[["ID", "Task", "Load_num", "AlphaPower", gaze_col]].dropna().copy()
    if sub.shape[0] < 20:
        continue

    # Centre continuous predictors
    sub["Gaze_c"] = sub[gaze_col] - sub[gaze_col].mean()
    sub["Load_c"] = sub["Load_num"] - sub["Load_num"].mean()

    # Full 3-way model
    f_full = 'AlphaPower ~ Gaze_c * Load_c * C(Task, Treatment("sternberg"))'
    f_no3  = 'AlphaPower ~ (Gaze_c + Load_c + C(Task, Treatment("sternberg")))' \
             ' ** 2'  # all 2-way, drop 3-way

    res_full = safe_fit(f_full, sub)
    res_no3  = safe_fit(f_no3, sub)

    if res_full is None or res_no3 is None:
        print(f"  {gaze_col}: model did not converge.")
        continue

    lrt_3way = drop1_lrt(res_full, res_no3)

    print(f"\n  Gaze predictor: {gaze_col}")
    print(f"    N = {sub.shape[0]} obs, {sub['ID'].nunique()} subjects")
    print(f"    3-way interaction LRT: χ²({lrt_3way['df_diff']}) = "
          f"{lrt_3way['LR']:.2f}, p = {fmt_p(lrt_3way['p'])}")
    print(f"    Full model fixed effects:")
    for term in res_full.params.index:
        b  = res_full.params[term]
        se = res_full.bse[term] if term in res_full.bse.index else np.nan
        p  = res_full.pvalues[term] if term in res_full.pvalues.index else np.nan
        print(f"      {term:55s}  β={b:+.6f}  SE={se:.6f}  p={fmt_p(p)}")
        log("3_JointTask", "both", f"Alpha~{gaze_col}*Load*Task", term, "beta", b)
        log("3_JointTask", "both", f"Alpha~{gaze_col}*Load*Task", term, "p", p)

    log("3_JointTask", "both", f"Alpha~{gaze_col}*Load*Task",
        "3way_LRT", "chi2", lrt_3way["LR"])
    log("3_JointTask", "both", f"Alpha~{gaze_col}*Load*Task",
        "3way_LRT", "p", lrt_3way["p"])


# ============================================================
# %% 4. PUPIL-SIZE COVARIATE  (Alpha ~ Gaze * Condition + Pupil + (1|ID))
# ============================================================
print("\n" + "=" * 70)
print("  4. PUPIL COVARIATE: Alpha ~ Gaze_c * Cond + PupilSize_c + (1|ID)")
print("=" * 70)

for task in ["sternberg", "nback"]:
    dat = load_long(task)
    ref = COND_ORDER[task][0]

    for gaze_col in GAZE_VARS:
        if gaze_col not in dat.columns or "PupilSize" not in dat.columns:
            continue
        sub = dat[["ID", "Condition", "AlphaPower", gaze_col, "PupilSize"]].dropna().copy()
        if sub.shape[0] < 10:
            continue

        sub["Gaze_c"]  = sub[gaze_col]  - sub[gaze_col].mean()
        sub["Pupil_c"] = sub["PupilSize"] - sub["PupilSize"].mean()

        # Model WITHOUT pupil
        f_nopup = f'AlphaPower ~ Gaze_c * C(Condition, Treatment("{ref}"))'
        # Model WITH pupil
        f_pup   = f'AlphaPower ~ Gaze_c * C(Condition, Treatment("{ref}")) + Pupil_c'

        res_nopup = safe_fit(f_nopup, sub)
        res_pup   = safe_fit(f_pup, sub)

        if res_nopup is None or res_pup is None:
            print(f"  [{task}] {gaze_col}: did not converge.")
            continue

        lrt = drop1_lrt(res_pup, res_nopup)

        print(f"\n  [{task.upper()}] {gaze_col} + PupilSize")
        print(f"    Pupil covariate LRT: χ²({lrt['df_diff']}) = {lrt['LR']:.2f}, "
              f"p = {fmt_p(lrt['p'])}")

        # Compare gaze beta with and without pupil
        gaze_terms = [t for t in res_nopup.params.index if "Gaze_c" in t]
        for term in gaze_terms:
            b_no = res_nopup.params.get(term, np.nan)
            b_pu = res_pup.params.get(term, np.nan)
            p_no = res_nopup.pvalues.get(term, np.nan)
            p_pu = res_pup.pvalues.get(term, np.nan)
            print(f"    {term:50s}")
            print(f"      Without pupil: β={b_no:+.5f}  p={fmt_p(p_no)}")
            print(f"      With    pupil: β={b_pu:+.5f}  p={fmt_p(p_pu)}")
            log("4_PupilCov", task, f"Alpha~{gaze_col}*Cond+Pupil",
                term, "beta_noPupil", b_no)
            log("4_PupilCov", task, f"Alpha~{gaze_col}*Cond+Pupil",
                term, "beta_withPupil", b_pu)

        # Report pupil term itself
        if "Pupil_c" in res_pup.params.index:
            b_pup = res_pup.params["Pupil_c"]
            p_pup = res_pup.pvalues["Pupil_c"]
            print(f"    Pupil_c  β={b_pup:+.5f}  p={fmt_p(p_pup)}")
            log("4_PupilCov", task, f"Alpha~{gaze_col}*Cond+Pupil",
                "Pupil_c", "beta", b_pup)
            log("4_PupilCov", task, f"Alpha~{gaze_col}*Cond+Pupil",
                "Pupil_c", "p", p_pup)


# ============================================================
# %% 5. WITHIN-SUBJECT TRIAL CORRELATIONS
#        Per-subject r(alpha, gaze), Fisher z, group-level t-test
# ============================================================
print("\n" + "=" * 70)
print("  5. WITHIN-SUBJECT TRIAL CORRELATIONS (per-subject r → Fisher z → t)")
print("=" * 70)

for task in ["sternberg", "nback"]:
    trl = load_trials(task)

    for alpha_col, alpha_label in ALPHA_VARS_TRIAL[task]:
        if alpha_col not in trl.columns:
            print(f"  [{task}] {alpha_col}: column not found. Skipping.")
            continue

        for gaze_col in GAZE_VARS_TRIAL[task]:
            sub = trl[["ID", "Condition", alpha_col, gaze_col]].dropna()
            if sub.shape[0] < 50:
                continue

            # --- A) Collapsed across conditions ---
            rs_all = []
            for sid, g in sub.groupby("ID"):
                if g.shape[0] < 10:
                    continue
                r, _ = spearmanr(g[alpha_col], g[gaze_col])
                if np.isfinite(r):
                    rs_all.append(r)

            if len(rs_all) < 5:
                continue

            rs_all = np.array(rs_all)
            zs_all = np.arctanh(rs_all)  # Fisher z
            t_all, p_all = sp_stats.ttest_1samp(zs_all, 0)
            mean_r = np.tanh(np.mean(zs_all))  # back-transform mean z to r

            print(f"\n  [{task.upper()}] {gaze_col} × {alpha_col}  "
                  f"({alpha_label}, collapsed)")
            print(f"    N subjects = {len(rs_all)}")
            print(f"    Mean r (back-transformed) = {mean_r:.4f}")
            print(f"    Fisher z one-sample t: t({len(rs_all)-1}) = {t_all:.3f}, "
                  f"p = {fmt_p(p_all)}")

            log("5_TrialCorr", task,
                f"r({alpha_label},{gaze_col})_collapsed",
                "mean_r", "r", mean_r)
            log("5_TrialCorr", task,
                f"r({alpha_label},{gaze_col})_collapsed",
                "Fisher_z_ttest", "t", t_all)
            log("5_TrialCorr", task,
                f"r({alpha_label},{gaze_col})_collapsed",
                "Fisher_z_ttest", "p", p_all)

            # --- B) Per condition ---
            for cond in COND_ORDER[task]:
                csub = sub[sub["Condition"] == cond]
                rs_cond = []
                for sid, g in csub.groupby("ID"):
                    if g.shape[0] < 6:
                        continue
                    r, _ = spearmanr(g[alpha_col], g[gaze_col])
                    if np.isfinite(r):
                        rs_cond.append(r)
                if len(rs_cond) < 5:
                    continue

                rs_cond = np.array(rs_cond)
                zs_cond = np.arctanh(rs_cond)
                t_c, p_c = sp_stats.ttest_1samp(zs_cond, 0)
                mean_r_c = np.tanh(np.mean(zs_cond))

                print(f"    [{cond}] mean r = {mean_r_c:.4f}, "
                      f"t({len(rs_cond)-1}) = {t_c:.3f}, p = {fmt_p(p_c)}")

                log("5_TrialCorr", task,
                    f"r({alpha_label},{gaze_col})_{cond}",
                    "mean_r", "r", mean_r_c)
                log("5_TrialCorr", task,
                    f"r({alpha_label},{gaze_col})_{cond}",
                    "Fisher_z_ttest", "p", p_c)


# ============================================================
# %% 6. FORMAL MEDIATION  (Load → Gaze → Alpha)
#        Baron-Kenny steps + Sobel test
# ============================================================
print("\n" + "=" * 70)
print("  6. MEDIATION: Load → Gaze → Alpha  (Baron-Kenny + Sobel)")
print("=" * 70)

for task in ["sternberg", "nback"]:
    dat = load_long(task)
    # Numeric load (1, 2, 3) for continuous predictor
    cond_to_num = {c: i + 1 for i, c in enumerate(COND_ORDER[task])}
    dat["Load_num"] = dat["Condition"].map(cond_to_num).astype(float)

    for gaze_col in GAZE_VARS:
        if gaze_col not in dat.columns:
            continue
        sub = dat[["ID", "Load_num", "AlphaPower", gaze_col]].dropna().copy()
        if sub.shape[0] < 10:
            continue

        # Centre
        sub["Load_c"] = sub["Load_num"] - sub["Load_num"].mean()
        sub["Gaze_c"] = sub[gaze_col]  - sub[gaze_col].mean()

        # Step 1 (c path): Alpha ~ Load
        res_c = safe_fit("AlphaPower ~ Load_c", sub)
        # Step 2 (a path): Gaze ~ Load
        res_a = safe_fit(f"{gaze_col} ~ Load_c", sub)
        # Step 3 (b + c' path): Alpha ~ Load + Gaze
        res_bc = safe_fit("AlphaPower ~ Load_c + Gaze_c", sub)

        if any(r is None for r in [res_c, res_a, res_bc]):
            print(f"  [{task}] {gaze_col}: mediation models did not converge.")
            continue

        c  = res_c.params.get("Load_c", np.nan)
        pc = res_c.pvalues.get("Load_c", np.nan)
        a  = res_a.params.get("Load_c", np.nan)
        pa = res_a.pvalues.get("Load_c", np.nan)
        b  = res_bc.params.get("Gaze_c", np.nan)
        pb = res_bc.pvalues.get("Gaze_c", np.nan)
        cp = res_bc.params.get("Load_c", np.nan)  # c prime
        pcp = res_bc.pvalues.get("Load_c", np.nan)

        # Sobel test
        se_a = res_a.bse.get("Load_c", np.nan)
        se_b = res_bc.bse.get("Gaze_c", np.nan)
        ab = a * b
        se_ab = np.sqrt(a**2 * se_b**2 + b**2 * se_a**2)
        z_sobel = ab / se_ab if se_ab > 0 else np.nan
        p_sobel = 2 * (1 - norm.cdf(abs(z_sobel))) if np.isfinite(z_sobel) else np.nan

        # Proportion mediated
        prop_med = ab / c if abs(c) > 1e-12 else np.nan

        print(f"\n  [{task.upper()}] Mediator: {gaze_col}")
        print(f"    c  (Load→Alpha):           β={c:+.5f}  p={fmt_p(pc)}")
        print(f"    a  (Load→Gaze):            β={a:+.5f}  p={fmt_p(pa)}")
        print(f"    b  (Gaze→Alpha|Load):      β={b:+.5f}  p={fmt_p(pb)}")
        print(f"    c' (Load→Alpha|Gaze):      β={cp:+.5f}  p={fmt_p(pcp)}")
        print(f"    Indirect effect (a×b):     {ab:+.5f}")
        print(f"    Sobel z = {z_sobel:.3f}, p = {fmt_p(p_sobel)}")
        print(f"    Proportion mediated:       {prop_med:.3f}" if np.isfinite(prop_med) else
              f"    Proportion mediated:       NA (c ≈ 0)")

        for lbl, val in [("c", c), ("a", a), ("b", b), ("c_prime", cp),
                         ("ab_indirect", ab), ("sobel_z", z_sobel),
                         ("sobel_p", p_sobel), ("prop_mediated", prop_med)]:
            log("6_Mediation", task, f"Load→{gaze_col}→Alpha", lbl, "value", val)


# ============================================================
# %% Save all results
# ============================================================
results_df = pd.DataFrame(all_results)
out_csv = os.path.join(output_dir, "AOC_extra_glmms_results.csv")
results_df.to_csv(out_csv, index=False)
print(f"\n{'=' * 70}")
print(f"  All results saved → {out_csv}")
print(f"  Total entries: {len(all_results)}")
print(f"{'=' * 70}")


# ============================================================
# %% 7. VISUALIZATIONS
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid", font_scale=1.05)

fig_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/stats/extra_glmms"
os.makedirs(fig_dir, exist_ok=True)

print("\n" + "=" * 70)
print("  7. CREATING FIGURES")
print("=" * 70)

# ---- Palettes & helpers ----
_C = {
    "sternberg": ["#4C72B0", "#55A868", "#C44E52"],
    "nback":     ["#4C72B0", "#55A868", "#C44E52"],
}
_P = {
    "sternberg": dict(zip(["WM2", "WM4", "WM6"], _C["sternberg"])),
    "nback":     dict(zip(["NB1", "NB2", "NB3"], _C["nback"])),
}
_GL = {
    "GazeDeviationFull": "Gaze Deviation",
    "ScanPathLengthFull": "Scan Path Length",
    "MSRateFull": "MS Rate",
    "GazeDeviation": "Gaze Deviation",
    "ScanPathLength": "Scan Path Length",
    "MSRate": "MS Rate",
}


def _rc(data, x, y, order, palette, ax, title="", ylabel=""):
    """Raincloud-style: violin + inner box + jittered strip."""
    sns.violinplot(data=data, x=x, y=y, order=order, palette=palette,
                   ax=ax, inner="box", linewidth=0.7, saturation=0.45, cut=0)
    sns.stripplot(data=data, x=x, y=y, order=order, palette=palette,
                  ax=ax, size=2.5, alpha=0.3, jitter=True)
    if title:
        ax.set_title(title, fontweight="bold", fontsize=11)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_xlabel("")


# ------------------------------------------------------------------
# Fig 1 — Trial-level binned Gaze → Alpha scatter
# ------------------------------------------------------------------
print("  Fig 1 — Trial-level gaze–alpha binned scatter")
for task in ["sternberg", "nback"]:
    trl = load_trials(task)
    for acol, alab in ALPHA_VARS_TRIAL[task]:
        if acol not in trl.columns:
            continue
        gvs = GAZE_VARS_TRIAL[task]
        ng = len(gvs)
        fig, axes = plt.subplots(1, ng, figsize=(5.5 * ng, 4.5), sharey=True)
        if ng == 1:
            axes = [axes]
        for ax, gc in zip(axes, gvs):
            sub = trl[["ID", "Condition", acol, gc]].dropna()
            for j, cond in enumerate(COND_ORDER[task]):
                cs = sub[sub["Condition"] == cond].copy()
                if len(cs) < 30:
                    continue
                cs["_b"] = pd.qcut(cs[gc], 15, duplicates="drop")
                bm = cs.groupby("_b", observed=True).agg(
                    xm=(gc, "mean"), ym=(acol, "mean"),
                    ye=(acol, "sem")).reset_index()
                ax.errorbar(bm["xm"], bm["ym"], yerr=bm["ye"], fmt="o-",
                            color=_C[task][j], label=cond,
                            ms=4, lw=1.2, capsize=2, alpha=0.85)
            ax.set_xlabel(_GL.get(gc, gc))
            ax.set_title(_GL.get(gc, gc), fontweight="bold")
            ax.legend(fontsize=8)
        axes[0].set_ylabel(f"Alpha Power\n({alab})")
        fig.suptitle(f"{task.upper()} — Trial-level Gaze → Alpha",
                     fontweight="bold", fontsize=13)
        fig.tight_layout()
        fn = f"fig1_gaze_alpha_{task}_{acol}.png"
        fig.savefig(os.path.join(fig_dir, fn), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"    → {fn}")


# ------------------------------------------------------------------
# Fig 2 — Raincloud of key DVs by Condition
# ------------------------------------------------------------------
print("  Fig 2 — Raincloud: DVs by Condition")
for task in ["sternberg", "nback"]:
    dat = load_long(task)
    dvs = [d for d in key_dvs if d in dat.columns]
    nd = len(dvs)
    fig, axes = plt.subplots(1, nd, figsize=(4.5 * nd, 5))
    if nd == 1:
        axes = [axes]
    for ax, dv in zip(axes, dvs):
        _rc(dat, "Condition", dv, COND_ORDER[task], _P[task], ax,
            title=dv.replace("AlphaPower", "Alpha Power"), ylabel=dv)
    fig.suptitle(f"{task.upper()} — DVs by Condition",
                 fontweight="bold", fontsize=13, y=1.02)
    fig.tight_layout()
    fn = f"fig2_rainclouds_{task}.png"
    fig.savefig(os.path.join(fig_dir, fn), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    → {fn}")


# ------------------------------------------------------------------
# Fig 3 — Alpha by Task (joint model context)
# ------------------------------------------------------------------
print("  Fig 3 — Alpha by Task (joint model)")
_jf = []
for t in ["sternberg", "nback"]:
    d = load_long(t)
    d["Task"] = "Sternberg" if t == "sternberg" else "N-back"
    _jf.append(d[["ID", "Task", "AlphaPower"]].dropna())
_jdf = pd.concat(_jf, ignore_index=True)
fig, ax = plt.subplots(figsize=(4, 5))
_rc(_jdf, "Task", "AlphaPower", ["Sternberg", "N-back"],
    {"Sternberg": "#4C72B0", "N-back": "#C44E52"}, ax,
    title="Alpha Power by Task", ylabel="Alpha Power (a.u.)")
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "fig3_alpha_by_task.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("    → fig3_alpha_by_task.png")


# ------------------------------------------------------------------
# Fig 4 — Within-subject Gaze–Alpha correlation rainclouds
# ------------------------------------------------------------------
print("  Fig 4 — Within-subject Gaze–Alpha correlations")
for task in ["sternberg", "nback"]:
    trl = load_trials(task)
    for acol, alab in ALPHA_VARS_TRIAL[task]:
        if acol not in trl.columns:
            continue
        gvs = GAZE_VARS_TRIAL[task]
        ng = len(gvs)
        fig, axes = plt.subplots(1, ng, figsize=(5 * ng, 5), sharey=True)
        if ng == 1:
            axes = [axes]
        for ax, gc in zip(axes, gvs):
            sub = trl[["ID", "Condition", acol, gc]].dropna()
            rows = []
            for cond in COND_ORDER[task]:
                for sid, g in sub[sub["Condition"] == cond].groupby("ID"):
                    if len(g) < 6:
                        continue
                    r, _ = spearmanr(g[acol], g[gc])
                    if np.isfinite(r):
                        rows.append({"Condition": cond,
                                     "Fisher_z": np.arctanh(r)})
            if not rows:
                continue
            cdf = pd.DataFrame(rows)
            cdf["Condition"] = pd.Categorical(
                cdf["Condition"], categories=COND_ORDER[task], ordered=True)
            _rc(cdf, "Condition", "Fisher_z", COND_ORDER[task],
                _P[task], ax, title=_GL.get(gc, gc), ylabel="Fisher z")
            ax.axhline(0, color="black", ls="--", lw=0.8, alpha=0.5)
        fig.suptitle(
            f"{task.upper()} — Within-Subj Gaze × {acol} ({alab})",
            fontweight="bold", fontsize=12)
        fig.tight_layout()
        fn = f"fig4_within_corr_{task}_{acol}.png"
        fig.savefig(os.path.join(fig_dir, fn), dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"    → {fn}")


# ------------------------------------------------------------------
# Fig 5 — Mediation path coefficients
# ------------------------------------------------------------------
print("  Fig 5 — Mediation path coefficients")
_med = [r for r in all_results if r["Section"] == "6_Mediation"]
if _med:
    _mdf = pd.DataFrame(_med)
    _paths = ["c", "a", "b", "c_prime", "ab_indirect"]
    _plbl = {
        "c": "c (total)", "a": "a (Load→Gaze)",
        "b": "b (Gaze→Alpha|Load)", "c_prime": "c\u2032 (direct)",
        "ab_indirect": "a×b (indirect)",
    }
    _pcol = {
        "c": "#4C72B0", "a": "#55A868", "b": "#C44E52",
        "c_prime": "#8172B2", "ab_indirect": "#CCB974",
    }
    for task in ["sternberg", "nback"]:
        td = _mdf[_mdf["Task"] == task]
        if td.empty:
            continue
        td = td.copy()
        td["Mediator"] = td["Model"].str.extract(r"Load→(.+?)→Alpha")[0]
        mediators = td["Mediator"].dropna().unique()
        nm = len(mediators)
        fig, axes = plt.subplots(1, nm, figsize=(5 * nm, 4), sharey=True)
        if nm == 1:
            axes = [axes]
        for ax, med in zip(axes, mediators):
            ms = td[td["Mediator"] == med]
            labels, vals, cols = [], [], []
            for p in _paths:
                row = ms[ms["Term"] == p]
                if row.empty:
                    continue
                labels.append(_plbl.get(p, p))
                vals.append(float(row["Value"].iloc[0]))
                cols.append(_pcol.get(p, "#999"))
            ypos = list(range(len(labels)))
            ax.barh(ypos, vals, color=cols,
                    edgecolor="black", lw=0.5, alpha=0.85)
            ax.set_yticks(ypos)
            ax.set_yticklabels(labels)
            ax.axvline(0, color="black", lw=0.8)
            ax.set_xlabel("β")
            ax.set_title(_GL.get(med, med), fontweight="bold")
        fig.suptitle(f"{task.upper()} — Mediation Paths",
                     fontweight="bold", fontsize=13)
        fig.tight_layout()
        fn = f"fig5_mediation_{task}.png"
        fig.savefig(os.path.join(fig_dir, fn), dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"    → {fn}")

print(f"\n  All figures saved to → {fig_dir}")
