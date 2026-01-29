# AOC Multiverse — Alpha ~ Gaze × Condition

## Overview

Multiverse analysis for Sternberg and N-back: **alpha ~ gaze_measure × condition + (1|subjectID)**.  
MATLAB builds long-format CSVs; R fits the LMM per universe and plots specification curve + panel (DALAS-style).

## Files

- **AOC_multiverse_prep.m** — Loads merged Sternberg/N-back data, builds decision grid, writes `multiverse_sternberg.csv` and `multiverse_nback.csv`. Compatible with Science Cloud (`ispc`) and Mac paths.
- **AOC_multiverse_sternberg.R** — Reads Sternberg CSV, fits LMM per universe, saves specification figure to AOC figures folder.
- **AOC_multiverse_nback.R** — Same for N-back.

## Decisions (specification space)

| Dimension    | Options |
|------------|----------|
| Electrodes | all, posterior, parietal, occipital |
| 1/f        | FOOOFed, non-FOOOFed |
| Latency    | 0–500 ms, 0–1000 ms, 0–2000 ms |
| Alpha      | canonical (8–14 Hz), IAF |
| Gaze       | gaze_density (Fixations), scan_path_length, gaze_velocity, microsaccades |

**Note:** The MATLAB script builds the table from per-subject power and gaze files and computes trial-level alpha (raw and FOOOF) and gaze (fixations, SPL, velocity, microsaccades) for all electrode/latency/alpha/gaze combinations, so each of the 192 universes has the correct alpha and gaze_value per trial.

## Model: alpha ~ gaze × condition + (1|subjectID)

- **Rationale:** Tests whether the gaze–alpha relationship differs by condition (WM load), with subject random intercepts.
- **Possible improvements:**
  1. **Random slopes:** `(gaze_value + Condition | subjectID)` or `(gaze_value | subjectID)` if you expect subject-specific gaze–alpha slopes.
  2. **Covariates:** Age, sex, or trial count if they might confound alpha or gaze.
  3. **Scale:** Consider scaling (e.g. z-score) `gaze_value` and/or `alpha` within subject or globally for interpretability and convergence.
  4. **Distribution:** If alpha is skewed, consider `log(alpha)` or a GLMM with appropriate family.
  5. **Multiple comparisons:** For multiverse, report proportion of significant universes and/or median estimate; avoid single-universe NHST as the main claim.
  6. **Convergence:** Use `lmerControl(optimizer = "bobyqa")` (or `optimizer = "nloptwrap"`) and check convergence per universe; flag or exclude failed fits.

## Running

1. **MATLAB (Science Cloud or Mac):**  
   Run `AOC_multiverse_prep.m` once. It builds the tables from per-subject files (no merged file needed). Paths: on Windows (Science Cloud) `base_data = W:\Students\Arne\AOC`, `base_features = W:\...\data\features`; on Mac the script uses `/Volumes/...`. CSVs are written to the script’s folder (`multiverse_sternberg.csv`, `multiverse_nback.csv`).

2. **R:**  
   Set working directory to this folder (or set env `AOC_MULTIVERSE_DIR` to this path), then:
   ```r
   source("AOC_multiverse_sternberg.R")
   source("AOC_multiverse_nback.R")
   ```
   Figure path is configurable: set env **`AOC_MULTIVERSE_FIGURES`** (e.g. on Science Cloud) or it defaults to  
   `'/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse'`. R creates the folder if needed.

## Science Cloud checklist (before running MATLAB)

- **Paths:** Script uses `ispc` → Windows paths `W:\Students\Arne\AOC` and `W:\...\data\features`. Ensure this drive/folder exists and contains subject folders.
- **Subject list:** From `setup('AOC')` if on path, else from folder names under `base_features`. At least one subject folder must exist.
- **Per subject, per task (Sternberg / N-back):**
  - **Required:** `base_features/<subjectID>/eeg/power_*_early_trials.mat`, `power_*_full_trials.mat`; `base_features/<subjectID>/gaze/dataET_sternberg.mat` or `dataET_nback.mat`. Subject is skipped if any of these are missing.
  - **0–500 ms alpha:** From precomputed `eeg/tfr_*_trials.mat` if present; otherwise from long time-domain EEG: script loads `dataEEG_TFR_*.mat` or `dataEEG_*.mat`, selects the 0–0.5 s epoch, and runs mtmfft to get trial-level power. So 0–500 ms is always computed when preprocessing output (dataEEG or dataEEG_TFR) exists.
- **Channel labels:** Taken from the first subject’s `power_stern_early_trials.mat` (field `powload2_early.label`). At least one subject must have Sternberg power so the script does not error at “Resolving channel sets”.
- **On path:** FieldTrip (for `ft_selectdata`, `ft_freqanalysis`). Optional: `ft_freqanalysis_Arne_FOOOF` for FOOOFed alpha; `detect_microsaccades` for microsaccade counts (otherwise NaNs in try/catch).
- **Output:** If no subject has both EEG and gaze for a task, the script errors with a clear “table is empty” message.

## Output figures

- **AOC_multiverse_sternberg.svg** / **.png** — Specification curve (estimates by universe, per term) + specification panel (decision options by universe).
- **AOC_multiverse_nback.svg** / **.png** — Same for N-back.

---

## Refined choices (from your answers)

1. **Electrode sets:** occipital = label contains O or I; parietal = contains P; posterior = contains PO, O, or I; all = all channels. Always bilateral.
2. **1/f:** raw power vs fitted model minus aperiodic (FOOOF – aperiodic).
3. **Latency:** aligned to stimulus; windows 0–500 ms, 0–1000 ms, 0–2000 ms after stimulus.
4. **Alpha:** canonical = 8–14 Hz; IAF as in AOC_eeg_fex_sternberg.m (peak in 8–14, power in IAF−4/+2 Hz with guards).
5. **Gaze density:** as in your pipeline (Fixations: count of L_fixation/R_fixation events per trial in analysis window).
6. **Gaze velocity:** same epoch as alpha; same epoch for EEG and ET always.
7. **Microsaccades:** per trial (rate or count per trial). In the **0–500 ms** window, microsaccades are often suppressed post-stimulus, so many trials may have NaN; the pipeline writes these as NaN, and R drops them via `complete.cases()` and skips universes with fewer than 10 valid rows—no error.
8. **Scan path length:** normalized (e.g. by duration or screen).
9. **Exclusion:** not necessary (preprocessed).
10. **Condition:** 3 WM loads (Sternberg: 2/4/6; N-back: 1/2/3); factor with default treatment coding.
11. **Covariates:** no Age/Sex.
12. **Level:** trial-level multiverse (MATLAB computes trial-level alpha and gaze per universe).
13. **N:** report N in panel.
14. **Figures:** same style as DALAS Specification_Channel_Group.R.
15. **Figure path:** configurable (env var).
16. **Ordering:** by estimate.
17. **Terms to plot:** see recommendation below.
18. **Robustness:** no.

---

## Recommendation for #19 (which terms to plot)

**Recommendation:** Plot **gaze_value** (main effect of gaze on alpha) and **gaze_value:Condition** (interaction: does the gaze–alpha slope differ by WM load?) in the specification curve. Optionally add **Condition** main effects (e.g. Condition4, Condition6 for Sternberg) in a separate small panel or same figure.

**Reasoning:** The main question is “does gaze predict alpha?” (gaze_value) and “does that relationship depend on load?” (interaction). Condition main effects are secondary (alpha by load) but standard to report. So: **primary = gaze_value; secondary = gaze_value:Condition; optional = Condition**. R scripts order universes by the **gaze_value** estimate and show the specification curve only for **gaze_value** and **gaze_value:Condition** (Condition main effects are not plotted).

---

## 20 questions (original list, for reference)

1. Electrode sets: occipital O/I, parietal P, posterior PO+O+I, all. Bilateral.
2. Posterior: single bilateral ROI.
3. 1/f: raw vs FOOOF minus aperiodic.
4. Latency: 0–500, 0–1000, 0–2000 ms after stimulus.
5. Canonical: 8–14 Hz; IAF as in your script.
6. IAF power: as in your script (IAF−4/+2 Hz).
7. Gaze density: Fixations per trial (as in your files).
8. Gaze velocity: same epoch as alpha.
9. Microsaccades: per trial.
10. SPL: normalized.
11. Exclusion: not necessary.
12. Condition: 3 WM loads.
13. No Age/Sex.
14. Trial-level.
15. Report N in panel.
16. DALAS-style figures.
17. Configurable figure path.
18. Order by estimate.
19. Plot gaze_value + interaction (see recommendation above).
20. No robustness.
