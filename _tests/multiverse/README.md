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

**Note:** Current MATLAB script uses the *existing* merged matrices (one electrode set, one latency, one alpha/FOOOF). Alpha and gaze values are repeated across electrode/latency/FOOOF/alpha-type universes; only the *gaze measure* varies the data (Fixations, ScanPathLength, GazeVelocity, MSRate). To populate all 192 universes with distinct alpha/gaze, run full feature extraction for each combination (or add that logic to this script).

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
   Open `AOC_multiverse_prep.m`, ensure paths point to your `merged_data_*_trials.mat` and per-subject `power_*_early_trials.mat` / `power_*_full_trials.mat` under `base_features/<subjectID>/eeg/`. Run the script. CSVs are written next to the script.

2. **R:**  
   Set working directory to this folder (or set env `AOC_MULTIVERSE_DIR` to this path), then:
   ```r
   source("AOC_multiverse_sternberg.R")
   source("AOC_multiverse_nback.R")
   ```
   Figure path is configurable: set env **`AOC_MULTIVERSE_FIGURES`** (e.g. on Science Cloud) or it defaults to  
   `'/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/tests/multiverse'`. R creates the folder if needed.

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
7. **Microsaccades:** per trial (rate or count per trial).
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
