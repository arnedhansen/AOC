# AOC Multiverse Analysis

## Overview

Multiverse analysis for Sternberg and N-back tasks, examining:

1. **Gaze → Alpha relationship:** `alpha ~ gaze_value * Condition + (1|subjectID)`
2. **Condition effect on Alpha:** `alpha ~ Condition + (1|subjectID)`
3. **Aperiodic ~ Gaze:** `aperiodic_{exponent,offset} ~ gaze_value + (1|subjectID)` — tests whether eye movements predict the aperiodic spectral component (Tröndle & Langer, 2026)
4. **Aperiodic ~ Condition:** `aperiodic_{exponent,offset} ~ Condition + (1|subjectID)` — tests whether WM load affects aperiodic parameters

MATLAB builds long-format trial-level CSVs; R fits LMMs per universe and produces specification curve figures (600 dpi, publication-ready).

## Files

- **AOC_multiverse_prep.m** — Loads per-subject EEG and gaze data, computes all multiverse dimensions (Hanning tapers) at trial level, writes `multiverse_sternberg.csv` and `multiverse_nback.csv`. Compatible with Science Cloud (`ispc`) and Mac paths.
- **AOC_multiverse_prep_subject.m** — Subject-level data preparation. Recomputes everything from scratch on trial-averaged data: power spectra are averaged across trials before alpha extraction, FOOOF is fitted to averaged spectra (not per-trial), gaze metrics are per-trial averaged to subject means. Writes `multiverse_sternberg_subject.csv` and `multiverse_nback_subject.csv`. Much faster than trial-level (~75× fewer FOOOF calls).
- **AOC_multiverse_sternberg_analysis.R** — Reads Sternberg CSV, fits LMMs via multiverse package (Sarma et al., 2021), saves 5 result CSVs.
- **AOC_multiverse_sternberg_visualize.R** — Loads Sternberg result CSVs, produces 6 specification curve / forest plot figures (600 dpi).
- **AOC_multiverse_nback_analysis.R** — Same analysis pipeline for N-back.
- **AOC_multiverse_nback_visualize.R** — Same visualization pipeline for N-back (6 figures).
- **AOC_multiverse_sternberg_analysis_subject.R** — Subject-level Sternberg analysis (loads pre-aggregated subject-level CSV from `AOC_multiverse_prep_subject.m`).
- **AOC_multiverse_sternberg_visualize_subject.R** — Subject-level Sternberg visualization (figures with `_subject` suffix).
- **AOC_multiverse_nback_analysis_subject.R** — Subject-level N-back analysis.
- **AOC_multiverse_nback_visualize_subject.R** — Subject-level N-back visualization.

## Decisions (specification space — 7 dimensions, 1440 universes per task)

Ordered by when the decision occurs in the processing pipeline:

| Dimension      | Options |
|----------------|---------|
| Latency        | 0–500 ms, 0–1000 ms, 0–2000 ms, 1000–2000 ms (4) |
| Electrodes     | posterior, occipital (2) |
| FOOOF          | FOOOFed, non-FOOOFed (2) |
| EEG baseline   | raw, dB, percentage change from `[-0.5, -0.25] s` (3) |
| Alpha band     | canonical (8–14 Hz), IAF (2) |
| Gaze measure   | gaze deviation, gaze velocity, scan path length, BCEA (4) |
| Gaze baseline  | raw, dB, percentage change from `[-0.5, -0.25] s` (3) |

**Total:** 4 × 2 × 2 × 2 × 2 × 4 × 2 = **512 universes per task**.

## Processing details

- **Spectral method:** Hanning tapers throughout (via `ft_freqanalysis`, `method = 'mtmfft'`).
- **Power sources (hierarchy):** Pre-computed Hanning files (0–1 s, 0–2 s) → precomputed TFRs (0–500 ms, 1–2 s) → time-domain EEG via `ft_freqanalysis` (fallback for any remaining gaps).
- **FOOOF (trial-level):** `ft_freqanalysis_Arne_FOOOF` on raw time-domain data per trial. FOOOF result is baseline-independent (same value written for all three EEG baseline options).
- **FOOOF (subject-level):** `ft_freqanalysis_Arne_FOOOF` on all condition trials at once with `keeptrials = 'no'` — internally averages the spectrum across trials before FOOOF fitting. More comparable to standard subject-level pipelines.
- **Late retention window:** 1000–2000 ms captures the late retention interval.
- **Gaze deviation:** Mean Euclidean distance from screen center (400, 300) px per time window. The main registered gaze metric.
- **BCEA:** Bivariate Contour Ellipse Area (95% confidence).
- **Gaze baseline (dB):** `10 * log10(task / baseline)` from `[-0.5, -0.25] s`.
- **Gaze baseline (% change):** `(task - base) / |base| × 100` from `[-0.5, -0.25] s`.
- **EEG baseline (dB):** `10 * log10(task_spectrum / mean_baseline_spectrum)` from `[-0.5, -0.25] s`.
- **EEG baseline (% change):** `(task_spectrum - mean_baseline_spectrum) / |mean_baseline_spectrum| × 100` from `[-0.5, -0.25] s`. Baseline mean uses all trials for a stable estimate.
- **Aperiodic parameters:** Offset and exponent are extracted from `fooofparams.aperiodic_params` during FOOOF fitting. Only available for FOOOFed universes (NaN otherwise). The aperiodic multiverse is a reduced dimension space: only latency and electrodes matter (alpha band, EEG baseline, and FOOOF/no-FOOOF dimensions are not applicable). Motivated by Tröndle & Langer (2026) showing that ocular contamination increases aperiodic offset and steepens the slope.
- **Robust z-scoring:** Median + MAD, winsorized at ±3, applied per universe before model fitting.
- **Unstable universe filtering:** Universes with SE > 95th percentile are dropped.

## Output figures (per task)

All figures are 600 dpi PNG. Y-axis limits are symmetric and derived from the full extent of confidence intervals.

| Figure | Filename suffix | Description |
|--------|----------------|-------------|
| 1 | `_estimate` | `alpha ~ gaze` — Specification curve sorted by FOOOF then estimate (highest WM load), with analysis decision panel below |
| 2 | `_grouped` | `alpha ~ gaze (grouped)` — Specification curve ordered by decision hierarchy (FOOOF → Latency → Electrodes → ...), with analysis decision panel below |
| 3 | `_condition_alpha` | `alpha ~ condition` — Condition effect on alpha (highest vs. reference), EEG-only universes, sorted by processing-stage hierarchy |
| 4 | `_interaction` | `alpha ~ gaze × condition` — Interaction term (gaze × highest condition), all 7 dimensions, sorted by processing-stage hierarchy |
| 5 | `_condition_gaze` | `gaze ~ condition` — Condition effect on gaze (highest vs. reference), gaze-only universes |
| 6 | `_aperiodic` | Aperiodic component — combined forest plot. Top: exponent + offset ~ gaze (grouped by gaze measure, 4 × 16 universes). Bottom: exponent + offset ~ condition (8 labeled universes). |

### Figure details

- **Figures 1 & 2** show the `gaze_value` slope at the **highest WM load** (Set size 6 for Sternberg, 3-back for N-back). Universes are split by FOOOF (FOOOF left, No FOOOF right).
- **Figure 3** tests whether alpha power differs at the highest WM load vs. reference (Condition as factor). Only EEG dimensions are relevant, so it uses 5 EEG-only dimensions.
- **Figure 4** shows the gaze × condition interaction term from the full model, using all 7 dimensions.
- **Figure 5** tests whether gaze differs at the highest WM load vs. reference, using 3 gaze-only dimensions.
- **Figure 6** is a combined forest plot for the aperiodic component. Top section: aperiodic exponent and offset ~ gaze, grouped by gaze measure (Gaze Deviation, Gaze Velocity, Scan Path Length, BCEA) with 16 universes per group. Bottom section: aperiodic ~ condition with 8 labeled universes (latency × electrodes). Motivated by Tröndle & Langer (2026): ocular contamination increases aperiodic offset and steepens the slope, justifying why spectral parameterization matters for the main multiverse results.
- **Panel strip labels** are left-aligned and ordered by processing stage. Y-axis is labeled "Analysis Decision".

## Output CSVs (per task)

| File | Contents |
|------|----------|
| `multiverse_{task}.csv` | Trial-level data (from `AOC_multiverse_prep.m`) |
| `multiverse_{task}_subject.csv` | Subject-level data (from `AOC_multiverse_prep_subject.m`) |
| `multiverse_{task}_results.csv` | Interaction model results (`gaze_value` and `gaze_value:Condition` terms) |
| `multiverse_{task}_conditions_results.csv` | Per-condition simple model results (`gaze_value` slope per WM load) |
| `multiverse_{task}_condition_results.csv` | Condition → alpha results (highest condition factor contrast, EEG-only) |
| `multiverse_{task}_interaction_results.csv` | Gaze × condition interaction term (highest condition contrast) |
| `multiverse_{task}_condition_gaze_results.csv` | Condition → gaze results (highest condition factor contrast, gaze-only) |
| `multiverse_{task}_aperiodic_gaze_results.csv` | Aperiodic (exponent + offset) ~ gaze (4D: latency × electrodes × gaze × gaze baseline) |
| `multiverse_{task}_aperiodic_condition_results.csv` | Aperiodic (exponent + offset) ~ condition (2D: latency × electrodes) |

## Quality control

All control scripts live in `_controls/multiverse/`. Output data and figures go to the shared server.

| Script | What it checks | Figures | Data |
|--------|---------------|---------|------|
| `AOC_multiverse_control_fooof.m` | FOOOF R² goodness-of-fit (threshold 0.90). Scatterplots per subject/latency/electrode, R² vs aperiodic params, histogram. | `.../figures/controls/multiverse/FOOOF/` | `.../data/controls/multiverse/` |
| `AOC_control_nan_rate.R` | NaN rate per dimension for alpha and gaze_value. Heatmap of gaze NaN by gaze_measure × latency, bar chart per dimension. | `.../figures/controls/multiverse/NaN_rate/` | `.../data/controls/multiverse/NaN_rate/` |
| `AOC_control_trial_count.R` | Valid trial counts per subject × condition (in a reference universe). Bar charts overall and faceted by gaze measure. | `.../figures/controls/multiverse/trial_count/` | `.../data/controls/multiverse/trial_count/` |
| `AOC_control_winsorization.R` | Robust z-score winsorization rate (% clipped at ±3 MADs per universe). Histogram of rates, breakdown by dimension, worst-case universes. | `.../figures/controls/multiverse/winsorization/` | `.../data/controls/multiverse/winsorization/` |
| `AOC_control_baseline_sanity.R` | Baseline correction sanity: distribution of alpha/gaze by baseline type, extreme value rates for pct_change, Inf/NaN counts. | `.../figures/controls/multiverse/baseline_sanity/` | `.../data/controls/multiverse/baseline_sanity/` |
| `AOC_control_concordance.R` | Trial-level vs subject-level concordance: correlates LMM estimates across matched universes. Scatterplots with r values. | `.../figures/controls/multiverse/concordance/` | `.../data/controls/multiverse/concordance/` |

## Paths

| What | Path |
|------|------|
| Scripts | `/Users/Arne/Documents/GitHub/AOC/multiverse/` |
| CSV data | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/` |
| Control data | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/controls/multiverse/` (subfolders per check) |
| Figures (trial-level) | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse/` |
| Figures (subject-level) | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse/subject-level/` |
| Control figures | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/controls/multiverse/` (subfolders per check) |
| Control scripts | `/Users/Arne/Documents/GitHub/AOC/_controls/multiverse/` |

## Running

1. **MATLAB — trial-level (Science Cloud or Mac):**
   Run `AOC_multiverse_prep.m` once. Paths auto-detect via `ispc`: Windows uses `W:\Students\Arne\AOC`, Mac uses `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC`. CSVs are written to `data/features/`.

1b. **MATLAB — subject-level (Science Cloud or Mac):**
   Run `AOC_multiverse_prep_subject.m`. Same path auto-detection. Computes everything from scratch on trial-averaged data (does NOT require the trial-level CSVs). Writes `multiverse_*_subject.csv` to `data/features/`.

2. **R (analysis — fits models, saves result CSVs):**
   ```r
   source("AOC_multiverse_sternberg_analysis.R")
   source("AOC_multiverse_nback_analysis.R")
   ```

3. **R (visualization — loads result CSVs, generates figures):**
   ```r
   source("AOC_multiverse_sternberg_visualize.R")
   source("AOC_multiverse_nback_visualize.R")
   ```

   Paths are configurable via environment variables:
   - `AOC_MULTIVERSE_DIR` — CSV directory (default: `.../data/features`)
   - `AOC_MULTIVERSE_FIGURES` — Figure output (default: `.../figures/multiverse`)

   The visualization scripts can be re-run independently (e.g., to tweak aesthetics) without re-running the analysis.

4. **R (subject-level — loads pre-aggregated subject-level CSVs):**
   ```r
   source("AOC_multiverse_sternberg_analysis_subject.R")
   source("AOC_multiverse_sternberg_visualize_subject.R")
   source("AOC_multiverse_nback_analysis_subject.R")
   source("AOC_multiverse_nback_visualize_subject.R")
   ```
   Requires `multiverse_*_subject.csv` from `AOC_multiverse_prep_subject.m` (step 1b). Robust z-scoring is applied to subject means (not individual trials). Output CSVs and figures use `_subject` suffix.

5. **R (quality controls — requires trial-level CSVs and, for concordance, result CSVs):**
   ```r
   source("../_controls/multiverse/AOC_control_nan_rate.R")
   source("../_controls/multiverse/AOC_control_trial_count.R")
   source("../_controls/multiverse/AOC_control_winsorization.R")
   source("../_controls/multiverse/AOC_control_baseline_sanity.R")
   source("../_controls/multiverse/AOC_control_concordance.R")     # needs result CSVs from steps 2+4
   ```

## Science Cloud checklist

- **Drive:** `W:\Students\Arne\AOC` must exist with subject folders under `data/features/`.
- **Per subject, per task:** Required files in `eeg/`: `power_*_early_trials.mat`, `power_*_full_trials.mat`. Required in `gaze/`: `dataET_sternberg.mat` or `dataET_nback.mat`. Subject is skipped if any are missing.
- **0–500 ms and 1–2 s alpha:** From `tfr_*_trials.mat` if present, otherwise from time-domain EEG (`dataEEG_TFR_*.mat` or `dataEEG_*.mat`).
- **Channel labels:** From first subject's power file.
- **On path:** FieldTrip, `ft_freqanalysis_Arne_FOOOF` (for FOOOF), `detect_microsaccades` (for microsaccades; NaN on failure).

## Electrode sets

| Set | Channels |
|-----|----------|
| Posterior | Union of occipital + parietal (contains PO, O, I, or P-without-F) |
| Occipital | Label contains O or I |

## Notes

- **Microsaccades** in the 0–500 ms window are often suppressed post-stimulus; many trials will have NaN. The pipeline writes NaN, R drops them via `complete.cases()` and skips universes with < 10 valid rows.
- **Condition coding:** Sternberg: Set size 2/4/6; N-back: 1/2/3-back.
- **CSV columns:** Trial-level CSVs contain `aperiodic_offset` and `aperiodic_exponent` columns (NaN for non-FOOOFed universes). Subject-level CSVs contain the same columns (extracted from trial-averaged FOOOF).
- **Trial-level analysis:** Each row in the CSV is a single trial × universe combination.
- **Subject-level analysis:** Uses `AOC_multiverse_prep_subject.m` to compute features from scratch on trial-averaged data. EEG spectra are averaged across trials before alpha extraction (including FOOOF). Gaze metrics are computed per-trial and then averaged. This is fundamentally different from simply aggregating trial-level CSV results: FOOOF fitted to averaged spectra produces different (typically cleaner) results than averaging per-trial FOOOF values. Output CSVs have ~360k rows (vs ~31M trial-level). R scripts apply robust z-scoring to subject means (not individual trials).
- **Error bars** in all plots are 95% confidence intervals from the LMM fits.
- **Plot titles** use model notation (e.g., `alpha ~ gaze`, `alpha ~ condition`).
- **Baseline options:** Both dB and percentage change are computed for EEG and gaze. dB compresses extreme ratios logarithmically; % change is linear and more interpretable. Extreme values from small baselines are handled by robust z-scoring (median + MAD, winsorize ±3) in R.
