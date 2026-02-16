# AOC Multiverse Analysis

## Overview

Multiverse analysis for Sternberg and N-back tasks, examining:

1. **Gaze → Alpha relationship:** `alpha ~ gaze_value * Condition + (1|subjectID)`
2. **Condition effect on Alpha:** `alpha ~ Condition_numeric + (1|subjectID)`

MATLAB builds long-format trial-level CSVs; R fits LMMs per universe and produces specification curve figures (600 dpi, publication-ready).

## Files

- **AOC_multiverse_prep.m** — Loads per-subject EEG and gaze data, computes all multiverse dimensions (Hanning tapers), writes `multiverse_sternberg.csv` and `multiverse_nback.csv`. Compatible with Science Cloud (`ispc`) and Mac paths.
- **AOC_multiverse_sternberg.R** — Reads Sternberg CSV, fits LMMs, produces 3 specification curve figures + result CSVs.
- **AOC_multiverse_nback.R** — Same for N-back.

## Decisions (specification space — 7 dimensions, 1536 universes per task)

Ordered by when the decision occurs in the processing pipeline:

| Dimension      | Options |
|----------------|---------|
| Latency        | 0–500 ms, 0–1000 ms, 0–2000 ms, 1000–2000 ms (4) |
| Electrodes     | posterior, parietal, occipital (3) — "all" is computed but excluded from plotting |
| FOOOF          | FOOOFed, non-FOOOFed (2) |
| EEG baseline   | raw, dB-baselined `[-0.5, -0.25] s` (2) |
| Alpha band     | canonical (8–14 Hz), IAF (2) |
| Gaze measure   | scan path length, gaze velocity, microsaccades, BCEA (4) |
| Gaze baseline  | raw, percentage change from `[-0.5, -0.25] s` (2) |

**Plotted universes:** 4 × 3 × 2 × 2 × 2 × 4 × 2 = **1536 per task** (2048 total in CSV, "all" electrodes filtered out in R).

## Processing details

- **Spectral method:** Hanning tapers throughout (via `ft_freqanalysis`, `method = 'mtmfft'`).
- **Power sources (hierarchy):** Pre-computed Hanning files (0–1 s, 0–2 s) → precomputed TFRs (0–500 ms, 1–2 s) → time-domain EEG via `ft_freqanalysis` (fallback for any remaining gaps).
- **FOOOF:** `ft_freqanalysis_Arne_FOOOF` on raw time-domain data per trial. FOOOF result is baseline-independent (same value written for both raw and dB baseline options).
- **Late retention window:** 1000–2000 ms captures the late retention interval.
- **BCEA:** Bivariate Contour Ellipse Area (95% confidence).
- **Gaze baseline:** Percentage change: `(task - base) / |base| × 100` from `[-0.5, -0.25] s`.
- **EEG dB baseline:** Pre-computed `_bl` files where available, or `10*log10(task / mean_baseline_spectrum)` from `[-0.5, -0.25] s`.
- **Robust z-scoring:** Median + MAD, winsorized at ±3, applied per universe before model fitting.
- **Unstable universe filtering:** Universes with SE > 95th percentile are dropped.

## Output figures (per task)

All figures are 600 dpi PNG. Y-axis limits are symmetric and derived from the full extent of confidence intervals.

| Figure | Filename suffix | Description |
|--------|----------------|-------------|
| 1 | `_estimate` | `alpha ~ gaze` — Specification curve sorted by FOOOF then estimate (highest WM load), with analysis decision panel below |
| 2 | `_grouped` | `alpha ~ gaze (grouped)` — Specification curve ordered by decision hierarchy (FOOOF → Latency → Electrodes → ...), with analysis decision panel below |
| 3 | `_condition` | `alpha ~ condition` — Condition effect on alpha, using EEG-only universes, sorted purely by estimate (no FOOOF split) |

### Figure details

- **Figures 1 & 2** show the `gaze_value` slope at the **highest WM load** (Set size 6 for Sternberg, 3-back for N-back). The interaction model results are saved to CSV but not plotted. Universes are split by FOOOF (FOOOF left, No FOOOF right).
- **Figure 3** tests whether alpha power changes with increasing WM load (linear trend). Only EEG dimensions are relevant (gaze dimensions don't affect alpha), so it uses deduplicated EEG-only universes sorted purely by estimate (no FOOOF-first split). Panel shows the 5 EEG dimensions ordered by processing stage.
- **Panel strip labels** are left-aligned and ordered by processing stage. Y-axis is labeled "Analysis Decision".

## Output CSVs (per task)

| File | Contents |
|------|----------|
| `multiverse_{task}.csv` | Trial-level data (from MATLAB) |
| `multiverse_{task}_results.csv` | Interaction model results (`gaze_value` and `gaze_value:Condition` terms) |
| `multiverse_{task}_conditions_results.csv` | Per-condition simple model results (`gaze_value` slope per WM load) |
| `multiverse_{task}_condition_results.csv` | Condition → alpha results (linear WM load slope per EEG universe) |

## Paths

| What | Path |
|------|------|
| Scripts | `/Users/Arne/Documents/GitHub/AOC/multiverse/` |
| CSV data | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/features/` |
| Figures | `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/multiverse/` |

## Running

1. **MATLAB (Science Cloud or Mac):**
   Run `AOC_multiverse_prep.m` once. Paths auto-detect via `ispc`: Windows uses `W:\Students\Arne\AOC`, Mac uses `/Volumes/g_psyplafor_methlab$/Students/Arne/AOC`. CSVs are written to `data/features/`.

2. **R:**
   ```r
   source("AOC_multiverse_sternberg.R")
   source("AOC_multiverse_nback.R")
   ```
   Paths are configurable via environment variables:
   - `AOC_MULTIVERSE_DIR` — CSV directory (default: `.../data/features`)
   - `AOC_MULTIVERSE_FIGURES` — Figure output (default: `.../figures/multiverse`)

## Science Cloud checklist

- **Drive:** `W:\Students\Arne\AOC` must exist with subject folders under `data/features/`.
- **Per subject, per task:** Required files in `eeg/`: `power_*_early_trials.mat`, `power_*_full_trials.mat`. Required in `gaze/`: `dataET_sternberg.mat` or `dataET_nback.mat`. Subject is skipped if any are missing.
- **0–500 ms and 1–2 s alpha:** From `tfr_*_trials.mat` if present, otherwise from time-domain EEG (`dataEEG_TFR_*.mat` or `dataEEG_*.mat`).
- **Channel labels:** From first subject's power file.
- **On path:** FieldTrip, `ft_freqanalysis_Arne_FOOOF` (for FOOOF), `detect_microsaccades` (for microsaccades; NaN on failure).

## Electrode sets

| Set | Channels |
|-----|----------|
| Parietal | Label contains P (excluding F) |
| Posterior | Label contains PO, O, or I |
| Occipital | Label contains O or I |

Note: "All channels" is still computed in the MATLAB CSV but filtered out in R before plotting.

## Notes

- **Microsaccades** in the 0–500 ms window are often suppressed post-stimulus; many trials will have NaN. The pipeline writes NaN, R drops them via `complete.cases()` and skips universes with < 10 valid rows.
- **Condition coding:** Sternberg: Set size 2/4/6; N-back: 1/2/3-back.
- **Trial-level analysis:** Each row in the CSV is a single trial × universe combination.
- **Error bars** in all plots are 95% confidence intervals from the LMM fits.
- **Plot titles** use model notation (e.g., `alpha ~ gaze`, `alpha ~ condition`).
