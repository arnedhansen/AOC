### README for AOC Study (Sternberg & N-back)

Sternberg and N-back tasks. Combined EEG and Eye-Tracking (ET) analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. [https://doi.org/10.22541/au.172466871.17083913/v1](https://doi.org/10.22541/au.172466871.17083913/v1)

The titles below correspond to the folder names. Apart from the Python and R scripts in ’4_stats’, all files are MATLAB scripts.

## paradigms

The experimental paradigms ‚AOC_NBack.m’ and ‚AOC_Sternberg.m‘ are executed using ‚master.m‘. The dependencies can be found in the paradigms folder as well.

## 1_preprocessing

### 1_cut

Raw EEG and ET files are cut by `AOC_Cutting.m` using `AOC_DataCuttingFunction.m`.

### 2_automagic

The cut EEG and ET files are preprocessed using Automagic (Pedroni et al., 2019). In the first step in Automagic, the bad channels are detected using the EEGLAB plugin clean_rawdata (Mullen et al., 2015). An electrode is defined as bad when recorded data from that electrode is correlated at less than .85 to an estimate based on neighboring electrodes. Furthermore, an electrode is defined as bad if it had more line noise relative to its signal than all other electrodes (4 standard deviations). Finally, if an electrode had a longer flat line than 5 s, it is considered bad. These bad channels will be removed from the original EEG data. The data will be filtered using a 0.1 Hz high pass filter (-6 dB cut off: 0.05 Hz) using the EEGLAB function pop_eegfiltnew (Widmann & Schröger, 2012). Line noise will be removed using a ZapLine method with a passband edge of 50 Hz (De Cheveigné, 2020), removing 7 power line components.

An independent component analysis (ICA) for ocular artifact correction will be applied with both the optimized ICA training (OPTICAT) function (Dimigen, 2020) and the pre-trained classifier ICLabel (Pion-Tonachini et al., 2019). OPTICAT enhances ocular artifact removal from EEG data, particularly saccadic spike activity. Initially, EEG data is high-pass filtered at 2 Hz, thereby preserving high-frequency components (> 40 Hz) that characterize saccadic spikes (Keren et al., 2010). Then, the contribution of saccadic spike activity in the EEG input to ICA is overweighted. This is achieved by extracting 30 ms long EEG segments around saccade onsets (−20 to +10 ms) identified via eye tracking, and appending these segments repeatedly to the EEG, resulting in EEG data double the length of the original data. An ICA is then trained on these overweighted data, and the resulting ICA labels are saved. Now, ICLabel is used on the original data by identifying artifact components and rating the probability of these artifacts being muscle, heart or eye activity, line noise or channel noise. All independent components receiving a probability rating by ICLabel of > 0.8 to be one of these non-brain artifacts are merged with the labels provided by OPTICAT. These components will be removed from the data and the remaining components will be back-projected on the original data. This is followed by interpolation of bad electrodes using the spherical interpolation method.

Afterwards, the quality of the data is automatically and objectively assessed in Automagic, thus increasing research reproducibility by having objective measures for data quality. The data of each individual block will be classified regarding data quality using the following exclusion criteria:

(1) The proportion of high-amplitude data points in the signal (> 30 μV) is larger than 0.2
(2) More than 20% of the time points show variance larger than 15 μV across electrodes
(3) 40% of the electrodes show high variance (15 μV)
(4) The proportion of bad electrodes is higher than 0.4

Any data file of a block exceeding any one of these criteria will be rated as bad and excluded from further analyses.

### 3_merge

The EEG files that were preprocessed by Automagic are subsequently merged with their respective ET files using `AOC_mergeData.m`.

### 4_preprocessing

`AOC_preprocessing_nback.m` and `AOC_preprocessing_sternberg.m` segment the merged EEG+ET data into epochs around stimulus onset, convert to FieldTrip format, and write separate EEG and ET files per condition. N-back: loads 1/2/3; Sternberg: loads 2/4/6. Epoch windows: N-back [-1.5 2.5] s, Sternberg [-2 3.5] s. Run after merge; outputs go to `data/features/<subject>/eeg` and `.../gaze`.

## 2_feature_extraction

`behavioral/` → accuracy, RT; `eeg/` → alpha power at IAF over occipital channels, IAF, lateralization; `gaze/` → gaze deviation, microsaccade rate, fixations, saccades, scan path length, pupil. `AOC_master_matrix_nback.m` and `AOC_master_matrix_sternberg.m` merge these into `merged_data_*_nback.mat/.csv` and `merged_data_*_sternberg.mat/.csv`. `AOC_demographics.m` attaches age, gender, handedness, ocular dominance from the VP table. **Run order:** behavioral → eeg → gaze → demographics → master_matrix. The CSVs are the input for the Python raincloud/ANOVA scripts.

## 3_visualization

**Behavioral:** Accuracy and RT by condition (`behavioral/`). **EEG:** Power spectra (`powspctrm/`), TFRs—raw, FOOOF, baseline-relative (`tfr/`), ERPs (`erp/`), topoplots incl. baseline and raster (`topos/`), lateralization (`lateralization/`). **Gaze:** Deviation over time (`deviation/`), heatmaps and split heatmaps (`heatmap/`), microsaccades (`microsaccades/`), scan path length (`scanPathLength/`). **Interactions:** Alpha×velocity crosscorr, pupil×condition, TFR×scan path. All read from `data/features/`; figures to `figures/`. Run after feature extraction.

## 4_stats

### Omnibus (MATLAB, FieldTrip)

`AOC_omnibus_prep.m` loads single-subject TFR (FOOOF), applies baseline [-.5 -.25] s, builds high–low load contrasts (N-back: 3−1; Sternberg: 6−2) and grand averages, saves `omnibus_data.mat`. `AOC_omnibus.m` loads that, extracts posterior alpha spectra, runs cluster-based permutation (F-stat for load, t for omnibus), and produces ROI raincloud/box plots with paired t-tests.

### Rainclouds (Python)

`AOC_stats_glmms_rainclouds.py` produces raincloud plots, repeated-measures ANOVA, and mixed models for all variables; input: `merged_data_*_nback.csv` and `merged_data_*_sternberg.csv`. Python helpers (`stats_helpers`, `rainclouds_plotting_helpers`, `mixedlm_helpers`, `export_model_table`) come from [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions). Adapt `base_dir` and input paths in the script to your setup.

## Additional Files

### AOC_MASTER_ANALYSIS.m

Runs the full MATLAB pipeline (merge → 4_preprocessing → _controls → feature extraction → visualization → omnibus_prep → omnibus) with try/catch and a log. Set `basePath` to your repo root. It does *not* run: 1_cut, 2_automagic, or the Python stats.

### _controls

Optional checks: ET calibration/validation, baseline effects, trial exclusions, missing data, paradigm durations, recording order. Run after 4_preprocessing; see `AOC_MASTER_ANALYSIS.m` for placement.

### Dependencies

`startup` and `setup('AOC')` (paths, subject list; `setup` is project-specific). For plots: `cbrewer` (colormaps), `layANThead` (ANT Neuro cap layout), `shadedErrorBar` (power spectra). Many scripts hardcode data roots (e.g. `/Volumes/...` or `W:\...`); change these to your `data/` location.

## multiverse

Multiverse analysis for Sternberg and N-back tasks, examining:

1. **Gaze → Alpha relationship:** `alpha ~ gaze_value * Condition + (1|subjectID)`
2. **Condition effect on Alpha:** `alpha ~ Condition + (1|subjectID)`
3. **Aperiodic ~ Gaze:** `aperiodic_{exponent,offset} ~ gaze_value + (1|subjectID)` — tests whether eye movements predict the aperiodic spectral component (Tröndle & Langer, 2026)
4. **Aperiodic ~ Condition:** `aperiodic_{exponent,offset} ~ Condition + (1|subjectID)` — tests whether WM load affects aperiodic parameters

MATLAB builds long-format trial-level CSVs; R fits LMMs per universe and produces specification curve figures (600 dpi, publication-ready).

### Files

- **AOC_multiverse_prep.m** — Loads per-subject EEG and gaze data, computes all multiverse dimensions (Hanning tapers) at trial level, writes `multiverse_sternberg.csv` and `multiverse_nback.csv`. Compatible with Science Cloud (`ispc`) and Mac paths.
- **AOC_multiverse_prep_fooof_only.m** — Fast trial-level FOOOF refresh that reuses existing trial-level multiverse CSVs and recomputes only FOOOF-derived fields (`alpha` for FOOOFed rows plus `aperiodic_offset`/`aperiodic_exponent`), then writes mode-tagged CSVs.
- **AOC_multiverse_prep_subject.m** — Subject-level data preparation. Recomputes everything from scratch on trial-averaged data: power spectra are averaged across trials before alpha extraction, FOOOF is fitted to averaged spectra (not per-trial), gaze metrics are per-trial averaged to subject means. Writes `multiverse_sternberg_subject.csv` and `multiverse_nback_subject.csv`. Much faster than trial-level (~75× fewer FOOOF calls).
- **AOC_multiverse_sternberg_analysis.R** — Reads Sternberg CSV, fits LMMs via multiverse package (Sarma et al., 2021), saves 5 result CSVs.
- **AOC_multiverse_sternberg_visualize.R** — Loads Sternberg result CSVs, produces 6 specification curve / forest plot figures (600 dpi).
- **AOC_multiverse_nback_analysis.R** — Same analysis pipeline for N-back.
- **AOC_multiverse_nback_visualize.R** — Same visualization pipeline for N-back (6 figures).
- **AOC_multiverse_sternberg_analysis_subject.R** — Subject-level Sternberg analysis (loads pre-aggregated subject-level CSV from `AOC_multiverse_prep_subject.m`).
- **AOC_multiverse_sternberg_visualize_subject.R** — Subject-level Sternberg visualization (figures with `_subject` suffix).
- **AOC_multiverse_nback_analysis_subject.R** — Subject-level N-back analysis.
- **AOC_multiverse_nback_visualize_subject.R** — Subject-level N-back visualization.

### Decisions (specification space — 7 dimensions, 1440 universes per task)

Ordered by when the decision occurs in the processing pipeline:


| Dimension     | Options                                                   |
| ------------- | --------------------------------------------------------- |
| Latency       | 0–500 ms, 0–1000 ms, 0–2000 ms, 1000–2000 ms (4)          |
| Electrodes    | posterior, occipital (2)                                  |
| FOOOF         | FOOOFed, non-FOOOFed (2)                                  |
| EEG baseline  | raw, dB, percentage change from `[-0.5, -0.25] s` (3)     |
| Alpha band    | canonical (8–14 Hz), IAF (2)                              |
| Gaze measure  | gaze deviation, gaze velocity, scan path length, BCEA (4) |
| Gaze baseline | raw, dB, percentage change from `[-0.5, -0.25] s` (3)     |


**Total:** 4 × 2 × 2 × 2 × 2 × 4 × 2 = **512 universes per task**.

### Processing details

- **Spectral method:** Hanning tapers throughout (via `ft_freqanalysis`, `method = 'mtmfft'`).
- **Power sources (hierarchy):** Pre-computed Hanning files (0–1 s, 0–2 s) → precomputed TFRs (0–500 ms, 1–2 s) → time-domain EEG via `ft_freqanalysis` (fallback for any remaining gaps).
- **FOOOF (trial-level):** `ft_freqanalysis_Arne_FOOOF` on raw time-domain data per trial. In `AOC_multiverse_prep.m` this now has a mode toggle: `singleFFT` (legacy, one spectrum from the full window), `welch500_50` (500 ms segments with 50% overlap, averaged before FOOOF; default), or `BOTH` (runs both modes in one script execution and writes mode-tagged outputs). FOOOF result remains baseline-independent (same value written for all three EEG baseline options).
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

### Output figures (per task)

All figures are 600 dpi PNG. Y-axis limits are symmetric and derived from the full extent of confidence intervals.


| Figure | Filename suffix    | Description                                                                                                                                                                         |
| ------ | ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1      | `_estimate`        | `alpha ~ gaze` — Specification curve sorted by FOOOF then estimate (highest WM load), with analysis decision panel below                                                            |
| 2      | `_grouped`         | `alpha ~ gaze (grouped)` — Specification curve ordered by decision hierarchy (FOOOF → Latency → Electrodes → ...), with analysis decision panel below                               |
| 3      | `_condition_alpha` | `alpha ~ condition` — Condition effect on alpha (highest vs. reference), EEG-only universes, sorted by processing-stage hierarchy                                                   |
| 4      | `_interaction`     | `alpha ~ gaze × condition` — Interaction term (gaze × highest condition), all 7 dimensions, sorted by processing-stage hierarchy                                                    |
| 5      | `_condition_gaze`  | `gaze ~ condition` — Condition effect on gaze (highest vs. reference), gaze-only universes                                                                                          |
| 6      | `_aperiodic`       | Aperiodic component — combined forest plot. Top: exponent + offset ~ gaze (grouped by gaze measure, 4 × 16 universes). Bottom: exponent + offset ~ condition (8 labeled universes). |


### Notes

- **Microsaccades** in the 0–500 ms window are often suppressed post-stimulus; many trials will have NaN. The pipeline writes NaN, R drops them via `complete.cases()` and skips universes with < 10 valid rows.
- **Condition coding:** Sternberg: Set size 2/4/6; N-back: 1/2/3-back.
- **CSV columns:** Trial-level CSVs contain `aperiodic_offset` and `aperiodic_exponent` columns (NaN for non-FOOOFed universes). Subject-level CSVs contain the same columns (extracted from trial-averaged FOOOF).
- **Trial-level analysis:** Each row in the CSV is a single trial × universe combination.
- **Subject-level analysis:** Uses `AOC_multiverse_prep_subject.m` to compute features from scratch on trial-averaged data. EEG spectra are averaged across trials before alpha extraction (including FOOOF). Gaze metrics are computed per-trial and then averaged. This is fundamentally different from simply aggregating trial-level CSV results: FOOOF fitted to averaged spectra produces different (typically cleaner) results than averaging per-trial FOOOF values. Output CSVs have ~360k rows (vs ~31M trial-level). R scripts apply robust z-scoring to subject means (not individual trials).
- **Error bars** in all plots are 95% confidence intervals from the LMM fits.
- **Plot titles** use model notation (e.g., `alpha ~ gaze`, `alpha ~ condition`).
- **Baseline options:** Both dB and percentage change are computed for EEG and gaze. dB compresses extreme ratios logarithmically; % change is linear and more interpretable. Extreme values from small baselines are handled by robust z-scoring (median + MAD, winsorize ±3) in R.

