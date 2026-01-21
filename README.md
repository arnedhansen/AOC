### README for AOC Study (Sternberg & N-back)

Sternberg and N-back tasks. Combined EEG and Eye-Tracking (ET) analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. https://doi.org/10.22541/au.172466871.17083913/v1

The following are short descriptions of what the files in the respective folders in this repository do. The titles below correspond to the folder names. Apart from the Python scripts in ’4_stats’, all files are MATLAB scripts.

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
`AOC_omnibus_prep.m` loads single-subject TFR (FOOOF), applies baseline [-.5 -.25] s, builds high–low load contrasts (N-back: 3−1; Sternberg: 6−2) and grand averages, saves `omnibus_data.mat`. `AOC_omnibus.m` loads that, extracts posterior alpha spectra, runs cluster-based permutation (F-stat for load, t for omnibus), and produces ROI raincloud/box plots with paired t-tests. **Run order:** omnibus_prep → omnibus.

### Rainclouds (Python)
`AOC_stats_glmms_rainclouds.py` produces raincloud plots, repeated-measures ANOVA, and mixed models for all variables; input: `merged_data_*_nback.csv` and `merged_data_*_sternberg.csv`. Python helpers (`stats_helpers`, `rainclouds_plotting_helpers`, `mixedlm_helpers`, `export_model_table`) come from [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions). Adapt `base_dir` and input paths in the script to your setup.

## Additional Files

### AOC_MASTER_ANALYSIS.m
Runs the full MATLAB pipeline (merge → 4_preprocessing → _controls → feature extraction → visualization → omnibus_prep → omnibus) with try/catch and a log. Set `basePath` to your repo root. It does *not* run: 1_cut, 2_automagic, or the Python stats.

### _controls
Optional checks: ET calibration/validation, baseline effects, trial exclusions, missing data, paradigm durations, recording order. Run after 4_preprocessing; see `AOC_MASTER_ANALYSIS.m` for placement.

### Dependencies
`startup` and `setup('AOC')` (paths, subject list; `setup` is project-specific). For plots: `cbrewer` (colormaps), `layANThead` (ANT Neuro cap layout), `shadedErrorBar` (power spectra). Many scripts hardcode data roots (e.g. `/Volumes/...` or `W:\...`); change these to your `data/` location. 
