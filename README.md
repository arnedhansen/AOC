### README for AOC Study (Sternberg and N back)

Sternberg and N back tasks. Combined EEG and Eye Tracking ET analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. [https://doi.org/10.22541/au.172466871.17083913/v1](https://doi.org/10.22541/au.172466871.17083913/v1) Raw EEG data are provided via OpenNeuro. [https://openneuro.org/](https://openneuro.org/) Study materials and analysis code are distributed via OSF. [https://osf.io/xfujb/](https://osf.io/xfujb/)

The titles below correspond to the folder names. Apart from the Python and R scripts in `4_stats`, all files are MATLAB scripts.

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

`behavioral/` produces accuracy and RT features. `eeg/` produces EEG features used in the manuscript analyses. `gaze/` produces gaze deviation, microsaccade, fixation, saccade, scan path length, and pupil features. `AOC_master_matrix_nback.m` and `AOC_master_matrix_sternberg.m` merge outputs into `merged_data_*_nback.mat/.csv` and `merged_data_*_sternberg.mat/.csv`. `AOC_demographics.m` adds age, gender, handedness, and ocular dominance from the VP table.

## 3_visualization

Behavioral plots are in `behavioral/`. EEG visualization scripts used in the run all pipeline include `tfr/` and `topos/`. Gaze visualization scripts include `heatmap/`, `deviation/`, and `microsaccades/`. All visualization scripts read from `data/features/` and write figures to `figures/`.

## 4_stats

### Rainclouds (Python)

`AOC_stats_rainclouds.py` produces raincloud plots, repeated measures ANOVA, and mixed models for all variables; input: `merged_data_*_nback.csv` and `merged_data_*_sternberg.csv`. Python helpers (`stats_helpers`, `rainclouds_plotting_helpers`, `mixedlm_helpers`, `export_model_table`) come from [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions). Adapt `base_dir` and input paths in the script to your setup.

## Additional Files

### AOC_Master_RUN_ALL.m

Runs the main MATLAB analysis pipeline with absolute `run(...)` calls for both Windows and macOS Linux style paths. Update these paths before execution.

## Execution Order

Primary execution order follows stagewise scripts in `1_preprocessing`, `2_feature_extraction`, and `3_visualization`.

Recommended workflow:
1. Run preprocessing scripts in order: `1_cut`, `2_automagic`, `3_merge`, then `4_preprocessing`.
2. Run feature extraction scripts in `2_feature_extraction`.
3. Build final analysis matrices.
4. Generate manuscript figures.
5. Run confirmatory statistical models from `4_stats`.

## Setup for New Machines

All scripts are designed to run from a user editable project root.

Required setup actions:
1. Clone or download this repository.
2. Set local root paths in the setup script.
3. Verify availability of required MATLAB toolboxes and external packages.
4. Ensure all helper code is referenced from `functions`.
5. Ensure third party software is available via `toolboxes` or documented system installation.

## Run Instructions

### 1) Initial Setup
1. Open MATLAB in the repository root.
2. Run:
   1. `startup`
   2. `setup('AOC')`
3. Set environment variable `AOC_OSF_ROOT` if the project root should be forced explicitly.

### 2) Data Layout
1. Put processed inputs under:
   1. `data/raw`
   2. `data/features`
2. Participant metadata are read from `2_feature_extraction/behavioral/AOC_VP_List_anonymized.csv`, generated from the full source workbook after exclusion filtering.
3. Exclusion lists are provided in `2_feature_extraction/AOC_exclusion_participants.rtf` and `2_feature_extraction/AOC_exclusion_participants_info.rtf`.
4. Ensure required merged CSV files are available in `data/features` for stats scripts.

### 3) MATLAB Pipeline
1. Run scripts stagewise: preprocessing folders first, then `2_feature_extraction`, then `3_visualization`, then `splits`.

### 4) Confirmatory Statistics
1. Edit the CSV placeholders in:
   1. `4_stats/AOC_stats_confirmatory_models_nback.R`
   2. `4_stats/AOC_stats_confirmatory_models_sternberg.R`
2. Run each script with `Rscript`.

### 5) Raincloud Figures
1. Run `python 4_stats/AOC_stats_rainclouds.py`.
2. Output figures are written to `data/figures/stats/rainclouds`.

### Dependencies

`startup` and `setup('AOC')` handle path setup and subject configuration. Many scripts use absolute data roots such as `/Volumes/...`, `W:\...`, or user specific local paths. Update these to your local `data/` location before running.
