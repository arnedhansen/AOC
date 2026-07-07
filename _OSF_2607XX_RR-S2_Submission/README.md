# AOC: Alpha Oculomotor Control

This repository contains the analysis code prepared for OSF release for the study:

Modulations of Posterior Alpha Power During Working Memory Co Vary With Task Dependent Eye Movement Patterns.

The project analyzes combined EEG and eye tracking data from Sternberg and N back paradigms to quantify behavioral performance, oculomotor measures, and posterior alpha activity during working memory.

## Study Context

The manuscript is a Stage 2 Registered Report in psychophysiology.  
Raw EEG data are provided via OpenNeuro.  
Study materials and analysis code are distributed via OSF.

OSF project: https://osf.io/xfujb/  
OpenNeuro: https://openneuro.org/

## Repository Scope

This OSF package contains the scripts required to reproduce the published analysis workflow from processed inputs onward.

Included:
1. MATLAB scripts for preprocessing outputs, feature extraction, visualization, and derived data assembly.
2. Statistical scripts for confirmatory model execution.
3. Local helper code under `functions`.
4. Required third party software references and local toolbox structure under `toolboxes`.

Not included:
1. Raw recording files.
2. Environment specific private paths.
3. Legacy or exploratory code not required for manuscript reproduction.

## Directory Structure

1. `1_preprocessing`  
   Final preprocessing scripts used after data cutting, Automagic processing, and merge steps.
2. `2_feature_extraction`  
   Behavioral, EEG, and gaze feature extraction plus matrix assembly.
3. `3_visualization`  
   Figure generation scripts for EEG and gaze outputs.
4. `4_stats`  
   Statistical model scripts for manuscript results.
5. `_controls`  
   Quality control and sensitivity scripts retained when required for interpretation.
6. `functions`  
   Project local helper functions used by MATLAB, Python, and R scripts.
7. `toolboxes`  
   Third party dependencies when redistribution is permitted, with install notes for external dependencies when redistribution is restricted.

## Execution Order

Primary execution order follows stagewise scripts in `1_preprocessing`, `2_feature_extraction`, and `3_visualization`.

Recommended workflow:
1. Run preprocessing scripts in order: `1_cut`, `2_automagic`, `3_merge`, then `4_preprocessing`.
2. Run feature extraction scripts in `2_feature_extraction`.
3. Build final analysis matrices.
4. Generate manuscript figures.
5. Run confirmatory statistical models from `4_stats`.

## Statistical Models

A minimal R script is provided in `4_stats` for confirmatory manuscript models.

Design principle:
1. Read the final CSV input file.
2. Fit predefined confirmatory models only.
3. Save model results in a transparent table format.

No exploratory analyses, plotting, or auxiliary processing are included in this minimal model runner.

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

### 5) Raincloud Figure Stage
1. Run `python 4_stats/AOC_stats_rainclouds.py`.
2. Output figures are written to `data/figures/stats/rainclouds`.

### 6) Notes
1. Scripts rely on vendored helpers from `functions`.
2. External packages not redistributed must be installed separately as documented in this README.

## Software Requirements

1. MATLAB with required signal processing and statistics functionality.
2. FieldTrip.
3. EEGLAB.
4. Python environment for selected analysis and figure scripts, where applicable.
5. R environment for confirmatory model scripts.

Dependency requirements are listed below in this README.


## Dependency Requirements

### MATLAB Runtime Dependencies
1. MATLAB release with support for table I O, mixed modeling, and plotting used by the scripts.
2. Statistics and Machine Learning Toolbox.
3. Signal Processing Toolbox.
4. FieldTrip toolbox. Redistribution is not included in this package.
5. EEGLAB toolbox for preprocessing conversion steps. Redistribution is not included in this package.

### Bundled MATLAB Helpers
1. Project helpers are vendored in `functions`.
2. Third party plotting helpers are bundled in `toolboxes`:
   1. `shadedErrorBar.m`
   2. `cbrewer.m`
   3. `layANThead.mat`

### Python Dependencies
1. Python 3.10 or newer recommended.
2. Required packages:
   1. `numpy`
   2. `pandas`
   3. `matplotlib`
   4. `seaborn`
   5. `scipy`
   6. `statsmodels`
3. Vendored Python helpers in `functions`:
   1. `stats_helpers.py`
   2. `rainclouds_plotting_helpers.py`

### R Dependencies
1. R 4.2 or newer recommended.
2. Required packages:
   1. `lme4`
   2. `lmerTest`
   3. `car`

### External Installation Notes
1. Install FieldTrip from [https://www.fieldtriptoolbox.org/download/](https://www.fieldtriptoolbox.org/download/).
2. Install EEGLAB from [https://sccn.ucsd.edu/eeglab/download.php](https://sccn.ucsd.edu/eeglab/download.php).
3. Place external toolboxes under `toolboxes` or add them to the MATLAB path before running.

## Reproducibility Notes

1. Scripts avoid machine specific hardcoded paths in the OSF release.
2. Output folders are created when absent.
3. Console messages use standardized task tags and subject progress reporting.
4. Script headers provide concise purpose and output declarations.

## Outputs

Main outputs include:
1. Processed feature tables for N back and Sternberg paradigms.
2. Publication figures for EEG and gaze analyses.
3. Confirmatory model summaries for manuscript reporting.

## Contact

Correspondence: arne.hansen@psychologie.uzh.ch
