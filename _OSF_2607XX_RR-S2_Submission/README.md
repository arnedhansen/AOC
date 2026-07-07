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

Primary execution order is defined by:

`AOC_Master_RUN_ALL.m`

Recommended workflow:
1. Run preprocessing outputs and feature extraction in MATLAB.
2. Build final analysis matrices.
3. Generate manuscript figures.
4. Run confirmatory statistical models from `4_stats`.

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
2. Ensure required merged CSV files are available in `data/features` for stats scripts.

### 3) MATLAB Pipeline
1. Run `AOC_Master_RUN_ALL.m` to execute preprocessing through visualization and split analyses in sequence.

### 4) Confirmatory Statistics
1. Edit the CSV placeholders in:
   1. `4_stats/AOC_stats_confirmatory_models_nback_minimal.R`
   2. `4_stats/AOC_stats_confirmatory_models_sternberg_minimal.R`
2. Run each script with `Rscript`.

### 5) Raincloud Figure Stage
1. Run `python 4_stats/AOC_stats_glmms_rainclouds.py`.
2. Output figures are written to `data/figures/stats/rainclouds`.

### 6) Notes
1. Scripts rely on vendored helpers from `functions`.
2. External packages not redistributed must be installed separately as listed in `DEPENDENCY_MANIFEST.md`.

## Software Requirements

1. MATLAB with required signal processing and statistics functionality.
2. FieldTrip.
3. EEGLAB.
4. Python environment for selected analysis and figure scripts, where applicable.
5. R environment for confirmatory model scripts.

Exact versions should be documented in the repository level dependency manifest in this OSF package.

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
