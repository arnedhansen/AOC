# AOC OSF Dependency Manifest

## MATLAB Runtime Dependencies
1. MATLAB release with support for table I O, mixed modeling, and plotting used by the scripts.
2. Statistics and Machine Learning Toolbox.
3. Signal Processing Toolbox.
4. FieldTrip toolbox. Redistribution is not included in this package.
5. EEGLAB toolbox for preprocessing conversion steps. Redistribution is not included in this package.

## Bundled MATLAB Helpers
1. Project helpers are vendored in `functions`.
2. Third party plotting helpers are bundled in `toolboxes`:
   1. `shadedErrorBar.m`
   2. `cbrewer.m`
   3. `layANThead.mat`

## Python Dependencies
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

## R Dependencies
1. R 4.2 or newer recommended.
2. Required packages:
   1. `lme4`
   2. `lmerTest`
   3. `car`

## External Installation Notes
1. Install FieldTrip from [https://www.fieldtriptoolbox.org/download/](https://www.fieldtriptoolbox.org/download/).
2. Install EEGLAB from [https://sccn.ucsd.edu/eeglab/download.php](https://sccn.ucsd.edu/eeglab/download.php).
3. Place external toolboxes under `toolboxes` or add them to the MATLAB path before running.
