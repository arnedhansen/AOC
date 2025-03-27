**README for AOC Study (Sternberg & N-back)**

Sternberg and N-back tasks. Combined EEG and Eye-Tracking (ET) analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. https://doi.org/10.22541/au.172466871.17083913/v1

The following are short descriptions of what the files in the respective folders in this repository do. The titles below correspond to the folder names. Apart from the R scripts in ’4_stats’, all files are MATLAB scripts.

**paradigms**

The experimental paradigms ‚AOC_NBack.m’ and ‚AOC_Sternberg.m‘ are executed using ‚master.m‘. The dependencies can be found in the paradigm folder as well.

**1_preprocessing**


1_cut
EEG and ET files are cut by ‚AOC_Cutting.m‘ using ‚AOC_DataCuttingFunction.m‘. The electrode layout from ‚standard_1005.elc‘ is used for this. 

2_automagic
The cut EEG and ET files are preprocessed using Automagic (Pedroni et al., 2019).

3_merge
The EEG files that were preprocessed by Automagic are subsequently merged with their respective ET files using ‚mergeData.m‘. 

4_preprocessing
‚AOC_preprocessing_nback.m‘ and ‚AOC_preprocessing_sternberg.m‘ segment the merged data and save distinct data files for EEG and ET data for each condition.

**2_feature_extraction**

The behavioral files extract accuracy and reaction times. The EEG files extract alpha power at IAF over occipital electrodes. The gaze files extract gaze deviation and microsaccade rate. All files save structure arrays for each subject and each condition, which are combined in ,AOC_master_matrix_nback.m’ and ,AOC_master_matrix_sternberg.m’.

**3_visualization**

EEG: Visualizes IAF and alpha power in power spectra, topoplots, and TFR.
Gaze: Visualizes gaze deviation in heatmaps. 

**4_stats**

Overview: Creates boxplots and percentage change barplots for all variables of interest. Visualizes associations.
Rainclouds: R scripts for raincloud plots for all variables.
GLMM: R scripts for the GLMMs.

**Additional Files**

Additional necessary files and functions: colorbrewer for colormaps, layANThead for the electrode layout of the ANT Neuro eeg cap, and shadedErrorBar for the power spectra plots. 
