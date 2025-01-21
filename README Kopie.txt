**README for AOC Study (Sternberg & N-back)**

Sternberg and N-back tasks. Combined EEG and ET analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. https://doi.org/10.22541/au.172466871.17083913/v1

The following are short descriptions of what the files in the respective folders in this repository do. The titles below correspond to the folder names. Apart from the R scripts in ’4_stats’, all files are MATLAB scripts.

**paradigms**

The experimental paradigms ‚AOC_nback.m’ and ‚AOC_Sternberg.m‘ are executed using ‚master.m‘. The dependencies can be found in this folder as well.

**1_preprocessing**

EEG and ET files are cut by ‚doCutting.m‘ using the ‚cutData‘ function. The electrode layout from ‚standard_1005.elc‘ is used for this. The cut EEG and ET files are  preprocessed using Automagic and subsequently merged using ‚mergeData.m‘. ‚AOC_preprocessing_nback.m‘ and ‚AOC_preprocessing_sternberg.m‘ segment the merged data and save distinct data files for EEG and ET data for each condition.

**2_feature_extraction**

The behavioral files extract accuracy and reaction times. The eeg files extract alpha power over occipital electrodes. The gaze files extract gaze deviation. All files save structure arrays for each subject and each condition, which are combined in ,AOC_master_matrix_nback.m’ and ,AOC_master_matrix_sternberg.m’. Note, however, that the alpha power structure array is only created with ‚AOC_eeg_alpha_power_nback.m‘ and  ‚AOC_eeg_alpha_power_sternberg.m‘, respectively.

**3_visualization**

Behavioral: Visualizes accuracy and reaction times. Gaze: Visualizes gaze deviation (in euclidean distances) and microsaccades. Alpha Power: Calculates IAF and alpha power, saves the structure array and visualizes as boxplots, power spectra and topoplots.

**4_stats**

R scripts for the GLMMs.

**Additional Files**

Additional necessary files and functions: colorbrewer for colormaps, layANThead for the electrode layout of the ANT Neuro eeg cap, and shadedErrorBar for the power spectra plots. 
