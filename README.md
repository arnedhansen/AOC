**README for AOC Study (Sternberg & N-back)**

Sternberg and N-back tasks. Combined EEG and Eye-Tracking (ET) analysis of neural signatures of oculomotor control in the alpha band. Published in Psychophysiology as a Registered Report. https://doi.org/10.22541/au.172466871.17083913/v1

The following are short descriptions of what the files in the respective folders in this repository do. The titles below correspond to the folder names. Apart from the R scripts in ’4_stats’, all files are MATLAB scripts.

**paradigms**

The experimental paradigms ‚AOC_NBack.m’ and ‚AOC_Sternberg.m‘ are executed using ‚master.m‘. The dependencies can be found in the paradigm folder as well.

**1_preprocessing**


1_cut
EEG and ET files are cut by ‚AOC_Cutting.m‘ using ‚AOC_DataCuttingFunction.m‘. 

2_automagic
The cut EEG and ET files are preprocessed using Automagic (Pedroni et al., 2019). In the first step in Automagic, the bad channels are detected using the EEGLAB plugin clean_rawdata (Mullen et al., 2015). An electrode is defined as bad when recorded data from that electrode is correlated at less than .85 to an estimate based on neighboring electrodes. Furthermore, an electrode is defined as bad if it had more line noise relative to its signal than all other electrodes (4 standard deviations). Finally, if an electrode had a longer flat line than 5 s, it is considered bad. These bad channels will be removed from the original EEG data. The data will be filtered using a 0.1 Hz high pass filter (-6 dB cut off: 0.05 Hz) using the EEGLAB function pop_eegfiltnew (Widmann & Schröger, 2012). Line noise will be removed using a ZapLine method with a passband edge of 50 Hz (De Cheveigné, 2020), removing 7 power line components.

An independent component analysis (ICA) for ocular artifact correction will be applied with both the optimized ICA training (OPTICAT) function (Dimigen, 2020) and the pre-trained classifier ICLabel (Pion-Tonachini et al., 2019). OPTICAT enhances ocular artifact removal from EEG data, particularly saccadic spike activity. Initially, EEG data is high-pass filtered at 2 Hz, thereby preserving high-frequency components (> 40 Hz) that characterize saccadic spikes (Keren et al., 2010). Then, the contribution of saccadic spike activity in the EEG input to ICA is overweighted. This is achieved by extracting 30 ms long EEG segments around saccade onsets (−20 to +10 ms) identified via eye tracking, and appending these segments repeatedly to the EEG, resulting in EEG data double the length of the original data. An ICA is then trained on these overweighted data, and the resulting ICA labels are saved. Now, ICLabel is used on the original data by identifying artifact components and rating the probability of these artifacts being muscle, heart or eye activity, line noise or channel noise. All independent components receiving a probability rating by ICLabel of > 0.8 to be one of these non-brain artifacts are merged with the labels provided by OPTICAT. These components will be removed from the data and the remaining components will be back-projected on the original data. This is followed by interpolation of bad electrodes using the spherical interpolation method.

Afterwards, the quality of the data is automatically and objectively assessed in Automagic, thus increasing research reproducibility by having objective measures for data quality. The data of each individual block will be classified regarding data quality using the following exclusion criteria:

(1) The proportion of high-amplitude data points in the signal (> 30 μV) is larger than 0.2
(2) More than 20% of the time points show variance larger than 15 μV across electrodes
(3) 40% of the electrodes show high variance (15 μV)
(4) The proportion of bad electrodes is higher than 0.4

Any data file of a block exceeding any one of these criteria will be rated as bad and excluded from further analyses.

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
Rainclouds: R scripts for raincloud plots (combination of boxplot, scatterplot and density plot) for all variables.
GLMM: R scripts for the GLMMs.

**Additional Files**

Additional necessary files and functions: colorbrewer for colormaps, layANThead for the electrode layout of the ANT Neuro EEG cap, and shadedErrorBar for the power spectra plots. 
