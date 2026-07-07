# AOC BIDS migration scripts

Scripts for building and populating `/Volumes/g_psyplafor_methlab_data$/OCC/AOC_BIDS`.

## Feature tables (done)

OSF submission CSVs were converted to BIDS TSV in `derivatives/features/group/`:

- `group_task-nback_desc-mastermatrix_features.tsv`
- `group_task-sternberg_desc-mastermatrix_features.tsv`

Source:

- `_OSF_2607XX_RR-S2_Submission/4_stats/AOC_merged_data_nback.csv`
- `_OSF_2607XX_RR-S2_Submission/4_stats/AOC_merged_data_sternberg.csv`

A `participant_id` column (`sub-XXXX`) was prepended. Original `ID` was retained.

## Preprocessed data copy

```bash
bash BIDS/AOC_BIDS_copy_preproc.sh
```

Copies merged EEG+ET `.mat` files from the analysis `data/merged` folder into `derivatives/preproc/sub-XXXX/eeg/`, excluding `RestingEO` files. Progress is logged to `AOC_BIDS/preproc_copy.log`.

## Raw data conversion

```matlab
startup
setup('AOC')
addpath('/Users/Arne/Documents/GitHub/AOC/BIDS')
AOC_BIDS_convert_raw()
```

Test on one subject first:

```matlab
AOC_BIDS_convert_raw(struct('subjects', {{'301'}}, 'overwrite', true))
```

### Raw source recommendation

Use task-segmented `.mat` files in each subject folder, not `archive/*.cnt`.

| Source | Use for BIDS? | Reason |
|---|---|---|
| `*_task_EEG.mat` | Yes | Task-segmented EEG with channel locations |
| `*_task_ET.mat` | Yes | Synchronized eye tracking segments |
| `*_task.mat` / `*_NBack_block*_task.mat` | Yes | Behavioral logs |
| `archive/*.cnt` | No | Continuous recording; users would need to re-cut |
| `archive/*.edf` | No | EyeLink exports, not EEG |
| Training / Resting / Photodiode | No | Excluded from BIDS export |

Raw conversion writes EEGLAB `.set` files plus BIDS sidecars (`_eeg.json`, `_channels.tsv`, `_events.tsv`).
