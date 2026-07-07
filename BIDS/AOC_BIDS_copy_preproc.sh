#!/usr/bin/env bash
# Copy merged EEG+ET preprocessed files into BIDS derivatives/preproc.
# Excludes RestingEO files. Renames to BIDS derivative convention.
set -euo pipefail

SRC="/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/merged"
DST="/Volumes/g_psyplafor_methlab_data$/OCC/AOC_BIDS/derivatives/preproc"
LOG="/Volumes/g_psyplafor_methlab_data$/OCC/AOC_BIDS/preproc_copy.log"

mkdir -p "$DST"
: > "$LOG"

count=0
skipped=0

for subdir in "$SRC"/*/; do
    id=$(basename "$subdir")
    [[ "$id" =~ ^[0-9]+$ ]] || continue

    bids_sub=$(printf "sub-%04d" "$id")
    outdir="$DST/$bids_sub/eeg"
    mkdir -p "$outdir"

    shopt -s nullglob
    for f in "$subdir"*_merged.mat; do
        base=$(basename "$f")

        if [[ "$base" == *RestingEO* ]]; then
            skipped=$((skipped + 1))
            continue
        fi

        if [[ "$base" =~ _EEG_ET_Nback_block([0-9]+)_merged\.mat$ ]]; then
            run=$(printf "%02d" "${BASH_REMATCH[1]}")
            out="${bids_sub}_task-nback_run-${run}_desc-preproc_eeg.mat"
        elif [[ "$base" =~ _EEG_ET_Sternberg_block([0-9]+)_merged\.mat$ ]]; then
            run=$(printf "%02d" "${BASH_REMATCH[1]}")
            out="${bids_sub}_task-sternberg_run-${run}_desc-preproc_eeg.mat"
        else
            echo "SKIP unrecognized: $base" | tee -a "$LOG"
            skipped=$((skipped + 1))
            continue
        fi

        if [[ -f "$outdir/$out" ]]; then
            echo "EXISTS $outdir/$out" >> "$LOG"
            continue
        fi

        cp "$f" "$outdir/$out"
        count=$((count + 1))
        echo "COPIED $base -> $outdir/$out" >> "$LOG"
    done
done

echo "Done. Copied $count files, skipped $skipped entries." | tee -a "$LOG"
