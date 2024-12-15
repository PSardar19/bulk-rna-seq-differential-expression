#!/bin/bash

# Define paths. Change the path of folders as required
BASE_DIR="/scratch_tmp/grp/msc_appbio/group8"
ALIGNMENT_DIR="$BASE_DIR/output"
COUNTS_DIR="$BASE_DIR/data/processed/counts_folder"

# Create counts_folder if it doesn't exist
mkdir -p "$COUNTS_DIR"

# Iterate over alignment directories
for dir in "$ALIGNMENT_DIR"/alignment_*; do
  if [ -d "$dir" ]; then
    SAMPLE_ID=$(basename "$dir" | sed 's/alignment_//')  # Extract sample ID
    COUNTS_FILE="$dir/counts.txt"
    if [ -f "$COUNTS_FILE" ]; then
      cp "$COUNTS_FILE" "$COUNTS_DIR/${SAMPLE_ID}_counts.txt"
      echo "Copied $COUNTS_FILE to $COUNTS_DIR/${SAMPLE_ID}_counts.txt"
    else
      echo "Warning: counts.txt not found in $dir"
    fi
  fi
done
