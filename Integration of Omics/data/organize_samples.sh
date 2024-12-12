#!/bin/bash

# Set the raw data directory where the original files are located
RAW_DATA_PATH="/work/Spontan_AF_project/raw-data"

# Loop through each unique sample name by extracting the base name without _R1 or _R2 suffix
for file in ${RAW_DATA_PATH}/*_R1.fq.qz; do
    # Extract sample name by removing the _R1 suffix
    SAMPLE_NAME=$(basename "$file" | sed 's/_R1.fq.qz//')
    
    # Create a subdirectory for each sample (if it doesn't exist)
    SAMPLE_DIR="${RAW_DATA_PATH}/${SAMPLE_NAME}"
    if [ ! -d "$SAMPLE_DIR" ]; then
        mkdir "$SAMPLE_DIR"
    fi
    
    # Move the R1 and R2 files into the corresponding subdirectory
    mv "${RAW_DATA_PATH}/${SAMPLE_NAME}_R1.fq.qz" "$SAMPLE_DIR/"
    mv "${RAW_DATA_PATH}/${SAMPLE_NAME}_R2.fq.qz" "$SAMPLE_DIR/"
    
    echo "Moved ${SAMPLE_NAME}_R1.fq.qz and ${SAMPLE_NAME}_R2.fq.qz to $SAMPLE_DIR"
done
