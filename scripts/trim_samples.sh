#!/bin/bash
#SBATCH --job-name=trimmomatic_job       # Job name
#SBATCH --output=../logs/slurm-%j.out    # Standard output and error log

# Set project and directory variables
PROJECT_DIR="$(dirname "$(pwd)")"
INPUT_DIR=${1:-$PROJECT_DIR/data/raw}                  # First argument: Input directory (default: data/raw)
OUTPUT_DIR=${2:-$PROJECT_DIR/data/processed}           # Second argument: Output directory (default: data/processed)
THREADS=${3:-4}                                        # Third argument: Number of threads (default: 4)
HEADCROP=${4:-15}                                      # Fourth argument: Number of bases to crop (default: 15)

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Load Trimmomatic module (if required)
module load trimmomatic/0.39-gcc-13.2.0

# Process all FASTQ files in the input directory
for INPUT_FILE in "$INPUT_DIR"/*.fastq; do
    if [ -f "$INPUT_FILE" ]; then
        # Get the base name of the file (without path or extension)
        BASE_NAME=$(basename "$INPUT_FILE" .fastq)

        # Define the output file name
        OUTPUT_FILE="$OUTPUT_DIR/${BASE_NAME}_trimmed.fastq"

        # Run Trimmomatic for single-end reads
        echo "Processing $INPUT_FILE with HEADCROP:$HEADCROP"
        trimmomatic SE -threads "$THREADS" -phred33 "$INPUT_FILE" "$OUTPUT_FILE" HEADCROP:"$HEADCROP"

        # Check if the command was successful
        if [ $? -eq 0 ]; then
            echo "$INPUT_FILE trimmed successfully. Output: $OUTPUT_FILE"
        else
            echo "Trimmomatic failed for $INPUT_FILE. Check logs for details."
        fi
    else
        echo "No FASTQ files found in $INPUT_DIR."
    fi
done
