#!/bin/bash

# Define paths
baseDirectory="$(dirname "$(pwd)")" # Root directory
genomeDirectory="$baseDirectory/data/raw/ncbi_dataset/data/GCF_000146045.2"
resultsDirectory="$baseDirectory/output/alignment_SRR1166443"

# Ensure the results directory exists
mkdir -p "$resultsDirectory"

# STAR alignment and generating the BAM file
singularity exec --bind "$resultsDirectory","$baseDirectory" "$baseDirectory/containers/star-alignment_latest.sif" \
    STAR --runThreadN 12 \
         --genomeDir "$genomeDirectory" \
         --readFilesCommand cat \
         --readFilesIn "$baseDirectory/data/processed/SRR1166443_trimmed.fastq" \
         --outFileNamePrefix "$resultsDirectory/theAlignment" \
         --outSAMtype BAM SortedByCoordinate

# Index the BAM file with samtools
singularity exec --bind "$resultsDirectory","$baseDirectory" "$baseDirectory/containers/samtools_latest.sif" \
    samtools index "$resultsDirectory/theAlignmentAligned.sortedByCoord.out.bam"

# Perform feature counting with subread
singularity exec --bind "$resultsDirectory","$baseDirectory" "$baseDirectory/containers/subread_latest.sif" \
    featureCounts -p -g gene_id -a "$genomeDirectory/genomic.gtf" -o "$resultsDirectory/counts.txt" "$resultsDirectory/theAlignmentAligned.sortedByCoord.out.bam"

