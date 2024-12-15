#!/bin/bash

# Define paths
baseDirectory="$(dirname "$(pwd)")"  # Root project directory
genomeDirectory="$baseDirectory/data/raw/ncbi_dataset/data/GCF_000146045.2"
fastaFile="$genomeDirectory/GCF_000146045.2_R64_genomic.fna"
gtfFile="$genomeDirectory/genomic.gtf"

# Create genome index using Singularity
singularity exec --bind "$genomeDirectory","$baseDirectory" "$baseDirectory/containers/star-alignment_latest.sif" \
    STAR --runThreadN 12 \
         --runMode genomeGenerate \
         --genomeDir "$genomeDirectory" \
         --genomeFastaFiles "$fastaFile" \
         --sjdbGTFfile "$gtfFile" \
         --sjdbGTFfeatureExon exon \
         --sjdbOverhang 100

echo "Genome indexing completed. Index files are stored in $genomeDirectory"
