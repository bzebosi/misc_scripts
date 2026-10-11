#!/bin/bash
# Copy this RNA-seq configuration file into your working directory
# and edit only the settings below.

#--------------------------------------------------------------------------------------------------------
# User settings
#--------------------------------------------------------------------------------------------------------

# Reference directory
ref_dir="/nfs5/BPP/Leiboff_Lab/Brian/ref_genome"

# Genomes to process
genomes=("B73v5")


# Genome FASTA URLs
declare -A genome_urls=(
    [B73v5]="https://download.maizegdb.org/Genomes/B73/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
)

# Annotation URLs
declare -A annotation_urls=(
    [B73v5]="https://download.maizegdb.org/Genomes/B73/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
)



declare -A project_dir=(
    [Ts3]="/nfs5/BPP/Leiboff_Lab/Brian/ts3_rnaseq"
    [Knauer]="/nfs5/BPP/Leiboff_Lab/Brian/rna_seq/Knauer"

)

# Projects to process
projects=(
    "Knauer"
)