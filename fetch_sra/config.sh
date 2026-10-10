#!/bin/bash

# Project directory
proj_dir="/nfs5/BPP/Leiboff_Lab/Brian/rna_seq"

# Projects
projects=( "GSE137715" )

# Number of threads
#threads=8
threads=${SLURM_CPUS_PER_TASK:-8}

# SRA accessions
declare -A sra_groups=(
    ["GSE137715"]="SRR10145220 SRR10145221 SRR10145222 SRR10145223 SRR10145224 SRR10145225"
)