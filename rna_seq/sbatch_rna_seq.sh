#!/bin/bash

#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64  # Number of cores
#SBATCH --job-name=rnaseq
#SBATCH --output="rnaseq_%j.out"
#SBATCH --error="rnaseq_%j.err"
#SBATCH -A leiboff_lab           # Lab account (priority)
#SBATCH -p leiboff_lab           # Lab partition
#SBATCH -w cerebro               # Target cerebro node (optional)

# Define variables
env_name="rnaseq_packages"
script=/nfs5/BPP/Leiboff_Lab/Brian/scripts/misc_scripts/rna_seq.sh


# Load Conda environment management
source /local/cqls/software/x86_64/miniforge3/etc/profile.d/conda.sh

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
logmsg() { echo "$(timestamp): $*"; }

# Create Conda environment if it doesn't exist
if ! conda env list | awk '{print $1}' | grep -Fxq "${env_name}"; then
    logmsg "Creating Conda environment '${env_name}'..."

    conda create \
        --name "${env_name}" \
        --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority \
        fastp fastqc multiqc star samtools subread jq \
        -y || {
            logmsg "ERROR: Failed to create Conda environment: ${env_name}"
            exit 1
        }
else
    logmsg "Conda environment '${env_name}' already exists."
fi

# Activate the environment
conda activate "${env_name}" || {
    logmsg "Error: Failed to activate Conda environment '${env_name}'."
    exit 1
}

# Make allocated cores available
export threads="${SLURM_CPUS_PER_TASK}"

# Export environment details for reproducibility
conda env export -n "${env_name}" > "${env_name}.yaml"

# Run processing script
if [ -x "${script}" ]; then
    logmsg "Running script: ${script}"
    "${script}" && logmsg "Successfully ran ${script}."
else
    logmsg "Error: Script not found or not executable: ${script}"
    conda deactivate
    exit 1
fi

# Deactivate and clean up environment
conda deactivate

logmsg "SLURM job completed."
