#!/bin/bash

# Load configuration
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# -------------------------------------------------------------------------------------------------------
# log functions for updates and status
# -------------------------------------------------------------------------------------------------------
logmsg() { echo "$(date '+%Y-%m-%d %H:%M:%S'): $*" ; }


# -------------------------------------------------------------------------------------------------------
# Create one or more directories
# -------------------------------------------------------------------------------------------------------
create_dir() {
    local dir

    # Loop through each directory provided
    for dir in "$@"; do
        
        # check if the directory already exists
        if [[ -d "${dir}" ]]; then
            logmsg "Directory already exists: ${dir}"

            # create the directory
            elif mkdir -p "${dir}"; then
                logmsg "Directory created: ${dir}"
            
            # Report an error if the directory failed to create
            else
                logmsg "Error: Failed to create directory: ${dir}"
                return 1
        fi
    done
}

fetch_sra(){
    local sra="$1"
    local project="$2"
    local reads="${proj_dir}/${project}/reads/${sra}"
    local tmp="${reads}/.tmp_${sra}"

    local fq1="${reads}/R1_fq.gz"
    local fq2="${reads}/R2_fq.gz"

    # create directories 
    create_dir "${reads}" "${tmp}" || return 1

    # skip if both files already exist
    if [[ -s "${fq1}" && -s "${fq2}" ]]; then
        logmsg "${sra} already exist in $(basename ${reads}); skip download."
        return 0
    else
        # Download and convert SRA to FASTQ
        logmsg "Downloading ${sra} for ${project}."

        if ! fasterq-dump "${sra}" \
            --split-files --threads "${threads}" \
            --outdir "${tmp}" --temp "${tmp}"; then
            logmsg "ERROR : Download failed for ${sra}."
            return 1
        fi
    fi


    # Compress and rename
    local mate fastq gz

    for mate in 1 2; do
        fastq="${tmp}/${sra}_${mate}.fastq"
        gz="${reads}/R${mate}.fq.gz"

        if [[ -s "$gz" ]] && gzip -t "${gz}"; then
            logmsg "$(basename "$gz") already exists; skipping."
            continue
        fi

        if [[ ! -s "$fastq" ]]; then
            logmsg "ERROR: Missing $(basename "$fastq")"
            return 1
        fi

        gzip -c "$fastq" > "${gz}.tmp" || return 1
        gzip -t "${gz}.tmp" || return 1
        mv "${gz}.tmp" "$gz" || return 1

        logmsg "Compressed $(basename "$fastq") to $(basename "$gz")"
    done

    # Verify both files
    if [[ -s "$fq1" && -s "$fq2" ]]; then
        logmsg "${sra} successfully downloaded and compressed."
        rm -rf -- "$tmp"
        return 0
    fi

    return 1

}


# Loop through projects
for project in "${projects[@]}"; do
    for sra in ${sra_groups[$project]}; do
        fetch_sra "$sra" "$project" || exit 1
    done
done