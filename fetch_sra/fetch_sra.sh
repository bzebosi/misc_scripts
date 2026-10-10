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
    local reads="${proj_dir}/${project}/01_raw_reads"
    local tmp="${reads}/.tmp_${sra}"

    local fq1="${reads}/R1.fq.gz"
    local fq2="${reads}/R2.fq.gz"

    local single="${tmp}/${sra}.fastq"
    local pair1="${tmp}/${sra}_1.fastq"
    local pair2="${tmp}/${sra}_2.fastq"

    # create directories 
    create_dir "${reads}" "${tmp}" || return 1

    # Skip if paired-end reads already exist
    if [[ -s "${fq1}" && -s "${fq2}" ]] && gzip -t "${fq1}" && gzip -t "${fq2}"; then
        logmsg "${sra} paired-end reads already exist; skipping."
        return 0
    fi


    # Download and convert SRA to FASTQ
    if [[ ! -s "${single}" && ! -s "${pair1}" && ! -s "${pair2}" ]]; then
        logmsg "Downloading ${sra} for ${project}."

        if ! fasterq-dump "${sra}" \
            --split-files --threads "${threads}" \
            --outdir "${tmp}" --temp "${tmp}"; then
            logmsg "ERROR : Download failed for ${sra}."
            return 1
        fi
    fi

    # Detect single-end or paired-end
    local mates=2

    if [[ -s "${pair1}" && -s "${pair2}" ]]; then
        logmsg "Paired-end reads detected for ${sra}."
    elif [[ -s "${single}" ]]; then
        mates=1
        logmsg "Single-end reads detected for ${sra}."
    else
        logmsg "ERROR: Missing FASTQ files for ${sra}."
        return 1
    fi

    # Compress and rename
    local mate fastq gz

    for ((mate=1; mate<=mates; mate++)); do

        if [[ "$mates" -eq 1 ]]; then
            fastq="${single}"
            gz="${fq1}"
        else
            fastq="${tmp}/${sra}_${mate}.fastq"
            gz="${reads}/${sra}_${mate}.fq.gz"
        fi

        if [[ -s "$gz" ]] && gzip -t "$gz"; then
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

    # Verify output files
    # Verify output files
    if [[ -s "${fq1}" ]] && gzip -t "${fq1}"; then
        if [[ "$mates" -eq 1 ]] || { [[ -s "${fq2}" ]] && gzip -t "${fq2}"; }; then
            logmsg "${sra} successfully downloaded and compressed."
            rm -rf -- "${tmp}"
            return 0
        fi
    fi

    return 1
}


# Loop through projects
for project in "${projects[@]}"; do
    for sra in ${sra_groups[$project]}; do
        fetch_sra "$sra" "$project" || exit 1
    done
done