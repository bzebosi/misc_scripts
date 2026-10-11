#!/bin/bash

# load your settings:
source "$(dirname "$0")/config.sh"

# Logging functions
#--------------------------------------------------------------------------------------------------------
timestamp() { date '+%Y-%s m-%d %H:%M:%S'; }
logmsg() { echo "$(timestamp): $*"; }

# Validate project and return its directory
#--------------------------------------------------------------------------------------------------------
get_project_dir() {
    local project="${1}"
    local dir

    # check if project was provided
    if [[ -z "${project}" ]]; then
        logmsg "ERROR: No project specified."
        return 1
    fi

    dir="${project_dir[$project]:-}"

    # check if project exits in the array
    if [[ -z "${dir}" ]]; then
        logmsg "ERROR: Unknown project ${project}."
        logmsg "Available projects: ${!project_dir[*]}"
        return 1
    fi

    # Check if the project directory exists
    if [[ ! -d "${dir}" ]]; then
        logmsg "ERROR: Project directory does not exist: ${dir}"
        return 1
    fi

    echo "${dir}"
}

# Create one or more directories
#--------------------------------------------------------------------------------------------------------
create_dir() {
    local dir

    for dir in "$@"; do
        if [[ -d "${dir}" ]]; then
            logmsg "Directory already exists: ${dir}"

        elif mkdir -p "${dir}"; then
            logmsg "Directory created: ${dir}"

        else
            logmsg "ERROR: Failed to create directory: ${dir}"
            return 1
        fi
    done
}



# Directory names
#--------------------------------------------------------------------------------------------------------
reads_dir="01_raw_reads"
fastp_dir="02_fastp"
fastqc_dir="03_fastqc"
trim_reads_dir="trim_reads"
reports_dir="reports"


trim_reads() {
    local project="${1}"
    local project_path
    project_path="$(get_project_dir "${project}")" || return 1
    local raw_dir="${project_path}/${reads_dir}"
    local fastp_path="${project_path}/${fastp_dir}"
    local trim_dir="${fastp_path}/${trim_reads_dir}"
    local report_dir="${fastp_path}/${reports_dir}"
    local summary="${report_dir}/fastp_summary.tsv"

    create_dir "${trim_dir}" "${report_dir}" || return 1

    logmsg "Project: ${project}"
    logmsg "Project directory: ${project_path}"
    logmsg "Raw reads directory: ${raw_dir}"
    logmsg "Trimmed reads directory: ${trim_dir}"
    logmsg "Fastp reports directory: ${report_dir}"

    # Initialize summary file 
    if [[ -s "${summary}" ]]; then
        if [[ "$(head -n 1 "${summary}" | cut -f2)" != "Read_Type" ]]; then
            logmsg "ERROR: Existing summary does not contain Read_Type."
            return 1
        fi
        logmsg "${summary} already exists."
    else
        # Add a header to the summary file
        echo -e \
        "Sample\tRead_Type\tFastp_Version\tTotal_Reads_Before\tTotal_Bases_Before\tQ20_Bases_Before\t\
        Q30_Bases_Before\tGC_Content_Before\tTotal_Reads_After\tTotal_Bases_After\tQ20_Bases_After\t\
        Q30_Bases_After\tGC_Content_After\tPassed_Filter_Reads\tLow_Quality_Reads\tToo_Many_N_Reads\t\
        Too_Short_Reads\tDuplication_Rate\tInsert_Size_Peak" > "${summary}"
        
        logmsg "summary file created at: $summary"
    fi

    shopt -s nullglob

    local R1

    for R1 in "${raw_dir}"/*.fq.gz; do
        local fname="$(basename "${R1}")"
        local sbase R2 O1 O2 readtype

        # Skip paired-end R2
        [[ "${fname}" == *_R2.fq.gz ]] && continue

        # Detect single-end or paired-end reads
        if [[ "${fname}" == *_R1.fq.gz ]]; then
            readtype="PE"
            sbase="${fname%_R1.fq.gz}"
            R2="${raw_dir}/${sbase}_R2.fq.gz"
            O1="${trim_dir}/${sbase}_trim_R1.fq.gz"
            O2="${trim_dir}/${sbase}_trim_R2.fq.gz"

            # Check matching R2 for paired-end reads
            if [[ ! -s "${R2}" ]]; then
                logmsg "WARN: Missing R2 file for ${sbase}: ${R2}"
                continue
            fi
        else
            readtype="SE"
            sbase="${fname%.fq.gz}"
            R2=""
            O1="${trim_dir}/${sbase}_trim.fq.gz"
            O2=""
        fi

        local ht="${report_dir}/${sbase}_fastp_report.html"
        local jt="${report_dir}/${sbase}_fastp_report.json"
        local metrics=""
        local row=""

        logmsg "Processing ${sbase} (${readtype})"

        # Skip if sample is already summarized and trimmed outputs exist
        if [[ -s "${O1}" && ( "${readtype}" == "SE" || -s "${O2}" ) && -s "${jt}" ]]; then
            # Check if the sample is already in the summary
            if awk -F '\t' -v sample="${sbase}" -v type="${readtype}" '
                # Ignore the header and find the row for the current sample
                NR > 1 && $1 == sample && $2 == type {
                    # Verify the row has the expected number of columns
                    if (NF == 19 ) {

                        # Assume the row is complete
                        found = 1
                    
                        # Check every column for missing values
                        for (i = 1; i <= NF; i++ ){

                            # If any column is empty, mark the row as incomplete
                            if ($i == "" ) {
                                found = 0
                                break
                            }
                        }
                    }
                }

                # Return success (0) if the row is complete,
                # otherwise return failure (1)
                END {
                    exit !found
                }

                ' "${summary}"; then
                logmsg "${sbase} is already trimmed and summarized. Skipping."
                continue
            fi
        fi

        # Run fastp if trimmed files do not exist
        if [[ ! -s "${O1}" || ! -s "${jt}" ]] || [[ "${readtype}" == "PE" && ! -s "${O2}" ]]; then
            if [[ "${readtype}" == "PE" ]]; then
                if fastp -i "${R1}" -I "${R2}" -o "${O1}" -O "${O2}" --detect_adapter_for_pe --thread "${threads:-8}" -h "${ht}" -j "${jt}"; then
                    logmsg "fastp complete for ${sbase}"
                else
                    logmsg "ERROR: fastp failed for ${sbase}"
                    continue
                fi
            else
                if fastp -i "${R1}" -o "${O1}" --thread "${threads:-8}" -h "${ht}" -j "${jt}"; then
                    logmsg "fastp complete for ${sbase}"
                else
                    logmsg "ERROR: fastp failed for ${sbase}"
                    continue
                fi
            fi
        else
            logmsg "${sbase} already trimmed. Skipping fastp..."
        fi

        # Extract fastp metrics
        if [[ -s "${jt}" ]]; then
            metrics="$(
                jq -r '
                    [.summary.fastp_version, .summary.before_filtering.total_reads, .summary.before_filtering.total_bases,
                     .summary.before_filtering.q20_bases, .summary.before_filtering.q30_bases, .summary.before_filtering.gc_content,
                     .summary.after_filtering.total_reads, .summary.after_filtering.total_bases, .summary.after_filtering.q20_bases,
                     .summary.after_filtering.q30_bases, .summary.after_filtering.gc_content, .filtering_result.passed_filter_reads,
                     .filtering_result.low_quality_reads, .filtering_result.too_many_N_reads, .filtering_result.too_short_reads,
                     .duplication.rate, (.insert_size.peak // "NA") ] | map(if . == null then "NA" else . end) | @tsv' "${jt}"
                )" || metrics=""
        fi

        # Build the exact row and append only if not already present
        if [[ -z "${metrics}" ]]; then
            logmsg "ERROR: Could not extract fastp metrics for ${sbase}"
            continue
        fi

        row="${sbase}\t${readtype}\t${metrics}"

        # Add statistics if the sample is not already present
        if awk -F '\t' -v sample="${sbase}" -v type="${readtype}" '
            NR > 1 && $1 == sample && $2 == type && NF == 19 {
                found = 1

                for (i = 1; i <= NF; i++ ) {
                    if ($i == ""){
                    found = 0 
                    break
                    }
                }
            }

            END {
                exit !found
            }
        ' "${summary}"; then

            logmsg "Summary already contains ${sbase}. Skipping summary update."
        else
            # remove any old or incomplete row for this sample
            local summary_tmp="${summary}.tmp.$$"
            awk -F '\t' -v sample="${sbase}" 'NR == 1 || $1 != sample' "${summary}" > "${summary_tmp}" || return 1

            mv "${summary_tmp}" "${summary}" || return 1
            
            # Add the new complete row
            echo -e "${row}" >> "${summary}"
            logmsg "Statistics added for ${sbase}"
        fi
    done

    shopt -u nullglob

    logmsg "Read trimming completed for project: ${project}"
}

# Fastqc
run_fastqc() {
    local project="${1}"
    local stage="${2}"
    local project_path
    project_path="$(get_project_dir "${project}")" || return 1
    local input_dir
    local output_dir
    local read_suffix

    # Set input and output directories
    if [[ "${stage}" == "raw" ]]; then
        input_dir="${project_path}/${reads_dir}"
        output_dir="${project_path}/${fastqc_dir}/raw"
        read_suffix=""
    elif [[ "${stage}" == "trimmed" ]]; then
        input_dir="${project_path}/${fastp_dir}/${trim_reads_dir}"
        output_dir="${project_path}/${fastqc_dir}/trimmed"
        read_suffix="_trim"
    else
        logmsg "ERROR: Unknown FastQC stage: ${stage}"
        return 1
    fi

    create_dir "${output_dir}" || return 1
    logmsg "Running FastQC on ${stage} reads for project: ${project}"

    shopt -s nullglob

    # Run FastQC on raw or trimmed reads
    local R1
    for R1 in "${input_dir}"/*.fq.gz; do
        local fname="$(basename "${R1}")"
        local sample R2 readtype
        local files=()

        # Skip paired-end R2
        [[ "${fname}" == *"${read_suffix}_R2.fq.gz" ]] && continue

        # Detect single-end or paired-end reads
        if [[ "${fname}" == *"${read_suffix}_R1.fq.gz" ]]; then
            readtype="PE"
            sample="${fname%${read_suffix}_R1.fq.gz}"
            R2="${input_dir}/${sample}${read_suffix}_R2.fq.gz"

            # Check matching R2 for paired-end reads
            if [[ ! -s "${R2}" ]]; then
                logmsg "WARN: Missing R2 for ${sample}. Skipping."
                continue
            fi

            files=("${R1}" "${R2}")

        elif [[ "${fname}" == *"${read_suffix}.fq.gz" ]]; then
            readtype="SE"
            sample="${fname%${read_suffix}.fq.gz}"
            files=("${R1}")
        else
            logmsg "WARN: Unexpected FASTQ filename: ${fname}. Skipping."
            continue
        fi

        logmsg "Processing ${sample} (${readtype})"
        
        # Skip if FastQC reports already exist
        local report_prefix="${output_dir}/${sample}${read_suffix}"
        if [[ "${readtype}" == "PE" ]]; then
            if [[ -s "${report_prefix}_R1_fastqc.html" && -s "${report_prefix}_R1_fastqc.zip" && \
                  -s "${report_prefix}_R2_fastqc.html" && -s "${report_prefix}_R2_fastqc.zip" ]]; then
                logmsg "FastQC already completed for ${sample}. Skipping."
                continue
            fi

        else
            if [[ -s "${report_prefix}_fastqc.html" &&  -s "${report_prefix}_fastqc.zip" ]]; then
                logmsg "FastQC already completed for ${sample}. Skipping."
                continue
            fi
        fi

        # Run FastQC
        if fastqc --threads "${threads:-8}" --outdir "${output_dir}" "${files[@]}"; then
            logmsg "FastQC completed for ${sample}."
        else
            logmsg "ERROR: FastQC failed for ${sample}."
            shopt -u nullglob
            return 1
        fi

    done

    shopt -u nullglob

    logmsg "FastQC completed for ${stage} reads: ${project}"
}

###############################################################################
# downlad_genome function
#   • Looks up URL in genome_urls[] using the tag (e.g. b73, mo17…)
#   • Downloads the compressed FASTA to $ref_dir if missing
#   • Logs every step
###############################################################################

# Download genome if it does not already exist
get_genome() {
    local gbase="${1}"
    local genome_dir="${ref_dir}/${gbase}"

    create_dir "${genome_dir}" || exit 1

    local file url

    for url in "${genome_urls[$gbase]}" "${annotation_urls[$gbase]}"; do

        if [[ -z "${url}" ]]; then
            logmsg "ERROR: Missing download URL for ${gbase}"
            exit 1
        fi

        file="${genome_dir}/$(basename "${url}")"

        if [[ -s "${file}" ]]; then
            logmsg "$(basename "${file}") already exists. Skipping."
        else
            logmsg "Downloading $(basename "${file}")..."
            
            if wget -O "${file}" "${url}"; then
                logmsg "Successfully downloaded $(basename "${file}")"
            else
                logmsg "ERROR: Failed to download $(basename "${file}")"
                rm -f "${file}"
                exit 1
            fi
        fi
    done
}


# Build STAR genome index
index_genome() {
    local gbase="${1}"
    local genome_dir="${ref_dir}/${gbase}"
    local genome_gz="${genome_dir}/$(basename "${genome_urls[$gbase]}")"
    local genome_fasta="${genome_gz%.gz}"
    local annotation_gff="${genome_dir}/$(basename "${annotation_urls[$gbase]}")"
    local star_index="${genome_dir}/STAR_index"

    # Check required files
    if [[ ! -s "${genome_gz}" ]]; then
        logmsg "ERROR: Genome file not found: ${genome_gz}"
        exit 1
    fi

    if [[ ! -s "${annotation_gff}" ]]; then
        logmsg "ERROR: Annotation file not found: ${annotation_gff}"
        exit 1
    fi

    create_dir "${star_index}" || exit 1

    # Decompress genome if needed
    if [[ ! -s "${genome_fasta}" ]]; then
        logmsg "Decompressing genome: ${gbase}"

        if gzip -dc "${genome_gz}" > "${genome_fasta}"; then
            logmsg "Genome decompressed successfully."
        else
            logmsg "ERROR: Failed to decompress genome."
            rm -f "${genome_fasta}"
            exit 1
        fi
    fi


    # Skip if STAR index already exists
    if [[ -s "${star_index}/SA" ]]; then
        logmsg "STAR index already exists for ${gbase}. Skipping."
    else
        logmsg "Building STAR index for ${gbase}..."

        if STAR \
            --runMode genomeGenerate \
            --runThreadN "${threads:-8}" \
            --genomeDir "${star_index}" \
            --genomeFastaFiles "${genome_fasta}" \
            --sjdbGTFfile "${annotation_gff}" \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang 149; then

            logmsg "STAR index successfully built for ${gbase}."
        else
            logmsg "ERROR: Failed to build STAR index for ${gbase}."
            exit 1
        fi
    fi
}


map_reads(){
    local project="${1}"
    local gbase="${2}"
    local project_path
    project_path="$(get_project_dir "${project}")" || return 1
    local genome_dir="${ref_dir}/${gbase}"
    local star_index="${genome_dir}/STAR_index"
    local trim_dir="${project_path}/${fastp_dir}/${trim_reads_dir}"
    local align_dir="${project_path}/04_alignment"
    local bam_dir="${align_dir}/bam"
    local log_dir="${align_dir}/logs"
    local sj_dir="${align_dir}/sj"

    create_dir "${bam_dir}" "${log_dir}" "${sj_dir}" || return 1

    if [[ ! -s "${star_index}/SA" ]]; then
        logmsg "ERROR: STAR index not found or incomplete: ${star_index}"
        return 1
    fi

    shopt -s nullglob

    local reads=( "${trim_dir}"/*.fq.gz )

    if [[ ${#reads[@]} -eq 0 ]]; then
        logmsg "ERROR: No trimmed reads found in ${trim_dir}"
        shopt -u nullglob
        return 1
    fi
    
    local R1
    for R1 in "${reads[@]}"; do
        local fname="$(basename "${R1}")"
        local sample R2 readtype
        local files=()

        # Skip paired-end R2
        [[ "${fname}" == *_trim_R2.fq.gz ]] && continue

        # Detect single-end or paired-end reads
        if [[ "${fname}" == *_trim_R1.fq.gz ]]; then
            readtype="PE"
            sample="${fname%_trim_R1.fq.gz}"
            R2="${trim_dir}/${sample}_trim_R2.fq.gz"

            # Check matching R2 for paired-end reads
            if [[ ! -s "${R2}" ]]; then
                logmsg "WARN: Missing R2 for ${sample}. Skipping."
                continue
            fi

            files=("${R1}" "${R2}")

        elif [[ "${fname}" == *_trim.fq.gz ]]; then
            readtype="SE"
            sample="${fname%_trim.fq.gz}"
            files=("${R1}")

        else
            logmsg "WARN: Unexpected FASTQ filename: ${fname}. Skipping."
            continue
        fi

        local prefix="${align_dir}/${sample}_"
        local bam="${bam_dir}/${sample}_Aligned.sortedByCoord.out.bam"

        # Skip if BAM already exists
        if [[ -s "${bam}" ]]; then
            logmsg "${sample} already mapped. Skipping."
        else
            logmsg "Mapping ${sample} (${readtype})..."
            if STAR \
                --runThreadN "${threads:-8}" --genomeDir "${star_index}" --readFilesIn "${files[@]}" \
                --readFilesCommand zcat --outFileNamePrefix "${prefix}" --outSAMtype BAM SortedByCoordinate; then
                
                # Organize STAR outputs
                if mv "${prefix}Aligned.sortedByCoord.out.bam" "${bam_dir}/" && \
                   mv "${prefix}"{Log.final.out,Log.out,Log.progress.out} "${log_dir}/" && \
                   mv "${prefix}SJ.out.tab" "${sj_dir}/"; then

                    rm -rf "${prefix}_STARtmp"
                    logmsg "Successfully mapped ${sample}."
                else
                    logmsg "ERROR: Failed to organize STAR outputs for ${sample}."
                    shopt -u nullglob
                    return 1
                fi
            else
                logmsg "ERROR: STAR mapping failed for ${sample}."
                shopt -u nullglob
                return 1
            fi
        fi

    done

    shopt -u nullglob
    logmsg "STAR mapping completed for project: ${project}"

}


count_reads() {
    local project="${1}"
    local gbase="${2}"
    local project_path
    project_path="$(get_project_dir "${project}")" || return 1
    local genome_dir="${ref_dir}/${gbase}"
    local annotation_gtf="${genome_dir}/$(basename "${annotation_urls[$gbase]}" .gff3).gtf"
    local bam_dir="${project_path}/04_alignment/bam"
    local count_dir="${project_path}/05_featurecounts"
    local counts="${count_dir}/counts.tsv"

    create_dir "${count_dir}" || return 1

    if [[ ! -s "${annotation_gtf}" ]]; then
        logmsg "ERROR: GTF annotation not found: ${annotation_gtf}"
        return 1
    fi

    shopt -s nullglob

    local bams=( "${bam_dir}"/*.bam )

    if [[ ${#bams[@]} -eq 0 ]]; then
        logmsg "ERROR: No BAM files found in ${bam_dir}"
        shopt -u nullglob
        return 1
    fi

    if [[ -s "${counts}" ]]; then
        logmsg "featureCounts already completed. Skipping."
    else
        logmsg "Counting reads..."

        local bam paired
        local readtype=""
        local options=()

        for bam in "${bams[@]}"; do
            paired="$(samtools view -c -f 1 "${bam}")" || return 1

            if [[ "${paired}" -gt 0 ]]; then
                if [[ "${readtype}" == "SE" ]]; then
                    logmsg "ERROR: Mixed SE and PE BAM files."
                    shopt -u nullglob
                    return 1
                fi
                readtype="PE"
            else
                if [[ "${readtype}" == "PE" ]]; then
                    logmsg "ERROR: Mixed SE and PE BAM files."
                    shopt -u nullglob
                    return 1
                fi
                readtype="SE"
            fi
        done

        if [[ "${readtype}" == "PE" ]]; then
            options=(-p --countReadPairs)
        fi


        if featureCounts \
            -T "${threads:-8}" \
            "${options[@]}" \
            -t exon \
            -g gene_id \
            -a "${annotation_gtf}" \
            -o "${counts}" \
            "${bams[@]}"; then

            logmsg "featureCounts completed successfully."
        else
            logmsg "ERROR: featureCounts failed."
            rm -f "${counts}" "${counts}.summary"
            shopt -u nullglob
            return 1
        fi
    fi

    shopt -u nullglob
}



# Run trimming once per project
for project in "${projects[@]}"; do
    trim_reads "${project}"
    run_fastqc "${project}" "raw"
    run_fastqc "${project}" "trimmed"
done


# Download reference genomes and build STAR indices
for gbase in "${genomes[@]}"; do
    get_genome "${gbase}"
    index_genome "${gbase}"
done

# Map reads
for project in "${projects[@]}"; do
    for gbase in "${genomes[@]}"; do
        map_reads "${project}" "${gbase}"
        count_reads "${project}" "${gbase}"
    done
done
