#!/bin/bash

#######################################
#                                     #
#      Bowtie2 alignment script       #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #
# Date: 11-11-24                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: Runs end-to-end Bowtie2 alignment for metagenomic             #
# sequence data, including alignment and post-processing with samtools.      #
#                                                                            #
#   - bowtie2 v2.5.1: For sequence alignment                                 #
#     Paper: https://doi.org/10.1038/nmeth.1923                              #
#     Tool: https://github.com/BenLangmead/bowtie2                           #
#   - samtools v1.17: For processing alignment results                       #
#     Tool: https://github.com/samtools/samtools                             #
#   - parallel: For parallel processing of multiple samples                  #
#     Tool: https://www.gnu.org/software/parallel/                           #
#                                                                            #
# Usage: ./run_bowtie2_alignment.sh [options] <r1_file> <r2_file> <db_path>  #
#        [output_dir]                                                        #
#                                                                            #
# Last updated: 20-01-25                                                     #
##############################################################################

# Strict error handling
# -e: Exit immediately if a command exits with non-zero status
# -u: Treat unset variables as an error
# -o pipefail: Return value of a pipeline is the status of the last command to exit with a non-zero status

set -euo pipefail

################################################################################
# Global Variables and Settings                                                #
################################################################################

# Script information
SCRIPT_NAME=$(basename "$0")
VER="1.0.0"

# Default settings
THREADS=16
MODE="PE"
SEED=141205
DRY_RUN=false
DEBUG=false

# Required tools
TOOLS=("bowtie2" "samtools" "parallel")

# source helper functions from bash-utils
# provides logging, debugging, dependency checks etc
if [ ! -f "utils/bash-utils" ]; then
    echo "Bash utility functions not found: utils/bash-utils"
    exit 1
fi

source utils/bash-utils

################################################################################
# Functions                                                                    #
################################################################################

# Function: help
# Description: Display usage information
# Arguments: None
help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [options] <r1_file> <r2_file> <db_path> [output_dir]

Required arguments:
    r1_file     Path to R1 FASTQ file (processed by fastp)
    r2_file     Path to R2 FASTQ file (processed by fastp)
    db_path     Path to Bowtie2 database (without .bt2 extension)
    output_dir  Path to output directory (optional). Default: ./bowtie-output

Options:
Options:
    -t, --threads INT    Number of threads to use (default: 16)
    -m, --mode STR      Sequencing mode: PE or SE (default: PE)
    -s, --seed INT      Random seed (default: 141205)
    -d, --dry-run      Show commands without executing them
    --debug            Enable debug logging
    -h, --help         Show this help message
    -v, --version      Show version information

Example:
    ${SCRIPT_NAME} -t 24 input_R1.fq.gz input_R2.fq.gz /path/to/bt2db/ref
EOF
}

# Function: run_alignment
# Description: Execute Bowtie2 alignment
# Arguments:
#   $1: R1 file
#   $2: R2 file
#   $3: Database path
#   $4: Output prefix
run_alignment() {
    local r1=$1
    local r2=$2
    local db=$3
    local output_prefix=$4
    
    if [[ "${MODE}" == "PE" ]]; then
        if [[ "$r1" == *.gz ]]; then
            execute bowtie2 \
                --seed "${SEED}" \
                --threads "${THREADS}" \
                --end-to-end \
                --very-sensitive \
                --time -q \
                -x "$db" \
                -1 <(gunzip -c "$r1") \
                -2 <(gunzip -c "$r2") \
                -S "${output_prefix}.sam"
        else
            execute bowtie2 \
                --seed "${SEED}" \
                --threads "${THREADS}" \
                --end-to-end \
                --very-sensitive \
                --time -q \
                -x "$db" \
                -1 "$r1" \
                -2 "$r2" \
                -S "${output_prefix}.sam"
        fi
    else
        # Single-end mode
        if [[ "$r1" == *.gz ]]; then
            execute bowtie2 \
                --seed "${SEED}" \
                --threads "${THREADS}" \
                --end-to-end \
                --very-sensitive \
                --time -q \
                -x "$db" \
                -U <(gunzip -c "$r1") \
                -S "${output_prefix}.sam"
        else
            execute bowtie2 \
                --seed "${SEED}" \
                --threads "${THREADS}" \
                --end-to-end \
                --very-sensitive \
                --time -q \
                -x "$db" \
                -U "$r1" \
                -S "${output_prefix}.sam"
        fi
    fi
}

# Function: process_sam
# Description: Process SAM file with samtools
# Arguments:
#   $1: Input SAM file
#   $2: Output prefix
process_sam() {
    local sam_file=$1
    local output_prefix=$2
    
    # Extract unmapped reads
    execute samtools fastq -@ "${THREADS}" \
        -f 12 -F 256 \
        "$sam_file" > "${output_prefix}_unmapped.fq"
        
    # Convert to BAM, sort, and index
    execute samtools view -@ "${THREADS}" --bam "$sam_file" > "${output_prefix}.bam"
    execute samtools sort -@ "${THREADS}" -O bam -o "${output_prefix}.bam" "$sam_file"
    execute samtools index -M -@ "${THREADS}" "${output_prefix}.bam"
    
    # Generate alignment statistics
    execute samtools stats -@ "${THREADS}" "${output_prefix}.bam" > "${output_prefix}.stats"

    execute samtools flagstats -@ "${THREADS}" "${output_prefix}.bam" > "${output_prefix}.flagstats"

}

################################################################################
# Main script body                                                              #
################################################################################

main() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                help
                exit 0
                ;;
            -v|--version)
                echo "${SCRIPT_NAME} version ${VER}"
                exit 0
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -m|--mode)
                MODE="$2"
                shift 2
                ;;
            -s|--seed)
                SEED="$2"
                shift 2
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift
                ;;
            --debug)
                DEBUG=true
                shift
                ;;
            *)
                break
                ;;
        esac
    done

    # Check required arguments
    if [ $# -lt 3 ]; then
        log_error "Error: Missing required arguments"
        help
        exit 1
    fi

    local r1_file=$1
    local r2_file=$2
    local db_path=$3
    local output_dir="${4:-./bowtie-output}"
    
    # Create output directory
    if [[ "$DRY_RUN" != "true" ]]; then
        mkdir -p "$output_dir"
    else
        echo "Would create directory: $output_dir"
    fi
    
    # Validate inputs
    verify_dependencies "${TOOLS[@]}"
    
    if [[ ! -f "$r1_file" ]]; then
        log_error "Error: R1 file not found: $r1_file"
        exit
    fi
    
    if [[ "${MODE}" == "PE" && ! -f "$r2_file" ]]; then
        log_error "Error: R2 file not found: $r2_file"
        exit 1
    fi

    # Generate output prefix
    local basename=$(basename "$r1_file" | sed 's/_R1.*$//')
    local output_prefix="${output_dir}/${basename}"
    
    # Run alignment
    log_info "Starting alignment for ${basename}..."
    run_alignment "$r1_file" "$r2_file" "$db_path" "$output_prefix"
    
    # Process alignment results
    log_info "Processing alignment results..."
    process_sam "${output_prefix}.sam" "$output_prefix"
    
    # Cleanup
    if [[ -f "${output_prefix}.sam" ]]; then
        execute rm "${output_prefix}.sam"
    fi
    
    log_info "Alignment complete for ${basename}"
}

# If being sourced, don't execute main
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
