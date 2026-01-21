#!/bin/bash

#######################################
#                                     #
#      BAM Coverage Calculator        #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #
# Date: 2024-02-13                                                           #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script calculates genomic coverage from BAM files using  #
#              bedtools. It processes multiple BAM files in parallel and     #
#              generates coverage statistics for each file.                  #
#                                                                            #
# Dependencies:                                                              #
#   - bedtools: For calculating genomic coverage                             #
#   - samtools: For checking if BAM files are sorted                         #
#   - parallel: For parallel processing of multiple files                    #
#                                                                            #
# Citations:                                                                 #
#   - BEDTools                                                               #
#      - Paper: https://doi.org/10.1093/bioinformatics/btq033                # 
#      - Available at: https://github.com/arq5x/bedtools2                    #
#   - Gnu Parallel                                                           #
#      - Paper: https://doi.org/10.5281/ZENODO.1146014                       #
#      - Available at: https://www.gnu.org/software/parallel/                #
#   - SAMtools                                                               #
#      - Paper: https://doi.org/10.1093/gigascience/giab008                  #
#      - Available at: https://github.com/samtools/samtools                  #
#                                                                            #
# Usage: ./calculate_coverage.sh -i /path/to/bam/dir -o /path/to/output      #
#                                                                            #
# Last updated: 2024-02-13                                                   #
##############################################################################

# Strict error handling
set -euo pipefail

# Source utility functions
# shellcheck disable=SC1091
source "utils/bash-utils"

# Required tools
TOOLS=("bedtools" "samtools" "parallel")

################################################################################
# Functions                                                                    #
################################################################################

# Function: check_bam_sorted
# Description: Verifies if a BAM file is sorted
# Arguments:
#   $1: Path to BAM file
# Usage: check_bam_sorted "sample.bam"
check_bam_sorted() {
    local bam_file=$1
    if ! samtools view -H "$bam_file" | grep -q "SO:coordinate"; then
        log_error "BAM file not sorted: $bam_file"
        return 1
    fi
    return 0
}

# Function: calculate_coverage
# Description: Calculates coverage for a single BAM file
# Arguments:
#   $1: Path to BAM file
#   $2: Output directory
# Usage: calculate_coverage "sample.bam" "/path/to/output"
calculate_coverage() {
    local bam_file=$1
    local output_dir=$2
    local basename=$(basename "${bam_file%.*}")
    
    log_info "Processing $basename"
    
    if ! check_bam_sorted "$bam_file"; then
        log_error "Skipping unsorted BAM file: $bam_file"
        return 1
    fi
    
    bedtools genomecov -ibam "$bam_file" > "${output_dir}/${basename}_coverage.txt"
}

################################################################################
# Main Script Logic                                                            #
################################################################################

main() {
    local input_dir=""
    local output_dir=""
    local threads=16
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                input_dir="$2"
                shift 2
                ;;
            -o|--output)
                output_dir="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -h|--help)
                help
                exit 0
                ;;
            *)
                log_error "Unknown parameter: $1"
                help
                exit 1
                ;;
        esac
    done
    
    # Validate input parameters
    if [[ -z "$input_dir" || -z "$output_dir" ]]; then
        log_error "Input and output directories are required"
        help
        exit 1
    fi
    
    # Create output directory if it doesn't exist
    create_directory "$output_dir"
    
    # Verify dependencies
    verify_dependencies "${TOOLS[@]}"
    
    # Find all BAM files
    mapfile -t bam_files < <(find "$input_dir" -maxdepth 1 -type f -name "*.bam")
    
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_error "No BAM files found in $input_dir"
        exit 1
    fi
    
    log_info "Found ${#bam_files[@]} BAM files to process"
    
    # Process BAM files in parallel
    export -f calculate_coverage
    export -f log_info
    export -f log_error
    export -f check_bam_sorted
    
    parallel --jobs "$threads" --eta \
        calculate_coverage {} "$output_dir" ::: "${bam_files[@]}"
    
    log_info "Coverage calculation complete"
}

# Execute main function
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi