#!/bin/bash
#
# Copyright (C): 2020 - Gert Hulselmans
#
# Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping only reads:
#            - which overlap with SNPs in the VCF file
#            - and which have a cell barcode contained in the cell barcode list
#          Keeping only relevant reads for dsc-pileup can speedup it up several hunderd times.



# Function to check if any of the programs in a pipe failed.
check_exit_codes () {
    local GET_PIPESTATUS="${PIPESTATUS[@]}";
    local exit_code;

    for exit_code in ${GET_PIPESTATUS} ; do
        if [ ${exit_code} -ne 0 ] ; then
             return ${exit_code};
        fi
    done

    return 0;
}



# Check if necessary programs are installed.
check_if_programs_exists () {
    local exit_code=0;

    # Check if bedtools is installed.
    if ! type bedtools > /dev/null 2>&1 ; then
        printf 'Error: "bedtools" could not be found in PATH.\n';
        exit_code=2;
    fi

    # Check if samtools is installed.
    if ! type samtools > /dev/null 2>&1 ; then
        printf 'Error: "samtools" could not be found in PATH.\n';
        exit_code=2;
    fi

    # Check if samtools 1.10 or higher is installed (needs to have "-D STR:FILE" option).
    if ! samtools view --help 2>&1 | grep -q -- '-D STR:FILE' ; then
        printf 'Error: The version of "samtools" should be 1.10 or higher.\n';
        exit_code=2;
    fi

    return ${exit_code};
}



filter_bam_file_for_popscle_dsc_pileup () {
    local input_bam_filename="${1}";
    local barcodes_tsv_filename="${2}";
    local vcf_filename="${3}";
    local output_bam_filename="${4}";
    local threads="${5}";

    local threads_regex="^[1-9][0-9]*"
    if [[ ! ${threads} =~ ${threads_regex} ]]; then
            threads=8 
    fi

    local exit_code=0;

    if [ ${#@} -lt 4 ] ; then
        printf 'Usage:   filter_bam_file_for_popscle_dsc_pileup input_bam_filename barcodes_tsv_filename vcf_filename output_bam_filename [threads]\n\n';
        printf 'Purpose: Filter BAM file for usage with popscle dsc-pileup by keeping reads:\n';
        printf '           - which overlap with SNPs in the VCF file\n';
        printf '           - and which have a cell barcode contained in the cell barcode list\n';
        printf '         Keeping only relevant reads for popscle dsc-pileup can speedup it up several hunderd times.\n\n';

        return 1;
    fi

    if [ ! -f  "${input_bam_filename}" ] ; then
        printf 'Error: Input (CellRanger) BAM file "%s" could not be found.\n' "${input_bam_filename}";
        return 2;
    fi

    if [ ! -f  "${barcodes_tsv_filename}" ] ; then
        printf 'Error: File with barcodes "%s" could not be found.\n' "${barcodes_tsv_filename}";
        return 2;
    fi

    if [ ! -f  "${vcf_filename}" ] ; then
        printf 'Error: File with unique SNPs per sample "%s" could not be found.\n' "${vcf_filename}";
        return 2;
    fi

    # Check if bedtools and samtools are in PATH.
    if ! check_if_programs_exists ; then
        return 2;
    fi

    # Create much smaller BAM file for dsc-pileup of popscle:
    #   - Convert VCF file with unique SNPs for each sample
    #     to a BED file and merge adjacent SNP regions to one.
    #   - Only include reads that contain a SNP position
    #     and which contain a cell barcode of interest.
    if [ "${barcodes_tsv_filename%.gz}".gz = "${barcodes_tsv_filename}" ] ; then
        # Barcodes file is compressed with gzip.
        bedtools merge -i "${vcf_filename}" \
          | samtools view\
                -@ "${threads}" \
                -L - \
                -D CB:<(zcat "${barcodes_tsv_filename}") \
                -o "${output_bam_filename}" \
                "${input_bam_filename}";

        # Check if any of the previous commands failed.
        check_exit_codes;

        exit_code=$?;
    else
        # Barcodes file is uncompressed.
        bedtools merge -i "${vcf_filename}" \
          | samtools view\
                -@ "${threads}" \
                -L - \
                -D CB:"${barcodes_tsv_filename}" \
                -o "${output_bam_filename}" \
                "${input_bam_filename}";

        # Check if any of the previous commands failed.
        check_exit_codes;

        exit_code=$?;
    fi


    return ${exit_code};
}



filter_bam_file_for_popscle_dsc_pileup "${@}";