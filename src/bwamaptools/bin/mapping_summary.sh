#!/usr/bin/env bash

nbr_threads=2;

sampleId="${1}";
bam="${2}";

if [ ${#@} -ne 2 ] ; then
    printf 'Usage: mapping_summary.sh sampleId bam_file\n' >&2;
    exit 1;
fi


# Run samtools stat:
samtools stat -@ "${nbr_threads}" "${bam}" > "${sampleId}.stat"

# Uniquely mapped reads (BWA):
UMR=$(samtools view -@ "${nbr_threads}" -F 0x4 -F 0x100 -F 0x800 "${bam}" | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l);

# Fraction of total read pairs mapped confidently to genome (>30 mapq):
CMR=$(samtools view -@ "${nbr_threads}" -c -F 0x4 -F 0x100 -F 0x800 -q 30 "${bam}");


# Output file:
printf "\t${sampleId}\n" > "${sampleId}.mapping_stats.tsv";
grep '^SN' "${sampleId}.stat" | cut -f 2,3 >> "${sampleId}.mapping_stats.tsv";
printf "Uniquely mapped reads:\t${UMR}\n" >> "${sampleId}.mapping_stats.tsv";
printf "Reads mapped with MAPQ>30:\t${CMR}\n" >> "${sampleId}.mapping_stats.tsv";

