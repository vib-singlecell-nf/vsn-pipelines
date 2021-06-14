#!/usr/bin/env bash

sampleId="${1}";
bam="${2}";

if [ ${#@} -ne 2 ] ; then
    printf 'Usage: mapping_summary.sh sampleId bam_file\n' >&2;
    exit 1;
fi


# Get mapping statistics from BAM file:
#   - Read BAM file and write uncompressed BAM.
#   - Uncompressed BAM file is written to each samtools command with tee (writes to each specified file and stdout).
#   - samtools commands:
#       - Get samtools statistics with:
#           samtools stat "${bam}" > "${sampleId}.stat"
#       - Uniquely mapped reads (BWA):
#           samtools view -c -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' "${bam}"
#       - Fraction of total read pairs mapped confidently to genome (>30 mapq):
#           samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 "${bam}"
#   - Only use threads for "samtools stat". Using it with any of the other samtools commands
#     makes everything slower than not using any threads at all.
samtools view -u "${bam}" \
  | tee \
        >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' - > "${sampleId}.uniquely_mapped_reads.txt") \
        >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 - > "${sampleId}.fraction_total_read_pairs.txt") \
  | samtools stat -@ 2 - > "${sampleId}.stat"


# Output file:
printf "\t${sampleId}\n" > "${sampleId}.mapping_stats.tsv";

grep '^SN' "${sampleId}.stat" | cut -f 2,3 >> "${sampleId}.mapping_stats.tsv";

printf "Uniquely mapped reads:\t" >> "${sampleId}.mapping_stats.tsv";
cat "${sampleId}.uniquely_mapped_reads.txt" >> "${sampleId}.mapping_stats.tsv";

printf "Reads mapped with MAPQ>30:\t" >> "${sampleId}.mapping_stats.tsv";
cat "${sampleId}.fraction_total_read_pairs.txt" >> "${sampleId}.mapping_stats.tsv";

rm "${sampleId}.uniquely_mapped_reads.txt" "${sampleId}.fraction_total_read_pairs.txt";
