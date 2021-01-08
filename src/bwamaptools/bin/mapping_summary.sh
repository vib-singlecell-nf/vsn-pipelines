#!/usr/bin/env bash

sampleId="$1"
bam="$2"


# run samtools stat:
samtools stat ${bam} > ${sampleId}.stat

# uniquely mapped reads (BWA):
UMR=$(samtools view -F 0x4 -F 0x100 -F 0x800 ${bam} | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l)

# Fraction of total read pairs mapped confidently to genome (>30 mapq)
CMR=$(samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 ${bam})


# output file:
printf "\t${sampleId}\n" > ${sampleId}.mapping_stats.tsv
grep ^SN ${sampleId}.stat | cut -f 2,3 >> ${sampleId}.mapping_stats.tsv
printf "Uniquely mapped reads:\t${UMR}\n" >> ${sampleId}.mapping_stats.tsv
printf "Reads mapped with MAPQ>30:\t${CMR}\n" >> ${sampleId}.mapping_stats.tsv

