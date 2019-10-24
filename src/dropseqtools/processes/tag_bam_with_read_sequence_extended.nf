
process DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_Cell.bam'), emit: bam
        tuple file('*.unaligned_tagged_Cellular.bam_summary.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_Cell.bam \
			SUMMARY=${sample}.unaligned_tagged_Cellular.bam_summary.txt \
			BASE_RANGE=${params.baseRange} \
			BASE_QUALITY=${params.baseQuality} \
			BARCODED_READ=${params.barcodedRead} \
			DISCARD_READ=${params.discardRead} \
			TAG_NAME=${params.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params.numBasesBelowQuality}        
        """
}

process DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_CellMolecular.bam'), emit: bam
        tuple file('*.unaligned_tagged_Molecular.bam_summary.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_CellMolecular.bam \
			SUMMARY=${sample}.unaligned_tagged_Molecular.bam_summary.txt \
			BASE_RANGE=${params.baseRange} \
			BASE_QUALITY=${params.baseQuality} \
			BARCODED_READ=${params.barcodedRead} \
			DISCARD_READ=${params.discardRead} \
			TAG_NAME=${params.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params.numBasesBelowQuality}        
        """
}
