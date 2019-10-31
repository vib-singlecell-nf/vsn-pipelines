
process SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE {

	container params.sc.dropseqtools.container
    publishDir "${params.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_Cell.bam'), emit: bam
        tuple file('*.unaligned_tagged_Cellular.bam_summary.txt'), emit: report
    script:
        """
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_Cell.bam \
			SUMMARY=${sample}.unaligned_tagged_Cellular.bam_summary.txt \
			BASE_RANGE=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.baseRange} \
			BASE_QUALITY=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.baseQuality} \
			BARCODED_READ=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.barcodedRead} \
			DISCARD_READ=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.discardRead} \
			TAG_NAME=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode.numBasesBelowQuality}        
        """
}

process SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR {
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
			BASE_RANGE=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.baseRange} \
			BASE_QUALITY=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.baseQuality} \
			BARCODED_READ=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.barcodedRead} \
			DISCARD_READ=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.discardRead} \
			TAG_NAME=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular.numBasesBelowQuality}        
        """
}
