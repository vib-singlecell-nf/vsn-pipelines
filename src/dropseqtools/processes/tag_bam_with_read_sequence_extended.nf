
process SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE {

	container params.sc.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path('*.unaligned_tagged_Cell.bam'), emit: bam
		tuple file('*.unaligned_tagged_Cellular.bam_summary.txt'), emit: report

	script:
		processParams = params.sc.dropseqtools.tag_unaligned_bam_with_cellbarcode
		"""
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_Cell.bam \
			SUMMARY=${sample}.unaligned_tagged_Cellular.bam_summary.txt \
			BASE_RANGE=${processParams.baseRange} \
			BASE_QUALITY=${processParams.baseQuality} \
			BARCODED_READ=${processParams.barcodedRead} \
			DISCARD_READ=${processParams.discardRead} \
			TAG_NAME=${processParams.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${processParams.numBasesBelowQuality}        
		"""

}

process SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR {

    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path('*.unaligned_tagged_CellMolecular.bam'), emit: bam
		tuple file('*.unaligned_tagged_Molecular.bam_summary.txt'), emit: report

	script:
		processParams = params.sc.dropseqtools.tag_unaligned_bam_with_cellmolecular
		"""
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_CellMolecular.bam \
			SUMMARY=${sample}.unaligned_tagged_Molecular.bam_summary.txt \
			BASE_RANGE=${processParams.baseRange} \
			BASE_QUALITY=${processParams.baseQuality} \
			BARCODED_READ=${processParams.barcodedRead} \
			DISCARD_READ=${processParams.discardRead} \
			TAG_NAME=${processParams.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${processParams.numBasesBelowQuality}        
		"""

}
