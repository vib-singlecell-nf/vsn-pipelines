
process SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE {

	container params.tools.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path('*.unaligned_tagged_Cell.bam'), emit: bam
		tuple file('*.unaligned_tagged_Cellular.bam_summary.txt'), emit: report

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.tools.dropseqtools.tag_unaligned_bam_with_cellbarcode)
		processParams = sampleParams.local
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
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path('*.unaligned_tagged_CellMolecular.bam'), emit: bam
		tuple file('*.unaligned_tagged_Molecular.bam_summary.txt'), emit: report

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.tools.dropseqtools.tag_unaligned_bam_with_cellmolecular)
		processParams = sampleParams.local
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
