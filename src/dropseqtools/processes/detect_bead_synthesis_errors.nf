nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS {

	container params.getToolParams("dropseqtools").container
	publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path("*.final_cleaned.bam"), emit: bam
		tuple file("*.synthesis_stats.txt"), emit: stats
		// tuple file("*.synthesis_stats.summary.txt"), emit: statsSummary

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("dropseqtools").detect_repair_barcode_synthesis_errors)
		processParams = sampleParams.local
		"""
		DetectBeadSynthesisErrors \
			I=${bam} \
			O=${sample}.final_cleaned.bam \
			OUTPUT_STATS=${sample}.synthesis_stats.txt \
			SUMMARY=${sample}.synthesis_stats.summary.txt \
			NUM_BARCODES=${processParams.numBarcodes * 2} \
			PRIMER_SEQUENCE=${processParams.primerSequence} \
			TMP_DIR=$DWMAX/tmp
		"""

}
