nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS {

	container params.sc.dropseqtools.container
	publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    	tuple val(sample), path(bam)

	output:
		tuple val(sample), path("*.final_cleaned.bam"), emit: bam
		tuple file("*.synthesis_stats.txt"), emit: stats
		// tuple file("*.synthesis_stats.summary.txt"), emit: statsSummary

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.dropseqtools.detect_repair_barcode_synthesis_errors)
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