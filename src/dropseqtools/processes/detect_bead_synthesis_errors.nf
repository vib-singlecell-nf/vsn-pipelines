nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS {
    
	container params.sc.dropseqtools.container
	publishDir "${params.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
		tuple val(sample), file("*.final_cleaned.bam"), emit: bam
		tuple file("*.synthesis_stats.txt"), emit: stats
		// tuple file("*.synthesis_stats.summary.txt"), emit: statsSummary
    script:
        """
		DetectBeadSynthesisErrors \
			I=${bam} \
			O=${sample}.final_cleaned.bam \
			OUTPUT_STATS=${sample}.synthesis_stats.txt \
			SUMMARY=${sample}.synthesis_stats.summary.txt \
			NUM_BARCODES=${params.sc.dropseqtools.detect_repair_barcode_synthesis_errors.numBarcodes * 2} \
			PRIMER_SEQUENCE=${params.sc.dropseqtools.detect_repair_barcode_synthesis_errors.primerSequence} \
            TMP_DIR=$DWMAX/tmp
        """    
}