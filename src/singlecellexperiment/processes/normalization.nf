#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__NORMALIZATION {
	//publishDir "${params.out_dir}/${samplename}", mode: 'copy'
	container params.sce.container
	input:
	tuple val(sampleId), file(scerds)
	output:
	tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__NORMALIZATION.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.sce.normalization)
  		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/normalization.R --sceObj ${scerds} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__NORMALIZATION.rds" \
		${(processParams.quickCluster == null) ? '':"--quickCluster"} \
		${(processParams.useSumFactors == null) ? '':"--useSumFactors"}
	"""
}
