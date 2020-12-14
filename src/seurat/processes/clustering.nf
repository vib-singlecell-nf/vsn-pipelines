#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__CLUSTERING {
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	val assay
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__CLUSTERING_${assay}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.clustering)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/clustering.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__CLUSTERING_${assay}.rds" \
		--assay ${assay} \
		${(processParams.dimsToUse == null) ? '' :'--dimsToUse ' +processParams.dimsToUse } \
		${(processParams.resToUse == null) ? '' :'--resToUse ' + processParams.resToUse} \
		${(processParams.perplexity == null) ? '' : '--perplexity ' + processParams.perplexity} \
		${(processParams.assayForCrossModalityGraphs == null) ? '' : '--assayForCrossModalityGraphs ' + processParams.assayForCrossModalityGraphs}
	"""
}
