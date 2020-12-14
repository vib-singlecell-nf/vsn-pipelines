#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__MARKER_GENES{
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__MARKER_GENES.rds")
	file("Plots/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.markerGenes)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/markerGenes.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__MARKER_GENES.rds" \
		--markerGenes "${processParams.markerGenes}"
	"""
}
