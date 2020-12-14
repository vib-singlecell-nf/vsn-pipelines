#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__HASHTAGS_SELECTION{
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern : "${sampleId}_logQC.txt"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__HASHTAG_SELECTION.rds")
	file("${sampleId}_logQC.txt")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.hashtagsSelection)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/hashtagsSelection.R --inputSeuratRds "${sobj}" \
		--output "${sampleId}.SEURAT__HASHTAG_SELECTION.rds" \
		--hashtags "${processParams.hashtags}"
	"""
}
