#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__SCTRANSFORM{
	//publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__SCTRANSFORM.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.scTransform)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/scTransform.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__SCTRANSFORM.rds" \
		${(processParams.regressSubsamples == null || processParams.regressSubsamples == "false") ? '': '--regressSubsamples'}
	"""
}
