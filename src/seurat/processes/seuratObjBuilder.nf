#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__SEURAT_OBJECT_BUILDER {
	//publishDir "${params.global.outDir}/${sampleId}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(sampleId), val(cmat)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__SEURAT_OBJECT_BUILDER.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.seuratObjBuilder)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/seuratObjBuilder.R --inputMatrix "${cmat}" \
		--output "${sampleId}.SEURAT__SEURAT_OBJECT_BUILDER.rds" \
		--sample $sampleId \
		${(!(processParams.minFeaturesGEX == null)) ? '--minFeaturesGEX ' + processParams.minFeaturesGEX: ''} \
		${(!(processParams.minCellsGEX == null)) ? '--minCellsGEX ' + processParams.minCellsGEX: ''}
	"""
}
