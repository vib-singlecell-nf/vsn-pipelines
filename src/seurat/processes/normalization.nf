#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__NORMALIZATION{
	//publishDir "${params.global.outDir}/${sampleId}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	val assayType
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__NORMALIZATION_${assayType}.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.normalization)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/normalization.R --inputSeuratRds "${sobj}" \
		--output "${sampleId}.SEURAT__NORMALIZATION_${assayType}.rds" \
		--assay "${assayType}" \
		${(!(processParams.containsKey == null)) ? '--margin ' + processParams.margin: ''} \
		${(!(processParams.scalefactor == null)) ? '--scalefactor ' + processParams.scalefactor: ''} \
		${(!(processParams.normalizationMethod == null)) ? '--normalizationMethod ' + processParams.normalizationMethod: ''}
	"""
}
