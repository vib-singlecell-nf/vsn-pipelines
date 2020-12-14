#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__FIND_VARIABLE_FEATURES{
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	val assayType
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__FIND_VARIABLE_FEATURES_${assayType}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.findVariableFeatures)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/findVariableFeatures.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__FIND_VARIABLE_FEATURES_${assayType}.rds" \
		--assay "$assayType" \
		${(!(processParams.selectionMethod == null)) ? '--selectionMethod ' + processParams.selectionMethod: ''} \
		${(!(processParams.loesSpan == null)) ? '--loesSpan ' + processParams.loesSpan: ''} \
		${(!(processParams.clipMax == null)) ? '--clipMax ' + processParams.clipMax: ''} \
		${(!(processParams.meanFunction == null)) ? '--meanFunction ' + processParams.meanFunction: ''} \
		${(!(processParams.dispersionFunction == null)) ? '--dispersionFunction ' + processParams.dispersionFunction: ''} \
		${(!(processParams.numBin == null)) ? '--numBin ' + processParams.numBin: ''} \
		${(!(processParams.binningMethod == null)) ? '--binningMethod ' + processParams.binningMethod: ''} \
		${(!(processParams.nfeatures == null)) ? '--nfeatures ' + processParams.nfeatures: ''} \
		${(!(processParams.meanCutoff == null)) ? '--meanCutoff ' + processParams.meanCutoff: ''} \
		${(!(processParams.dispersionCutoff == null)) ? '--dispersionCutoff ' + processParams.dispersionCutoff: ''} \
		${(!(processParams.verbose == null)) ? '--verbose ' + processParams.verbose: ''}
	"""
}
