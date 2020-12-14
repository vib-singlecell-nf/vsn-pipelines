#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__THRESHOLDFILTERING{
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/RNA/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__THRESHOLDFILTERING.rds")
	file("Plots/RNA/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.thresholdFiltering)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/thresholdFiltering.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__THRESHOLDFILTERING.rds" \
		--scriptFunctions ${scriptDir}/script_functions.R \
		${(!(processParams.nmad_low_feature == null)) ? '--nmad_low_feature ' + processParams.nmad_low_feature: ''} \
		${(processParams.nmad_high_feature == null) ? '' : '--nmad_high_feature ' + processParams.nmad_high_feature} \
		${(processParams.nmad_low_UMI == null) ? '' : '--nmad_low_UMI ' + processParams.nmad_low_UMI} \
		${(processParams.nmad_high_UMI == null) ? '' : '--nmad_high_UMI ' + processParams.nmad_high_UMI} \
		${(processParams.nmad_high_mito == null) ? '' : '--nmad_high_mito ' + processParams.nmad_high_mito} \

	"""
}
