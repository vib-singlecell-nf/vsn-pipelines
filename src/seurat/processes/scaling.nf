#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__SCALING{
	//publishDir "${params.global.outDir}/${sampleId}", mode: 'symlink'
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	val assayType
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__SCALING_${assayType}.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.scaling)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/scaling.R --inputSeuratRds ${sobj} \
		--output "${sampleId}.SEURAT__SCALING_${assayType}.rds" \
		--assay $assayType \
		${(!(processParams.features == null)) ? '--features ' + processParams.features: ''} \
		${(!(processParams.varsToRegress == null)) ? '--varsToRegress ' + processParams.varsToRegress: ''} \
		${(!(processParams.splitBy == null)) ? '--splitBy ' + processParams.splitBy: ''} \
		${(!(processParams.modelUse == null)) ? '--modelUse ' + processParams.modelUse: ''} \
		${(!(processParams.useUmi == null)) ? '--useUmi ' + processParams.useUmi: ''} \
		${(!(processParams.doScale == null)) ? '--doScale ' + processParams.doScale: ''} \
		${(!(processParams.doCenter == null)) ? '--doCenter ' + processParams.doCenter: ''} \
		${(!(processParams.scaleMax == null)) ? '--scaleMax ' + processParams.scaleMax: ''} \
		${(!(processParams.blockSize == null)) ? '--blockSize ' + processParams.blockSize: ''} \
		${(!(processParams.minCellsToBlock == null)) ? '--minCellsToBlock ' + processParams.minCellsToBlock: ''} \
		${(!(processParams.verbose == null)) ? '--verbose ' + processParams.verbose: ''}
	"""
}
