#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__DIMENSIONALITY_REDUCTION_PCA {
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(seuratobj)
	val assay
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__DIMENSIONALITY_REDUCTION_PCA_${assay}.rds")
	file("Plots/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.dimensionalityReduction.pca)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/pca.R --seuratObj ${seuratobj} \
		--output "${sampleId}.SEURAT__DIMENSIONALITY_REDUCTION_PCA_${assay}.rds" \
		--assay $assay \
		${(processParams.nPcs == null) ? '' : '--nPcs ' + processParams.nPcs} \
		${(processParams.nPlotedPcs == null) ? '' : '--nPlotedPcs ' + processParams.nPlotedPcs} \
		${(processParams.scaleData == null || processParams.scaleData != "true") ? '' :'--scaleData'} \
		${(processParams.removeScaledData == null || processParams.removeScaledData != "true") ? '' :'--removeScaledData'} \
		${(processParams.excludePatternHVG == null) ? '' : '--excludePatternHVG ' + processParams.excludePatternHVG} \
		${(processParams.diagnosticPlots == null) ? '' : '--diagnosticPlots'}
	"""
}

process SEURAT__DIMENSIONALITY_REDUCTION_TSNE {
	//publishDir "${params.out_dir}/${sampleId}", mode: 'copy'
	container params.Seurat.container
	input:
	tuple val(sampleId), file(seuratobj)
	val assayType
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__DIMENSIONALITY_REDUCTION_TSNE_${assayType}.rds")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.dimensionalityReduction.tsne)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/tsne.R --inputSeuratRds ${seuratobj} \
		--output "${sampleId}.SEURAT__DIMENSIONALITY_REDUCTION_TSNE_${assayType}.rds" \
		--assay "${assayType}" \
		${(!(processParams.dims == null)) ? '--dims ' + processParams.dims: ''} \
		${(!(processParams.reduction == null)) ? '--reduction ' + processParams.reduction: ''} \
		${(!(processParams.seedUse == null)) ? '--seedUse ' + processParams.seedUse: ''} \
		${(!(processParams.tsneMethod == null)) ? '--tsneMethod ' + processParams.tsneMethod: ''} \
		${(!(processParams.addIter == null)) ? '--addIter ' + processParams.addIter: ''} \
		${(!(processParams.dimEmbed == null)) ? '--dimEmbed ' + processParams.dimEmbed: ''} \
		${(!(processParams.perplexity == null)) ? '--perplexity ' + processParams.perplexity: ''}
	"""
}

/*
process SEURAT__DIMANSIONALITY_REDUCTION_UMAP {
	//publishDir "${params.out_dir}/${sampleId}", mode: 'copy'
	input:
	tuple val(sampleId), file(seuratobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__DIMANSIONALITY_REDUCTION_UMAP.rds")
	script:
	"""
	Rscript ${params.pdir}/modules/dimensionalityReduction/umap.R --inputSeuratRds ${seuratobj} \
		--output "${sampleId}.SEURAT__DIMANSIONALITY_REDUCTION_UMAP.rds"
	"""
}
*/
