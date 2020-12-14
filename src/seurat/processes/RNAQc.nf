#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process SEURAT__RNA_QC{
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/RNA/**"
	container params.Seurat.container
	input:
	tuple val(sampleId), file(sobj)
	output:
	tuple val(sampleId), file("${sampleId}.SEURAT__RNA_QC.rds")
	file("Plots/RNA/**")
	script:
	def sampleParams = params.parseConfig(sampleId, params.global,params.Seurat.RNAQc)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/RNAQc.R --seuratObj "${sobj}" \
		--output "${sampleId}.SEURAT__RNA_QC.rds" \
		${(processParams.mitoGenes == null || processParams.mitoGenes == "false") ? '' : '--mitoGenes'} \
	"""
}


//${(processParams.covidGenes == null || processParams.covidGenes == "false") ? '' : '--covidGenes'} \
//${(processParams.rbcGenes == null || processParams.rbcGenes == "false") ? '' : '--rbcGenes'} \
//${(processParams.genomeName1 == null) ? '' : '--genomeName1 ' + processParams.genomeName1} \
//${(processParams.genomeName2 == null) ? '' : '--genomeName2 ' + processParams.genomeName2}
