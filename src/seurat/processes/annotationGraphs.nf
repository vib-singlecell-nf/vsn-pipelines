#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = (params.global.standAlone != true) ? "${workflow.projectDir}/src/seurat/bin" : "${workflow.projectDir}/bin"

process annotation_graphs {
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "Plots/**"
	publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern: "allClusters_${sampleId}.xlsx"
	container params.Seurat.container
  input:
	tuple val(sampleId), file(seuratobj)
	file(markersfile)
	val(assay)
  output:
	file("Plots/***")
	file("allClusters_${sampleId}.xlsx")
  script:
  	def realmarkersfile = markersfile.name != 'NO_FILE' ? "--markersFile $markersfile" : ''
	"""
	Rscript ${scriptDir}/annotationGraphs.R --seuratObj ${seuratobj} \
	--assay $assay \
	$realmarkersfile
	"""
}

workflow SEURAT__ANNOTATION_GRAPHS {
	take:
		input
		assayname
	main:

		if(params.Seurat.annotationGraphs.containsKey(assayname)){
			assayParams = params.Seurat.annotationGraphs."${assayname}"
		} else {
			assayParams = params.Seurat.annotationGraphs
		}

		annotation_graphs(input,file(assayParams.markersFile),assayname)
}
