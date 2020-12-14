#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__REPORT_GENERATOR {
	publishDir "${params.outdir}/${sampleId}", mode: 'move', pattern: "${sampleId}_report.html"
	container params.sce.container
	input:
	tuple val(sampleId), file(rdsfile)
	output:
	file("${sampleId}_report.html")
	script:
	"""
	cp -L ${scriptDir}/reportGenerator.Rmd report.Rmd
	cp -L ${rdsfile} final.rds
	R -e 'rmarkdown::render("report.Rmd",output_file = "${sampleId}_report.html",params = list(input="final.rds",id="${sampleid}"))'
	rm report.Rmd final.rds
	"""
}
