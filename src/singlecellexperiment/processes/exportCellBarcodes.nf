#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"


process SINGLE_CELL_EXPERIMENT__EXPORT_CELLBARCODES {
  //publishDir "${params.ddir}/${samplename}", mode: 'copy'
  container params.sce.container
  input:
  tuple val(sampleId), file(cmat)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__EXPORT_CELLBARCODES.csv")
  script:
  """
  Rscript ${scriptDir}/exportCellBarcodes.R --rdsFile ${cmat} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__EXPORT_CELLBARCODES.csv"
  """
}
