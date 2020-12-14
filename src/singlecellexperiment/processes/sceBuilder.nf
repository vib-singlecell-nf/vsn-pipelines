#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__SCE_BUILDER {
  //publishDir "${params.ddir}/${samplename}", mode: 'copy'
  container params.sce.container
  input:
  tuple val(sampleId), val(cmat)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__SCE_BUILDER.rds")
  script:
  """
  Rscript ${scriptDir}/sceBuilder.R --inputMatrix ${cmat} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__SCE_BUILDER.rds"
  """
}
