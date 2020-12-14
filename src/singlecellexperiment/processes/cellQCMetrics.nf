#!/usr/bin/env nextflow
nextflow.preview.dsl=2

scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__CELL_QC_METRICS {
  //publishDir "${params.out_dir}/${sampleId}", mode: 'copy'
  container params.sce.container
  input:
  tuple val(sampleId), file(scerds)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__CELL_QC_METRICS.rds")
  script:
  def sampleParams = params.parseConfig(sampleId, params.global, params.sce.cellQCMetrics)
		processParams = sampleParams.local
  """
  Rscript ${scriptDir}/cellQCMetrics.R --rdsFile ${scerds} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__CELL_QC_METRICS.rds" \
		${(!(processParams.useMitochondrialGenes == null)) ? '--mito': ''}
  """
}
