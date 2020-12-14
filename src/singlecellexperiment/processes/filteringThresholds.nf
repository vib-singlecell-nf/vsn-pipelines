#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__FILTERING_THRESHOLDS {
  //publishDir "${params.ddir}/${samplename}", mode: 'copy'
  container params.sce.container
  input:
  tuple val(sampleId), file(scerds)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__FILTERING_THRESHOLDS.rds")
  script:
  def sampleParams = params.parseConfig(sampleId, params.global, params.sce.filteringThresholds)
		processParams = sampleParams.local
  """
  Rscript ${scriptDir}/filteringThresholds.R --rdsFile ${scerds} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__FILTERING_THRESHOLDS.rds" \
		${(!(processParams.nbMADs == null)) ? '--nmads ' + processParams.nbMADs: ''} \
		${(!(processParams.useMitochondrialGenes == null)) ? '--mito': ''}
  """
}
