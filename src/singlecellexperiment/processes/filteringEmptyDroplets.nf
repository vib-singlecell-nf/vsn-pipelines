#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__FILTERING_EMPTYDROPLETS {
  //publishDir "${params.ddir}/${samplename}", mode: 'copy'
  container params.sce.container
  input:
  tuple val(sampleId), file(cmat)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__FILTERING_EMPTYDROPLETS.rds")
  script:
  def sampleParams = params.parseConfig(sampleId, params.global, params.sce.filteringEmptyDroplets)
		processParams = sampleParams.local
  """
  Rscript ${scriptDir}/filteringEmptyDroplets.R --rdsFile ${cmat} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__FILTERING_EMPTYDROPLETS.rds" \
		${(!(processParams.fdrThreshold == null)) ? '--fdr-threshold ' + processParams.fdrThreshold: ''} \
		${(!(processParams.nbUMILowerBound == null)) ? '--lower ' + processParams.nbUMILowerBound: ''}
  """
}
