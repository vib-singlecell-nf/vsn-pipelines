#!/usr/bin/env nextflow
nextflow.preview.dsl=2
scriptDir = "${workflow.projectDir}/src/singlecellexperiment/bin"

process SINGLE_CELL_EXPERIMENT__PCA_FILTERING {
  publishDir "${params.global.outdir}/${sampleId}", mode: 'symlink', pattern:"Plots/RNA/**"
  container params.sce.container
  input:
  tuple val(sampleId), file(scerds)
  output:
  tuple val(sampleId), file("${sampleId}.SINGLE_CELL_EXPERIMENT__PCA_FILTERING.rds")
  file("Plots/RNA/**")
  script:
  def sampleParams = params.parseConfig(sampleId, params.global, params.sce.pcaFiltering)
		processParams = sampleParams.local
  """
  Rscript ${scriptDir}/pcaFiltering.R --sceObj ${scerds} \
		--output "${sampleId}.SINGLE_CELL_EXPERIMENT__PCA_FILTERING.rds" \
		--scriptFunctions ${scriptDir}/script_functions_BioIT.R \
		${(!(processParams.removeOutliers == null)) ? '--removeOutliers': ''}
  """
}
