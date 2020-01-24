nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

process SC__SCANPY__FEATURE_SELECTION {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__FEATURE_SELECTION.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.feature_selection)
		processParams = sampleParams.local
		"""
		${binDir}feature_selection/sc_select_variable_genes.py \
			--method ${processParams.featureSelectionMethod} \
			${(processParams.containsKey('featureSelectionMinMean')) ? '--min-mean ' + processParams.featureSelectionMinMean : ''} \
			${(processParams.containsKey('featureSelectionMaxMean')) ? '--max-mean ' + processParams.featureSelectionMaxMean : ''} \
			${(processParams.containsKey('featureSelectionMinDisp')) ? '--min-disp ' + processParams.featureSelectionMinDisp : ''} \
			${(processParams.containsKey('featureSelectionMaxDisp')) ? '--max-disp ' + processParams.featureSelectionMaxDisp : ''} \
			$f \
			"${sampleId}.SC__SCANPY__FEATURE_SELECTION.${processParams.off}"
		"""

}
