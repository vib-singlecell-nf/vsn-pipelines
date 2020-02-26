nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

include '../../utils/processes/utils.nf'

/**
 * STATIC VERSION OF SCANPY MARKER GENES
 */
process SC__SCANPY__MARKER_GENES {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
  
  	input:
		// Expects 
		// - normalizedTransformedData to be an AnnData containing log normalised data
    	tuple \
			val(sampleId), \
			path(normalizedTransformedData), \
			path(clusteredData)
  
  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__MARKER_GENES.${processParams.off}")
  
  	script:
    	def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.marker_genes)
		processParams = sampleParams.local
		"""
		${binDir}cluster/sc_marker_genes.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('groupby')) ? '--groupby ' + processParams.groupby : ''} \
			${(processParams.containsKey('ngenes')) ? '--ngenes ' + processParams.ngenes : ''} \
			$normalizedTransformedData \
			$clusteredData \
			"${sampleId}.SC__SCANPY__MARKER_GENES.${processParams.off}"
		"""

}

/**
 * BENCHMARK VERSION OF SCANPY MARKER GENES
 */
process SC__SCANPY__PARAM_EXPLORE_MARKER_GENES {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate/markers/${isParamNull(clusteringMethod) ? "default": clusteringMethod.toLowerCase()}/${isParamNull(clusteringResolution) ? "res_": clusteringResolution}", mode: 'symlink', overwrite: true
  
  	input:
		// Expects 
		// - normalizedTransformedData to be an AnnData containing log normalised data
    	tuple \
			val(sampleId), \
			path(normalizedTransformedData), \
			path(clusteredData), \
			val(clusteringMethod), \
			val(clusteringResolution)
  
  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__PARAM_EXPLORE_MARKER_GENES.${uuid}.${processParams.off}"), \
			val(clusteringMethod), \
			val(clusteringResolution)
  
  	script:
    	def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.marker_genes)
		processParams = sampleParams.local
		// In parameter exploration mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		uuid = UUID.randomUUID().toString().substring(0,8)
		"""
		${binDir}cluster/sc_marker_genes.py \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${!isParamNull(clusteringMethod) ? '--groupby ' + clusteringMethod : ''} \
			${(processParams.containsKey('ngenes')) ? '--ngenes ' + processParams.ngenes : ''} \
			$normalizedTransformedData \
			$clusteredData \
			"${sampleId}.SC__SCANPY__PARAM_EXPLORE_MARKER_GENES.${uuid}.${processParams.off}"
		"""

}