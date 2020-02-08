nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast

@TupleConstructor()
class SC__SCANPY__CLUSTERING_PARAMS {
	String iff = null;
	String off = null;
	String clusteringMethod = null; ArrayList<String> clusteringMethods = null
    Float resolution = null; ArrayList<Float> resolutions = null
	String report_ipynb = null;

    // ArrayTuple get() {
	//    return *tuple(method,resolution)
    // }

	void displayMessage(tag) {
		Channel.from('').view {
			"""
------------------------------------------------------------------
\u001B[32m Benchmarking SC__SCANPY__CLUSTERING step... \u001B[0m
\u001B[32m Tag: ${tag} \u001B[0m
\u001B[32m Parameters tested: \u001B[0m
\u001B[32m - method: \u001B[0m \u001B[33m     ${clusteringMethods instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${clusteringMethods} \u001B[0m
\u001B[32m - resolution: \u001B[0m \u001B[33m ${resolutions instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${resolutions} \u001B[0m
------------------------------------------------------------------
            """
        }
	}

	// Define a function to check if the current process is running in benchmark mode
	boolean isBenchmarkMode() {
		return (clusteringMethods instanceof List
				|| resolutions instanceof List
		)
	}

	DataflowBroadcast $(tag) {
		// Prepare argument stream
		def $method = Channel.from(clusteringMethods == null ? "NULL" : clusteringMethods)
		def $resolution = Channel.from(resolutions == null ? "NULL" : resolutions)
		displayMessage(tag)
		return $method.combine($resolution)
	}
}

def SC__SCANPY__CLUSTERING_PARAMS(params) {
	return (new SC__SCANPY__CLUSTERING_PARAMS(params))
}

/**
 * DEFAULT VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__CLUSTERING {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.clustering)
		processParams = sampleParams.local
		"""
		${binDir}cluster/sc_clustering.py \
			${(processParams.containsKey('clusteringMethod')) ? '--method ' + processParams.clusteringMethod : ''} \
			${(processParams.containsKey('resolution')) ? '--resolution ' + processParams.resolution : ''} \
			$f \
			"${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}"
		"""

}

/**
 * BENCHMARK VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__BENCHMARK_CLUSTERING {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate/clustering/${method == null ? "default": method.toLowerCase()}/${resolution == null ? "default" : "res_" + resolution}", mode: 'symlink', overwrite: true

  	input:
    	tuple \
			val(sampleId), \
			path(f), \
			val(inertParams), \
			val(method), \
			val(resolution)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"), val(method), val(resolution)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.clustering)
		processParams = sampleParams.local
		"""
		${binDir}cluster/sc_clustering.py \
			${method != null ? '--method ' + method : ''} \
			${resolution != null ? '--resolution ' + resolution : ''} \
			$f \
			"${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"
		"""

}
