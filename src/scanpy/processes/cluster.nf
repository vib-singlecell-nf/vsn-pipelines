nextflow.enable.dsl=2

import java.nio.file.Paths
import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCANPY__CLUSTERING_PARAMS {

	Script env = null;
	Map params = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	String report_ipynb = null;
	boolean preflight_checks = null;
	// Parameters benchmarkable
	String method = null; ArrayList<String> methods = null;
    Float resolution = null; ArrayList<Float> resolutions = null;

	void setEnv(env) {
		this.env = env
	}

	void setParams(params) {
		this.params = params
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	void displayMessage(tag) {
		if(!this.params?.quiet) {
        	Channel.from('').view {
			"""
------------------------------------------------------------------
\u001B[32m Parameter exploration of SC__SCANPY__CLUSTERING step... \u001B[0m
\u001B[32m Tag: ${tag} \u001B[0m
\u001B[32m Parameters tested: \u001B[0m
\u001B[32m - method: \u001B[0m \u001B[33m     ${methods instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${methods} \u001B[0m
\u001B[32m - resolution: \u001B[0m \u001B[33m ${resolutions instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${resolutions} \u001B[0m
------------------------------------------------------------------
            """
        	}
		}
	}

	String getMethodAsArgument(method) {
		return !this.env.isParamNull(method) ? '--method ' + method : ''
	}

	String getResolutionAsArgument(resolution) {
		return !this.env.isParamNull(resolution) ? '--resolution ' + resolution : ''
	}

	int numParamsBenchmarked() {
		def paramsBenchmarked = [ 
			methods instanceof List,
			resolutions instanceof List
		]
		def sum = { result, i -> result + (i ? 1 : 0) }
		return paramsBenchmarked.inject(0, sum)
	}

	int numParams() {
		return 2 // Total number of parameters implemented for benchmarking
	}

	// Define a function to check if the current process is running in parameter exploration mode
	boolean isParameterExplorationModeOn() {
		return (methods instanceof List
				|| resolutions instanceof List
		)
	}

	DataflowBroadcast $(tag) {
		// Prepare argument stream
		def _method = method == null ? "NULL" : method
		def _resolution = resolution == null ? "NULL" : resolution
		def $method = Channel.from(methods == null ? _method : methods)
		def $resolution = Channel.from(resolutions == null ? _resolution : resolutions)
		displayMessage(tag)
		return $method.combine($resolution)
	}

}

process SC__SCANPY__CLUSTERING_PREFLIGHT_CHECKS {

	container params.getToolParams("scanpy").container
	label 'compute_resources__mem'

	input:
    	tuple \
			val(sampleId), \
			path(f)

	output:
		tuple \
			val(sampleId), \
			path(f)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("scanpy").clustering)
		processParams = sampleParams.local
		methodAsArguments = processParams?.methods ? processParams.methods.collect({ '--method' + ' ' + it }).join(' ') : '--method ' + processParams.method
		resolutionAsArguments = processParams?.resolutions ? processParams?.resolutions.collect({ '--resolution' + ' ' + it }).join(' ') : '--resolution ' + processParams.resolution
		"""
		${binDir}/cluster/sc_clustering_preflight_checks.py \
			--seed ${params.global.seed} \
			${methodAsArguments} \
			${resolutionAsArguments} \
			$f \
		"""

}

def SC__SCANPY__CLUSTERING_PARAMS(params) {
	return (new SC__SCANPY__CLUSTERING_PARAMS(params))
}

/**
 * DEFAULT VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__CLUSTERING {

  	container params.getToolParams("scanpy").container
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("scanpy").clustering)
		processParams = sampleParams.local
		"""
		${binDir}/cluster/sc_clustering.py \
			--seed ${params.global.seed} \
			${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
			${(processParams.containsKey('resolution')) ? '--resolution ' + processParams.resolution : ''} \
			$f \
			"${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}"
		"""

}

/**
 * BENCHMARK VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__PARAM_EXPLORE_CLUSTERING {

  	container params.getToolParams("scanpy").container
  	publishDir "${params.global.outdir}/data/intermediate/clustering/${isParamNull(method) ? "default": method.toLowerCase()}/${isParamNull(resolution) ? "default" : "res_" + resolution}", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

  	input:
    	tuple \
			val(sampleId), \
			path(f), \
			val(stashedParams), \
			val(method), \
			val(resolution)

  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__PARAM_EXPLORE_CLUSTERING.${processParams.off}"), \
			val(method), \
			val(resolution)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("scanpy").clustering)
		processParams = sampleParams.local
		def _processParams = new SC__SCANPY__CLUSTERING_PARAMS()
		_processParams.setEnv(this)
		_processParams.setParams(params)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}/cluster/sc_clustering.py \
			--seed ${params.global.seed} \
			${_processParams.getMethodAsArgument(method)} \
			${_processParams.getResolutionAsArgument(resolution)} \
			$f \
			"${sampleId}.SC__SCANPY__PARAM_EXPLORE_CLUSTERING.${processParams.off}"
		"""

}
