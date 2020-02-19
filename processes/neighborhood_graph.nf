nextflow.preview.dsl=2

import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.util.ArrayTuple
import nextflow.script.ScriptBinding

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

include '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCANPY__NEIGHBORHOOD_GRAPH_PARAMS {

	Script env = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	// Parameters benchmarkable
    Integer nComps = null; ArrayList<Integer> nCompss = null;
	Integer nNeighbors = null; ArrayList<Integer> nNeighborss = null;

	void setEnv(env) {
		this.env = env
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	void displayMessage(tag = null) {
		Channel.from('').view {
			"""
------------------------------------------------------------------
\u001B[32m Benchmarking SC__SCANPY__NEIGHBORHOOD_GRAPH step... \u001B[0m
\u001B[32m Tag: ${tag} \u001B[0m
\u001B[32m Parameters tested: \u001B[0m
\u001B[32m - nComps: \u001B[0m \u001B[33m     ${nComps instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${nComps} \u001B[0m
------------------------------------------------------------------
            """
        }
	}

	String getNCompsAsArgument(nComps) {
		// Check if nComps is both dynamically and if statically set
		if(!this.env.isParamNull(nComps) && this.configParams.containsKey('nComps'))
			throw new Exception("SC__SCANPY__NEIGHBORHOOD_GRAPH: nComps is both statically (" + nComps + ") and dynamically (" + this.configParams["nComps"] + ") set. Choose one.")
		if(!this.env.isParamNull(nComps))
			return '--n-comps ' + nComps.replaceAll("\n","")
		return this.configParams.containsKey('nComps') ? '--n-comps ' + this.configParams.nComps: ''
	}

	// Define a function to check if the current process is running in parameter exploration mode
	boolean isParameterExplorationModeOn() {
		return (nComps instanceof List)
	}

	DataflowBroadcast $(tag = null) {
		// Prepare argument stream
		def $nComps = Channel.from("NULL")
		if(isParameterExplorationModeOn()) {
			displayMessage(tag)
			$nComps = Channel.from(nComps)
		}
		return $nComps
	}

	ArrayTuple asTuple() {
	   	return tuple(nComps)
    }

}

process SC__SCANPY__NEIGHBORHOOD_GRAPH {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

  	input:
        tuple \
            val(sampleId), \
            path(f), \
            val(stashedParams), \
			val(nComps)

	output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SCANPY__NEIGHBORHOOD_GRAPH.${processParams.off}"), \
            val(stashedParams),
			val(nComps)

	script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.filter)
		processParams = sampleParams.local
        // In parameter exploration mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		// Output file will only be tagged with UUID if in parameter exploration mode
		uuid = UUID.randomUUID().toString().substring(0,8)
		// Cannot call constructor with parameter if nComps is not provided (aka NULL), type do not match
		def _processParams = new SC__SCANPY__NEIGHBORHOOD_GRAPH_PARAMS()
		_processParams.setEnv(this)
		_processParams.setConfigParams(processParams)
        """
        ${binDir}nn/sc_neighborhood_graph.py \
            $f \
            ${sampleId}.SC__SCANPY__NEIGHBORHOOD_GRAPH.${processParams.off} \
			${'--seed ' + (params.global.containsKey('seed') ? params.global.seed: params.seed)} \
            ${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${_processParams.getNCompsAsArgument(nComps)}
        """

}
