nextflow.enable.dsl=2

import java.nio.file.Paths
import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.util.ArrayTuple
import nextflow.script.ScriptBinding

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCANPY__NEIGHBORHOOD_GRAPH_PARAMS {

	Script env = null;
	Map params = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	// Parameters benchmarkable
    Integer nPcs = null; ArrayList<Integer> nPcss = null;
	Integer nNeighbors = null; ArrayList<Integer> nNeighborss = null;

	void setEnv(env) {
		this.env = env
	}

	void setParams(params) {
		this.params = params
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	void displayMessage(tag = null) {
		if(!this.params?.quiet) {
			Channel.from('').view {
				"""
	------------------------------------------------------------------
	\u001B[32m Parameter exploration of SC__SCANPY__NEIGHBORHOOD_GRAPH step... \u001B[0m
	\u001B[32m Tag: ${tag} \u001B[0m
	\u001B[32m Parameters tested: \u001B[0m
	\u001B[32m - nPcs: \u001B[0m \u001B[33m     ${nPcss instanceof List} \u001B[0m
	\u001B[32m   - values: \u001B[0m \u001B[33m   ${nPcss} \u001B[0m
	------------------------------------------------------------------
				"""
			}
		}
	}

	String getNPcsAsArgument(nPcs) {
		// Check if nPcs is both dynamically and if statically set
		if(!this.env.isParamNull(nPcs) && this.configParams.containsKey('nPcs'))
			throw new Exception("SC__SCANPY__NEIGHBORHOOD_GRAPH: nPcs is both statically (" + this.configParams["nPcs"] + ") and dynamically (" + nPcs + ") set. Choose one.")
		if(!this.env.isParamNull(nPcs))
			return '--n-pcs ' + nPcs.replaceAll("\n","")
		return this.configParams.containsKey('nPcs') ? '--n-pcs ' + this.configParams.nPcs: ''
	}

	// Define a function to check if the current process is running in parameter exploration mode
	boolean isParameterExplorationModeOn() {
		return (nPcs instanceof List)
	}

	DataflowBroadcast $(tag = null) {
		// Prepare argument stream
		def $nPcs = Channel.from("NULL")
		if(isParameterExplorationModeOn()) {
			displayMessage(tag)
			$nPcs = Channel.from(nPcs)
		}
		return $nPcs
	}

	ArrayTuple asTuple() {
	   	return tuple(nPcs)
    }

}

process SC__SCANPY__NEIGHBORHOOD_GRAPH {

  	container params.getToolParams("scanpy").container
    label 'compute_resources__mem'

  	input:
        tuple \
            val(sampleId), \
            path(f), \
            val(stashedParams), \
			val(nPcs)

	output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SCANPY__NEIGHBORHOOD_GRAPH.${processParams.off}"), \
            val(stashedParams),
			val(nPcs)

	script:
        def sampleParams = params.parseConfig(
			sampleId,
			params.global,
			params.getToolParams("scanpy").neighborhood_graph
		)
		processParams = sampleParams.local
        // In parameter exploration mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		if(!isParamNull(stashedParams))
			uuid = stashedParams.findAll { it != 'NULL' }.join('_')
		// Cannot call constructor with parameter if nPcs is not provided (aka NULL), type do not match
		def _processParams = new SC__SCANPY__NEIGHBORHOOD_GRAPH_PARAMS()
		_processParams.setEnv(this)
		_processParams.setParams(params)
		_processParams.setConfigParams(processParams)
        """
        ${binDir}/nn/sc_neighborhood_graph.py \
            $f \
            ${sampleId}.SC__SCANPY__NEIGHBORHOOD_GRAPH.${processParams.off} \
			--seed ${params.global.seed} \
            ${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${_processParams.getNPcsAsArgument(nPcs)}
        """

}
