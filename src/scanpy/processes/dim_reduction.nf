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
class SC__SCANPY__DIM_REDUCTION_PARAMS {

	Script env = null;
	Map params = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	String method = null;
	// Parameters benchmarkable
    Integer nComps = null; ArrayList<Integer> nCompss = null;

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
\u001B[32m Parameter exploration of SC__SCANPY__DIM_REDUCTION step... \u001B[0m
\u001B[32m Tag: ${tag} \u001B[0m
\u001B[32m Parameters tested: \u001B[0m
\u001B[32m - nComps: \u001B[0m \u001B[33m     ${nComps instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${nComps} \u001B[0m
------------------------------------------------------------------
            """
        	}
		}
	}

	String getNCompsAsArgument(nComps) {
		// Check if nComps is both dynamically and if statically set
		if(!this.env.isParamNull(nComps) && this.configParams.containsKey('nComps'))
			throw new Exception("SC__SCANPY__DIM_REDUCTION: nComps is both statically (" + this.configParams["nComps"] + ") and dynamically (" + nComps + ") set. Choose one.")
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

def SC__SCANPY__DIM_REDUCTION_PARAMS(params) {
	return (new SC__SCANPY__DIM_REDUCTION_PARAMS(params))
}

process SC__SCANPY__DIM_REDUCTION {

	container params.sc.scanpy.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

	input:
		tuple \
			val(sampleId), \
			path(data), \
			val(stashedParams), \
			val(nComps)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${!isParamNull(stashedParams) ? uuid + '.' : ''}${processParams.off}"), \
			val(stashedParams), \
			val(nComps)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.dim_reduction.get(params.method))
		processParams = sampleParams.local
		// In parameter exploration mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		if(!isParamNull(stashedParams))
			uuid = stashedParams.findAll { it != 'NULL' }.join('_')
		method = processParams.method.replaceAll('-','').toUpperCase()
		// Cannot call constructor with parameter if nComps is not provided (aka NULL), type do not match
		def _processParams = new SC__SCANPY__DIM_REDUCTION_PARAMS()
		_processParams.setEnv(this)
		_processParams.setParams(params)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}/dim_reduction/sc_dim_reduction.py \
			--seed ${params.global.seed} \
			--method ${processParams.method} \
			${(processParams.containsKey('svdSolver')) ? '--svd-solver ' + processParams.svdSolver : ''} \
			${(processParams.containsKey('perplexity')) ? '--perplexity ' + processParams.perplexity : ''} \
			${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${_processParams.getNCompsAsArgument(nComps)} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
            --n-jobs ${task.cpus} \
			${(processParams.containsKey('useFastTsne')) ? '--use-fast-tsne ' + processParams.useFastTsne : ''} \
			$data \
			"${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${!isParamNull(stashedParams) ? uuid + '.' : ''}${processParams.off}"
		"""

}
