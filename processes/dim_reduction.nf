nextflow.preview.dsl=2

import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.util.ArrayTuple
import nextflow.script.ScriptBinding

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

include '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCANPY__DIM_REDUCTION_PARAMS {

	Script env = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	String dimReductionMethod = null;
	// Parameters benchmarkable
    Integer nComps = null; ArrayList<Integer> nCompss = null;

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
\u001B[32m Benchmarking SC__SCANPY__DIM_REDUCTION step... \u001B[0m
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
			throw new Exception("SC__SCANPY__DIM_REDUCTION: nComps is both statically (" + nComps + ") and dynamically (" + this.configParams["nComps"] + ") set. Choose one.")
		if(!this.env.isParamNull(nComps))
			return '--n-comps ' + nComps.replaceAll("\n","")
		return this.configParams.containsKey('nComps') ? '--n-comps ' + this.configParams.nComps: ''
	}

	// Define a function to check if the current process is running in benchmark mode
	boolean isBenchmarkMode() {
		return (nComps instanceof List)
	}

	DataflowBroadcast $(tag = null) {
		// Prepare argument stream
		def $nComps = Channel.from("NULL")
		if(isBenchmarkMode()) {
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
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple \
			val(sampleId), \
			path(data), \
			val(inertParams), \
			val(nComps)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${!isParamNull(inertParams) ? uuid + '.' : ''}${processParams.off}"), \
			val(inertParams)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.dim_reduction.get(params.method))
		processParams = sampleParams.local
		// In benchmark mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		// Output file will only be tagged with UUID if in benchmark mode
		uuid = UUID.randomUUID().toString().substring(0,8)
		method = processParams.dimReductionMethod.replaceAll('-','').toUpperCase()
		// Cannot call constructor with parameter if nComps is not provided (aka NULL), type do not match
		def _processParams = new SC__SCANPY__DIM_REDUCTION_PARAMS()
		_processParams.setEnv(this)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}dim_reduction/sc_dim_reduction.py \
			--method ${processParams.dimReductionMethod} \
			${(processParams.containsKey('svdSolver')) ? '--svd-solver ' + processParams.svdSolver : ''} \
			${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${_processParams.getNCompsAsArgument(nComps)} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
			${(processParams.containsKey('nJobs')) ? '--n-jobs ' + processParams.nJobs : ''} \
			${(processParams.containsKey('useFastTsne') && processParams.useFastTsne) ? '--use-fast-tsne' : ''} \
			$data \
			"${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${!isParamNull(inertParams) ? uuid + '.' : ''}${processParams.off}"
		"""

}
