import groovy.transform.TupleConstructor
import nextflow.util.ArrayTuple

nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

@TupleConstructor()
class SC__SCANPY__DIM_REDUCTION_PARAMS {
    Integer nComps = null

	static String setNComps(nComps, processParams) {
		// Check if nComps is both dynamically and if statically set
		if(nComps != null && processParams.containsKey('nComps'))
			throw new Exception("SC__SCANPY__DIM_REDUCTION: nComps is both statically and dynamically set. Choose one.")
		if(nComps)
			return '--n-comps ' + nComps.replaceAll("\n","")
		return processParams.containsKey('nComps') ? '--n-comps ' + processParams.nComps: ''
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
		tuple val(sampleId), \
			path(data), \
			val(inertParams), \
			val(nComps)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${uuid}.${processParams.off}"), \
			val(inertParams)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.dim_reduction.get(params.method))
		processParams = sampleParams.local
		// In benchmark mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		uuid = UUID.randomUUID().toString().substring(0,8)
		method = processParams.dimReductionMethod.replaceAll('-','').toUpperCase()
		"""
		${binDir}dim_reduction/sc_dim_reduction.py \
			--method ${processParams.dimReductionMethod} \
			${(processParams.containsKey('svdSolver')) ? '--svd-solver ' + processParams.svdSolver : ''} \
			${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${SC__SCANPY__DIM_REDUCTION_PARAMS.setNComps(nComps, processParams)} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
			${(processParams.containsKey('nJobs')) ? '--n-jobs ' + processParams.nJobs : ''} \
			${(processParams.containsKey('useFastTsne') && processParams.useFastTsne) ? '--use-fast-tsne' : ''} \
			$data \
			"${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${uuid}.${processParams.off}"
		"""

}
