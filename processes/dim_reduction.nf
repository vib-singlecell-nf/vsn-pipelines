nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  	binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  	binDir = ""
}

process SC__SCANPY__DIM_REDUCTION {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), path(f), val(nComps)

	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${processParams.off}")

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.dim_reduction.get(params.method))
		processParams = sampleParams.local
		method = processParams.dimReductionMethod.replaceAll('-','').toUpperCase()
		// Check if nComps is both dynamically and if statically set
		if(nComps && processParams.containsKey('nComps'))
			throw new Exception("SC__SCANPY__DIM_REDUCTION: nComps is both statically and dynamically set. Choose one.")
		_nComps = processParams.containsKey('nComps') ? processParams.nComps: ''
		if(nComps)
			_nComps = nComps
		
		"""
		${binDir}dim_reduction/sc_dim_reduction.py \
			--method ${processParams.dimReductionMethod} \
			${(processParams.containsKey('svdSolver')) ? '--svd-solver ' + processParams.svdSolver : ''} \
			${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			--n-comps ${_nComps} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
			${(processParams.containsKey('nJobs')) ? '--n-jobs ' + processParams.nJobs : ''} \
			${(processParams.containsKey('useFastTsne') && !processParams.useFastTsne) ? '' : '--use-fast-tsne'} \
			$f \
			"${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${processParams.off}"
		"""

}
