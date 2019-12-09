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
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${processParams.off}")

	script:
		processParams = params.sc.scanpy.dim_reduction.get(params.method)
		method = processParams.dimReductionMethod.replaceAll('-','').toUpperCase()
		"""
		${binDir}dim_reduction/sc_dim_reduction.py \
			--method ${processParams.dimReductionMethod} \
			${(processParams.containsKey('svdSolver')) ? '--svd-solver ' + processParams.svdSolver : ''} \
			${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''} \
			${(processParams.containsKey('nComps')) ? '--n-comps ' + processParams.nComps : ''} \
			${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
			${(processParams.containsKey('nJobs')) ? '--n-jobs ' + processParams.nJobs : ''} \
			${(processParams.containsKey('useFastTsne') && !processParams.useFastTsne) ? '' : '--use-fast-tsne'} \
			$f \
			"${sampleId}.SC__SCANPY__DIM_REDUCTION_${method}.${processParams.off}"
		"""

}
