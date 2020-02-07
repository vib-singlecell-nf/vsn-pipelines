nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

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
  	publishDir "${params.global.outdir}/data/intermediate/clustering/${method == "NULL" ? "default": method.toLowerCase()}/${resolution == "NULL" ? "res_": resolution}", mode: 'symlink', overwrite: true

  	input:
    	tuple \
			val(sampleId), \
			path(f), \
			val(method), \
			val(resolution)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"), val(method), val(resolution)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.clustering)
		processParams = sampleParams.local
		"""
		${binDir}cluster/sc_clustering.py \
			${method != "NULL" ? '--method ' + method : ''} \
			${resolution != "NULL" ? '--resolution ' + resolution : ''} \
			$f \
			"${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"
		"""

}
