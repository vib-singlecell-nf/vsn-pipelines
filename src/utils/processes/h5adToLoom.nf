nextflow.preview.dsl=2
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin/" : ""

process SC__H5AD_TO_LOOM {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

	input:
		// Expects:
		// - rawFilteredData to be h5ad file containing the raw filtered (gene + cell filtered) data
		// - data to be the h5ad file containing the final results to be stored in the loom
		tuple val(sampleId), path(rawFilteredData), path(data)

	output:
		tuple val(sampleId), path("${sampleId}.SC__H5AD_TO_LOOM.loom")

	script:
		"""
		${binDir}h5ad_to_loom.py \
			$rawFilteredData \
			$data \
			${(params.sc.containsKey('scope') && params.sc.scope.genome.length() > 0) ? '--nomenclature "' + params.sc.scope.genome + '"' : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_1.length() > 0 ) ? '--scope-tree-level-1 ' + params.sc.scope.tree.level_1 : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_2.length() > 0 ) ? '--scope-tree-level-2 ' + params.sc.scope.tree.level_2 : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_3.length() > 0 ) ? '--scope-tree-level-3 ' + params.sc.scope.tree.level_3 : ''} \
			"${sampleId}.SC__H5AD_TO_LOOM.loom"
		"""

}

process SC__H5AD_TO_FILTERED_LOOM {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

	input:
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.filtered.loom")

	script:
		"""
		${binDir}h5ad_to_filtered_loom.py \
			$f \
			"${sampleId}.filtered.loom"
		"""

}
