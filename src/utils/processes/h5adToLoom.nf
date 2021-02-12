nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__H5AD_TO_LOOM {

	container params.getToolParams("scanpy").container
    publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true, saveAs: { filename -> "${sampleId}.SCope_output.loom" }
    label 'compute_resources__mem'

	input:
		// Expects:
		// - rawFilteredData to be h5ad file containing the raw filtered (gene + cell filtered) data
		// - data to be one or more h5ad files containing the final results to be stored in the loom
		tuple \
			val(sampleId), \
			path(rawFilteredData), \
			path(data)

	output:
		tuple val(sampleId), \
		path("${sampleId}.SC__H5AD_TO_LOOM.loom")

	script:
		"""
		${binDir}/h5ad_to_loom.py \
			${(params.sc.containsKey('scope') && params.sc.scope.genome.length() > 0) ? '--nomenclature "' + params.sc.scope.genome + '"' : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_1.length() > 0 ) ? '--scope-tree-level-1 "' + params.sc.scope.tree.level_1 + '"'  : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_2.length() > 0 ) ? '--scope-tree-level-2 "' + params.sc.scope.tree.level_2 + '"'  : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.tree.level_3.length() > 0 ) ? '--scope-tree-level-3 "' + params.sc.scope.tree.level_3 + '"'  : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.containsKey('markers') && params.sc.scope.markers.log_fc_threshold.length() > 0 ) ? '--markers-log-fc-threshold ' + params.sc.scope.markers.log_fc_threshold : ''} \
			${(params.sc.containsKey('scope') && params.sc.scope.containsKey('markers') && params.sc.scope.markers.fdr_threshold.length() > 0 ) ? '--markers-fdr-threshold ' + params.sc.scope.markers.fdr_threshold : ''} \
			$data \
			$rawFilteredData \
			"${sampleId}.SC__H5AD_TO_LOOM.loom"
		"""

}

process SC__H5AD_TO_FILTERED_LOOM {

	container params.getToolParams("scanpy").container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple val(sampleId), path(f)

	output:
		tuple val(sampleId), path("${sampleId}.filtered.loom")

	script:
		"""
		${binDir}/h5ad_to_filtered_loom.py \
			$f \
			"${sampleId}.filtered.loom"
		"""

}
