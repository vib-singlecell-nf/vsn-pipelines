nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__H5AD_TO_LOOM {

	container params.tools.scanpy.container
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
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__H5AD_TO_LOOM.loom")

	script:
		"""
		${binDir}/h5ad_to_loom.py \
			${(params.utils?.scope.genome.length() > 0) ? '--nomenclature "' + params.utils?.scope.genome + '"' : ''} \
			${(params.utils?.scope.tree.level_1.length() > 0 ) ? '--scope-tree-level-1 "' + params.utils.scope.tree.level_1 + '"'  : ''} \
			${(params.utils?.scope.tree.level_2.length() > 0 ) ? '--scope-tree-level-2 "' + params.utils.scope.tree.level_2 + '"'  : ''} \
			${(params.utils?.scope.tree.level_3.length() > 0 ) ? '--scope-tree-level-3 "' + params.utils.scope.tree.level_3 + '"'  : ''} \
			${(params.utils?.scope.?.markers?.log_fc_threshold) ? '--markers-log-fc-threshold ' + params.utils.scope.markers.log_fc_threshold : ''} \
			${(params.utils?.scope.?.markers?.fdr_threshold) ? '--markers-fdr-threshold ' + params.utils.scope.markers.fdr_threshold : ''} \
			$data \
			$rawFilteredData \
			"${sampleId}.SC__H5AD_TO_LOOM.loom"
		"""

}

process SC__H5AD_TO_FILTERED_LOOM {

	container params.tools.scanpy.container
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
