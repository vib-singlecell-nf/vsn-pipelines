nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__H5AD_MERGE {

	container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		// Expects:
		// - data to be multiple h5ad files containing the final results to be merged
		tuple \
            val(sampleId), \
			path(data)

	output:
		tuple \
            val(sampleId), \
		    path("${sampleId}.SC__H5AD_MERGE.h5ad")

	script:
		"""
        ${binDir}/sc_h5ad_merge.py \
            * \
            "${sampleId}.SC__H5AD_MERGE.h5ad"
		"""

}
