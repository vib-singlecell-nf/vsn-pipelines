nextflow.enable.dsl=2

import java.nio.file.Paths
import static groovy.json.JsonOutput.*

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__H5AD_UPDATE_X_PCA {

	container params.getToolParams("scanpy").container
    label 'compute_resources__mem'

	input:
		tuple \
            val(sampleId), \
		    path(data), \
            path(xPca)

	output:
		tuple \
            val(sampleId), \
		    path("${sampleId}.SC__H5AD_UPDATE_X_PCA.h5ad")

	script:
		"""
		${binDir}/sc_h5ad_update.py \
			--x-pca ${xPca} \
			$data \
			"${sampleId}.SC__H5AD_UPDATE_X_PCA.h5ad"
		"""

}

process SC__H5AD_CLEAN {

	container params.getToolParams("scanpy").container
    label 'compute_resources__mem'

	input:
		tuple \
            val(sampleId), \
		    path(data), \
			val(stashedParams)

	output:
		tuple \
            val(sampleId), \
		    path("${sampleId}.SC__H5AD_CLEAN.h5ad"), \
			val(stashedParams)

	script:
		"""
		${binDir}/sc_h5ad_update.py \
			--empty-x \
			$data \
			"${sampleId}.SC__H5AD_CLEAN.h5ad"
		"""

}

process SC__H5AD_BEAUTIFY {

	container params.getToolParams("scanpy").container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple \
            val(sampleId), \
		    path(data), \
			val(stashedParams)

	output:
		tuple \
            val(sampleId), \
		    path("${sampleId}.SC__H5AD_BEAUTIFY.h5ad"), \
			val(stashedParams)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getUtilsParams("file_cleaner"))
        processParams = sampleParams.local

		obsColumnsToRemoveAsArgument = processParams.containsKey("obsColumnsToRemove") ? 
			processParams.obsColumnsToRemove.collect({ '--obs-column-to-remove' + ' ' + it }).join(' ') : 
			''
		"""
		${binDir}/sc_h5ad_update.py \
			${obsColumnsToRemoveAsArgument} \
			${processParams.containsKey("obsColumnMapper") ? "--obs-column-mapper '" + toJson(processParams.obsColumnMapper) + "'": ''} \
			${processParams.containsKey("obsColumnValueMapper") ? "--obs-column-value-mapper '" + toJson(processParams.obsColumnValueMapper) + "'": ''} \
			$data \
			"${sampleId}.SC__H5AD_BEAUTIFY.h5ad"
		"""

}

