nextflow.preview.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__H5AD_UPDATE_X_PCA {

	container params.sc.scanpy.container
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

	container params.sc.scanpy.container
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
