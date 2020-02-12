nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin/" : ""

process SC__H5AD_UPDATE_X_PCA {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

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
		${binDir}sc_h5ad_update.py \
			--x-pca ${xPca} \
			$data \
			"${sampleId}.SC__H5AD_UPDATE_X_PCA.h5ad"
		"""

}