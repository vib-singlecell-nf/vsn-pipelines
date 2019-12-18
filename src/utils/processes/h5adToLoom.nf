nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
	binDir = ""
}


process SC__H5AD_TO_LOOM {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

	input:
		// Expects:
		// - rawFilteredData to be h5ad file containing the raw filtered (gene + cell filtered) data
		// - data to be the h5ad file containing the final results to be stored in the loom
		tuple val(sampleId), file(data), file(rawFilteredData)

	output:
		tuple val(sampleId), path("${sampleId}.SC__H5AD_TO_LOOM.loom")
	
	when:
		sampleId != 'EMPTY'

	script:
		"""
		${binDir}h5ad_to_loom.py \
			$rawFilteredData \
			$data \
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
