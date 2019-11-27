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
	tuple val(sampleId), file(f)
	
	output:
	tuple val(sampleId), file("${sampleId}.SC__H5AD_TO_LOOM.loom")
	
	script:
	"""
	${binDir}h5ad_to_loom.py \
		$f \
		"${sampleId}.SC__H5AD_TO_LOOM.loom"
	"""

}

process SC__H5AD_TO_FILTERED_LOOM {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

	input:
	tuple val(sampleId), file(f)
	
	output:
	tuple val(sampleId), file("${sampleId}.filtered.loom")
	
	script:
	"""
	${binDir}h5ad_to_filtered_loom.py \
		$f \
		"${sampleId}.filtered.loom"
	"""

}
