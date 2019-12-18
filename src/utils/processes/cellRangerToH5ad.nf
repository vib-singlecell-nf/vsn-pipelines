nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
	binDir = ""
}

process SC__CELLRANGER_TO_H5AD {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple val(sampleId), file(f)

	output:
		tuple val(sampleId), path("${sampleId}.SC__CELLRANGER_TO_H5AD.h5ad")

	when:
		sampleId != 'EMPTY'

	script:
		processParams = params.sc.file_converter
		"""
		${binDir}sc_file_converter.py \
			${(processParams.containsKey('tagCellWithSampleId')) ? '--sample-id ' + sampleId : ''} \
			--input-format "10x_mtx" \
			--output-format "h5ad" \
            ${f} \
            "${sampleId}.SC__CELLRANGER_TO_H5AD.h5ad"
		"""

}
