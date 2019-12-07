nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${reportTitle}
*/

toolParams = params.sc.scenic

process GENERATE_REPORT {

	container toolParams.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${toolParams.scenicoutdir}/${sampleId}/notebooks", mode: 'link', overwrite: true

	input:
		file ipynb
		tuple val(sampleId), path(loom)
		val reportTitle

	output:
		tuple val(sampleId), path("${reportTitle}.ipynb")

	script:
		"""
		papermill ${ipynb} \
			${reportTitle}.ipynb \
			-p FILE $loom
		"""

}

process REPORT_TO_HTML {

	container toolParams.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${toolParams.scenicoutdir}/${sampleId}/notebooks", mode: 'link', overwrite: true

	input:
		tuple val(sampleId), path(ipynb)

	output:
		tuple val(sampleId), path("*.html")

	script:
		"""
		jupyter nbconvert ${ipynb} --to html
		"""

}
