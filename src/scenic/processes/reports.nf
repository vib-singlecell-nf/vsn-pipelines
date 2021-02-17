nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${reportTitle}
*/

toolParams = params.tools.scenic

process GENERATE_REPORT {

	container toolParams.container
	publishDir "${toolParams.scenicoutdir}/${sampleId}/notebooks", mode: 'link', overwrite: true
    label 'compute_resources__report'

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
	publishDir "${toolParams.scenicoutdir}/${sampleId}/notebooks", mode: 'link', overwrite: true
    label 'compute_resources__report'

	input:
		tuple val(sampleId), path(ipynb)

	output:
		tuple val(sampleId), path("*.html")

	script:
		"""
		jupyter nbconvert ${ipynb} --to html
		"""

}

