nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*

process UTILS__GENERATE_WORKFLOW_CONFIG_REPORT {

  	container params.utils.container
  	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

    input:
        path(ipynb)

	output:
		path("workflow_configuration_report.ipynb")

	script:
		"""
		papermill ${ipynb} \
			workflow_configuration_report.ipynb \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}

process UTILS__REPORT_TO_HTML {

	container params.utils.container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notbooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true
	label 'compute_resources__report'

	input:
		tuple \
			val(sampleId), \
			path(ipynb)

	output:
		file("*.html")

	script:
		"""
		jupyter nbconvert ${ipynb} --to html
		"""

}
