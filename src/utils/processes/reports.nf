nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

process UTILS__GENERATE_WORKFLOW_CONFIG_REPORT {

  	container params.utils.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true

    input:
        path(ipynb)

	output:
		path("workflow_configuration_report.ipynb")

	script:
		"""
		papermill ${ipynb} \
			workflow_configuration_report.ipynb \
            -p WORKFLOW_MANIFEST '${toJson(workflow.manifest)}' \
			-p WORKFLOW_PARAMETERS '${toJson(params)}'
		"""

}
