nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${reportTitle}
*/
process SC__SCANPY__GENERATE_REPORT {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true

	input:
		file ipynb
		tuple val(sampleId), path(adata)
		val(reportTitle)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.ipynb")

	script:
		def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.ipynb \
			-p FILE $adata \
			-p WORKFLOW_MANIFEST '${toJson(workflow.manifest)}' \
			-p WORKFLOW_PARAMETERS '${toJson(paramsCopy)}'
		"""

}

// QC report takes two inputs, so needs it own process
process SC__SCANPY__GENERATE_DUAL_INPUT_REPORT {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true

  	input:
		file(ipynb)
		tuple val(sampleId), file(data1), file(data2)
		val reportTitle

  	output:
    	tuple val(sampleId), file("${sampleId}.${reportTitle}.ipynb")

  	script:
	  	def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.ipynb \
			-p FILE1 $data1 -p FILE2 $data2 \
			-p WORKFLOW_MANIFEST '${toJson(workflow.manifest)}' \
			-p WORKFLOW_PARAMETERS '${toJson(paramsCopy)}'
		"""

}

process SC__SCANPY__REPORT_TO_HTML {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notbooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true

	input:
		tuple val(sampleId), path(ipynb)

	output:
		file("*.html")

	script:
		"""
		jupyter nbconvert ${ipynb} --to html
		"""

}

process SC__SCANPY__MERGE_REPORTS {

	container params.sc.scanpy.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notbooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true

	input:
		tuple val(sampleId), path(ipynbs)
		val(reportTitle)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.ipynb")

	script:
		"""
		nbmerge ${ipynbs} -o "${sampleId}.${reportTitle}.ipynb"
		"""

}
