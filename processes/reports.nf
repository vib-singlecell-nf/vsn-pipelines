nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

/* 
 * STATIC VERSION GENERATE REPORT
 * 
 * General reporting function: 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
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

/* 
 * BENCHMARK VERSION OF SCANPY CLUSTERING GENERATE REPORT
 * 
 * General reporting function: 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
 */
process SC__SCANPY__BENCHMARK_CLUSTERING_GENERATE_REPORT {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/notebooks/intermediate/clustering/${method == "NULL" ? "default": method.toLowerCase()}/${resolution == "NULL" ? "res_": resolution}", mode: 'symlink', overwrite: true

	input:
		file ipynb
		tuple \
			val(sampleId), \
			path(adata), \
			val(method), \
			val(resolution)
		val(reportTitle)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.${uuid}.ipynb")

	script:
		def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
		// File output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		uuid = UUID.randomUUID().toString().substring(0,8)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.${uuid}.ipynb \
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
    	tuple val(sampleId), file("${sampleId}.${reportTitle}.${uuid}.ipynb")

  	script:
	  	def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
		uuid = UUID.randomUUID().toString().substring(0,8)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.${uuid}.ipynb \
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
	// copy final "merged_report" to notebooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true

	input:
		tuple val(sampleId), path(ipynbs)
		val(reportTitle)
		val(isBenchmarkMode)

	output:
		tuple val(sampleId), path("${sampleId}.${isBenchmarkMode ? uuid + "." : ""}${reportTitle}.ipynb")

	script:
		uuid = UUID.randomUUID().toString().substring(0,8)
		"""
		nbmerge ${ipynbs} -o "${sampleId}.${isBenchmarkMode ? uuid + "." : ""}${reportTitle}.ipynb"
		"""

}
