nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*
import org.yaml.snakeyaml.Yaml

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

/* 
 * STATIC VERSION GENERATE REPORT
 * 
 * General reporting function: 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
 */
process SC__SCANPY__GENERATE_REPORT {

  	container params.getToolParams("scanpy").container
  	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

	input:
		file ipynb
		tuple val(sampleId), path(adata)
		val(reportTitle)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.ipynb")

	script:
		def reportParams = new Yaml().dump(annotations_to_plot: params.getToolParams("scanpy").report.annotations_to_plot)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.ipynb \
			-p FILE $adata \
			-y "${reportParams}" \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}

/* 
 * PARAMETER EXPLORATION VERSION OF SCANPY CLUSTERING GENERATE REPORT
 * 
 * General reporting function: 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
 */
process SC__SCANPY__PARAM_EXPLORE_CLUSTERING_GENERATE_REPORT {

  	container params.getToolParams("scanpy").container
  	publishDir "${params.global.outdir}/notebooks/intermediate/clustering/${isParamNull(method) ? "default": method.toLowerCase()}/${isParamNull(resolution) ? "res_": resolution}", mode: 'symlink', overwrite: true
    label 'compute_resources__report'

	input:
		file ipynb
		tuple \
			val(sampleId), \
			path(adata), \
			val(method), \
			val(resolution)
		val(reportTitle)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.${reportTitle}.${uuid}.ipynb"), \
			val(method), \
			val(resolution)

	script:
		// In parameter exploration mode, file output needs to be tagged with a unique identitifer because of:
		// - https://github.com/nextflow-io/nextflow/issues/470
		stashedParams = [method, resolution]
		if(!isParamNull(stashedParams))
			uuid = stashedParams.findAll { it != 'NULL' }.join('_')
		def reportParams = new Yaml().dump(annotations_to_plot: params.getToolParams("scanpy").report.annotations_to_plot)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.${uuid}.ipynb \
			-p FILE $adata \
			-y "${reportParams}" \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}

// QC report takes two inputs, so needs it own process
process SC__SCANPY__GENERATE_DUAL_INPUT_REPORT {

	container params.getToolParams("scanpy").container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

  	input:
		file(ipynb)
		tuple \
			val(sampleId), \
			file(data1), \
			file(data2), \
			val(stashedParams)
		val(reportTitle)
		val(isParameterExplorationModeOn)

  	output:
    	tuple \
			val(sampleId), \
			file("${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + "." : ''}ipynb"), \
			val(stashedParams)

  	script:
		if(!isParamNull(stashedParams))
			uuid = stashedParams.findAll { it != 'NULL' }.join('_')
		def reportParams = new Yaml().dump(annotations_to_plot: params.getToolParams("scanpy").report.annotations_to_plot)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + "." : ''}ipynb \
			-p FILE1 $data1 -p FILE2 $data2 \
			-y "${reportParams}" \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}

process SC__SCANPY__REPORT_TO_HTML {

	container params.getToolParams("scanpy").container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notbooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true
    label 'compute_resources__report'

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

	container params.getToolParams("scanpy").container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notebooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true
    label 'compute_resources__report'

	input:
		tuple \
			val(sampleId), \
			path(ipynbs), \
			val(stashedParams)
		val(reportTitle)
		val(isParameterExplorationModeOn)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + '.' : ''}ipynb")

	script:
		if(!isParamNull(stashedParams))
			uuid = stashedParams.findAll { it != 'NULL' }.join('_')
		"""
		nbmerge \
			${ipynbs} \
			-o "${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + '.' : ''}ipynb"
		"""

}
