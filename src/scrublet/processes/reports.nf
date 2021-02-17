nextflow.enable.dsl=2

process SC__SCRUBLET__DOUBLET_DETECTION_REPORT {

	container params.getToolParams("scrublet").container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	label 'compute_resources__report'

  	input:
		file(ipynb)
		tuple \
			val(sampleId), \
            file(scrubletObjectFile), \
			file(adataWithScrubletInfo), \
			file(adataWithDimRed)

  	output:
    	tuple \
			val(sampleId), \
			file("${sampleId}.SC_Scrublet_doublet_detection_report.ipynb")

  	script:
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.SC_Scrublet_doublet_detection_report.ipynb \
			-p SCRUBLET_OBJECT_FILE $scrubletObjectFile \
            -p H5AD_WITH_SCRUBLET_INFO $adataWithScrubletInfo \
            -p H5AD_WITH_DIM_RED $adataWithDimRed \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}
