nextflow.enable.dsl=2

import java.nio.file.Paths
import static groovy.json.JsonOutput.*

toolParams = params.tools.bap

process GENERATE_REPORT {

    container toolParams.container
    publishDir "${params.global.outdir}/notebooks/", mode: params.utils.publish.mode
    label 'compute_resources__report'

    input:
        path(ipynb)
        tuple val(sampleId),
              path(bap_final)
        val(reportTitle)

    output:
        tuple val(sampleId),
              path("${sampleId}.${reportTitle}.ipynb")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_multiplet)
        processParams = sampleParams.local
        bap_params = toJson(processParams)
        """
        papermill ${ipynb} \
            --report-mode \
            ${sampleId}.${reportTitle}.ipynb \
            -p SAMPLE ${sampleId} \
            -p WORKFLOW_PARAMETERS '${bap_params}' \
        """
}


process REPORT_TO_HTML {

    container toolParams.container
    publishDir "${params.global.outdir}/notebooks/", mode: params.utils.publish.mode
    label 'compute_resources__report'

    input:
        tuple val(sampleId),
              path(ipynb)

    output:
        file("*.html")

    script:
        """
        jupyter nbconvert ${ipynb} --to html
        """
}

