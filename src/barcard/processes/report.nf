nextflow.enable.dsl=2

import java.nio.file.Paths
import static groovy.json.JsonOutput.*

toolParams = params.tools.barcard

process GENERATE_REPORT {
    container "vibsinglecellnf/bap:2021-04-27-3b48f4b"
    publishDir "${params.global.outdir}/notebooks/", mode: params.utils.publish.mode
    label 'compute_resources__report'

    input:
        path(ipynb)
        tuple val(sampleId),
              path("${sampleId}.barcard.overlap.tsv")
        //val(reportTitle)

    output:
        tuple val(sampleId),
              path("${sampleId}.barcard_otsu.ipynb"),
              path("${sampleId}.barcard_kneeplot.png"),
              path("${sampleId}.barcard.overlap.otsu_filtered.tsv")

    script:
        //def sampleParams = params.parseConfig(sampleId)
        //processParams = sampleParams
        //barcardParams = toJson(processParams)

        //def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_multiplet)
        //processParams = sampleParams.local
        //barcard_params = toJson(processParams)
        """
        mkdir .cache/
        mkdir .cache/black/
        mkdir .cache/black/21.4b1/

        papermill ${ipynb} \
            ${sampleId}.barcard_otsu.ipynb \
            --report-mode \
            -p SAMPLE ${sampleId} \
            -p BARCARD_OVERLAP_TSV '${sampleId}.barcard.overlap.tsv'
        """
}


process REPORT_TO_HTML {
    container "vibsinglecellnf/bap:2021-04-27-3b48f4b"
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

