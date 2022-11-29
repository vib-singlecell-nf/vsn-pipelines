nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""


process DETECT_BARCODE_MULTIPLETS {
    container params.tools.barcard.container
    label 'compute_resources__barcard__detect_barcode_multiplets'
    publishDir "${params.global.outdir}/data/reports/barcard/", mode: 'copy'

    input:
        tuple val(sampleId),
              path(fragments)

    output:
        tuple val(sampleId),
              //path(fragments),
              path("${sampleId}.barcard.overlap.tsv")

    script:
        //def sampleParams = params.parseConfig(sampleId, params.global)
        //processParams = sampleParams.local
        """
        set -euo pipefail

        chromosome_regex='^(chr)?([0-9]+|[XY])\$'
        calculate_jaccard_index_cbs.py -i ${fragments} -o ${sampleId}.barcard.overlap.tsv -t 1000 -c \${chromosome_regex}
        """
}
