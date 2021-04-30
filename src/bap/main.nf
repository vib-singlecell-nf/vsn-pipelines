nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include { SC__BAP__BARCODE_MULTIPLET_PIPELINE; } from './processes/barcode_multiplet.nf' params(params)
include {
    BAM_TO_FRAGMENTS as BAP_BAM_TO_FRAGMENTS;
} from './../../src/sinto/main.nf' addParams(tools_sinto_fragments_barcodetag: params.tools.bap.barcode_multiplet.drop_tag)

include {
    PUBLISH as PUBLISH_BAP_FRAGMENTS;
    PUBLISH as PUBLISH_BAP_FRAGMENTS_INDEX;
} from "../../src/utils/workflows/utils.nf" params(params)

include {
    GENERATE_REPORT;
    REPORT_TO_HTML;
} from './processes/report.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow get_bam {

    take:
        bam
        fragments

    emit:
        bam
}


workflow BAP__BARCODE_MULTIPLET_PIPELINE {

    take:
        bam

    main:

        bap = SC__BAP__BARCODE_MULTIPLET_PIPELINE(bam.map { it -> tuple(it[0], it[1], it[2]) })

        GENERATE_REPORT(
            file(workflow.projectDir + params.tools.bap.barcode_multiplet.report_ipynb),
            bap.map { it -> tuple(it[0], it[3]) },
            "BAP_multiplet_report"
        ) |
        REPORT_TO_HTML

        // generate a fragments file:
        fragments = BAP_BAM_TO_FRAGMENTS(bap.map {it -> tuple(it[0], it[1], it[2])})

        // publish fragments output:
        PUBLISH_BAP_FRAGMENTS(fragments, 'bap.sinto.fragments.tsv', 'gz', 'bap/fragments_sinto', false)
        PUBLISH_BAP_FRAGMENTS_INDEX(fragments.map{ it -> tuple(it[0], it[2]) }, 'bap.sinto.fragments.tsv.gz', 'tbi', 'bap/fragments_sinto', false)

    emit:
        fragments

}

