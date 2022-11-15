nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    CREATE_FRAGMENTS_FROM_BAM as BARCARD__CREATE_FRAGMENTS_FROM_BAM
} from './processes/create_fragments_from_bam.nf' params(params)

include {
    DETECT_BARCODE_MULTIPLETS as BARCARD__DETECT_BARCODE_MULTIPLETS;
} from './processes/detect_barcode_multiplets.nf' params(params)

include {
    MERGE_BARCODE_MULTIPLETS as BARCARD__MERGE_BARCODE_MULTIPLETS;
} from './processes/merge_barcode_multiplets.nf' params(params)

include {
    GENERATE_REPORT;
    REPORT_TO_HTML;
} from './processes/report.nf' params(params)

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    BWAMAPTOOLS__INDEX_BED;
} from './../../src/bwamaptools/processes/index.nf' params(params)
include {
    PUBLISH as PUBLISH_FRAGMENTS;
    PUBLISH as PUBLISH_FRAGMENTS_INDEX;
} from "../utils/workflows/utils.nf" params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow BAM_TO_FRAGMENTS {

    take:
        bam

    main:

        // sampleID, frag, frag idx
        fragments = BARCARD__CREATE_FRAGMENTS_FROM_BAM(bam)

        //fragments_sort = SINTO__SORT_FRAGMENTS(fragments)
        //index = BWAMAPTOOLS__INDEX_BED(fragments_sort)

        // join bed index into the fragments channel:
        //fragments_out = fragments_sort.join(index)

    emit:
        fragments
        //fragments_out

}


workflow DETECT_BARCODE_MULTIPLETS {

    take:
        fragments

    main:

//        barcard_multiplets = BARCARD__DETECT_BARCODE_MULTIPLETS(fragments.map { it -> tuple(it[0], it[1][0], it[1][1]) })
        barcard_multiplets = BARCARD__DETECT_BARCODE_MULTIPLETS(fragments.map { it -> tuple(it[0], it[1]) })

        //GENERATE_REPORT(
        //    file(workflow.projectDir + params.tools.barcard.barcode_multiplet.report_ipynb),
        //    barcard_multiplets.map { it -> tuple(it[0], it[3]) },
        //    "BARCARD__multiplet_report"
        //) |
        //REPORT_TO_HTML

        GENERATE_REPORT(
            file(workflow.projectDir + params.tools.barcard.barcode_multiplet.report_ipynb),
            barcard_multiplets.map { it -> tuple(it[0], it[1]) },
            //"BARCARD__otsu_filtering_report"
        ) |
        REPORT_TO_HTML

    emit:
        barcard_multiplets

}
