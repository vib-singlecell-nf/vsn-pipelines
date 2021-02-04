nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include { 
    SC__SINTO__FRAGMENTS;
    SC__SINTO__SORT_FRAGMENTS;
} from './processes/fragments.nf' params(params)
include { SC__BWAMAPTOOLS__INDEX_BED; } from './../../src/bwamaptools/processes/index.nf' params(params)
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

        //println("${params.tools.sinto.fragments.barcodetag}")
        fragments = SC__SINTO__FRAGMENTS(bam)
        fragments_sort = SC__SINTO__SORT_FRAGMENTS(fragments)
        index = SC__BWAMAPTOOLS__INDEX_BED(fragments_sort)

        // join bed index into the fragments channel:
        fragments_out = fragments_sort.join(index)

    emit:
        fragments_out

}

