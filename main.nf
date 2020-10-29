nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include { 
    SC__SINTO__FRAGMENTS;
    SC__SINTO__SORT_FRAGMENTS;
} from './processes/fragments.nf' params(params)
include { SC__BWAMAPTOOLS__INDEX_BED; } from './../../src/bwamaptools/processes/index.nf' params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow BAM_TO_FRAGMENTS {

    take:
        bam

    main:

        fragments = SC__SINTO__FRAGMENTS(bam)
        fragments_sort = SC__SINTO__SORT_FRAGMENTS(fragments)
        index = SC__BWAMAPTOOLS__INDEX_BED(fragments_sort)

        // join bam index into the bam channel:
        fragments_out = fragments_sort.join(index)

    emit:
        fragments_out

}

