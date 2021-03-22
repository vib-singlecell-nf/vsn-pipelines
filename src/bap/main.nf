nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    BAP__BARCODE_MULTIPLET_PIPELINE as BARCODE_MULTIPLET_PIPELINE;
} from './processes/barcode_multiplet.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow get_bam {

    take:
        bam
        fragments

    emit:
        bam
}

workflow BAP__BARCODE_MULTIPLET_WF {

    take:
        bam

    main:

        bap = BARCODE_MULTIPLET_PIPELINE(bam.map { it -> tuple(it[0], it[1], it[2]) })

    emit:
        bap

}

