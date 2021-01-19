nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    DSC_PILEUP_FILTERED
} from './dsc_pileup.nf' params(params)
include {
    SC__POPSCLE__FREEMUXLET;
    SC__POPSCLE__DEMUXLET;
} from '../processes/demuxlet.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow cellranger_output_to_bam_barcodes {

    take:
        data // standard data channel [ sampleId, path, type, format]

    main:
        out = data.map{ it -> [it[0],
                               file(it[1]+"/possorted*bam.bam"),
                               file(it[1]+"/filtered_*_bc_matrix/barcodes.tsv*")
                               ] }
                               .view()

    emit:
        out

}

workflow FREEMUXLET {

    take:
        data

    main:
        DSC_PILEUP_FILTERED(data) |
        SC__POPSCLE__FREEMUXLET
    
    emit:
        SC__POPSCLE__FREEMUXLET.out
}

workflow DEMUXLET {

    take:
        data

    main:
        vcf = file(params.sc.popscle.vcf)
        DSC_PILEUP_FILTERED(data)
        SC__POPSCLE__DEMUXLET(DSC_PILEUP_FILTERED.out, vcf)
    
    emit:
        SC__POPSCLE__DEMUXLET.out
}

