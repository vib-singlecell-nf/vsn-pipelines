nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include DSC_PILEUP_FILTERED from './dsc_pileup.nf' params(params)
include SC__POPSCLE__FREEMUXLET from '../processes/demuxlet.nf' params(params)
include SC__POPSCLE__DEMUXLET from '../processes/demuxlet.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

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
        DSC_PILEUP_FILTERED(data) |
        SC__POPSCLE__DEMUXLET(vcf)
    
    emit:
        SC__POPSCLE__DEMUXLET.out
}

