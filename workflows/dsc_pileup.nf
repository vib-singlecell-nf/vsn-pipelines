nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include SC__FILE_CONVERTER from '../../utils/processes/utils.nf' params(params)

include SC__POPSCLE__DSC_PILEUP from '../processes/dsc_pileup.nf' params(params)
include SC__POPSCLE__PREFILTER_DSC_PILEUP from '../processes/dsc_pileup.nf' params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow DSC_PILEUP_FILTERED {

    take:
        data

    main:
        SC__POPSCLE__PREFILTER_DSC_PILEUP(data) |
        SC__POPSCLE__DSC_PILEUP

}

