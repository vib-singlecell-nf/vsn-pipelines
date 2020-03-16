nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include '../src/utils/processes/files.nf' params(params)
include '../src/utils/processes/utils.nf' params(params)
include '../src/utils/workflows/utils.nf' params(params)

include DSC_PILEUP_FILTERED from '../src/popscle/workflows/dsc_pileup.nf' params(params)

workflow popscle {

    take:
        data

    main:
        // run the pipeline
        data = data.map {
                it -> tuple(it[0], it[1])
            }
        print(data)
        out = DSC_PILEUP_FILTERED( data )

    // emit:

}
