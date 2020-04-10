nextflow.preview.dsl=2

// Utils
include '../src/utils/processes/utils.nf' params(params)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

// Pipeline
include SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE from '../src/scanpy/workflows/single_sample.nf' params(params)

workflow single_sample {

    take:
        data

    main:
        // run the pipeline
        SC__FILE_CONVERTER( data )
        SCANPY__SINGLE_SAMPLE( SC__FILE_CONVERTER.out )

    emit:
        finalh5ad = SCANPY__SINGLE_SAMPLE.out.final_processed_data
        filteredloom = SCANPY__SINGLE_SAMPLE.out.filtered_loom
        scopeloom = SCANPY__SINGLE_SAMPLE.out.final_processed_scope_loom

}
