nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

// Utils
include '../src/utils/processes/utils.nf' params(params)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

// Pipeline
include single_sample as SCANPY__SINGLE_SAMPLE from '../src/scanpy/main.nf' params(params)

workflow single_sample {

    take:
        data

    main:
        // run the pipeline
        SC__FILE_CONVERTER( data )
        SCANPY__SINGLE_SAMPLE( SC__FILE_CONVERTER.out )

        // Conversion
        // Convert h5ad to X (here we choose: loom format)
        if(params.sc.scanpy.containsKey("filter")) {
            filteredloom = SC__H5AD_TO_FILTERED_LOOM( SCANPY__SINGLE_SAMPLE.out.filtered_data )
            // In parameter exploration mode, this automatically merge all the results into the resulting loom
            scopeloom = FILE_CONVERTER(
                SCANPY__SINGLE_SAMPLE.out.final_processed_data.groupTuple(),
                'loom',
                SCANPY__SINGLE_SAMPLE.out.filtered_data
            )
        } else {
            filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONVERTER.out )
            scopeloom = FILE_CONVERTER(
                SCANPY__SINGLE_SAMPLE.out.final_processed_data.groupTuple(),
                'loom',
                SC__FILE_CONVERTER.out
            )
        }

        // Publishing
        finalH5ad = SC__PUBLISH_H5AD(
            SCANPY__SINGLE_SAMPLE.out.final_processed_data.map { 
                it -> tuple(it[0], it[1], null)
            },
            params.global.project_name+".single_sample.output"
        )

    emit:
        finalH5ad
        filteredloom
        scopeloom

}
