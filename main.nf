nextflow.preview.dsl=2

// Should be set in case this pipeline is run with other pipelines (e.g.: single_sample)

include {
    getDataChannel;
} from './../channels/channels.nf' addParams(off : "sce_rds")
include {
    SC__FILE_CONVERTER;
} from './../utils/processes/utils' params(params)
include {
    PUBLISH;
} from "./../utils/workflows/utils" params(params)

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    DECONTX_FILTER
} from "./workflows/decontXFilter" params(params)

//////////////////////////////////////////////////////
// Define the workflow

// run decontx
workflow decontx {

    main:
        getDataChannel \
            | SC__FILE_CONVERTER \
            | DECONTX_FILTER

        if(params.utils.containsKey("publish")) {
            PUBLISH(
                DECONTX_FILTER.out.decontx_filtered,
                "CELDA_DECONTX",
                "h5ad",
                null,
                false
            )
        }
        
    emit:
        decontx_filtered = DECONTX_FILTER.out.decontx_filtered
        outlier_table = DECONTX_FILTER.out.outlier_table

}