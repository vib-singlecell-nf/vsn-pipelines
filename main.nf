nextflow.preview.dsl=2

include {
    getDataChannel;
} from './../channels/channels.nf' params(params)
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
} from "./workflows/decontXFilter"

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
                DECONTX_FILTER.out,
                "CELDA_DECONTX",
                "h5ad",
                null,
                false
            )
        
    emit:
        DECONTX_FILTER.out

}