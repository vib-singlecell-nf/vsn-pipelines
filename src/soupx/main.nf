nextflow.enable.dsl=2

// Should be set in case this pipeline is run with other pipelines (e.g.: single_sample)

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
    SOUPX_CORRECT;
} from "./workflows/soupX" params(params)

//////////////////////////////////////////////////////
// Define the workflow

// run decontx
workflow soupx {

    main:

        if(params.data.size() != 1 || !params.data?.tenx?.cellranger_mex) {
            throw new Exception("VSN ERROR: SoupX currenlty works for 'cellranger_mex' data format.")
        }

        data = getDataChannel()
        processed = SOUPX_CORRECT( data )


        if(params.utils?.publish) {
            PUBLISH(
                processed,
                "SOUPX_CORRECT",
                "h5ad",
                null,
                false
            )
        }
        
    emit:
        soupx_processed = processed

}
