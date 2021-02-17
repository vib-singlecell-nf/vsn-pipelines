nextflow.enable.dsl=2

// Utils
include {
    clean;
    SC__FILE_CONVERTER;
} from '../src/utils/processes/utils.nf' params(params)

// Pipeline
include {
    SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
} from '../src/scanpy/workflows/single_sample.nf' params(params)
include {
    SC__SCANPY__CLUSTERING_PARAMS;
} from '../src/scanpy/processes/cluster.nf' params(params)
include {
    SC__DIRECTS__SELECT_DEFAULT_CLUSTERING
} from '../src/directs/processes/selectDefaultClustering.nf'

workflow single_sample {

    take:
        data

    main:
        /*******************************************
        * Run the pipeline
        */
        SC__FILE_CONVERTER( data )
        SCANPY__SINGLE_SAMPLE( SC__FILE_CONVERTER.out )

        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.tools.scanpy.clustering) )

        // Select a default clustering when in parameter exploration mode
        if(params.tools?.directs && clusteringParams.isParameterExplorationModeOn()) {
            scopeloom = SC__DIRECTS__SELECT_DEFAULT_CLUSTERING(
                SCANPY__SINGLE_SAMPLE.out.final_processed_scope_loom
            )
        } else {
            scopeloom = SCANPY__SINGLE_SAMPLE.out.final_processed_scope_loom
        }

    emit:
        filteredloom = SCANPY__SINGLE_SAMPLE.out.filtered_loom
        scanpyh5ad = SCANPY__SINGLE_SAMPLE.out.final_processed_scanpy_h5ad
        scopeloom = scopeloom

}
