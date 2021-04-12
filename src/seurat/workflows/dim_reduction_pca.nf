nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__DIM_REDUCTION as SC__SEURAT__DIM_REDUCTION__PCA;
} from '../processes/dim_reduction.nf' params(params + [method: "pca"])
include {
    PCACV__FIND_OPTIMAL_NPCS;
} from './../../pcacv/processes/runPCACV' params(params)

workflow DIM_REDUCTION_PCA {
    take:
        data

    main:
        if (params.containsKey("pcacv")) {
            PCACV__FIND_OPTIMAL_NPCS( data )
            out = SC__SEURAT__DIM_REDUCTION__PCA(
                data.join(
                    PCACV__FIND_OPTIMAL_NPCS.out.optimalNumberPC.map {
                        it -> tuple(it[0], it[1])
                    }
                )
            )
        } else {
            out = SC__SEURAT__DIM_REDUCTION__PCA( 
                data.map { it -> tuple(it[0], it[1], null) } 
            )
        }

    emit:
        out
}