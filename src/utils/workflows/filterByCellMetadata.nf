//
// Version: 
// Test: 
// Command: 
//
/*
 * QC workflow 
 * Source:
 * 
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include SC__PREPARE_OBS_FILTER from './../processes/h5adSubset' params(params)
include SC__APPLY_OBS_FILTER from './../processes/h5adSubset' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow FILTER_BY_CELL_META_DATA {
    get:
        data
    main:
        Channel
            .from(params.sc.cell_filter.filters)
            .set{ filters }
        SC__PREPARE_OBS_FILTER( 
            data,
            filters
        )
        out = SC__APPLY_OBS_FILTER( 
            data,
            SC__PREPARE_OBS_FILTER.out.map{ it -> it[1] }.collect()
        )
    emit:
        out
}
