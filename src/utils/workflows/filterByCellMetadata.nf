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

workflow FILTER_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad)
        data

    main:
        Channel
            .from(params.sc.cell_filter.filters)
            .set{ filters }
        SC__PREPARE_OBS_FILTER( data.combine(filters) )
        out = SC__APPLY_OBS_FILTER( SC__PREPARE_OBS_FILTER.out.groupTuple() )

    emit:
        out

}
