//
// Version:
// Test:
// Command: 
// 
//
/*
 * Remote run test
 * Source:
 * 
 * Steps considered: 

 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

// print all parameters:
// println(prettyPrint(toJson( params )))


//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include CELLRANGER from '../src/cellranger/main.nf' params(params)
include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params + params.global)
include CLUSTER_IDENTIFICATION from '../src/scanpy/workflows/cluster_identification.nf' params(params + params.global)
include SC__H5AD_TO_LOOM from '../src/utils/processes/h5ad_to_loom.nf' params(params + params.global)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5ad_to_loom.nf' params(params + params.global)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params)

workflow single_sample {
    
    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data )
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out )
    NORMALIZE_TRANSFORM( QC_FILTER.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    DIM_REDUCTION( HVG_SELECTION.out )
    CLUSTER_IDENTIFICATION( DIM_REDUCTION.out )
    scopeloom = SC__H5AD_TO_LOOM( CLUSTER_IDENTIFICATION.out )
    // SC__SCANPY__MERGE_REPORTS(  , "merged_report" )
    emit:
        filteredloom
        scopeloom
    
}

