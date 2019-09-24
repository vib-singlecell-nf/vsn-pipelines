//
// Version: 0.1.0
// Test: passed
// Command: 
//  nextflow run src/singlecelltxbenchmark/pipelines/bec__bbknn -profile singularity --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny
//
/*
 * BEC__BBKNN workflow 
 * Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb
 * 
 * Steps considered: 
 * - filter (cell, gene) + qc report
 * - normalize
 * - concatenate the batches
 * - feature selection
 * - log transform
 * - feature scaling
 * - dimensionality reduction (PCA)
 * - batch effect correction using python package bbknn (Park et al. (2018), Fast Batch Alignment of Single Cell Transcriptomes Unifies Multiple Mouse Cell Atlases into an Integrated Landscape)
 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

include groupParams from '../../utils/utils.nf'

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces
PARAMS_GROUPED = groupParams( params )

println(prettyPrint(toJson(PARAMS_GROUPED)))

include SC__FILE_CONVERTER_HELP from '../../processes/utils/utils'

include SC__FILE_CONVERTER from '../../processes/utils/utils' params(PARAMS_GROUPED['SC__FILE_CONVERTER'])
include SC__FILE_ANNOTATOR from '../../processes/utils/utils' params(PARAMS_GROUPED['SC__FILE_ANNOTATOR'])
include '../../processes/scanpy/filter' params(PARAMS_GROUPED['SC__SCANPY_FILTER'])
include SC__FILE_CONCATENATOR from '../../processes/utils/utils' params(PARAMS_GROUPED['SC__FILE_CONCATENATOR'])
include SC__SCANPY__DATA_TRANSFORMATION from '../../processes/scanpy/transform' params(PARAMS_GROUPED['SC__SCANPY__DATA_TRANSFORMATION'])
include SC__SCANPY__NORMALIZATION from '../../processes/scanpy/transform' params(PARAMS_GROUPED['SC__SCANPY__NORMALIZATION'])
include '../../processes/scanpy/feature_selection' params(PARAMS_GROUPED['SC__SCANPY__FEATURE_SELECTION'])
include SC__SCANPY__FEATURE_SCALING from '../../processes/scanpy/transform' params(PARAMS_GROUPED['SC__SCANPY__FEATURE_SCALING'])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../../processes/scanpy/dim_reduction' params(PARAMS_GROUPED['SC__SCANPY__DIM_REDUCTION__PCA'])
include '../../processes/scanpy/batch_effect_correct' params(PARAMS_GROUPED['SC__SCANPY__BATCH_EFFECT_CORRECT'])
include SC__SCANPY__CLUSTERING from '../../processes/scanpy/cluster' params(PARAMS_GROUPED['SC__SCANPY__CLUSTERING'])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../../processes/scanpy/dim_reduction' params(PARAMS_GROUPED['SC__SCANPY__DIM_REDUCTION__UMAP'])
include SC__H5AD_TO_LOOM from '../../processes/utils/h5ad_to_loom'

//////////////////////////////////////////////////////
//  Get the data
include getChannel as getTenXChannel from '../../channels/tenx.nf'

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow bbknn {
    get:
        data
    main:
        // SC__FILE_CONVERTER_HELP()
        // SC__FILE_CONVERTER_HELP.out.subscribe { println it }
        SC__FILE_CONVERTER( data )
        SC__FILE_ANNOTATOR( SC__FILE_CONVERTER.out )
        SC__SCANPY__GENE_FILTER( SC__FILE_ANNOTATOR.out )
        SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        SC__SCANPY__FILTER_QC_REPORT(file(params.template_ipynb), SC__SCANPY__CELL_FILTER.out )
        SC__FILE_CONCATENATOR( SC__SCANPY__CELL_FILTER.out.collect() )
        SC__SCANPY__NORMALIZATION( SC__FILE_CONCATENATOR.out )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )
        SC__SCANPY__FEATURE_SELECTION( SC__SCANPY__DATA_TRANSFORMATION.out )
        SC__SCANPY__FEATURE_SCALING( SC__SCANPY__FEATURE_SELECTION.out )
        SC__SCANPY__DIM_REDUCTION__PCA( SC__SCANPY__FEATURE_SCALING.out )
        SC__SCANPY__BATCH_EFFECT_CORRECTION( SC__SCANPY__DIM_REDUCTION__PCA.out )
        SC__SCANPY__CLUSTERING( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__CLUSTERING.out )
        SC__H5AD_TO_LOOM(SC__SCANPY__DIM_REDUCTION__UMAP.out )
        // Not using t-SNE as it does not use the neighbour graph (which BBKNN alters) when constructing its dimensionality reduction
    emit:
        SC__H5AD_TO_LOOM.out
}

// Uncomment to test
workflow {
    main:
        bbknn( getTenXChannel( params.tenx_folder ) )
}