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

nextflow.preview.dsl=2

include groupParams from '../../utils/utils.nf'

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces


include SC__CELLRANGER__MKFASTQ from '../../processes/cellranger/mkfastq' params(params.sc.cellranger.mkfastq + params.global + params)
include SC__CELLRANGER__COUNT from '../../processes/cellranger/count' params(params.sc.cellranger.mkfastq + params.global)
include SC__FILE_CONVERTER from '../../processes/utils/utils' params(params.sc.file_converter + params.global)
include SC__FILE_ANNOTATOR from '../../processes/utils/utils' params(params.sc.file_annotator + params.global)
include '../../processes/scanpy/filter' params(params.sc.scanpy.filter + params.global)
include SC__FILE_CONCATENATOR from '../../processes/utils/utils' params(params.sc.file_concatenator + params.global)
include SC__SCANPY__DATA_TRANSFORMATION from '../../processes/scanpy/transform' params(params.sc.scanpy.data_transformation + params.global)
include SC__SCANPY__NORMALIZATION from '../../processes/scanpy/transform' params(params.sc.scanpy.normalization + params.global)
include '../../processes/scanpy/feature_selection' params(params.sc.scanpy.feature_selection + params.global)
include SC__SCANPY__FEATURE_SCALING from '../../processes/scanpy/transform' params(params.sc.scanpy.feature_scaling + params.global)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../../processes/scanpy/dim_reduction' params(params.sc.scanpy.dim_reduction.pca + params.global)
include '../../processes/scanpy/batch_effect_correct' params(params.sc.scanpy.batch_effect_correct + params.global)
include SC__SCANPY__CLUSTERING from '../../processes/scanpy/cluster' params(params.sc.scanpy.clustering + params.global)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../../processes/scanpy/dim_reduction' params(params.sc.scanpy.dim_reduction.umap + params.global)
include SC__H5AD_TO_LOOM from '../../processes/utils/h5ad_to_loom' params(params.global)

//////////////////////////////////////////////////////
//  Get the data
include getChannel as getTenXChannel from '../../channels/tenx.nf' params(params.global)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow {
    main:
        SC__CELLRANGER__MKFASTQ(file(params.sc.cellranger.mkfastq.csv), file(params.sc.cellranger.mkfastq.runFolder))
        SC__CELLRANGER__COUNT(file(params.sc.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
        SC__FILE_CONVERTER( SC__CELLRANGER__COUNT.out )
        SC__FILE_ANNOTATOR( SC__FILE_CONVERTER.out )
        SC__SCANPY__GENE_FILTER( SC__FILE_ANNOTATOR.out )
        SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        SC__SCANPY__FILTER_QC_REPORT(file(params.global.template_ipynb), SC__SCANPY__CELL_FILTER.out )
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
// workflow {
//     main:
//         bbknn( getTenXChannel( params.tenx_folder ) )
// }