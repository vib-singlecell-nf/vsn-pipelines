//
// Test: passed
// Command: 
//  nextflow src/pipelines/bec__bbknn.nf --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny -resume
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

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces
include SC__FILE_CONVERTER from '../processes/utils/utils' params(
    iff: '10x_mtx', 
    off: 'h5ad', 
    project_name: params.project_name
)

include '../processes/scanpy/filter' params(
    cellFilterMinNGenes: 20, // 200
    cellFilterMaxNGenes: -1, // 
    cellFilterMaxPercentMito: 0, //
    geneFilterMinNCells: 1, // 3
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__FILE_CONCATENATOR from '../processes/utils/utils' params(
    project_name: params.project_name,
    join: 'outer',
    off: 'h5ad'
)

include SC__SCANPY__DATA_TRANSFORMATION from '../processes/scanpy/transform' params(
    dataTransformationMethod: 'log1p',
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__NORMALIZATION from '../processes/scanpy/transform' params(
    normalizationMethod: 'cpx',
    normalizationNumberReadsPerCellFactor: 10000,
    iff: '10x_mtx',
    off: 'h5ad'
)

include '../processes/scanpy/feature_selection' params(
    featureSelectionMethod: 'mean_disp_plot',
    featureSelectionMinMean: 0.0125, // 0.125
    featureSelectionMaxMean: 3, // 2.5
    featureSelectionMinDisp: 0.5, //0.7
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__FEATURE_SCALING from '../processes/scanpy/transform' params(
    featureScalingMthod: 'zscore_scale',
    featureScalingMaxSD: 10,
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/scanpy/dim_reduction' params(
    dimReductionMethod: 'PCA',
    nComps: 25,
    iff: '10x_mtx',
    off: 'h5ad'
)

include '../processes/scanpy/batch_effect_correct' params(
    batchEffectCorrectionMethod: 'bbknn',
    neighbors_within_batch: 5,
    nPcs: 10,
    trim: 0,
    project_name: params.project_name,
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/scanpy/dim_reduction' params(
    dimReductionMethod: 'UMAP',
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/scanpy/dim_reduction' params(
    dimReductionMethod: 't-SNE',
    nJobs: 10,
    iff: '10x_mtx',
    off: 'h5ad'
)

//////////////////////////////////////////////////////
//  Get the data
include getChannel as getTenXChannel from '../channels/tenx.nf'

// def extractSample(path) {
//     (full, parentDir, id, genome) = (path =~ /(.+)\/(.+)\/outs\/filtered_gene_bc_matrices\/(.*)$/)[0]
//     return id
// }
// Channel
//     .from(
//         ['[tenx-path-1]',
//         '[tenx-path-2'])
//     .map { path -> tuple(extractSample( "${path}" ), file("${path}")) }
//     .set { dataChannel }
// dataChannel.view()

/*
 * Make the workflow 
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow test__from_raw__aggr_batch( f ) {
    main:
        SC__FILE_CONVERTER( f )
        SC__SCANPY__GENE_FILTER( SC__FILE_CONVERTER.out )
        SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        SC__SCANPY__FILTER_QC_REPORT( SC__SCANPY__CELL_FILTER.out )
        SC__FILE_CONCATENATOR( SC__SCANPY__CELL_FILTER.out.collect() )
        SC__SCANPY__NORMALIZATION( SC__FILE_CONCATENATOR.out )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )
        SC__SCANPY__FEATURE_SELECTION( SC__SCANPY__DATA_TRANSFORMATION.out )
        SC__SCANPY__FEATURE_SCALING( SC__SCANPY__FEATURE_SELECTION.out )
        SC__SCANPY__DIM_REDUCTION__PCA( SC__SCANPY__FEATURE_SCALING.out )
        SC__SCANPY__BATCH_EFFECT_CORRECTION( SC__SCANPY__DIM_REDUCTION__PCA.out )
    emit:
        SC__SCANPY__BATCH_EFFECT_CORRECTION.out
}

// Uncomment to test
// workflow {
//     main:
//         test__from_raw__aggr_batch( getTenXChannel( params.tenx_folder ) )   
// }
// workflow {
//     main:
//         test__from_raw__aggr_batch( dataChannel )
// }