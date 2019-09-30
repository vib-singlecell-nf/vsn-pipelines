//
// Tested: passed
// Command: 
//  nextflow src/processes/scanpy/transform_1.test.nf --tenx_folder data/01.count/**/filtered_feature_bc_matrix -resume
//

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces
include '../utils/utils' params(iff: '10x_mtx', off: 'h5ad')

include 'filter' params(
    cellFilterMinNGenes: 10, // 200
    cellFilterMaxNGenes: -1, // 4000
    cellFilterMaxPercentMito: 0, // 0.15
    geneFilterMinNCells: 1, // 3
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__NORMALIZATION from 'transform' params(
    normalizationMethod: 'cpx', // 200
    normalizationNumberReadsPerCellFactor: 10000,
    iff: '10x_mtx',
    off: 'h5ad'
)

include SC__SCANPY__DATA_TRANSFORMATION from 'transform' params(
    dataTransformationMethod: 'log1p',
    iff: '10x_mtx',
    off: 'h5ad'
)

//////////////////////////////////////////////////////
//  Get the data
include getChannel as getTenXChannel from '../../channels/tenx.nf'


/*
 * Make the test workflow 
 * Run the simple workflow for each 10xGenomics CellRanger output folders specified.
 * Steps considered: dimensionality reduction
 */ 
workflow test {
    get:
        data
    main:
        SC__FILE_CONVERTER( data )
        SC__SCANPY__NORMALIZATION( SC__SCANPY__CELL_FILTER.out )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )
    emit:
        SC__SCANPY__DATA_TRANSFORMATION.out
}

// Uncomment to test
// workflow {
//     main:
//         test( getTenXChannel( params.tenx_folder ) )   
// }

/*
 * Make the test workflow 
 * Run the workflow for each 10xGenomics CellRanger output folders specified 
 */ 
workflow test__from_raw {
    get:
        data
    main:
        SC__FILE_CONVERTER( data )
        SC__SCANPY__GENE_FILTER( SC__FILE_CONVERTER.out )
        SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        SC__SCANPY__FILTER_QC_REPORT( SC__SCANPY__CELL_FILTER.out )
        SC__SCANPY__NORMALIZATION( SC__SCANPY__CELL_FILTER.out )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )
    emit:
        SC__SCANPY__DATA_TRANSFORMATION.out
}

// Uncomment to test
// workflow {
//     main:
//         test__from_raw( getTenXChannel( params.tenx_folder ) )   
// }