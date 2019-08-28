nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes
include '../processes/scanpy/filter' params(
    cellFilter: {
        minNGenes: 200,
        maxNGenes: 4000,
        maxPercentMito: 0.15
    }, 
    geneFilter: {
        minNCells: 3
    }
)

// Transform
include '../processes/scanpy/transform' params(
    normalization: {
        method: 'cpx',
        numberReadsPerCellFactor: 10000
    }, 
    dataTransformation: {
        method: 'log1p'
    },
    featureScaling: {
        method: 'zscore_scale',
        maxSD: 10
    }
)

// Regress out
include '../processes/scanpy/adjust' params(
    normalization: {
        method: 'linear_regression',
        variablesToRegressOut: ['n_counts','percent_mito']
    }
)

// Feature selection
include '../processes/scanpy/feature_selection' params(
    featureSelection: {
        method: 'mean_disp_plot',
        minMean: 0.0125,
        maxMean: 3,
        minDisp: 0.5
    }
)

// Convert 
include SC_FILE_CONVERTER from 'utils' params(iff: '10x_mtx', off: 'h5ad')

///////////////////////////////////////////
//  Put the processes together
// Convert 
SC_FILE_CONVERTER( file(params.tenx_folder) )

// Filter
SC__SCANPY__GENE_FILTER( tenx_data_dirs )
SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
SC__SCANPY__FILTER_QC_REPORT( SC__SCANPY__CELL_FILTER.out )

// Transform
SC__SCANPY__NORMALIZATION( SC__SCANPY__GENE_FILTER.out )
SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )

// Feature selection
SC__SCANPY__FEATURE_SELECTION( transform.SC__SCANPY__DATA_TRANSFORMATION.out )

// Regress out
SC__SCANPY__ADJUSTMENT( SC__SCANPY__FEATURE_SELECTION.out )

// Scale
SC__SCANPY__FEATURE_SCALING( SC__SCANPY__ADJUSTMENT.out )
