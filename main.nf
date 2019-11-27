//
// Version: 0.1.0
// Test: passed
// Command: 
//  nextflow run src/singlecelltxbenchmark/pipelines/bec__bbknn -profile singularity --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny
//
/*
 */ 

nextflow.preview.dsl=2

// include groupParams from '../../utils/utils.nf'

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include SC__CELLRANGER__MKFASTQ from './processes/mkfastq'  params(params)
include SC__CELLRANGER__COUNT   from './processes/count'    params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow CELLRANGER {
    main:
        SC__CELLRANGER__MKFASTQ(file(params.sc.cellranger.mkfastq.csv), path(params.sc.cellranger.mkfastq.runFolder))
        SC__CELLRANGER__COUNT(file(params.sc.cellranger.count.transcriptome), SC__CELLRANGER__MKFASTQ.out.flatten())
    emit:
        SC__CELLRANGER__COUNT.out
}

