nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    SC__CELLRANGER__MKFASTQ;
} from './../processes/mkfastq' params(params)
include {
    SC__CELLRANGER__PREFLIGHT;
} from './../processes/preflight' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 

workflow MKFASTQ {

    take:
        mkfastq_csv
        runFolder
    main:
        SC__CELLRANGER__PREFLIGHT()
        SC__CELLRANGER__MKFASTQ(mkfastq_csv, runFolder)
        .flatMap()
        .map { fastq ->
            sample = file(fastq).getParent()
            tuple(
                sample.name,
                file(sample)
            )
        }
        .unique()
        .set { data }
        
    emit:
        data

}
