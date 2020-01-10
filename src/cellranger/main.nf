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

workflow MKFASTQ {

    take:
        mkfastq_csv
        runFolder
    main:
        SC__CELLRANGER__MKFASTQ(mkfastq_csv, runFolder)
        SC__CELLRANGER__MKFASTQ.out.view()
        fastqs = SC__CELLRANGER__MKFASTQ.out.map {
            fastqDirPath -> (full, parentDir, sampleId) = (fastqDirPath =~ /(.+)\/(.+)_fastqOut/)[0]
            return tuple(sampleId, fastqDirPath)
        }
    emit:
        fastqs

}


workflow CELLRANGER {

    take:
        mkfastq_csv
        runFolder
        transcriptome
    main:
        fastqs = MKFASTQ(mkfastq_csv, runFolder)
        SC__CELLRANGER__COUNT(transcriptome, fastqs)
    emit:
        SC__CELLRANGER__COUNT.out

}

