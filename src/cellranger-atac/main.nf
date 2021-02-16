nextflow.enable.dsl=2

// include groupParams from '../../utils/utils.nf'

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    SC__CELLRANGER_ATAC__MKFASTQ;
} from './processes/mkfastq' params(params)
include {
    SC__CELLRANGER_ATAC__COUNT;
} from './processes/count' params(params)
include {
    CELLRANGER_ATAC_COUNT_WITH_METADATA;
} from './workflows/cellRangerCountWithMetadata' params(params)


//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 

workflow MKFASTQ_ATAC {

    take:
        mkfastq_csv
        runFolder
    main:
        SC__CELLRANGER_ATAC__MKFASTQ(mkfastq_csv, runFolder)
        if(!params.containsKey('quiet')) SC__CELLRANGER_ATAC__MKFASTQ.out.view()
        fastqs = SC__CELLRANGER_ATAC__MKFASTQ.out.map {
            fastqDirPath -> (full, parentDir, sampleId) = (fastqDirPath =~ /(.+)\/(.+)_fastqOut/)[0]
            return tuple(sampleId, fastqDirPath)
        }
    emit:
        fastqs

}


workflow CELLRANGER_ATAC {

    take:
        mkfastq_csv
        runFolder
        reference
    main:
        fastqs = MKFASTQ_ATAC(mkfastq_csv, runFolder)
        SC__CELLRANGER_ATAC__COUNT(reference, fastqs)
    emit:
        SC__CELLRANGER_ATAC__COUNT.out

}

workflow CELLRANGER_ATAC_WITH_METADATA {

    main:
        CELLRANGER_ATAC_COUNT_WITH_METADATA(file("metadata_test.tsv"))

    emit:
        CELLRANGER_ATAC_COUNT_WITH_METADATA.out
}

