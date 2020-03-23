nextflow.preview.dsl=2

// include groupParams from '../../utils/utils.nf'

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include SC__CELLRANGER__MKFASTQ             from './processes/mkfastq'  params(params)
include SC__CELLRANGER__COUNT               from './processes/count'    params(params)
include CELLRANGER_COUNT_WITH_METADATA      from './workflows/cellRangerCountWithMetadata'    params(params)


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


workflow CELLRANGER {

    take:
        mkfastq_csv
        runFolder
        transcriptome
    main:
        data = MKFASTQ(mkfastq_csv, runFolder)

        // Allow to combine old demultiplexed data with new data
        if (params.sc.cellranger.count.fastqs instanceof Map) {
            // Remove default key
            Channel.from(params.sc.cellranger.count.fastqs.findAll {
                it.key != 'default' 
            }.collect { k, v -> 
                // Split possible multiple file paths
                if(v.contains(',')) {
                    v = Arrays.asList(v.split(',')).collect { fqs -> file(fqs) }
                } else {
                    v = file(v)
                }
                tuple(k, v) 
            })
            .concat(data)
            .groupTuple()
            .set { data }               
        }
        SC__CELLRANGER__COUNT(transcriptome, data)
    emit:
        SC__CELLRANGER__COUNT.out

}

workflow CELLRANGER_WITH_METADATA {

    main:
        CELLRANGER_COUNT_WITH_METADATA(file("metadata_test.tsv"))

    emit:
        CELLRANGER_COUNT_WITH_METADATA.out
}

