nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    SC__CELLRANGER__MKFASTQ;
} from './processes/mkfastq' params(params)
include {
    SC__CELLRANGER__PREFLIGHT;
} from './processes/preflight' params(params)
include {
    SC__CELLRANGER__COUNT;
} from './processes/count' params(params)
include {
    CELLRANGER_COUNT_WITH_METADATA;
} from './workflows/cellRangerCountWithMetadata' params(params)
include {
    MKFASTQ;
} from './workflows/mkfastq' params(params)


workflow CELLRANGER {

    take:
        mkfastq_csv
        runFolder
        transcriptome
    main:
        SC__CELLRANGER__PREFLIGHT()
        data = MKFASTQ(mkfastq_csv, runFolder)

        // Allow to combine old demultiplexed data with new data
        if (params.tools.cellranger.count.fastqs instanceof Map) {
            // Remove default key
            Channel.from(params.tools.cellranger.count.fastqs.findAll {
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
            // Group fastqs per sample
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
        SC__CELLRANGER__PREFLIGHT()
        CELLRANGER_COUNT_WITH_METADATA(file("metadata_test.tsv"))

    emit:
        CELLRANGER_COUNT_WITH_METADATA.out
}

