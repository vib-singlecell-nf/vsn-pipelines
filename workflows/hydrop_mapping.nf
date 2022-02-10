nextflow.enable.dsl=2

include {
    FASTP__CLEAN_AND_FASTQC;
} from '../src/fastp/processes/clean_and_fastqc.nf' params(params)

include {
    SC__STAR__SOLO_MAP_COUNT;
} from '../src/star/processes/solo_map_count.nf' params(params)



workflow hydrop_mapping {

    /*
    * Create a channel for input read files
    */
    Channel
        .fromFilePairs( params.tools.star.map_count.fastqs, size: 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.tools.star.map_count.fastqs}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard! E.g. \"/path/to/HYD__b69d3e__S1-HyD*_R{1,2}.fastq.gz\" " }
        // .map { println it[0] }
        .map { 
            if ( it[0].matches(/.*_L00[0-9]$/) ) {
                tuple( it[0].replaceFirst(/_L00[0-9]$/, ""), it[1][0], it[1][1])
            }
        }
        .groupTuple(sort: true)
        .map { it.size() == 3 ? tuple( it[0], it[1] + it[2]) :  tuple( it[0], it[1]) }
        .set { data }
    data.subscribe { println it }

    FASTP__CLEAN_AND_FASTQC( data )
    SC__STAR__SOLO_MAP_COUNT (
        FASTP__CLEAN_AND_FASTQC.out.fastq,
        Channel.fromList(params.tools.star.map_count.soloCBwhitelist).collect(),
    )

    emit:
        solo_outputs = SC__STAR__SOLO_MAP_COUNT.out.solo_outputs
        bam = SC__STAR__SOLO_MAP_COUNT.out.bam
}
