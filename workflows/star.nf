nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    SC__STAR__LOAD_GENOME;
} from '../src/star/processes/load_genome'  params(params)
include {
    SC__STAR__MAP_COUNT;
} from '../src/star/processes/map_count'  params(params)
include {
    SC__STAR__UNLOAD_GENOME;
} from '../src/star/processes/unload_genome'  params(params)
include {
    SC__STAR_CONCATENATOR;
} from '../src/utils/processes/utils.nf' params(params)

include {
    getChannel;
} as getSingleEndChannel from '../src/channels/singleend.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow star {

    main:
        SC__STAR__LOAD_GENOME( file(params.getToolParams("star").map_count.index) )
        SC__STAR__MAP_COUNT( 
            file(params.getToolParams("star").map_count.index),
            SC__STAR__LOAD_GENOME.out,
            getSingleEndChannel(params.getToolParams("star").map_count.fastqs)
        )
        SC__STAR__UNLOAD_GENOME(
            file(params.getToolParams("star").map_count.index),
            SC__STAR__MAP_COUNT.out.isDone.collect()
        )
        SC__STAR_CONCATENATOR( SC__STAR__MAP_COUNT.out.counts.map { it[1] }.collect() )

    emit:
        SC__STAR_CONCATENATOR.out

}
