nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    SC__STAR__LOAD_GENOME;
} from './processes/load_genome' params(params)
include {
    SC__STAR__MAP_COUNT;
} from './processes/map_count' params(params)
include {
    SC__STAR__UNLOAD_GENOME;
} from './processes/unload_genome' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * Run the workflow for each 10xGenomics CellRanger output folders specified.
 */ 
workflow star {

    main:
        SC__STAR__LOAD_GENOME( file(params.getToolParams("star").map_count.transcriptome) )
        SC__STAR__MAP_COUNT( file(params.getToolParams("star").map_count.transcriptome), SC__STAR__LOAD_GENOME.out, path(params.getToolParams("star").map_count.fastqs) )
        SC__STAR__UNLOAD_GENOME( file(params.getToolParams("star").map_count.transcriptome), SC__STAR__MAP_COUNT.out[0] )

    emit:
        SC__STAR__MAP_COUNT.out

}
