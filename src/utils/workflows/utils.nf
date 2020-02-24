nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// Imports

//////////////////////////////////////////////////////
//  Define the workflow 

workflow COMBINE_BY_PARAMS {

    take:
        // Expects (sampleId, data, unstashedParams)
        A
        // Expects (sampleId, data, [stashedParams])
        B
        params

    main:
        if(params != null && params.isParameterExplorationModeOn()) {
            out = A.concat(
                B.map { it -> tuple(it[0], it[1], *it[2]) } // Unstash params
            ).map {
                it -> tuple(it[2..(it.size()-1)], it[0], it[1]) // Stash params
            }.groupTuple(
                by: [0, params.numParams()-1]
            ).map { 
                it -> tuple(it[1], *it[2], it[0]) 
            }
        } else {
            out = A.join(B)
        }

    emit:
        out

}
