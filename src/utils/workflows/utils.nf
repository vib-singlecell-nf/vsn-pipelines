nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// Imports

//////////////////////////////////////////////////////
//  Define the workflow 

workflow COMBINE_BY_PARAMS {

    take:
        // Expects (sampleId, data, stashedParams, *params)
        A
        B
        params

    main:
        if(params != null && params.isBenchmarkMode()) {
            out = A.concat(
                B.map { it -> tuple(it[0], it[1], *it[2]) } // // Unstash params
            ).map {
                it -> tuple(it[2..(it.size()-1)], it[0], it[1])
            }.groupTuple(
                by: [0, params.numParams()-1]
            ).map { 
                it -> tuple(it[1], *it[2]) 
            }
        } else {
            out = A.join(B)
        }

    emit:
        out

}
