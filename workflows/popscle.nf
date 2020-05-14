nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    FREEMUXLET;
    DEMUXLET;
} from '../src/popscle/workflows/demuxlet.nf' params(params)

workflow freemuxlet {

    take:
        data

    main:
        // run the pipeline
        data = data.map {
                it -> tuple(it[0], it[1])
            }
        out = FREEMUXLET( data )

    // emit:

}

workflow demuxlet {

    take:
        data

    main:
        // run the pipeline
        data = data.map {
                it -> tuple(it[0], it[1])
            }
        out = DEMUXLET( data )

    // emit:

}
