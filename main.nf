import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

include '../utils/workflows/utils.nf' params(params)
INIT()
include '../utils/processes/utils.nf' params(params)
include './src/channels/channels' params(params)

// run multi-sample with bbknn, output a scope loom file
workflow bbknn {

    include bbknn as BBKNN from './workflows/bbknn' params(params)
    getDataChannel | BBKNN

}

// run multi-sample with mnncorrect, output a scope loom file
workflow mnncorrect {

    include mnncorrect as MNNCORRECT from './workflows/mnncorrect' params(params)
    getDataChannel | MNNCORRECT

}

// run multi-sample with bbknn, output a scope loom file
workflow harmony {

    include harmony as HARMONY from './workflows/harmony' params(params)
    getDataChannel | HARMONY

}

// run multi-sample with bbknn, then scenic from the filtered output:
workflow bbknn_scenic {

    include bbknn as BBKNN from './workflows/bbknn' params(params)
    include scenic_append as SCENIC_APPEND from './src/scenic/main.nf' params(params)
    getDataChannel | BBKNN
    SCENIC_APPEND(
        BBKNN.out.filteredloom,
        BBKNN.out.scopeloom
    )

}

// run multi-sample with harmony, then scenic from the filtered output:
workflow harmony_scenic {

    include harmony as HARMONY from './workflows/harmony' params(params)
    include scenic_append as SCENIC_APPEND from './src/scenic/main.nf' params(params)
    getDataChannel | HARMONY
    SCENIC_APPEND( 
        HARMONY.out.filteredloom,
        HARMONY.out.scopeloom 
    )

}


// run single_sample, output a scope loom file
workflow single_sample {

    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    getDataChannel | SINGLE_SAMPLE

}

workflow multi_sample {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    getDataChannel | MULTI_SAMPLE

}

// run single_sample, then scenic from the filtered output:
workflow single_sample_scenic {

    include scenic_append as SCENIC_APPEND from './src/scenic/main.nf' params(params)
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    getDataChannel | SINGLE_SAMPLE
    SCENIC_APPEND(
        SINGLE_SAMPLE.out.filteredloom,
        SINGLE_SAMPLE.out.scopeloom
    )

}

// run single_sample, then scenic from the previous input (not standalone):
workflow pipe_single_sample_scenic {

    take:
        data
    main:
        include scenic_append as SCENIC_APPEND from './src/scenic/main.nf' params(params)
        include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
        data | SINGLE_SAMPLE
        SCENIC_APPEND(
            SINGLE_SAMPLE.out.filteredloom,
            SINGLE_SAMPLE.out.scopeloom
        )

}


// run scenic directly from an existing loom file:
workflow scenic {

    include scenic as SCENIC from './src/scenic/main.nf' params(params)
    SCENIC( Channel.of( tuple("foobar", file(params.sc.scenic.filteredLoom))) )

}


// runs mkfastq, then CellRanger count:
workflow cellranger {

    include CELLRANGER from './src/cellranger/main.nf' params(params)
    CELLRANGER(
        file(params.sc.cellranger.mkfastq.csv),
        file(params.sc.cellranger.mkfastq.runFolder),
        file(params.sc.cellranger.count.transcriptome)
    )

    emit:
        CELLRANGER.out
}

workflow cellranger_libraries {

    include CELLRANGER_LIBRARIES from './src/cellranger/workflows/cellranger_libraries.nf' params(params)
    CELLRANGER_LIBRARIES(
        file(params.sc.cellranger.mkfastq.csv),
        file(params.sc.cellranger.mkfastq.runFolder),
        file(params.sc.cellranger.count.transcriptome),
        file(params.sc.cellranger.count.featureRef)
    )

    emit:
        CELLRANGER_LIBRARIES.out

}

workflow cellranger_metadata {

    include CELLRANGER_COUNT_WITH_METADATA from './src/cellranger/workflows/cellRangerCountWithMetadata' params(params)
    CELLRANGER_COUNT_WITH_METADATA(
        file(params.sc.cellranger.count.transcriptome),
        file(params.sc.cellranger.count.metadata)
    )
    emit:
        CELLRANGER_COUNT_WITH_METADATA.out

}

workflow cellranger_metadata_single_sample_scenic {

    cellranger_metadata | \
        map {
            it -> tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        } | \
        pipe_single_sample_scenic

}

workflow cellranger_count_libraries {

    include CELLRANGER_COUNT_WITH_LIBRARIES from './src/cellranger/workflows/cellRangerCountWithLibraries' params(params)
    CELLRANGER_COUNT_WITH_LIBRARIES(
        file(params.sc.cellranger.count.transcriptome),
        file(params.sc.cellranger.count.featureRef),
        params.sc.cellranger.count.libraries
    )

    emit:
        CELLRANGER_COUNT_WITH_LIBRARIES.out

}

workflow freemuxlet {
    include freemuxlet as FREEMUXLET from './workflows/popscle' params(params)
    getDataChannel | FREEMUXLET
}

workflow demuxlet {
    include demuxlet as DEMUXLET from './workflows/popscle' params(params)
    getDataChannel | DEMUXLET
}

// runs mkfastq, CellRanger count, then single_sample:
workflow single_sample_cellranger {

    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    data = cellranger()
    SINGLE_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )

}

workflow cellranger_multi_sample {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    data = cellranger()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )

}

workflow cellranger_multi_sample_demuxlet {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    include demuxlet as DEMUXLET from './workflows/popscle' params(params)
    data = cellranger()
    MULTI_SAMPLE(        
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        }
    )
    DEMUXLET(data)

}

workflow cellranger_libraries_multi_sample {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    data = cellranger_libraries()
    MULTI_SAMPLE(        
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        }
    )
}

workflow cellranger_libraries_freemuxlet_multi_sample {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    include freemuxlet as FREEMUXLET from './workflows/popscle' params(params)
    data = cellranger_libraries()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )
    FREEMUXLET(data)

}

workflow cellranger_libraries_demuxlet_multi_sample {

    include multi_sample as MULTI_SAMPLE from './workflows/multi_sample' params(params)
    include demuxlet as DEMUXLET from './workflows/popscle' params(params)
    data = cellranger_libraries()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )
    DEMUXLET(data)
}

workflow star {

    include star as STAR from './workflows/star' params(params)
    STAR()

}


workflow single_sample_star {

    include single_sample_star as SINGLE_SAMPLE_STAR from './workflows/single_sample_star' params(params)
    SINGLE_SAMPLE_STAR()

}

workflow nemesh {

    include nemesh as NEMESH from './workflows/nemesh' params(params)
    NEMESH()

}

workflow sra_cellranger_bbknn {

    main: 
        include getChannel as getSRAChannel from './src/channels/sra' params(params)
        include DOWNLOAD_FROM_SRA from './src/utils/workflows/downloadFromSRA' params(params)
        include SC__CELLRANGER__PREPARE_FOLDER from './src/cellranger/processes/utils.nf' params(params)
        include SC__CELLRANGER__COUNT from './src/cellranger/processes/count' params(params)
        include bbknn as BBKNN from './workflows/bbknn' params(params)
        
        // Run 
        DOWNLOAD_FROM_SRA( getSRAChannel( params.data.sra ) )
        SC__CELLRANGER__PREPARE_FOLDER( DOWNLOAD_FROM_SRA.out.groupTuple() )
        SC__CELLRANGER__COUNT(
            file(params.sc.cellranger.count.transcriptome),
            SC__CELLRANGER__PREPARE_FOLDER.out
        )
        BBKNN( 
            SC__CELLRANGER__COUNT.out.map {
                it -> tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
        )

    emit:
        filteredLoom = BBKNN.out.filteredloom
        scopeLoom = BBKNN.out.scopeloom

}

workflow sra_cellranger_bbknn_scenic {

    include scenic_append as SCENIC_APPEND from './src/scenic/main.nf' params(params)
    sra_cellranger_bbknn()
    SCENIC_APPEND(
        sra_cellranger_bbknn.out.filteredLoom,
        sra_cellranger_bbknn.out.scopeLoom
    )

}

