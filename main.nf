import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

if(!params.global.containsKey('seed')) {
    params.seed = workflow.manifest.version.replaceAll("\\.","").toInteger()

    Channel.from('').view {
            """
------------------------------------------------------------------
\u001B[32m No seed detected in the config \u001B[0m
\u001B[32m To ensure reproducibility the seed has been set to ${params.seed} \u001B[0m
------------------------------------------------------------------
            """
    }
}

// run multi-sample with bbknn, output a scope loom file
workflow bbknn {

    include bbknn_standalone as BBKNN from './workflows/bbknn' params(params)
    BBKNN()

}

// run multi-sample with mnncorrect, output a scope loom file
workflow mnncorrect {

    include mnncorrect as MNNCORRECT from './workflows/mnncorrect' params(params)
    MNNCORRECT()

}

// run multi-sample with bbknn, output a scope loom file
workflow harmony {

    include harmony_standalone as HARMONY from './workflows/harmony' params(params)
    HARMONY()

}

// run multi-sample with bbknn, then scenic from the filtered output:
workflow bbknn_scenic {

    include bbknn_standalone as BBKNN from './workflows/bbknn' params(params)
    include SCENIC_append from './src/scenic/main.nf' params(params)
    BBKNN()
    SCENIC_append( BBKNN.out.filteredloom, BBKNN.out.scopeloom )

}

// run multi-sample with harmony, then scenic from the filtered output:
workflow harmony_scenic {

    include harmony_standalone as HARMONY from './workflows/harmony' params(params)
    include SCENIC_append from './src/scenic/main.nf' params(params)
    HARMONY()
    SCENIC_append( HARMONY.out.filteredloom, HARMONY.out.scopeloom )

}


// run single_sample, output a scope loom file
workflow single_sample {

    include single_sample_standalone as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    SINGLE_SAMPLE()

}


// run single_sample, then scenic from the filtered output:
workflow single_sample_scenic {

    include SCENIC_append from './src/scenic/main.nf' params(params)
    include single_sample_standalone as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    SINGLE_SAMPLE()
    SCENIC_append( SINGLE_SAMPLE.out.filteredloom, SINGLE_SAMPLE.out.scopeloom )

}


// run scenic directly from an existing loom file:
workflow scenic {

    include SCENIC as SCENIC_WF from './src/scenic/main.nf' params(params)
    SCENIC_WF( Channel.of( tuple("foobar", file(params.sc.scenic.filteredLoom))) )

}


// runs mkfastq, then CellRanger count:
workflow cellranger {

    include CELLRANGER from './src/cellranger/main.nf' params(params)
    CELLRANGER(
        file(params.sc.cellranger.mkfastq.csv),
        file(params.sc.cellranger.mkfastq.runFolder),
        file(params.sc.cellranger.count.transcriptome)
    )

}

workflow cellranger_metadata {

    include CELLRANGER_COUNT_WITH_METADATA      from './src/cellranger/workflows/cellRangerCountWithMetadata'    params(params)
    CELLRANGER_COUNT_WITH_METADATA(file(params.sc.cellranger.count.metadata))

}


// runs mkfastq, CellRanger count, then single_sample:
workflow single_sample_cellranger {

    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    cellranger | SINGLE_SAMPLE

}

workflow h5ad_single_sample {

    include getChannel as getH5ADChannel from './src/channels/file' params(params)
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    data = getH5ADChannel( 
        params.data.h5ad.file_paths,
        params.data.h5ad.suffix
    ).view() | SINGLE_SAMPLE

}

workflow tsv_single_sample {

    include getChannel as getTSVChannel from './src/channels/file' params(params)
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    data = getTSVChannel( 
        params.data.tsv.file_paths,
        params.data.tsv.suffix
    ).view() | SINGLE_SAMPLE

}

workflow csv_single_sample {

    include getChannel as getCSVChannel from './src/channels/file' params(params)
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    data = getCSVChannel( 
        params.data.csv.file_paths,
        params.data.csv.suffix
    ).view() | SINGLE_SAMPLE

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
        BBKNN( SC__CELLRANGER__COUNT.out )

    emit:
        filteredLoom = BBKNN.out.filteredloom
        scopeLoom = BBKNN.out.scopeloom

}

workflow sra_cellranger_bbknn_scenic {

    include SCENIC_append from './src/scenic/main.nf' params(params)
    sra_cellranger_bbknn()
    SCENIC_append(
        sra_cellranger_bbknn.out.filteredLoom,
        sra_cellranger_bbknn.out.scopeLoom
    )

}

