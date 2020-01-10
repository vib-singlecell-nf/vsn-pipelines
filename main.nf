import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2


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

// run multi-sample with bbknn, then scenic from the filtered output:
workflow bbknn_scenic {

    include bbknn_standalone as BBKNN from './workflows/bbknn' params(params)
    include SCENIC_append from './src/scenic/main.nf' params(params)
    BBKNN()
    SCENIC_append( BBKNN.out.filteredloom, BBKNN.out.scopeloom )

}


// run single_sample, output a scope loom file
workflow single_sample {

    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    SINGLE_SAMPLE()

}


// run single_sample, then scenic from the filtered output:
workflow single_sample_scenic {

    include SCENIC_append from './src/scenic/main.nf' params(params)
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
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


// runs mkfastq, CellRanger count, then single_sample:
workflow single_sample_cellranger {

    cellranger | single_sample

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

    include getChannel as getSRAChannel from './src/channels/sra' params(params)
    include DOWNLOAD_FROM_SRA from './src/utils/workflows/downloadFromSRA' params(params)
    include SC__CELLRANGER__PREPARE_FOLDER from './src/cellranger/processes/utils.nf'
    include SC__CELLRANGER__COUNT from './src/cellranger/processes/count' params(params)
    include bbknn as BBKNN from './workflows/bbknn' params(params)
    
    // Run 
    DOWNLOAD_FROM_SRA( getSRAChannel( params.data.sra ) ).view()
    SC__CELLRANGER__PREPARE_FOLDER( DOWNLOAD_FROM_SRA.out.groupTuple() ).view()
    SC__CELLRANGER__COUNT(
        file(params.sc.cellranger.count.transcriptome),
        SC__CELLRANGER__PREPARE_FOLDER.out
    )
    BBKNN( SC__CELLRANGER__COUNT.out )

}

