import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2


// run multi-sample with bbknn, output a scope loom file
workflow bbknn {
    include bbknn as BBKNN from './workflows/bbknn' params(params)
    BBKNN()
}


// run multi-sample with bbknn, then scenic from the filtered output:
workflow bbknn_scenic {
    include bbknn as BBKNN from './workflows/bbknn' params(params)
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
    SCENIC_WF( file(params.sc.scenic.filteredloom) )
}


workflow star {
    include star as STAR from './workflows/star' params(params)
    STAR()
}

workflow scater {
    include scater as SCATER from './workflows/scater' params(params)
    SCATER()
}


workflow single_sample_star {
    include single_sample_star as SINGLE_SAMPLE_STAR from './workflows/single_sample_star' params(params)
    SINGLE_SAMPLE_STAR()
}

workflow nemesh {
    include nemesh as NEMESH from './workflows/nemesh' params(params)
    NEMESH()  
}

workflow r_single_sample {
    include r_single_sample as R_SINGLE_SAMPLE from './workflows/r_single_sample' params(params)
    R_SINGLE_SAMPLE()
}
