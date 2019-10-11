import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

workflow bbknn_scenic {
    include bbknn_scenic as BBKNN_SCENIC from './workflows/bbknn_scenic' params(params)
    BBKNN_SCENIC()

}

workflow single_sample {
    include single_sample as SINGLE_SAMPLE from './workflows/single_sample' params(params)
    SINGLE_SAMPLE()  
}

workflow star {
    include star as STAR from './workflows/star' params(params)
    STAR()
}

workflow single_sample_star {
    include single_sample_star as SINGLE_SAMPLE_STAR from './workflows/single_sample_star' params(params)
    SINGLE_SAMPLE_STAR()
}