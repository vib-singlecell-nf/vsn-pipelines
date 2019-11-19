import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

include QC_FILTER from '../src/scater/workflows/qc_filter.nf' params(params)

workflow scater  {
    
    data = getTenXChannel( params.global.tenx_folder )
    QC_FILTER( data )

    emit:
    QC_FILTER.out
}
