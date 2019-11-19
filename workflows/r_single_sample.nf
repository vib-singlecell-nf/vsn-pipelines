nextflow.preview.dsl=2

include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

include SC__FILE_CONVERTER from '../src/utils/processes/utils.nf' params(params.sc.file_converter + params.global + params)
include SC__DROPLET_UTILS__EMPTY_DROPS from '../src/dropletutils/processes/emptyDrops.nf' params(params)
include SC__SCATER__CELL_FILTER from '../src/scater/processes/filter.nf' params(params.sc.scater.filter + params.global + params)

workflow r_single_sample  {

    data = getTenXChannel( params.global.tenx_folder )
    SC__FILE_CONVERTER( data )
    SC__DROPLET_UTILS__EMPTY_DROPS( SC__FILE_CONVERTER.out )
    SC__SCATER__CELL_FILTER( SC__DROPLET_UTILS__EMPTY_DROPS.out )
    
    emit:
        SC__SCATER__CELL_FILTER.out

}
