nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// utils:
include SC__FILE_CONVERTER from '../../utils/processes/utils.nf' params(params.sc.file_converter + params.global + params)

// scanpy:
include SC__SCATER__CELL_FILTER from '../processes/filter.nf' params(params.sc.scater.filter + params.global + params)

workflow QC_FILTER {
    get:
        data
    main:
        data = SC__FILE_CONVERTER(data)
        SC__SCATER__CELL_FILTER(data)
    emit:
        SC__SCATER__CELL_FILTER.out
}