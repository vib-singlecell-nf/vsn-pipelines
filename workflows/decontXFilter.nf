nextflow.preview.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    SC__FILE_CONVERTER_FROM_SCE;
 } from '../../utils/processes/utils.nf' params(params)
include {
    FILTER_BY_CELL_METADATA
} from '../../utils/workflows/filterByCellMetadata' params(params)
include {
    PUBLISH as PUBLISH_H5AD_DECONTX_FILTERED;
 } from '../../utils/workflows/utils.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__CELDA__DECONTX
} from "../processes/runDecontX" params(params)


workflow DECONTX_FILTER {

    take:
        // Expects (sampleId, sceRds)
        data

    main:
        SC__CELDA__DECONTX(
            data
        )

        SC__FILE_CONVERTER_FROM_SCE(
            SC__CELDA__DECONTX.out.main,
            "h5ad",
            "counts"
        )

        FILTER_BY_CELL_METADATA(
            SC__FILE_CONVERTER_FROM_SCE.out,
            "celda.decontx"
        )

    emit:
        decontx_filtered = FILTER_BY_CELL_METADATA.out
        outlier_table = SC__CELDA__DECONTX.out.outlier_table

}