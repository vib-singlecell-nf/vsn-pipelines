nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    SC__FILE_CONVERTER_FROM_SCE;
} from '../../utils/processes/utils.nf' params(params)
include {
    getBaseName;
} from '../../utils/processes/files.nf' params(params)
include {
    FILTER_BY_CELL_METADATA;
} from '../../utils/workflows/filterByCellMetadata' params(params)
include {
    PUBLISH as PUBLISH_H5AD_DECONTX_FILTERED;
 } from '../../utils/workflows/utils.nf' params(params)
////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__CELDA__DECONTX;
} from "../processes/runDecontX" params(params)
include { 
    SC__CELDA__DECONTX_MERGE_OUTLIER_TABLES;
} from "../processes/utils" params(params)


workflow DECONTX {

    take:
        // Expects (sampleId, sceRds)
        data

    main:
        SC__CELDA__DECONTX(
            data
        )

        // Concat all outlier tables into 1
        SC__CELDA__DECONTX_MERGE_OUTLIER_TABLES(
            SC__CELDA__DECONTX.out.outlier_table.map {
                it -> it[1]
            }.toSortedList( 
                { a, b -> getBaseName(a, "CELDA__DECONTX") <=> getBaseName(b, "CELDA__DECONTX") }
            )
        )
    emit:
        main = SC__CELDA__DECONTX.out.main
        outlier_table = SC__CELDA__DECONTX.out.outlier_table

}

workflow DECONTX_FILTER {

    take:
        // Expects (sampleId, sceRds)
        data

    main:
        decontx = DECONTX( data )

        SC__FILE_CONVERTER_FROM_SCE(
            // Add no group
            decontx.main.map { it -> tuple(it[0], it[1], "NULL") },
            "h5ad",
            "counts"
        )

        FILTER_BY_CELL_METADATA(
            SC__FILE_CONVERTER_FROM_SCE.out,
            "celda.decontx"
        )

    emit:
        decontx_filtered = FILTER_BY_CELL_METADATA.out
        outlier_table = decontx.outlier_table

}

workflow DECONTX_CORRECT {

    take:
        // Expects (sampleId, sceRds)
        data

    main:
        decontx = DECONTX( data )

        SC__FILE_CONVERTER_FROM_SCE(
            decontx.main.map { it -> tuple(it[0], it[1], "NULL") },
            "h5ad",
            "decontXcounts"
        )

    emit:
        decontx_corrected = SC__FILE_CONVERTER_FROM_SCE.out
        outlier_table = decontx.outlier_table

}