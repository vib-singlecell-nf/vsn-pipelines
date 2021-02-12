nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    SC__FILE_CONVERTER_FROM_SCE;
} from '../../utils/processes/utils.nf' params(params)
////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SOUPX;
} from "../processes/runSoupX" params(params)


workflow SOUPX_CORRECT {

    take:
        // Expects (sampleId, CellRangerOuts|CellRangerMEX)
        data

    main:
        soupx = SC__SOUPX( data.map { it -> tuple(it[0], it[1]) } )

        SC__FILE_CONVERTER_FROM_SCE(
            soupx.main.map { it -> tuple(it[0], it[1], "NULL") },
            "h5ad",
            "soupXcounts"
        )

    emit:
        soupx_corrected = SC__FILE_CONVERTER_FROM_SCE.out

}