nextflow.enable.dsl=2

include {
    SC__H5AD_TO_FILTERED_LOOM
} from './../processes/h5adToLoom.nf' params(params)
include {
    FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    FILE_CONVERTER as FILE_CONVERTER_TO_SCANPY;
} from "./fileConverter"

// Convert to 
// - SCope-ready
// - Scanpy-ready files
workflow FINALIZE {

    take:
        rawFilteredData
        finalProcessedData
        fileOutputSuffix

    main:
        // Conversion
        // Convert h5ad to X (here we choose: loom format)
        filteredloom = SC__H5AD_TO_FILTERED_LOOM( rawFilteredData )
        FILE_CONVERTER_TO_SCOPE(
            finalProcessedData.groupTuple(),
            fileOutputSuffix,
            'mergeToSCopeLoom',
            rawFilteredData
        )
        FILE_CONVERTER_TO_SCANPY(
            finalProcessedData.groupTuple(),
            fileOutputSuffix,
            'mergeToScanpyH5ad',
            rawFilteredData
        )

    emit:
        filteredloom
        scopeloom = FILE_CONVERTER_TO_SCOPE.out
        scanpyh5ad = FILE_CONVERTER_TO_SCANPY.out

}