import nextflow.util.ArrayBag

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__H5AD_TO_LOOM;
} from './../processes/h5adToLoom.nf' params(params)
include {
    SC__H5AD_MERGE
} from "./../processes/h5adMerge.nf" params(params)
include {
    SC__SEURAT_RDS_TO_LOOM
} from "./../processes/seuratRdsToLoom.nf" params(params)
include {
    isParamNull;
    PUBLISH;
} from "./utils.nf" params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

inputFormatsAllowed = ['h5ad']
outputFormatsAllowed = ['loom', 'h5ad']

workflow FILE_CONVERTER {

    take:
        // Expects (sampleId, data[])
        data
        // Expects outputSuffix: string
        outputSuffix
        // Expects outputFormat: string
        outputFormat
        // Expects (sampleId, rawFilteredData)
        rawFilteredData

    main:
        out = Channel.empty()

        if(outputFormat == "mergeToSCopeLoom") {
            if(isParamNull(rawFilteredData)) {
                throw new Exception("VSN ERROR: Expecting rawFilteredData not to be null when outputFormat is "+ outputFormat)
            }
            out = SC__H5AD_TO_LOOM(
                rawFilteredData.combine(
                    data.map {
                        it -> tuple(it[0], it[1]) 
                    }, 
                    by: 0
                )
            )
        } else if(outputFormat == "mergeToScanpyH5ad") {
            out = SC__H5AD_MERGE(
                data.map {
                    it -> tuple(it[0], it[1])
                }
            )
        }  else if(outputFormat == "seuratRdsToSCopeLoom") {
            out = SC__SEURAT_RDS_TO_LOOM(
                data.map {
                    it -> tuple(it[0], it[1])
                }
            )
        } else {
            throw new Exception("VSN ERROR: Output format "+ outputFormat +"not supported")
        }

    emit:
        out

}

