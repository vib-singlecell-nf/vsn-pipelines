/*
 * Conversion workflow 
 * Source:
 * 
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include SC__H5AD_TO_LOOM from './../processes/h5adToLoom.nf' params(params)
include COMPRESS_HDF5 from './../processes/utils.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

inputFormatsAllowed = ['h5ad']
outputFormatsAllowed = ['loom']

workflow FILE_CONVERTER {

    get:
        // Expects (sampleId, rawFilteredData)
        rawFilteredData
        // Expects (sampleId, data)
        data
        // Expects (sampleId, outputFormat)
        outputFormat

    main:
        data.view()
            .branch {
                h5adToLoom: it[1].extension.toLowerCase() == 'h5ad' && outputFormat.toLowerCase() == 'loom'
                none: !inputFormatsAllowed.contains(it[1].extension.toLowerCase()) || !outputFormatsAllowed.contains(outputFormat.toLowerCase())
            }
            .set { convert }

        convert.h5adToLoom.view { 
            "Converting ${it[1].baseName} to ${it[1].baseName}.loom (w/ additional compression)..." 
        }
        SC__H5AD_TO_LOOM(
            rawFilteredData.join(convert.h5adToLoom).ifEmpty('Channel empty: no h5ad files were converted to the loom format.').view()
        )
        out = COMPRESS_HDF5(
            SC__H5AD_TO_LOOM.out
        )
        convert.none.view { 
            "Aborting conversion of ${it[1]} to ${it[1].baseName}.loom (not implemented)" 
        }

    emit:
        out

}
