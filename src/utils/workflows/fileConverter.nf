/*
 * Conversion workflow 
 * Source:
 * 
 */

 import nextflow.util.ArrayBag

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

    take:
        // Expects (sampleId, data[])
        data
        // Expects (sampleId, outputFormat)
        outputFormat
        // Expects (sampleId, rawFilteredData)
        rawFilteredData

    main:
        data.branch {
            h5adToLoom: (item = it[1] instanceof ArrayBag ? it[1][0] : it[1]) && item.extension.toLowerCase() == 'h5ad' && outputFormat.toLowerCase() == 'loom'
            none: (item = it[1] instanceof ArrayBag ? it[1][0] : it[1]) && !inputFormatsAllowed.contains(item.extension.toLowerCase()) || !outputFormatsAllowed.contains(outputFormat.toLowerCase())
        }
        .set { convert }

        convert.h5adToLoom.view {
            if(it[1].size() > 1) {
                """
------------------------------------------------------------------
\u001B[32m Aggregating multiple .h5ad files to ${it[1][0].baseName}.loom 
(w/ additional compression)...\u001B[0m
------------------------------------------------------------------
                """
            } else {
"""
------------------------------------------------------------------
\u001B[32m Converting ${it[1][0].baseName}.h5ad to ${it[1][0].baseName}.loom
(w/ additional compression)...\u001B[0m
------------------------------------------------------------------
"""
            }
        }
        SC__H5AD_TO_LOOM(
            rawFilteredData.combine(convert.h5adToLoom.map { it -> tuple(it[0], it[1]) }, by: 0).ifEmpty('Channel empty: no h5ad files were converted to the loom format.')
        )
        out = COMPRESS_HDF5(
            SC__H5AD_TO_LOOM.out
        )
        convert.none.view { 
"""
------------------------------------------------------------------
\u001B[31m Aborting conversion of ${it[1]} to ${it[1].baseName}.loom 
(not implemented) \u001B[0m
------------------------------------------------------------------
"""
        }

    emit:
        out

}
