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

    take:
        // Expects (sampleId, data[])
        data
        // Expects (sampleId, outputFormat)
        outputFormat
        // Expects (sampleId, rawFilteredData)
        rawFilteredData

    main:
        data.view()
            .branch {
                h5adToLoom: it[1][0].extension.toLowerCase() == 'h5ad' && outputFormat.toLowerCase() == 'loom'
                none: !inputFormatsAllowed.contains(it[1][0].extension.toLowerCase()) || !outputFormatsAllowed.contains(outputFormat.toLowerCase())
            }
            .set { convert }

        convert.h5adToLoom.view {
            if(it[1].size() > 1) {
                "Aggregating multiple .h5ad files to ${it[1][0].baseName}.loom (w/ additional compression)..."
            } else {
                "Converting ${it[1].baseName}.h5ad to ${it[1].baseName}.loom (w/ additional compression)..."
            }
        }
        SC__H5AD_TO_LOOM(
            rawFilteredData.combine(convert.h5adToLoom, by: 0).ifEmpty('Channel empty: no h5ad files were converted to the loom format.').view()
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
