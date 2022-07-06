nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Process imports:
include {
    isParamNull;
    getToolParams;
} from './../processes/utils.nf' params(params)
include {
    getChannel;
} from './../../channels/file' params(params)
include {
    SC__ANNOTATE_BY_CELL_METADATA;
} from './../processes/h5adAnnotate.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow ANNOTATE_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad, stashedParams?) : Channel
        data
        // Expects (sampleId, tsv) : (Channel || null)
        metadata
        // Describes: name of tool
        // Expects tool: (string || null)
        // Values
        // - tool != null:
        //   - The given tool is performing itself a cell-based annotation
        //   - params.tools[tool] should exist
        // - tool == null:
        //   - params.utils.cell_annotate should exist
        tool

    main:
        def workflowParams = isParamNull(tool) ?
            params.utils.cell_annotate :
            getToolParams(params.tools, tool)["cell_annotate"]
        def method = workflowParams.method
        if(method == 'aio') {
            out = SC__ANNOTATE_BY_CELL_METADATA( 
                data.map {
                    // Add NULL as stashedParams if needed
                    it -> it.size() == 2 ? tuple(it[0], it[1], 'NULL', file(workflowParams.cellMetaDataFilePath)) : tuple(it[0], it[1], it[2], file(workflowParams.cellMetaDataFilePath))
                },
                isParamNull(tool) ? 'NULL' : tool
            )
        } else if(method == 'obo') {
            if(metadata == null) {
                metadata = getChannel(
                    workflowParams.cellMetaDataFilePath,
                    workflowParams.sampleSuffixWithExtension,
                    'NULL'
                )
            }
            out = SC__ANNOTATE_BY_CELL_METADATA(
                data.map {
                    // Add NULL as stashedParams if needed
                    it -> it.size() == 2 ? tuple(it[0], it[1], 'NULL') : it
                }.combine(metadata, by: 0),
                isParamNull(tool) ? 'NULL' : tool
            )
        } else {
            throw new Exception("The given method '" + method + "' is not valid for cell_annotate.")
        }

    emit:
        out

}

workflow ANNOTATE_BY_CELL_METADATA_BY_PAIR {
    take:
        one
        two
        tool
    main:
        ANNOTATE_BY_CELL_METADATA(
            one.map {
                it -> tuple(it[0], it[1])
            },
            two.map {
                it -> tuple(it[0], it[1])
            },
            tool
        )
    emit:
        ANNOTATE_BY_CELL_METADATA.out
}

workflow STATIC__ANNOTATE_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad)
        data
        // Expects name of tool ([string] || null)
        tool

    main:
        out = ANNOTATE_BY_CELL_METADATA(
            data,
            null,
            tool
        )

    emit:
        out

}

