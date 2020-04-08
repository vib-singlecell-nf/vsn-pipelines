nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include './../../channels/file' params(params)
include './../processes/h5adAnnotate.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow ANNOTATE_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad) [Channel]
        data
        // Expects (sampleId, tsv) [Channel || null]
        metadata
        // Expects name of tool ([string] || null)
        tool

    main:
        def workflowParams = tool == null ? params.sc.cell_annotate : params.sc[tool].cell_annotate
        def method = workflowParams.method
        if(method == 'aio') {
            out = SC__ANNOTATE_BY_CELL_METADATA( 
                data.map {
                    it -> tuple(it[0], it[1], file(workflowParams.cellMetaDataFilePath))
                },
                tool
            )
        } else if(method == 'obo') {
            if(metadata == null) {
                metadata = getChannel(
                    workflowParams.cellMetaDataFilePath,
                    workflowParams.sampleSuffixWithExtension
                )
            }
            out = SC__ANNOTATE_BY_CELL_METADATA(
                data.join(metadata),
                tool
            )
        } else {
            throw new Exception("The given method '" + method + "' is not valid for cell_annotate.")
        }

    emit:
        out

}

workflow STATIC__ANNOTATE_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad)
        data

    main:
        out = ANNOTATE_BY_CELL_METADATA(
            data,
            null,
            null
        )

    emit:
        out

}

