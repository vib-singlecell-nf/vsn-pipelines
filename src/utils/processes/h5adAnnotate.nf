nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")

include {
    isParamNull;
    getToolParams;
} from './utils.nf' params(params)

def getPublishDir = { outDir, toolName ->
    if(isParamNull(toolName))
        return "${outDir}/data/intermediate"
    return "${outDir}/data/${toolName.toLowerCase()}"
}

def getMode = { toolName ->
    if(isParamNull(toolName))
        return 'symlink'
    return 'link'
}


process SC__ANNOTATE_BY_CELL_METADATA {

    container params.tools.scanpy.container
    publishDir "${getPublishDir(params.global.outdir,tool)}", mode: "${getMode(tool)}", overwrite: true
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            path(metadata)
        // Expects tool name [string || null]
        val(tool)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.${toolTag}SC__ANNOTATE_BY_CELL_METADATA.h5ad")

    script:
        def sampleParams = params.parseConfig(
            sampleId,
            params.global,
            isParamNull(tool) ? params.utils.cell_annotate : getToolParams(params.tools, tool)["cell_annotate"]
        )
		processParams = sampleParams.local
        toolTag = isParamNull(tool) ? '' : tool.toUpperCase() + '.'
        annotationColumnNamesAsArguments = processParams.containsKey("annotationColumnNames") ?
            processParams.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ')
            : ''
        """
        ${binDir}/sc_h5ad_annotate_by_cell_metadata.py \
            ${processParams.containsKey('method') ? '--method ' + processParams.method : ''} \
            --index-column-name ${processParams.indexColumnName} \
            --sample-id ${sampleId} \
            ${processParams.containsKey('sampleColumnName') ? '--sample-column-name ' + processParams.sampleColumnName : ''} \
            ${annotationColumnNamesAsArguments} \
            $f \
            ${metadata} \
            --output "${sampleId}.${toolTag}SC__ANNOTATE_BY_CELL_METADATA.h5ad"
        """

}

def getMetadataFilePath(processParams) {
    metadataFilePathAsArgument = ''
    if(processParams.containsKey("by")) {
        metadataFilePathAsArgument = processParams.by.containsKey('metadataFilePath') ? processParams.by.metadataFilePath : ''
    } else {
        // make it backward compatible (see sample_annotate_v1.config)
        metadataFilePathAsArgument = processParams.containsKey('metadataFilePath') ? processParams.metadataFilePath : metadataFilePathAsArgument
    }
    return metadataFilePathAsArgument
}

def hasMetadataFilePath(processParams) { 
    return getMetadataFilePath(processParams) != ''
}

process SC__ANNOTATE_BY_SAMPLE_METADATA {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__ANNOTATE_BY_SAMPLE_METADATA.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.sample_annotate)
        processParams = sampleParams.local

        // method / type param
        methodAsArgument = ''
        if(processParams.containsKey("by")) {
            methodAsArgument = processParams.by.containsKey('method') ? processParams.by.method : ''
        } else {
            // make it backward compatible (see sample_annotate_old_v1.config)
            methodAsArgument = processParams.containsKey('type') ? processParams.type : methodAsArgument
        }

        // metadataFilePath param
        metadataFilePathAsArgument = getMetadataFilePath(processParams)

        compIndexColumnNamesFromAdataAsArguments = ''
        compIndexColumnNamesFromMetadataAsArguments = ''
        annotationColumnNamesAsArguments = ''
        if(processParams.containsKey("by")) {
            compIndexColumnNamesFromAdataAsArguments = processParams.by.containsKey('compIndexColumnNames') ?
                processParams.by.compIndexColumnNames.collect { key, value -> return key }.collect({ '--adata-comp-index-column-name ' + ' ' + it }).join(' ') :
                ''
            compIndexColumnNamesFromMetadataAsArguments = processParams.by.containsKey('compIndexColumnNames') ?
                processParams.by.compIndexColumnNames.collect { key, value -> return value }.collect({ '--metadata-comp-index-column-name ' + ' ' + it }).join(' ') :
                ''
            annotationColumnNamesAsArguments = processParams.by.containsKey('annotationColumnNames') ?
                processParams.by.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ') :
                ''
        }

        //  samplecolumnName
        sampleColumnName = ''
        if(processParams.containsKey("by")) {
            if(!processParams.by.containsKey("sampleColumnName")) {
                throw new Exception("VSN ERROR: Missing sampleColumnName param in params.utils.sample_annotate.by.")
            }
            sampleColumnName = processParams.by.sampleColumnName
        } else {
            if(!processParams.containsKey("sampleColumnName")) {
                throw new Exception("VSN ERROR: Missing sampleColumnName param in params.utils.sample_annotate.")
            }
            // make it backward compatible (see sample_annotate_old_v1.config)
            sampleColumnName = processParams.sampleColumnName
        }

        """
        ${binDir}/sc_h5ad_annotate_by_sample_metadata.py \
            --sample-id ${sampleId} \
            ${methodAsArgument != '' ? '--method ' + methodAsArgument : '' } \
            ${metadataFilePathAsArgument != '' ? '--metadata-file-path ' + metadataFilePathAsArgument : '' } \
            ${'--sample-column-name ' + sampleColumnName} \
            ${compIndexColumnNamesFromAdataAsArguments} \
            ${compIndexColumnNamesFromMetadataAsArguments} \
            ${annotationColumnNamesAsArguments} \
            $f \
            "${sampleId}.SC__ANNOTATE_BY_SAMPLE_METADATA.${processParams.off}"
        """

}
