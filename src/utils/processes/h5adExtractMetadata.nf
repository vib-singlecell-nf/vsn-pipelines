nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__UTILS__EXTRACT_FEATURE_METADATA {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__UTILS__EXTRACT_FEATURE_METADATA.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.extract_feature_metadata)
		processParams = sampleParams.local
        columnNamesAsArguments = processParams.columnNames.collect({ '--column-name' + ' ' + it }).join(' ')
        """
        ${binDir}/sc_h5ad_extract_metadata.py \
            --axis feature \
            ${columnNamesAsArguments} \
            $f \
            "${sampleId}.SC__UTILS__EXTRACT_FEATURE_METADATA.tsv"
        """

}
