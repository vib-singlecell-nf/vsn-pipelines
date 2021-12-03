nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__UTILS__UPDATE_FEATURE_METADATA_INDEX {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f), path(additionalMetadata)

    output:
        tuple val(sampleId), path("${sampleId}.SC__UTILS__UPDATE_FEATURE_METADATA_INDEX.h5ad")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.update_feature_metadata_index)
		processParams = sampleParams.local
        """
        ${binDir}/sc_h5ad_update_metadata.py \
            --additional-metadata ${additionalMetadata} \
            --axis feature \
            --index-column-name ${processParams.indexColumnName} \
            --join-key ${processParams.joinKey} \
            $f \
            "${sampleId}.SC__UTILS__UPDATE_FEATURE_METADATA_INDEX.h5ad"
        """

}
