nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")

process SC__SEURAT_RDS_TO_LOOM {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true, saveAs: { filename -> "${sampleId}.SCope_output.loom" }
    label 'compute_resources__mem'

    input:
        // Expects:
        // - f to be a seurat rds file containing all data that should be made available in loom/SCope
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT_RDS_TO_LOOM.loom")

    script:
        """
        ${binDir}/seurat_rds_to_loom.R \
            --input $f \
            --output "${sampleId}.SC__SEURAT_RDS_TO_LOOM.loom"
		"""
}
