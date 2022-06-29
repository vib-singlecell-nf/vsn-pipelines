nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__MARKER_GENES_TO_XLSX {
    
    container params.tools.seurat.container
    publishDir \
        "${params.global.outdir}/data/seurat", \
        mode: "${params.utils.publish?.mode ? params.utils.publish.mode: 'link'}", \
        saveAs: { filename -> "${sampleId}.marker_genes.xlsx" }
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__MARKER_GENES.xlsx")

    script:
        """
        ${binDir}/utils/sc_marker_genes_to_xlsx.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__MARKER_GENES.xlsx
        """
}