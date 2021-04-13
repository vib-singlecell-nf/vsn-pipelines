nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SEURAT__MERGE {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
label 'compute_resources__default'

    input:
        file('*')

    output:
        tuple val(params.global.project_name), path("${params.global.project_name}.SC__SEURAT__MERGE.${off}")

    script:
        off = params.tools.seurat.merge.off
        """
        ${binDir}/merge/sc_merge.R \
            --output ${params.global.project_name}.SC__SEURAT__MERGE.${off} \
            *
        """

}