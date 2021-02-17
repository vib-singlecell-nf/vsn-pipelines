nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/archr/bin/" : ""

toolParams = params.getToolParams("archr")

process SC__ARCHR__CREATE_ARROW_UNFILTERED {

    container toolParams.container
    publishDir "${params.global.outdir}/archr/", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(fragments),
              path(fragments_index),
              val(filetype)

    output:
        tuple val(sampleId),
              path("${sampleId}.arrow")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.createArrow)
        processParams = sampleParams.local
        """
        ${binDir}createArrow_unfiltered.R \
        --fragments_file ${fragments} \
        --sample_name ${sampleId} \
        --filter_frags ${processParams.filter_frags} \
        --filter_tss ${processParams.filter_tss} \
        --min_frags ${processParams.min_frags} \
        --seed ${params.global.seed} \
        --threads ${task.cpus} \
        --genome ${toolParams.genome}
        """
}

