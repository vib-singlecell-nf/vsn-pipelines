nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/archr/bin/" : ""

toolParams = params.tools.archr

process SC__ARCHR__CELL_CALLING {

    container toolParams.container
    publishDir "${params.global.outdir}/archr/", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(arrow)

    output:
        tuple val(sampleId),
              path("${sampleId}-TSSEnrichment_vs_nFrags.pdf"),
              path("${sampleId}-qc_stats.txt")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.cell_calling)
		processParams = sampleParams.local
        """
        ${binDir}cell_calling.R \
        --arrow_file ${arrow} \
        --output_dir . \
        --filter_frags ${processParams.filter_frags} \
        --filter_tss ${processParams.filter_tss} \
        --seed ${params.global.seed} \
        --threads ${task.cpus} \
        --genome ${toolParams.genome}
        """
}

        //--sample_name ${sampleId} \
        //--min_frags ${processParams.min_frags} \
