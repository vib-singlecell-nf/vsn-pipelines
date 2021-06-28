nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SCTK__SATURATION {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(fragments),
              path(fragments_index)
        file(bc_whitelists)
        val(optional)

    output:
        tuple val(sampleId),
              path("${sampleId}.sampling_stats.tsv"),
              path("${sampleId}.saturation.png")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        def bc_wl_param = optional == 'RUN' ? '-w selected_barcodes/' + sampleId + '.cell_barcodes.txt' : ''
        def polars_max_threads = (task.cpus > 6) ? 6 : task.cpus
        """
        # Max threads polars is allowed to use (else will uses all cores).
        export POLARS_MAX_THREADS=${polars_max_threads};
        # Max threads pyarrow is allowed to use (else will uses all cores) (used to read the fragments file in the beginning).
        export OMP_NUM_THREADS=${polars_max_threads};
        calculate_saturation_from_fragments.py \
            -i ${fragments} \
            -o ${sampleId} \
            -p ${toolParams.saturation.percentages} \
            -m ${toolParams.saturation.min_frags_per_cb} \
            -s ${toolParams.saturation.sampling_fractions} \
            ${bc_wl_param}
        """
}

