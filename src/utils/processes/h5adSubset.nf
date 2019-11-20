nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SC__PREPARE_OBS_FILTER {
    input:
        tuple val(id), file(f)
        val(filterConfig)
    output:
        tuple val(id), file("${id}.SC__PREPARE_OBS_FILTER.${filterConfig.id}.txt")
    script:
        valuesToKeepFromFilterColumnAsArguments = filterConfig.valuesToKeepFromFilterColumn.collect({ '--value-to-keep-from-filter-column' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_prepare_obs_filter.py \
            --sample-id ${id} \
            --sample-column-name ${filterConfig.sampleColumnName} \
            --barcode-column-name ${filterConfig.barcodeColumnName} \
            --filter-column-name ${filterConfig.filterColumnName} \
            ${valuesToKeepFromFilterColumnAsArguments} \
            ${filterConfig.cellMetaDataFilePath} \
            "${id}.SC__PREPARE_OBS_FILTER.${filterConfig.id}.txt"
        """
}

process SC__APPLY_OBS_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(id), file(f)
        file(filters)
    output:
        tuple val(id), file("${id}.SC__APPLY_OBS_FILTER.${params.sc.cell_filter.off}")
    script:
        filtersAsArguments = filters.collect({ '--filter-file-path' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_apply_obs_filter.py \
            $f \
            --output "${id}.SC__APPLY_OBS_FILTER.${params.sc.cell_filter.off}" \
            $filtersAsArguments
        """
}
