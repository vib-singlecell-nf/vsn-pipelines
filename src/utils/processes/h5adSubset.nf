nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SC__PREPARE_OBS_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f), val(filterConfig)

    output:
        tuple val(sampleId), path(f), path("${sampleId}.SC__PREPARE_OBS_FILTER.${filterConfig.id}.txt")

    script:
        valuesToKeepFromFilterColumnAsArguments = filterConfig.valuesToKeepFromFilterColumn.collect({ '--value-to-keep-from-filter-column' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_prepare_obs_filter.py \
            --sample-id ${sampleId} \
            --sample-column-name ${filterConfig.sampleColumnName} \
            --barcode-column-name ${filterConfig.barcodeColumnName} \
            --filter-column-name ${filterConfig.filterColumnName} \
            ${valuesToKeepFromFilterColumnAsArguments} \
            ${filterConfig.cellMetaDataFilePath} \
            "${sampleId}.SC__PREPARE_OBS_FILTER.${filterConfig.id}.txt"
        """

}

process SC__APPLY_OBS_FILTER {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f), path(filters)

    output:
        tuple val(sampleId), path("${sampleId}.SC__APPLY_OBS_FILTER.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.cell_filter)
		processParams = sampleParams.local
        filtersAsArguments = filters.collect({ '--filter-file-path' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_apply_obs_filter.py \
            $f \
            --output "${sampleId}.SC__APPLY_OBS_FILTER.${processParams.off}" \
            $filtersAsArguments
        """

}
