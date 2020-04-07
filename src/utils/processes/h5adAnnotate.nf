nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin/" : ""

process SC__ANNOTATE_BY_CELL_METADATA {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple \
            val(sampleId), \
            path(f), \
            path(metadata)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__ANNOTATE_BY_CELL_METADATA.h5ad")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.cell_annotate)
		processParams = sampleParams.local
        annotationColumnNamesAsArguments = processParams.containsKey("annotationColumnNames") ?
            processParams.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ')
            : ''
        """
        ${binDir}sc_h5ad_annotate_by_cell_metadata.py \
            ${processParams.containsKey('method') ? '--method ' + processParams.method : ''} \
            --index-column-name ${processParams.indexColumnName} \
            ${annotationColumnNamesAsArguments} \
            $f \
            ${metadata} \
            --output "${sampleId}.SC__ANNOTATE_BY_CELL_METADATA.h5ad"
        """

}

process SC__ANNOTATE_BY_SAMPLE_METADATA() {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__ANNOTATE_BY_SAMPLE_METADATA.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.sample_annotate)
		processParams = sampleParams.local
        """
        ${binDir}sc_h5ad_annotate_by_sample_metadata.py \
            ${(processParams.containsKey('type')) ? '--type ' + processParams.type : ''} \
            ${(processParams.containsKey('metaDataFilePath')) ? '--meta-data-file-path ' + processParams.metaDataFilePath : ''} \
            $f \
            "${sampleId}.SC__ANNOTATE_BY_SAMPLE_METADATA.${processParams.off}"
        """

}
