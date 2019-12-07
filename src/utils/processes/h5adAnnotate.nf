nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SC__ANNOTATE_BY_CELL_METADATA {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__ANNOTATE_BY_CELL_METADATA.h5ad")

    script:
        processParams = params.sc.cell_annotate
        annotationColumnNamesAsArguments = processParams.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_annotate_by_cell_metadata.py \
            --index-column-name ${processParams.indexColumnName} \
            ${annotationColumnNamesAsArguments} \
            $f \
            ${processParams.cellMetaDataFilePath} \
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
        processParams = params.sc.sample_annotate
        """
        ${binDir}sc_h5ad_annotate_by_sample_metadata.py \
            ${(processParams.containsKey('type')) ? '--type ' + processParams.type : ''} \
            ${(processParams.containsKey('metaDataFilePath')) ? '--meta-data-file-path ' + processParams.metaDataFilePath : ''} \
            $f \
            "${sampleId}.SC__ANNOTATE_BY_SAMPLE_METADATA.${processParams.off}"
        """

}
