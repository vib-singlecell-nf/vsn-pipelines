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
        tuple val(id), file(f)
    output:
        tuple val(id), file("${id}.SC__ANNOTATE_BY_CELL_METADATA.h5ad")
    script:
        annotationColumnNamesAsArguments = params.sc.cell_annotate.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_annotate_by_cell_metadata.py \
            --index-column-name ${params.sc.cell_annotate.indexColumnName} \
            ${annotationColumnNamesAsArguments} \
            $f \
            ${params.sc.cell_annotate.cellMetaDataFilePath} \
            --output "${id}.SC__ANNOTATE_BY_CELL_METADATA.h5ad"
        """
}

process SC__ANNOTATE_BY_SAMPLE_METADATA() {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(id), file(f)
        file(metaDataFilePath)
    output:
        tuple val(id), file("${id}.SC__ANNOTATE_BY_SAMPLE_METADATA.${params.sc.sample_annotate.off}")
    script:
        """
        ${binDir}sc_h5ad_annotate_by_sample_metadata.py \
            ${(params.containsKey('type')) ? '--type ' + params.sc.sample_annotate.type : ''} \
            ${(params.containsKey('metaDataFilePath')) ? '--meta-data-file-path ' + metaDataFilePath.getName() : ''} \
            $f \
            "${id}.SC__ANNOTATE_BY_SAMPLE_METADATA.${params.sc.sample_annotate.off}"
    """
}
