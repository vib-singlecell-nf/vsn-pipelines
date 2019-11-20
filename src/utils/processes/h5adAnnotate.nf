nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

process SC__ANNOTATE_BY_CELL_META_DATA {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(id), file(f)
    output:
        tuple val(id), file("${id}.SC__ANNOTATE_BY_CELL_META_DATA.h5ad")
    script:
        annotationColumnNamesAsArguments = params.sc.cell_annotate.annotationColumnNames.collect({ '--annotation-column-name' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_annotate_by_cell_meta_data.py \
            --index-column-name ${params.sc.cell_annotate.indexColumnName} \
            ${annotationColumnNamesAsArguments} \
            $f \
            ${params.sc.cell_annotate.cellMetaDataFilePath} \
            --output "${id}.SC__ANNOTATE_BY_CELL_META_DATA.h5ad"
        """
}
