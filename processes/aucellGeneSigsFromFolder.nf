nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_aucell/", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file exprMat
    file regulonsFolder
    val type

    output:
    file "multi_runs_regulons_auc_${type}.tsv"
    // file params.output

    """
    ${binDir}aucell_genesigs_from_folder.py \
        $exprMat \
        $regulonsFolder \
        -o "multi_runs_regulons_auc_${type}.tsv" \
        --regulon-type ${type} \
        --auc-threshold ${params.sc.scenic.aucell.auc_threshold} \
        ${params.sc.scenic.aucell.containsKey('percentile_threshold') ? "--percentile-threshold " + params.sc.scenic.aucell.percentile_threshold : ""} \
        --gene-occurence-threshold ${params.sc.scenic.aucell.gene_occurence_threshold} \
        --num-workers ${params.sc.scenic.numWorkers} \
        --cell-id-attribute ${params.sc.scenic.cell_id_attribute} \
        --gene-attribute ${params.sc.scenic.gene_attribute}
    """
}
