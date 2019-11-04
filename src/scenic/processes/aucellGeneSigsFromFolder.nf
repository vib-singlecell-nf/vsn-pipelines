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
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_aucell/", mode: 'copy'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file exprMat
    file multiRunsAggrRegulonsFolder
    val type

    output:
    file "multi_runs_regulons_auc_${type}.tsv"
    // file params.output

    """
    ${binDir}aucell_genesigs_from_folder.py \
        $exprMat \
        $multiRunsAggrRegulonsFolder \
        -o "multi_runs_regulons_auc_${type}.tsv" \
        --min-genes ${params.sc.scenic.aucell.min_genes_regulon} \
        --auc-threshold ${params.sc.scenic.aucell.auc_threshold} \
        ${params.sc.scenic.aucell.containsKey('percentile_threshold') ? "--percentile-threshold " + params.sc.scenic.aucell.percentile_threshold : ""} \
        --min-regulon-gene-occurrence ${params.sc.scenic.aucell.min_regulon_gene_occurrence} \
        --num-workers ${params.sc.scenic.numWorkers} \
        --cell-id-attribute ${params.sc.scenic.cell_id_attribute} \
        --gene-attribute ${params.sc.scenic.gene_attribute}
    """
}
