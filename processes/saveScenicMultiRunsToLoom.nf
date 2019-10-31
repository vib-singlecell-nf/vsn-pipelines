nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_looms/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file exprMat
    file multiRunsAggrMotifEnrichmentTable
    file multiRunsAggrRegulonsFolder
    file multiRunsAggrRegulonsAUC
    val type

    output:
    file "multi_runs_regulons_auc_${type}.loom"
    // file params.output

    """
    ${binDir}save_SCENIC_multi_runs_to_loom.py \
        $exprMat \
        $multiRunsAggrMotifEnrichmentTable \
        $multiRunsAggrRegulonsFolder \
        $multiRunsAggrRegulonsAUC \
        -o "multi_runs_regulons_auc_${type}.loom" \
        --min-genes-regulon ${params.sc.scenic.aucell.min_genes_regulon} \
        --gene-occurence-threshold ${params.sc.scenic.aucell.gene_occurence_threshold} \
        --cell-id-attribute ${params.sc.scenic.cell_id_attribute} \
        --gene-attribute ${params.sc.scenic.gene_attribute} \
        --title "${params.global.project_name} - pySCENIC (${type})" \
        --nomenclature "${params.sc.scope.genome}" \
        --scope-tree-level-1 "${params.sc.scope.tree.level_1}" \
        --scope-tree-level-2 "${params.sc.scope.tree.level_2}" \
        --scope-tree-level-3 "${params.sc.scope.tree.level_3}"
    """
}

