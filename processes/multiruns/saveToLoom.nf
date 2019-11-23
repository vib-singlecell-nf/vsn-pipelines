nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.save_to_loom

process SAVE_MULTI_RUNS_TO_LOOM {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_looms/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=${processParams.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sampleId), file(filteredLoom), file(multiRunsAggrRegulons), file(multiRunsAggrRegulonsAUC)
    val type

    output:
    tuple val(sampleId), file("multi_runs_regulons_auc_${type}.loom")

    """
    ${binDir}save_multi_runs_to_loom.py \
        $filteredLoom \
        $multiRunsAggrRegulons \
        $multiRunsAggrRegulonsAUC \
        -o "multi_runs_regulons_auc_${type}.loom" \
        --min-genes-regulon ${toolParams.aucell.min_genes_regulon} \
        --min-regulon-gene-occurrence ${toolParams.aucell.min_regulon_gene_occurrence} \
        --cell-id-attribute ${toolParams.cell_id_attribute} \
        --gene-attribute ${toolParams.gene_attribute} \
        --title "${sampleId} - pySCENIC (${type})" \
        --nomenclature "${params.sc.scope.genome}" \
        --scope-tree-level-1 "${params.sc.scope.tree.level_1}" \
        --scope-tree-level-2 "${params.sc.scope.tree.level_2}" \
        --scope-tree-level-3 "${params.sc.scope.tree.level_3}"
    """

}
