nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SAVE_MULTI_RUNS_TO_LOOM {

    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}/multi_runs_looms/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=${params.sc.scenic.save_to_loom.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

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
        --min-genes-regulon ${params.sc.scenic.aucell.min_genes_regulon} \
        --min-regulon-gene-occurrence ${params.sc.scenic.aucell.min_regulon_gene_occurrence} \
        --cell-id-attribute ${params.sc.scenic.cell_id_attribute} \
        --gene-attribute ${params.sc.scenic.gene_attribute} \
        --title "${sampleId} - pySCENIC (${type})" \
        --nomenclature "${params.sc.scope.genome}" \
        --scope-tree-level-1 "${params.sc.scope.tree.level_1}" \
        --scope-tree-level-2 "${params.sc.scope.tree.level_2}" \
        --scope-tree-level-3 "${params.sc.scope.tree.level_3}"
    """
    
}
