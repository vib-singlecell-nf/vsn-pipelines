nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

toolParams = params.tools.scenic

process SAVE_MULTI_RUNS_TO_LOOM {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_looms/", mode: 'link', overwrite: true
    label 'compute_resources__scenic_multiruns'

    input:
		tuple val(sampleId), path(filteredLoom), path(multiRunsAggrRegulons), path(multiRunsAggrRegulonsAUC)
		val type

    output:
    	tuple val(sampleId), path("multi_runs_regulons_auc_${type}.loom")

	script:
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
