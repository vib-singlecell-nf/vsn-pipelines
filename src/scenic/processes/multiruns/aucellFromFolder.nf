nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic
def processParams = toolParams.aucell

process AUCELL_FROM_FOLDER {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_aucell/", mode: 'link', overwrite: true
    label 'compute_resources__scenic_multiruns'

    input:
        tuple val(sampleId), path(filteredLoom), path(multiRunsAggrRegulonsFolder)
        val type

    output:
        tuple val(sampleId), path("multi_runs_regulons_auc_${type}.tsv")

    script:
        """
        ${binDir}aucell_from_folder.py \
            $filteredLoom \
            $multiRunsAggrRegulonsFolder \
            -o "multi_runs_regulons_auc_${type}.tsv" \
            --min-genes ${processParams.min_genes_regulon} \
            --auc-threshold ${processParams.auc_threshold} \
            ${processParams.containsKey('percentile_threshold') ? "--percentile-threshold " + processParams.percentile_threshold : ""} \
            --min-regulon-gene-occurrence ${processParams.min_regulon_gene_occurrence} \
            --num-workers ${task.cpus} \
            --cell-id-attribute ${toolParams.cell_id_attribute} \
            --gene-attribute ${toolParams.gene_attribute}
        """

}
