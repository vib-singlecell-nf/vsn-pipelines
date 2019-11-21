nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process AUCELL_FROM_FOLDER {

    // Process will be submitted as job if params.sc.scenic.labels.processExecutor = 'qsub' (default)
    label params.sc.scenic.labels.processExecutor
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}/multi_runs_aucell/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sampleId), file(filteredLoom), file(multiRunsAggrRegulonsFolder)
    val type

    output:
    tuple val(sampleId), file("multi_runs_regulons_auc_${type}.tsv")

    """
    ${binDir}aucell_from_folder.py \
        $filteredLoom \
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
