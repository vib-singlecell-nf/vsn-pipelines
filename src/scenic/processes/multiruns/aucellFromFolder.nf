nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.aucell

process AUCELL_FROM_FOLDER {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label toolParams.labels.processExecutor
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_aucell/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

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
        --min-genes ${processParams.min_genes_regulon} \
        --auc-threshold ${processParams.auc_threshold} \
        ${processParams.containsKey('percentile_threshold') ? "--percentile-threshold " + processParams.percentile_threshold : ""} \
        --min-regulon-gene-occurrence ${processParams.min_regulon_gene_occurrence} \
        --num-workers ${toolParams.numWorkers} \
        --cell-id-attribute ${toolParams.cell_id_attribute} \
        --gene-attribute ${toolParams.gene_attribute}
    """

}
