nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

toolParams = params.sc.scenic
processParams = params.sc.scenic.aucell

process AUCELL_FROM_FOLDER {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/multi_runs_aucell/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

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
            --num-workers ${toolParams.numWorkers} \
            --cell-id-attribute ${toolParams.cell_id_attribute} \
            --gene-attribute ${toolParams.gene_attribute}
        """

}
