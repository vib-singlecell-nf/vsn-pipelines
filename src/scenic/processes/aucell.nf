nextflow.preview.dsl=2

process AUCELL {

    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}/aucell/${params.sc.scenic.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=${params.sc.scenic.aucell.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks params.sc.scenic.aucell.maxForks

    input:
    tuple val(sampleId), file(filteredLoom), file(regulons), val(runId)
    val type

    output:
    tuple val(sampleId), file(filteredLoom), file("${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__auc_" + type + ".loom": sampleId + "__auc_" + type + ".loom"}"), val(runId)

    """
    pyscenic aucell \
        $filteredLoom \
        $regulons \
        -o "${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__auc_" + type + ".loom": sampleId + "__auc_" + type + ".loom"}" \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute} \
        --num_workers ${params.sc.scenic.numWorkers}
    """

}
