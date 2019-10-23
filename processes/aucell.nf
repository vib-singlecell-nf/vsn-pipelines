nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

process SC__SCENIC__AUCELL {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/aucell/${params.sc.scenic.numRuns > 1 ? "run_" + runId : ""}", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks params.sc.scenic.maxForks

    input:
    val runId
    file exprMat
    file regulons
    val type

    output:
    file "${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__auc_" + type + ".loom": "auc_" + type + ".loom"}"
    // file params.output

    """
    pyscenic aucell \
        $exprMat \
        $regulons \
        -o "${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__auc_" + type + ".loom": "auc_" + type + ".loom"}" \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute} \
        --num_workers ${params.sc.scenic.numWorkers}
    """
}
