nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process GRNBOOST2_WITHOUT_DASK {

    // Process will be submitted as job if params.sc.scenic.labels.processExecutor = 'qsub' (default)
    label params.sc.scenic.labels.processExecutor
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}/grnboost2withoutDask/${params.sc.scenic.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=${params.sc.scenic.grn.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks params.sc.scenic.grn.maxForks
    
    input:
    tuple val(sampleId), file(filteredLoom), val(runId)
    file tfs

    output:
    tuple val(sampleId), file(filteredLoom), file("${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"}"), val(runId)

    script:
    """
    ${binDir}grnboost2_without_dask.py \
        $filteredLoom \
        $tfs \
        --output ${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"} \
        --num_workers ${params.sc.scenic.numWorkers} \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute}
    """

}

/* options to implement:
        --seed ${params.grn.seed} \

flag parameters not yet implemented:
        --transpose
*/
