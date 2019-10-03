nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__GRNBOOST2WITHOUTDASK {
    cache 'deep'
    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    val runId
    file filteredloom
    file tfs

    output:
    file "${params.sc.scenic.numRuns > 1 ? "run_" + runId +"_adj.tsv" : "adj.tsv"}"

    script:
    """
    ${binDir}grnboost2_without_dask.py \
        $filteredloom \
        $tfs \
        --output ${params.sc.scenic.numRuns > 1 ? "run_" + runId +"_adj.tsv" : "adj.tsv"} \
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
