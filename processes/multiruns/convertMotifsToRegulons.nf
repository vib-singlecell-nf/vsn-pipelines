nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_cistarget/", mode: 'link', overwrite: true
    // This process requires a large amount of memory especially for big datasets (force to use bigmem node)
    // This process is quite slow (could take more than 1h for big datasets, so keep 24h for now)
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=${params.sc.scenic.motifs_to_regulons.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file multiRunsAggrMotifEnrichmentTable
    file multiRunsAggrRegulonsFolder
    val type

    output:
    file "multi_runs_regulons_${type}.pkl.gz"

    """
    ${binDir}convert_multi_runs_features_to_regulons.py \
        $multiRunsAggrMotifEnrichmentTable \
        $multiRunsAggrRegulonsFolder \
        -o "multi_runs_regulons_${type}.pkl.gz"
        # --min-genes-regulon ${params.sc.scenic.aucell.min_genes_regulon} \
        # --min-regulon-gene-occurrence ${params.sc.scenic.aucell.min_regulon_gene_occurrence}
    """
}

