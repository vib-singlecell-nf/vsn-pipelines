nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/flybaser/bin/"
} else {
    binDir = "/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/GitHub/SingleCellTxBenchmark/src/flybaser/bin/"
}

process FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL {
    
    container params.flybaser.container
    publishDir "intermediate/data", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f)
    
    output:
        tuple val(sampleId), path("${sampleId}.FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.flybaser.convert_fbgn_to_gene_symbol)
		processParams = sampleParams.local
        """
        ${binDir}convertFBgnToGeneSymbol.R \
            ${f} \
            ${processParams.columnName} \
            "${sampleId}.FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.tsv"
        """

}
