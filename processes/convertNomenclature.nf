nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/flybaser/bin/"
} else {
    binDir = ""
}

process FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL {
    
    container params.flybaser.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    maxForks 1

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
