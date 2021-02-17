nextflow.enable.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/flybaser/bin/"
} else {
    binDir = ""
}

process FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL {
    
    container params.getToolParams("flybaser").container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)
    
    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("flybaser").convert_fbgn_to_gene_symbol)
		processParams = sampleParams.local
        """
        ${binDir}convertFBgnToGeneSymbol.R \
            ${f} \
            ${processParams.columnName} \
            "${sampleId}.FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.tsv"
        """

}
