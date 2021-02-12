nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT {
    
    container params.getToolParams("dropseqtools").container
    publishDir "${params.global.outdir}/00.refdata", mode: 'symlink'
    label 'compute_resources__default'

    input:
        file(annotation)
        file(seqdict)
    
    output:
        file("${seqdict.baseName}.refFlat")
    
    script:
        """
        ConvertToRefFlat \
            ANNOTATIONS_FILE=${annotation} \
            SEQUENCE_DICTIONARY=${seqdict} \
            OUTPUT=${seqdict.baseName}.refFlat
        """

}
