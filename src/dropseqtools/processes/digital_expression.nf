nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__DIGITAL_EXPRESSION {

    container params.getToolParams("dropseqtools").container
    publishDir "03.count", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sample), path(bam), val(tag), path(selectedBarcodes)

    output:
        tuple file("*.${tag}.cells_dge.txt.gz"), emit: dgem

    shell:
        """
        DigitalExpression \
            I=${bam} \
            O=${sample}.${tag}.cells_dge.txt.gz \
            SUMMARY=${sample}.${tag}.cells_dge.summary.txt \
            CELL_BC_FILE=${selectedBarcodes} \
            TMP_DIR=$DWMAX/tmp
        """

}
