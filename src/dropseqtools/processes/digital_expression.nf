nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__DIGITAL_EXPRESSION {

    container params.sc.dropseqtools.container
    publishDir "03.count", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

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