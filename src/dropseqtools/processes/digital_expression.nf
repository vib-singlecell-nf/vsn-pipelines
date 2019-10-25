nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__DIGITAL_EXPRESSION {
    publishDir "03.count", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam), val(tag), file(selectedBarcodes)
    output:
        tuple file("*.${tag}.cells_dge.txt.gz"), emit: dgem
    shell:
        """
        source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
        DigitalExpression \
            I=${bam} \
            O=${sample}.${tag}.cells_dge.txt.gz \
            SUMMARY=${sample}.${tag}.cells_dge.summary.txt \
            CELL_BC_FILE=${selectedBarcodes} \
            TMP_DIR=$DWMAX/tmp
        """
}