nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT {
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
        file(seqdict)
    output:
        file("${seqdict.baseName}.refFlat")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		ConvertToRefFlat \
			ANNOTATIONS_FILE=${annotation} \
            SEQUENCE_DICTIONARY=${seqdict} \
            OUTPUT=${seqdict.baseName}.refFlat
        """
}