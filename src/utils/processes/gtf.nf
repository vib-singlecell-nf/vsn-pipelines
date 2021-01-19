nextflow.enable.dsl=2

process FORMAT_GTF {

    publishDir "${params.global.outdir}/00.refdata", mode: 'symlink'
    label 'compute_resources__default'

    input:
        file(annotation)

    output:
        file "*.formatted.gtf"

    script:
        """
        sed -r 's/(.*); transcript_id (.*); (.*); gene_name (.*); \$/\\1; transcript_id \\2; \\3; gene_name \\4; transcript_name \\2;/' \
            ${annotation} \
            > ${annotation.baseName}.formatted.gtf
        """

}

process FORMAT_GTF_IGENOMES {

    publishDir "${params.global.outdir}/00.refdata", mode: 'symlink'
    label 'compute_resources__default'

    input:
        file(annotation)

    output:
        file "*.formatted.gtf"

    script:
        """
        sed -r 's/(.*); gene_name (.*); transcript_id (.*); (.*);\$/\\1; gene_name \\2; transcript_id \\3; \\4; transcript_name \\3;/' \
            ${annotation} \
            > ${annotation.baseName}.formatted.gtf
        """

}
