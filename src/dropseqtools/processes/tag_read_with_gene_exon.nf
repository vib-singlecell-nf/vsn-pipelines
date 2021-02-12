nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON {

    publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(bam)
        file(annotation)

    output:
        tuple val(sample), path("*.merged_gene-exon-tagged.bam")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("dropseqtools").tag_read_with_gene_exon)
		processParams = sampleParams.local
        """
        source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
        software load drop-seq_tools/1.12
        TagReadWithGeneExon \
            I=${bam} \
            O=${sample}.merged_gene-exon-tagged.bam \
            ANNOTATIONS_FILE=${annotation} \
            TAG=${processParams.tag}
        """

} 
