nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
        file(annotation)
    output:
        tuple val(sample), file("*.merged_gene-exon-tagged.bam")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagReadWithGeneExon \
			I=${bam} \
			O=${sample}.merged_gene-exon-tagged.bam \
			ANNOTATIONS_FILE=${annotation} \
			TAG=${params.dropseqtools.tag_read_with_gene_exon.tag}
        """    
} 
