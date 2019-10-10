nextflow.preview.dsl=2

process SC__STAR__UNLOAD_GENOME {

  container params.sc.star.container

  input:
    file(transcriptome)
    val all_done

  script:
    """

    echo "--genomeDir ${transcriptome}"

    STAR \
      --genomeLoad Remove \
      --genomeDir ${transcriptome}

    """
}
