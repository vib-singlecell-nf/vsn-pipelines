nextflow.preview.dsl=2

process SC__STAR__LOAD_GENOME {
  container params.sc.star.container

  input:
    file(transcriptome)

  output:
    val genomeLoaded 

  script:
    genomeLoaded = true
    """
    STAR \
      --genomeLoad LoadAndExit \
      --genomeDir ${transcriptome}
    """
}
