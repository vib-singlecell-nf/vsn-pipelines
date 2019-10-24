nextflow.preview.dsl=2

process SC__STAR__MAP_COUNT {

  container params.sc.star.container
  maxForks 2

  input:
    file(starIndex)
    val starIndexLoaded
    tuple val(sample), file(fastqs)

  output:
    val success
    file file("*ReadsPerGene.out.tab") optional params.sc.star.map_count.containsKey('quantMode') && params.sc.star.map_count.quantMode == "GeneCounts" ? true: false 
    tuple val(sample), file("*.STAR_Aligned.sortedByCoord.out.bam")

  script:
    success = true

    """
    STAR \
      --genomeLoad LoadAndKeep \
      --genomeDir ${starIndex} \
      ${(params.sc.star.map_count.containsKey('runThreadN')) ? '--runThreadN ' + params.sc.star.map_count.runThreadN: ''} \
      ${(params.sc.star.map_count.containsKey('limitBAMsortRAM')) ? '--limitBAMsortRAM ' + params.sc.star.map_count.limitBAMsortRAM: ''} \
      ${(params.sc.star.map_count.containsKey('outSAMtype')) ? '--outSAMtype ' + params.sc.star.map_count.outSAMtype: ''} \
      ${(params.sc.star.map_count.containsKey('quantMode')) ? '--quantMode ' + params.sc.star.map_count.quantMode: ''} \
      ${(params.sc.star.map_count.containsKey('outReadsUnmapped')) ? '--outReadsUnmapped ' + params.sc.star.map_count.outReadsUnmapped: ''} \
      --readFilesIn ${fastqs} \
      ${(fastqs.name.endsWith(".gz")) ? '--readFilesCommand zcat' : ''} \
      --outFileNamePrefix ${sample}.STAR_
    """
}
