nextflow.enable.dsl=2
process SC__STAR__SOLO_MAP_COUNT {

  container params.tools.star.container
  label 'compute_resources__star_map_count'

  input:
    file(transcriptome)
    val genome_loaded
    file(fastqs)

  output:
    val success
    // file '*ReadsPerGene.out.tab'

  script:
    sample = fastqs.getName()
    _sampleName = sample
    success = true

    """
    STAR \
      --genomeLoad LoadAndKeep \
      --soloType Droplet \
      --genomeDir ${transcriptome} \
      --runThreadN ${task.cpus} \
      ${(params.tools.star.map_count.containsKey('limitBAMsortRAM')) ? '--limitBAMsortRAM ' + params.tools.star.map_count.limitBAMsortRAM: ''} \
      ${(params.tools.star.map_count.containsKey('outSAMtype')) ? '--outSAMtype ' + params.tools.star.map_count.outSAMtype: ''} \
      ${(params.tools.star.map_count.containsKey('quantMode')) ? '--quantMode ' + params.tools.star.map_count.quantMode: ''} \
      ${(params.tools.star.map_count.containsKey('outReadsUnmapped')) ? '--outReadsUnmapped ' + params.tools.star.map_count.outReadsUnmapped: ''} \
      --readFilesIn ${fastqs} \
      ${(fastqs.name.endsWith(".gz")) ? '--readFilesCommand zcat' : ''} \
      --outFileNamePrefix ${_sampleName}
    """
}
