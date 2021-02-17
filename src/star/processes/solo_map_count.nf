nextflow.enable.dsl=2
process SC__STAR__SOLO_MAP_COUNT {

  container params.getToolParams("star").container
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
      ${(params.getToolParams("star").map_count.containsKey('limitBAMsortRAM')) ? '--limitBAMsortRAM ' + params.getToolParams("star").map_count.limitBAMsortRAM: ''} \
      ${(params.getToolParams("star").map_count.containsKey('outSAMtype')) ? '--outSAMtype ' + params.getToolParams("star").map_count.outSAMtype: ''} \
      ${(params.getToolParams("star").map_count.containsKey('quantMode')) ? '--quantMode ' + params.getToolParams("star").map_count.quantMode: ''} \
      ${(params.getToolParams("star").map_count.containsKey('outReadsUnmapped')) ? '--outReadsUnmapped ' + params.getToolParams("star").map_count.outReadsUnmapped: ''} \
      --readFilesIn ${fastqs} \
      ${(fastqs.name.endsWith(".gz")) ? '--readFilesCommand zcat' : ''} \
      --outFileNamePrefix ${_sampleName}
    """
}
