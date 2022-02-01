nextflow.enable.dsl=2
process SC__STAR__SOLO_MAP_COUNT {

  container params.tools.star.container
  publishDir "${params.global.outdir}/STARSolo/", mode: 'copy'
  label 'compute_resources__solo_map_count'

  input:
      tuple val(sample), path(fastqs)
      path(soloCBwhitelist)


  output:
    tuple val(sample), path("${sample}_Solo.out"), emit: solo_outputs
    tuple val(sample), path("${sample}_Aligned.sortedByCoord.out.bam"), emit: bam

  script:

    """
    STAR \
      --genomeDir ${params.tools.star.map_count.index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${fastqs} \
      --soloCBwhitelist ${soloCBwhitelist.toList().join(' ')} \
      ${(fastqs[0].name.endsWith(".gz")) ? '--readFilesCommand zcat' : ''} \
      ${(params.tools.star.map_count.containsKey('genomeLoad')) ? '--genomeLoad ' + params.tools.star.map_count.genomeLoad: ''} \
      ${(params.tools.star.map_count.containsKey('limitBAMsortRAM')) ? '--limitBAMsortRAM ' + params.tools.star.map_count.limitBAMsortRAM: ''} \
      ${(params.tools.star.map_count.containsKey('outSAMtype')) ? '--outSAMtype ' + params.tools.star.map_count.outSAMtype: ''} \
      ${(params.tools.star.map_count.containsKey('quantMode')) ? '--quantMode ' + params.tools.star.map_count.quantMode: ''} \
      ${(params.tools.star.map_count.containsKey('outReadsUnmapped')) ? '--outReadsUnmapped ' + params.tools.star.map_count.outReadsUnmapped: ''} \
      ${(params.tools.star.map_count.containsKey('soloType')) ? '--soloType ' + params.tools.star.map_count.soloType: ''} \
      ${(params.tools.star.map_count.containsKey('soloCBposition')) ? '--soloCBposition ' + params.tools.star.map_count.soloCBposition: ''} \
      ${(params.tools.star.map_count.containsKey('soloUMIposition')) ? '--soloUMIposition ' + params.tools.star.map_count.soloUMIposition: ''} \
      ${(params.tools.star.map_count.containsKey('soloCellFilter')) ? '--soloCellFilter ' + params.tools.star.map_count.soloCellFilter: ''} \
      ${(params.tools.star.map_count.containsKey('soloCBmatchWLtype')) ? '--soloCBmatchWLtype ' + params.tools.star.map_count.soloCBmatchWLtype: ''} \
      ${(params.tools.star.map_count.containsKey('outFilterMultimapNmax')) ? '--outFilterMultimapNmax ' + params.tools.star.map_count.outFilterMultimapNmax: ''} \
      ${(params.tools.star.map_count.containsKey('outSAMattributes')) ? '--outSAMattributes ' + params.tools.star.map_count.outSAMattributes: ''} \
      ${(params.tools.star.map_count.containsKey('bamRemoveDuplicatesType')) ? '--bamRemoveDuplicatesType ' + params.tools.star.map_count.bamRemoveDuplicatesType: ''} \
      ${(params.tools.star.map_count.containsKey('soloFeatures')) ? '--soloFeatures ' + params.tools.star.map_count.soloFeatures: ''} \
      --outFileNamePrefix ${sample}_

    gzip ${sample}_Solo.out/Gene/*/*
    """
}
