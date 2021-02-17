nextflow.enable.dsl=2

process SC__STAR__MAP_COUNT {

	container params.getToolParams("star").container
    label 'compute_resources__star_map_count'

	input:
		file(starIndex)
		val starIndexLoaded
		tuple val(sample), path(fastqs)

	output:
		val success, emit: isDone
		tuple val(sample), path("*ReadsPerGene.out.tab"), emit: counts optional processParams.containsKey('quantMode') && processParams.quantMode == "GeneCounts" ? true: false
		tuple val(sample), path("*.STAR_Aligned.sortedByCoord.out.bam"), emit: bam

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("star").map_count)
		processParams = sampleParams.local
		success = true
		"""
		STAR \
			--genomeLoad LoadAndKeep \
			--genomeDir ${starIndex} \
            --runThreadN ${task.cpus} \
			${(processParams.containsKey('limitBAMsortRAM')) ? '--limitBAMsortRAM ' + processParams.limitBAMsortRAM: ''} \
			${(processParams.containsKey('outSAMtype')) ? '--outSAMtype ' + processParams.outSAMtype: ''} \
			${(processParams.containsKey('quantMode')) ? '--quantMode ' + processParams.quantMode: ''} \
			${(processParams.containsKey('outReadsUnmapped')) ? '--outReadsUnmapped ' + processParams.outReadsUnmapped: ''} \
			--readFilesIn ${fastqs} \
			${(fastqs.name.endsWith(".gz")) ? '--readFilesCommand zcat' : ''} \
			--outFileNamePrefix ${sample}.STAR_
		"""

}
