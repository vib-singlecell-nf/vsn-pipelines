nextflow.preview.dsl=2

toolParams = params.sc.cellranger

def runCellRangerCount = {
	id,
	sample,
	fastqs,
	processParams,
	expectCells = null ->
	return (
		"""
		# --id: CellRanger will create a directory with this name in cellranger_parent_output_dir.
		# --sample: Start of FASTQ filenames that identifies a sample uniquely (multiple prefixes separated by ",").
		cellranger count \
			--id=${id} \
			--sample=${sample} \
			--fastqs=${fastqs.join(",")} \
			--transcriptome=${processParams.transcriptome} \
			${(processParams.containsKey('libraries')) ? '--libraries ' + processParams.libraries: ''} \
			${(processParams.containsKey('featureRef')) ? '--feature-ref ' + processParams.featureRef: ''} \
			${(processParams.containsKey('expectCells')) ? '--expect-cells ' + processParams.expectCells: ''} \
			${(expectCells && !processParams.containsKey('expectCells')) ? '--expect-cells ' + expectCells: ''} \
			${(processParams.containsKey('forceCells')) ? '--force-cells ' + processParams.forceCells: ''} \
			${(processParams.containsKey('nosecondary')) ? '--nosecondary ' + processParams.nosecondary: ''} \
			${(processParams.containsKey('noLibraries')) ? '--no-libraries ' + processParams.noLibraries: ''} \
			${(processParams.containsKey('chemistry')) ? '--chemistry ' + processParams.chemistry: ''} \
			${(processParams.containsKey('r1Length')) ? '--r1-length ' + processParams.r1Length: ''} \
			${(processParams.containsKey('r2Length')) ? '--r2-length ' + processParams.r2Length: ''} \
			${(processParams.containsKey('lanes')) ? '--lanes ' + processParams.lanes: ''} \
			${(processParams.containsKey('localCores')) ? '--localcores ' + processParams.localCores: ''} \
			${(processParams.containsKey('localMem')) ? '--localmem ' + processParams.localMem: ''} \
			${(processParams.containsKey('indicies')) ? '--indicies ' + processParams.indicies: ''} 
		"""
	)
}

process SC__CELLRANGER__COUNT {

	  label toolParams.labels.processExecutor
	  cache 'deep'
	  container toolParams.container
	  publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	  clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount}"
	  maxForks = toolParams.count.maxForks

    input:
		file(transcriptome)
		tuple val(sampleId), file(fastqs)

  	output:
    	tuple val(sampleId), file("${sampleId}/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local
		if(processParams.sample == '') {
			throw new Exception("Regards params.sc.cellranger.count: sample parameter cannot be empty")
		}
		runCellRangerCount(
			sampleId,
			sampleId,
			fastqs,
			processParams
		)

}

process SC__CELLRANGER__COUNT_WITH_METADATA {

	  label toolParams.labels.processExecutor
	  cache 'deep'
	  container toolParams.container
	  publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	  clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount}"
	  maxForks = toolParams.count.maxForks

    input:
		tuple \
			val(sampleId), \
			val(samplePrefix), \
			file(fastqs), \
			val(expectCells)

  	output:
    	tuple \
			val(sampleId), \
			file("${sampleId}/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local
		runCellRangerCount(
			sampleId,
			samplePrefix,
			fastqs,
			processParams,
			expectCells
		)

}