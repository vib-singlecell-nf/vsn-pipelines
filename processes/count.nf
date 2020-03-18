nextflow.preview.dsl=2

toolParams = params.sc.cellranger

def generateCellRangerCountCommandDefaults = {
	processParams,
	transcriptome,
	expectCells ->
	return (
		"""
		cellranger count \
			--transcriptome=${transcriptome} \
			${(processParams.containsKey('expectCells')) ? '--expect-cells ' + processParams.expectCells: ''} \
			${(expectCells && !processParams.containsKey('expectCells')) ? '--expect-cells ' + expectCells: ''} \
			${(processParams.containsKey('forceCells')) ? '--force-cells ' + processParams.forceCells: ''} \
			${(processParams.containsKey('nosecondary')) ? '--nosecondary ' + processParams.nosecondary: ''} \
			${(processParams.containsKey('chemistry')) ? '--chemistry ' + processParams.chemistry: ''} \
			${(processParams.containsKey('r1Length')) ? '--r1-length ' + processParams.r1Length: ''} \
			${(processParams.containsKey('r2Length')) ? '--r2-length ' + processParams.r2Length: ''} \
			${(processParams.containsKey('lanes')) ? '--lanes ' + processParams.lanes: ''} \
			${(processParams.containsKey('localCores')) ? '--localcores ' + processParams.localCores: ''} \
			${(processParams.containsKey('localMem')) ? '--localmem ' + processParams.localMem: ''} \
		"""
	)
}

def runCellRangerCount = {
	processParams,
	transcriptome,
	id,
	sample,
	fastqs,
	expectCells = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, expectCells) + \
		"""	\
		--id=${id} \
		--sample=${sample} \
		--fastqs=${fastqs.join(",")}
		"""
	)
}

def runCellRangerCountLibraries = {
	processParams,
	transcriptome,
	id,
	featureRef,
	libraries = null,
	expectCells = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, expectCells) + \
		""" \
		--id ${id} \
		--libraries ${libraries} \
		--feature-ref ${featureRef}
		"""
	)
}

process SC__CELLRANGER__COUNT {

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
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
			processParams,
			transcriptome,
			sampleId,
			sampleId,
			fastqs
		)

}

process SC__CELLRANGER__COUNT_WITH_LIBRARIES {

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

    input:
		file(transcriptome)
		tuple \
			val(sampleId), \
			file(featureRef), \
			file(libraries)

  	output:
    	tuple val(sampleId), file("${sampleId}/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local
		if(processParams.sample == '') {
			throw new Exception("Regards params.sc.cellranger.count: sample parameter cannot be empty")
		}
		runCellRangerCountLibraries(
			processParams,
			transcriptome,
			sampleId,
			featureRef,
			libraries
		)

}

process SC__CELLRANGER__COUNT_WITH_METADATA {

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

    input:
		file(transcriptome)
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
			processParams,
			transcriptome,
			sampleId,
			samplePrefix,
			fastqs,
			expectCells
		)

}