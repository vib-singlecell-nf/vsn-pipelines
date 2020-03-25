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
		--id=${id}_out \
		--sample=${sample} \
		--fastqs=${fastqs}
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
		--id ${id}_out \
		--libraries ${libraries} \
		--feature-ref ${featureRef}
		"""
	)
}

process SC__CELLRANGER__COUNT {

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

    input:
		path(transcriptome)
		tuple \
			val(sampleId), 
			path(fastqs, stageAs: "fastqs_??/*")

  	output:
    	tuple val(sampleId), path("${sampleId}_out/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local
		if(processParams.sample == '') {
			throw new Exception("Regards params.sc.cellranger.count: sample parameter cannot be empty")
		}
		fastqs = fastqs instanceof List ? fastqs.join(',') : fastqs
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
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

    input:
		path(transcriptome)
		path(featureRef)
		tuple \
			val(sampleId), \
			path(fastqs, stageAs: "fastqs_??/*"),
			val(sampleNames),
			val(assays)

  	output:
    	tuple val(sampleId), path("${sampleId}_out/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local

		if(processParams.sample == '') {
			throw new Exception("Regards params.sc.cellranger.count: sample parameter cannot be empty")
		}

		// We need to create the libraries.csv file here because it needs absolute paths

		csvData = "fastqs,sample,library_type\n"
		fastqs.eachWithIndex { fastq, ix -> 
			csvData += "\$PWD/${fastq},${sampleNames[ix]},${assays[ix]}\n"
		}

		"""
		echo "${csvData}" > libraries.csv
		""" + runCellRangerCountLibraries(
			processParams,
			transcriptome,
			sampleId,
			featureRef,
			"libraries.csv"
		)

}

process SC__CELLRANGER__COUNT_WITH_METADATA {

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

    input:
		path(transcriptome)
		tuple \
			val(sampleId), \
			val(samplePrefix), \
			path(fastqs), \
			val(expectCells)

  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}/outs")

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