nextflow.preview.dsl=2

toolParams = params.sc.cellranger

def generateCellRangerCountCommandDefaults = {
	processParams,
	transcriptome,
    task,
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
            --localcores=${task.cpus} \
            --localmem=${task.memory.toGiga()} \
		"""
	)
}

def runCellRangerCount = {
	processParams,
	transcriptome,
    task,
	id,
	sample,
	fastqs,
	expectCells = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, task, expectCells) + \
		"""	\
		--id=${id} \
		--sample=${sample} \
		--fastqs=${fastqs}
		"""
	)
}

def runCellRangerCountLibraries = {
	processParams,
	transcriptome,
    task,
	id,
	featureRef,
	libraries = null,
	expectCells = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, task, expectCells) + \
		""" \
		--id ${id} \
		--libraries ${libraries} \
		--feature-ref ${featureRef}
		"""
	)
}

process SC__CELLRANGER__COUNT {

	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
    label 'compute_resources__cellranger_count'

    input:
		path(transcriptome)
		tuple \
			val(sampleId), 
			path(fastqs, stageAs: "fastqs_??/*")

  	output:
    	tuple val(sampleId), path("${sampleId}/outs")

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
            task,
			sampleId,
			sampleId,
			fastqs
		)

}

process SC__CELLRANGER__COUNT_WITH_LIBRARIES {

	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
    label 'compute_resources__cellranger'

    input:
		path(transcriptome)
		path(featureRef)
		tuple \
			val(sampleId), \
			path(fastqs, stageAs: "fastqs_??/*"),
			val(sampleNames),
			val(assays)

  	output:
    	tuple val(sampleId), path("${sampleId}/outs")

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
            task,
			sampleId,
			featureRef,
			"libraries.csv"
		)

}

process SC__CELLRANGER__COUNT_WITH_METADATA {

	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
    label 'compute_resources__cellranger'

    input:
		path(transcriptome)
		tuple \
			val(sampleId), \
			val(samplePrefix), \
			path(fastqs, stageAs: "fastqs_??/*"), \
			val(expectCells)

  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, toolParams.count)
		processParams = sampleParams.local
		fastqs = fastqs instanceof List ? fastqs.join(',') : fastqs
		runCellRangerCount(
			processParams,
			transcriptome,
            task,
			sampleId,
			samplePrefix,
			fastqs,
			expectCells
		)

}
