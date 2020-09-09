nextflow.preview.dsl=2

include {
    isParamNull;
} from './../../utils/processes/utils.nf' params(params)

toolParams = params.sc.cellranger


def generateCellRangerCountCommandDefaults = {
	processParams,
	transcriptome,
	expectCells,
	chemistry ->
	_expectCells = null
	// --expect-cells argument
	if(!isParamNull(expectCells)) {
		_expectCells = expectCells
	} else if (processParams.containsKey('expectCells')) {
		_expectCells = processParams.expectCells
	}
	_chemistry = null
	// --chemistry argument
	if(!isParamNull(chemistry)) {
		_chemistry = chemistry
	} else if (processParams.containsKey('chemistry')) {
		_chemistry = processParams.chemistry
	}
	return (
		"""
		cellranger count \
			--transcriptome=${transcriptome} \
			${_expectCells ? '--expect-cells ' + _expectCells: ''} \
			${_chemistry ? '--chemistry ' + _chemistry: ''} \
			${(processParams.containsKey('forceCells')) ? '--force-cells ' + processParams.forceCells: ''} \
			${(processParams.containsKey('nosecondary')) ? '--nosecondary ' + processParams.nosecondary: ''} \
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
	expectCells = null,
	chemistry = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, expectCells, chemistry) + \
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
	id,
	featureRef,
	libraries = null,
	expectCells = null,
	chemistry = null ->
	return (
		generateCellRangerCountCommandDefaults(processParams, transcriptome, expectCells, chemistry) + \
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
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${toolParams.count.ppn} -l pmem=${toolParams.count.pmem} -l walltime=${toolParams.count.walltime} -A ${params.global.qsubaccount} -m abe -M ${params.global.qsubemail}"
	maxForks = toolParams.count.maxForks

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
		// Check if the current sample has multiple sequencing runs
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
			sampleId,
			featureRef,
			"libraries.csv"
		)

}

process SC__CELLRANGER__COUNT_WITH_METADATA {
	
	def getClusterOptions = { globalParams, processParams ->
		clusterOptions = ""
		clusterOptions += "-l nodes=${processParams.containsKey("nodes") ? processParams.nodes : 1}:ppn=${processParams.containsKey("ppn") ? processParams.ppn : 1} "
		clusterOptions += "-l pmem=${processParams.containsKey("pmem") ? processParams.pmem : 1} "
		clusterOptions += "-l walltime=${processParams.containsKey("walltime") ? processParams.walltime : '1:00:00'} "
		if(!globalParams.containsKey("qsubaccount"))
			throw new Exception("The param `qsubaccount` is not defined in the enclosing params.global config.")
		clusterOptions += "-A ${globalParams.qsubaccount} "
		clusterOptions += "-m abe "
		clusterOptions += "${globalParams.containsKey("qsubemail") && globalParams.qsubemail != '' ? '-M '+ globalParams.qsubemail : ''}"
		return clusterOptions
	}

	label toolParams.labels.processExecutor
	cache 'deep'
	container toolParams.container
	publishDir "${params.global.outdir}/counts", saveAs: {"${sampleId}/outs"}, mode: 'link', overwrite: true
	clusterOptions "${getClusterOptions(params.global,toolParams.count)}"
	maxForks = toolParams.count.maxForks

    input:
		path(transcriptome)
		tuple \
			val(sampleId), \
			val(samplePrefix), \
			path(fastqs, stageAs: "fastqs_??/*"), \
			val(expectCells), \
			val(chemistry)

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
			sampleId,
			samplePrefix,
			fastqs,
			expectCells,
			chemistry
		)

}