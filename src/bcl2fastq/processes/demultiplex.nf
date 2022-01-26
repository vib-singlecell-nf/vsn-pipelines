process BCL2FASTQ__DEMULTIPLEX {

	container params.tools.bcl2fastq.container

	input:
    	tuple \
			val(runName), \
			path(sampleSheet)

	output:
		tuple \
			path(fastqs), \
			path(logs)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.tools.bcl2fastq.demultiplex)
		processParams = sampleParams.local

        // Actually run the demultiplexing