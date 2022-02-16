process BCL2FASTQ__DEMULTIPLEX {

	container params.tools.bcl2fastq.container
    publishDir "${params.global.outdir}/fastqs-bcl2fastq"

	input:
    	// tuple \
		// 	val(runName), \
		// 	path(sampleSheet)
    	tuple \
			path(runFolder), \
			path(sampleSheet)
            
	output:
		tuple \
			path("**.fastq.gz"), \
			path("Stats/DemultiplexingStats.xml")

  	script:
		// def sampleParams = params.parseConfig(sampleId, params.global, params.tools.bcl2fastq.demultiplex)
		// processParams = sampleParams.local
        processParams = params.tools.bcl2fastq.demultiplex

        // Actually run the demultiplexing
		def mask_reads_param = processParams.mask_short_adapter_reads? "--mask-short-adapter-reads=1" : ""
        def tiles_param = processParams.lanes? "--tiles=s_" + processParams.lanes.join(',s_') : ""
        def lane_splitting_param = processParams.split_lanes? "": "--no-lane-splitting"
        def create_fastq_for_index_reads_param = processParams.create_fastq_for_index_reads? "--create-fastq-for-index-reads": ""

        // Actually run the demultiplexing
        """
        bcl2fastq \
            --runfolder-dir=${runFolder} \
            --sample-sheet=${sampleSheet} \
            --output-dir=./ \
            --processing-threads=${task.cpus} \
            ${mask_reads_param} \
            ${tiles_param} \
            ${lane_splitting_param} \
            ${create_fastq_for_index_reads_param} \
            --barcode-mismatches=${processParams.barcode_mismatches}
            --sample-sheet=${sampleSheet}
        """
}