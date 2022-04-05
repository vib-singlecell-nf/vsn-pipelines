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
		path("**.fastq.gz"), emit: fastqs
		path("Stats/DemultiplexingStats.xml"), emit: stats

  	script:
        processParams = params.tools.bcl2fastq.demultiplex

        // Actually run the demultiplexing
        """
        bcl2fastq \
            --runfolder-dir=${runFolder} \
            --sample-sheet=${sampleSheet} \
            --output-dir=./ \
            --processing-threads=${task.cpus} \
            ${processParams.mask_short_adapter_reads? "--mask-short-adapter-reads=1" : ""} \
            ${processParams.lanes? "--tiles=s_" + processParams.lanes.join(',s_') : ""} \
            ${processParams.split_lanes? "": "--no-lane-splitting"} \
            ${processParams.create_fastq_for_index_reads? "--create-fastq-for-index-reads": ""} \
            ${processParams.barcode_mismatches? "--barcode-mismatches=" + processParams.barcode_mismatches : ""} \
            --ignore-missing-filter
        """
}