nextflow.preview.dsl=2

process SC__CELLRANGER__COUNT {

	label params.sc.cellranger.labels.processExecutor
	cache 'deep'
	container params.sc.cellranger.container
	publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
	clusterOptions "-l nodes=1:ppn=${params.sc.cellranger.count.ppn} -l pmem=${params.sc.cellranger.count.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
	maxForks = params.sc.cellranger.count.maxForks

  	input:
		file(transcriptome)
		tuple val(sampleId), file(fastqs)

  	output:
    	tuple val(sampleId), file("${sampleId}/outs")

  	script:
	  	def sampleParams = params.parseConfig(sampleId, params.global, params.sc.cellranger.count)
		processParams = sampleParams.local
		"""
		cellranger count \
			--id=${sampleId} \
			--sample=${sampleId} \
			--fastqs=${fastqs.join(",")} \
			--transcriptome=${transcriptome} \
			${(processParams.containsKey('libraries')) ? '--libraries ' + processParams.libraries: ''} \
			${(processParams.containsKey('featureRef')) ? '--feature-ref ' + processParams.featureRef: ''} \
			${(processParams.containsKey('expectCells')) ? '--expect-cells ' + processParams.expectCells: ''} \
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

}
