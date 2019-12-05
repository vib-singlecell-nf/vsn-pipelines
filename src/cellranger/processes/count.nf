nextflow.preview.dsl=2

process SC__CELLRANGER__COUNT {

  	publishDir "${params.global.outdir}/counts", mode: 'link', overwrite: true
  	container params.sc.cellranger.container

  	input:
    file(transcriptome)
    tuple val(sampleId), file(fastqs)

  	output:
    tuple(_sampleName, file("${_sampleName}/outs") )

  	script:
	"""
	cellranger count \
		--id=${sampleId} \
		--sample=${sampleId} \
		--fastqs=${fastqs.join(",")} \
		--transcriptome=${transcriptome} \
		${(params.sc.cellranger.count.containsKey('libraries')) ? '--libraries ' + params.sc.cellranger.count.libraries: ''} \
		${(params.sc.cellranger.count.containsKey('featureRef')) ? '--feature-ref ' + params.sc.cellranger.count.featureRef: ''} \
		${(params.sc.cellranger.count.containsKey('expectCells')) ? '--expect-cells ' + params.sc.cellranger.count.expectCells: ''} \
		${(params.sc.cellranger.count.containsKey('forceCells')) ? '--force-cells ' + params.sc.cellranger.count.forceCells: ''} \
		${(params.sc.cellranger.count.containsKey('nosecondary')) ? '--nosecondary ' + params.sc.cellranger.count.nosecondary: ''} \
		${(params.sc.cellranger.count.containsKey('noLibraries')) ? '--no-libraries ' + params.sc.cellranger.count.noLibraries: ''} \
		${(params.sc.cellranger.count.containsKey('chemistry')) ? '--chemistry ' + params.sc.cellranger.count.chemistry: ''} \
		${(params.sc.cellranger.count.containsKey('r1Length')) ? '--r1-length ' + params.sc.cellranger.count.r1Length: ''} \
		${(params.sc.cellranger.count.containsKey('r2Length')) ? '--r2-length ' + params.sc.cellranger.count.r2Length: ''} \
		${(params.sc.cellranger.count.containsKey('lanes')) ? '--lanes ' + params.sc.cellranger.count.lanes: ''} \
		${(params.sc.cellranger.count.containsKey('localCores')) ? '--localcores ' + params.sc.cellranger.count.localCores: ''} \
		${(params.sc.cellranger.count.containsKey('localMem')) ? '--localmem ' + params.sc.cellranger.count.localMem: ''} \
		${(params.sc.cellranger.count.containsKey('indicies')) ? '--indicies ' + params.sc.cellranger.count.indicies: ''} 
	"""
}

