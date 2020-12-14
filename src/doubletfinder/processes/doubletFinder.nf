nextflow.preview.dsl=2

scriptDir = "${workflow.projectDir}/src/DoubletFinder/bin"

process DOUBLETFINDER__RUN {
	publishDir "${params.global.outdir}/${samplename}", mode: 'symlink', pattern: "Plots/**"
	container params.DoubletFinder.container
	input:
	tuple val(samplename), file(sobj)
	output:
	tuple val(samplename), file("${samplename}.DOUBLETFINDER.rds")
	file("Plots/**")
	script:
	def sampleParams = params.parseConfig(samplename, params.global, params.DoubletFinder)
		processParams = sampleParams.local
	"""
	Rscript ${scriptDir}/doubletFinder.R --seuratObj ${sobj} \
	--output ${samplename}.DOUBLETFINDER.rds \
	${(processParams.minPCT == null) ? '' : '--minPCT ' + processParams.minPCT} \
	${(processParams.maxPCT == null) ? '' : '--maxPCT ' + processParams.maxPCT} \
	${(processParams.pN == null) ? '' : '--pN ' + processParams.pN} \
	${(processParams.dimsToUse == null) ? '' : '--dimsToUse ' + processParams.dimsToUse} \
	${(processParams.cores == null) ? '' : '--cores ' + processParams.cores}
	"""
}
