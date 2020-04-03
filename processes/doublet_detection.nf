nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scrublet/bin/" : ""

include '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCRUBLET__DOUBLET_DETECTION_PARAMS {

	Script env = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String off = null;

	void setEnv(env) {
		this.env = env
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	String getNPrinCompsAsArgument(nPrinComps) {
		// Check if nComps is both dynamically and if statically set
		if(!this.env.isParamNull(nPrinComps) && this.configParams.containsKey('nPrinComps'))
			throw new Exception("SC__SCRUBLET__DOUBLET_DETECTION: nPrinComps is both statically (" + this.configParams["nPrinComps"] + ") and dynamically (" + nPrinComps + ") set. Choose one.")
		if(!this.env.isParamNull(nPrinComps))
			return '--n-prin-comps ' + nPrinComps.replaceAll("\n","")
		return this.configParams.containsKey('nPrinComps') ? '--n-prin-comps ' + this.configParams.nPrinComps: ''
	}

}

def SC__SCRUBLET__DOUBLET_DETECTION_PARAMS(params) {
	return (new SC__SCRUBLET__DOUBLET_DETECTION_PARAMS(params))
}

process SC__SCRUBLET__DOUBLET_DETECTION {

	container params.sc.scrublet.container
	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

	input:
		tuple \
			val(sampleId), \
			path(adataRaw), \
            path(adataHvg), \
			val(stashedParams), \
			val(nPrinComps)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCRUBLET__DOUBLET_DETECTION.${processParams.off}"), \
			val(stashedParams), \
			val(nPrinComps)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scrublet.scrub_doublets)
		processParams = sampleParams.local
		def _processParams = new SC__SCANPY__DIM_REDUCTION_PARAMS()
		_processParams.setEnv(this)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}sc_doublet_detection.py \
            ${(processParams.containsKey('syntheticDoubletUmiSubsampling')) ? '--synthetic-doublet-umi-subsampling ' + processParams.syntheticDoubletUmiSubsampling : ''} \
            ${(processParams.containsKey('minCounts')) ? '--min-counts ' + processParams.minCounts : ''} \
            ${(processParams.containsKey('minCells')) ? '--min-cells ' + processParams.minCells : ''} \
            ${(processParams.containsKey('minGeneVariabilityPctl')) ? '--min-gene-variability-pctl ' + processParams.minGeneVariabilityPctl : ''} \
            ${(processParams.containsKey('logTransform')) ? '--log-transform ' + processParams.logTransform : ''} \
            ${(processParams.containsKey('meanCenter')) ? '--mean-center ' + processParams.meanCenter : ''} \
            ${(processParams.containsKey('normalizeVariance')) ? '--normalize-variance ' + processParams.normalizeVariance : ''} \
            ${_processParams.getNPrinCompsAsArgument(nPrinComps)} \
            ${(processParams.containsKey('technology')) ? '--technology ' + processParams.technology : ''} \
			$adataRaw \
            $dataHvg \
			"${sampleId}.SC__SCRUBLET__DOUBLET_DETECTION.${processParams.off}"
		"""

}
